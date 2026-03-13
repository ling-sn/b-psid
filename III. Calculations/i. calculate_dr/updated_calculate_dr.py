## Use RNA-STAR conda environment
from pathlib import Path
import traceback
import argparse
import pandas as pd
import numpy as np
import pysam  # type: ignore
from pysam import FastaFile  # type: ignore
import pysamstats # type: ignore
import re
from multiprocessing import Pool

## Disable .loc indexing warning
pd.options.mode.chained_assignment = None

class BaseDelCounter:
   def create_tsv(self, df_calc: pd.DataFrame, key: dict, 
                  folder_name: str, output_tsv_name: str):
      """
      PURPOSE:
      Apply filter conditions based on filename
      ---
      WT:
      * BS files must have DeletionRate values of >= 0.3
      * NBS files must have DeletionRate values of <= 0.3
      ---
      Mutation (PUS7KO):
      * BS files must have DeletionRate values of <= 0.3
      """
      dr_pattern = key["DeletionRate"]
      cov_pattern = key["TotalCoverage"]
      # rr_pattern = key["RealRate"]

      ## Keep only RealRate >= 0.3
      # kept_rr = df_draft[df_draft[rr_pattern].ge(0.3)]

      ## Keep only rows where coverage >= 10)
      df_final = (
                     df_calc[df_calc[cov_pattern].ge(10)]
                     .sort_values(by = dr_pattern, ascending = False)
                 )

      ## Only filter files if WT or 7KO
      if re.search(fr"(WT|7KO).*", str(folder_name)):
         if re.search(fr"WT.*", str(folder_name)):
            if "_BS" in dr_pattern:
               df_final = df_final[df_final[dr_pattern].ge(0.3)]
            else: 
               df_final = df_final[df_final[dr_pattern].le(0.3)]

         if re.search(fr"7KO.*", str(folder_name)):
            if "_BS" in dr_pattern:
               df_final = df_final[df_final[dr_pattern].le(0.3)]
      
      df_final.to_csv(output_tsv_name, sep = "\t", index = False)

   def calc_rate(self, df_original: pd.DataFrame, df_count: pd.DataFrame, 
                 key: dict) -> pd.DataFrame:
      ## Calculate observed deletion rates
      coverage_list = [col for col in df_count.columns
                       if re.match("(A|C|G|T|Deletions)$", col)]
      df_count["TotalCoverage"] = df_count[coverage_list].sum(axis = 1)
      df_count["DeletionRate"] = df_count["Deletions"] / df_count["TotalCoverage"]
      
      ## Calculate real deletion rates
      df_calc = pd.merge(df_original, df_count, how = "left", 
                         on = ["Chrom", "GenomicModBase"]).dropna()
      num = df_calc["fit_b"] - df_calc["DeletionRate"]
      denom = (df_calc["fit_c"] * (df_calc["fit_b"] + df_calc["fit_s"] -
               df_calc["fit_s"] * df_calc["DeletionRate"] - 1))
      df_calc["RealRate"] = num/denom
      
      ## Rename columns
      df_calc = df_calc.rename(columns = {"A": key["A"], 
                                          "C": key["C"], 
                                          "G": key["G"], 
                                          "T": key["T"],
                                          "Deletions": key["Deletions"],
                                          "TotalCoverage": key["TotalCoverage"],
                                          "DeletionRate": key["DeletionRate"], 
                                          "RealRate": key["RealRate"]})
      return df_calc

   def count_base(self, unuar_dict: dict, input_bam_name: Path,
                  fasta_dir: Path, results: list):
      """
      PURPOSE:
      Counts number of bases and deletions for each UNUAR site 
      in each chromosome.
      ---
      VARIABLES:
      * chrom: Value in "Chrom" column (e.g., NW_018654708.1)
      * mod_base: Value in "GenomicModBase" column (e.g., 373)
      * **base_ct: Unpacks base_ct dict
      ---
      PYSAMSTATS:
      Description
      * pysamstats.stat_variation() returns an iterator for dict objects,
      where each dict holds data for different sites
      Arguments
      * start, end: 1-based genomic coordinates (half-open interval)
      * truncate: Ensures only positions within range are returned
      * no_dup: Don't count reads marked as duplicates
      """
      bamfile = pysam.AlignmentFile(input_bam_name, "rb")
      fastafile = FastaFile(fasta_dir) 

      try:
         for chrom in unuar_dict:
            ## Values in unuar_sites already stored in set (faster)
            positions = unuar_dict[chrom]
            base_keys = ["A", "C", "G", "T"]

            for stats in pysamstats.stat_variation(
               bamfile,
               fastafile,
               chrom,
               start = min(unuar_dict[chrom]),
               end = max(unuar_dict[chrom]) + 1,
               one_based = True,
               truncate = False,
               no_dup = True
            ):
               pos = stats["pos"]

               if pos not in positions:
                  continue

               results.append(
                  {
                     "Chrom": chrom,
                     "GenomicModBase": stats["pos"],
                     **{base: stats[base] for base in base_keys},
                     "Deletions": stats["deletions"]
                  }
               )
            
         bamfile.close()
      except Exception as e:
         print(f"Failed to count bases/deletions in UNUAR sites"
               f"within chromosome {chrom}: {e}")
         traceback.print_exc()
         raise

   def process_bam(self, bam: Path, unuar_dict: dict, fasta_dir: Path) -> pd.DataFrame:
      ## Turn string from list back into filepath
      input_bam_name = Path(bam)
      
      ## Count A, C, G, T and deletions @ each UNUAR site
      results = []
      self.count_base(unuar_dict, input_bam_name, fasta_dir, results)

      ## Convert results dict -> df
      counts = pd.DataFrame(results)
      return counts

class PrepData:
   def make_key(self, folder_name: str, base_key: str) -> str:
      """
      PURPOSE:
      Modifies names of dictionary keys based on the Rep # (detected via RegEx)
      and Sample Type (BS, NBS) in a given folder name.
      ---
      EXAMPLE:
      Inputs:
         folder_name = KEH-Rep1-7KO-HEK293T-Cyto-BS
         base_key = Deletions
      Output: 
         Rep1_Deletions_BS
      """
      rep = re.search(r"Rep\d+", str(folder_name))
      prefix = rep.group(0) + "_"
      
      for sample in ["BS", "NBS"]:
         if f"-{sample}" in str(folder_name):
            suffix = "_" + sample
            break
      
      return prefix + base_key + suffix

   def group_by_chrom(self, df: pd.DataFrame) -> dict:
      genome_coord = df[["Chrom", "GenomicModBase"]]
      grouped = genome_coord.groupby("Chrom")
      unuar_sites = grouped["GenomicModBase"].agg(set).to_dict()
      return unuar_sites

   def get_sample_group(self, folder_name: str) -> str:
      """
      PURPOSE:
      Given input folder names, extract the group name
      by returning the first capture group in RegEx.
      ---
      EXAMPLE: 
      'KEH-Rep1-7KO-HEK293T-Cyto-BS' -> '7KO-Cyto'
      """
      try:
         match = re.match(r"^.*(7KO|7LKO|WT)(?:.*)(-.*)(?:-.*)", folder_name)
      except Exception as e:
         print(f"Failed to RegEx match input folder to group: {e}")
         traceback.print_exc()
         raise
      return match.group(1) + match.group(2)

def main(folder_name: str, fasta: str):
   """
   PURPOSE: 
   Opens .bam in folder and runs calculations
   """
   current_path = Path.cwd()
   input_dir = current_path/"realignments"/folder_name
   fasta_dir = Path(fasta).expanduser()
   
   fwd = (
            pd.read_csv(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools"
                        "/B-PsiD_UNUAR_motif_sites_mRNA_hg38_fwd.tsv").expanduser(), 
                        sep = "\t")
            .drop_duplicates(subset = ["Chrom", "GenomicModBase"], keep = "first")
         )
   rev = (
            pd.read_csv(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools"
                        "/B-PsiD_UNUAR_motif_sites_mRNA_hg38_rev.tsv").expanduser(), 
                        sep = "\t")
            .drop_duplicates(subset = ["Chrom", "GenomicModBase"], keep = "first")
         )
   
   ## Initialize classes
   prep = PrepData()
   counter = BaseDelCounter()
   
   ## Define full dataframe containing sites on both +/- strand
   df_original = pd.concat([fwd, rev], ignore_index = True)

   ## Group sites by chromosome, then convert to dict
   fwd_unuar = prep.group_by_chrom(fwd)
   rev_unuar = prep.group_by_chrom(rev)
   
   try:
      if input_dir.is_dir():
         group_name = prep.get_sample_group(folder_name)
         processed_folder = current_path/"calculations"/group_name
         processed_folder.mkdir(exist_ok = True, parents = True)
         
         key = {base_key: prep.make_key(folder_name, base_key)
                          for base_key in ["A", "C", "G", "T",
                                           "Deletions",
                                           "TotalCoverage",
                                           "DeletionRate",
                                           "RealRate"]}
         
         types = ["fwd", "rev"]
         unuar_dicts = [fwd_unuar, rev_unuar]
         all_counts = []
         
         ## Obtain base-del counts for fwd/rev separately
         for type, unuar_dict in zip(types, unuar_dicts):
            bam = input_dir/f"{type}.bam"
            counts = counter.process_bam(bam, unuar_dict, fasta_dir)
            all_counts.append(counts)
         
         ## Concat fwd/rev count dataframes
         df_count = pd.concat(all_counts, ignore_index = True)
         
         """
         Process dataframe:
         1. Calculate DeletionRate and RealRate
         2. Filter out sites based on certain criteria
            and output final dataframe as TSV
         """
         ## Calculate rates and rename columns
         df_calc = counter.calc_rate(df_original, df_count, key)
         
         ## Filter out sites
         output_tsv_name = processed_folder/f"{folder_name}.tsv"
         counter.create_tsv(df_calc, key, folder_name, output_tsv_name)

   except Exception as e:
      print("Failed to calculate observed & real deletion rates in "
            f"{folder_name} and save as .tsv: {e}")
      traceback.print_exc()
      raise

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description = "Calculates observed and real deletion rates" 
                                                  "for every UNUAR site in a BAM file.")
   parser.add_argument("--folder_name", help = "Name of realignments folder", required = True)
   parser.add_argument("--fasta", help = "Directory to FASTA file", required = True)
   args = parser.parse_args()

   print("Calculating deletion rates...")
   main(args.folder_name, args.fasta)
   print("Process finished.")