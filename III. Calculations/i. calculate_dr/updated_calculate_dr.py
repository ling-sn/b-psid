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

def create_tsv(df: pd.Dataframe, counts: pd.DataFrame, key: dict, folder_name: str, output_tsv_name: str):
   ## Calculate observed deletion rates
   coverage_list = [col for col in counts.columns
                  if re.search("(A_|C_|G_|T_|Deletions_).*", col)]
   counts["TotalCoverage"] = counts[coverage_list].sum(axis = 1)
   total_del = counts[["Deletions_fwd", "Deletions_rev"]].sum(axis = 1)
   counts["DeletionRate"] = total_del / counts["TotalCoverage"]
   
   ## Calculate real deletion rates
   df_draft = pd.merge(df, counts, how = "left", 
                     on = ["Chrom", "GenomicModBase"]).dropna()
   num = df_draft["fit_b"] - df_draft["DeletionRate"]
   denom = (df_draft["fit_c"] * (df_draft["fit_b"] + df_draft["fit_s"] -
            df_draft["fit_s"] * df_draft["DeletionRate"] - 1))
   df_draft["RealRate"] = num/denom
   
   ## Rename columns
   df_draft = df_draft.rename(columns = {"A": key["A"], 
                                       "C": key["C"], 
                                       "G": key["G"], 
                                       "T": key["T"],
                                       "Deletions": key["Deletions"],
                                       "TotalCoverage": key["TotalCoverage"],
                                       "DeletionRate": key["DeletionRate"], 
                                       "RealRate": key["RealRate"]})

   ## Apply filter conditions based on filename
   """
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
   df_final = df_draft[df_draft[cov_pattern].ge(10)]

   ## Only output files if WT or 7KO
   if re.match(fr"(WT|7KO).*", str(folder_name)):
      if re.match(fr"WT.*", str(folder_name)):
         if "_BS" in dr_pattern:
               df_final = df_final[df_final[dr_pattern].ge(0.3)]
         else: 
               df_final = df_final[df_final[dr_pattern].le(0.3)]

      if re.match(fr"7KO.*", str(folder_name)):
         if "_BS" in dr_pattern:
               df_final = df_final[df_final[dr_pattern].le(0.3)]

      ## Save as .tsv output
      df_final = df_final.drop_duplicates().sort_values(by = dr_pattern, ascending = False)
      df_final.to_csv(output_tsv_name, sep = "\t", index = False)

def count_base(unuar_sites, input_bam_name, fasta_dir, results):
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
      for chrom in unuar_sites:
         ## Values in unuar_sites already stored in set (faster)
         positions = unuar_sites[chrom]
         base_keys = [
               "A", "C_fwd", "G_fwd", "T_fwd",
               "A_rev", "C_rev", "G_rev", "T_rev"
         ]

         for stats in pysamstats.stat_variation(
               bamfile,
               fastafile,
               chrom,
               start = min(unuar_sites[chrom]),
               end = max(unuar_sites[chrom]) + 1,
               one_based = True,
               truncate = True,
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

def process_bam(bam, processed_folder: str, fwd_unuar, rev_unuar,
               fasta_dir, df, key, folder_name):
   ## Turn string from list back into filepath
   input_bam_name = Path(bam)
   output_tsv_name = processed_folder/f"{input_bam_name.stem}.tsv"
   
   ## Count A, C, G, T and deletions @ each UNUAR site
   results = []
   for unuar_sites in [fwd_unuar, rev_unuar]:
      count_base(unuar_sites, input_bam_name, fasta_dir, results)

   ## Convert results dict -> df
   counts = pd.DataFrame(results)
   return counts

def make_key(folder_name: str, base_key: str):
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
   prefix = rep + "_"
   
   for sample in ["BS", "NBS"]:
      if f"-{sample}" in str(folder_name):
         suffix = "_" + sample
         break
   
   return prefix + base_key + suffix

def group_by_chrom(df: pd.DataFrame) -> dict:
   genome_coord = df[["Chrom", "GenomicModBase"]]
   grouped = genome_coord.groupby("Chrom")
   unuar_sites = grouped["GenomicModBase"].agg(set).to_dict()
   return unuar_sites

def get_sample_group(folder_name: str):
   """
   PURPOSE:
   Given input folder names, extract the group name
   by returning the first capture group in RegEx.
   ---
   EXAMPLE: 
   'KEH-Rep1-7KO-HEK293T-Cyto-BS' -> '7KO-Cyto-BS'
   """
   try:
      match = re.match(r"^.*(7KO|7LKO|WT)(?:.*)(-.*-.*)", folder_name)
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
   group_name = get_sample_group(folder_name)
   
   fwd = pd.read_csv(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools"
                     "/B-PsiD_UNUAR_motif_sites_mRNA_hg38_fwd.tsv").expanduser(), 
                     sep = "\t")
   rev = pd.read_csv(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools"
                     "/B-PsiD_UNUAR_motif_sites_mRNA_hg38_rev.tsv").expanduser(), 
                     sep = "\t")

   ## Group sites by chromosome, then convert to dict
   fwd_unuar = group_by_chrom(fwd)
   rev_unuar = group_by_chrom(rev)
   
   try:
      if input_dir.is_dir():
         processed_folder = current_path/"calculations"/group_name
         processed_folder.mkdir(exist_ok = True, parents = True)
         
         key = {base_key: make_key(folder_name, base_key) 
                  for base_key in ["A", "C", "G", "T",
                                 "Deletions",
                                 "TotalCoverage",
                                 "DeletionRate", 
                                 "RealRate"]}
         
         for bam_type in ["fwd", "rev"]:
               bam = input_dir/f"{bam_type}.bam"
               process_bam(bam, processed_folder, fwd_unuar, rev_unuar, 
                           fasta_dir, df, key, folder_name)

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