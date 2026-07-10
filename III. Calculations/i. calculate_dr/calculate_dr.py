"""
@author: ling-sn
Last updated: April 2026
Use B-PSID conda environment
"""
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
      ## Initialize class
      prep = PrepData()
      
      ## Collect all associated transcripts for each (Chrom, GenomicModBase) pair
      all_transcripts = prep.create_agg_dict(df_calc, 
                                             ["Chrom", "GenomicModBase"], 
                                             "TranscriptID")
      updated_transcripts = [(key[0], key[1], ", ".join(str(x) for x in value)) 
                             for key, value in all_transcripts.items()]
      df_transcripts = pd.DataFrame(updated_transcripts, 
                                    columns = ["Chrom", "GenomicModBase", 
                                               "AllAssocTranscripts"])
      
      ## Keep only RealRate >= 0.3
      # rr_pattern = key["RealRate"]
      # kept_rr = df_draft[df_draft[rr_pattern].ge(0.3)]

      ## Sort rows by deletion rate
      dr_pattern = key["DeletionRate"]
      df_draft = df_calc.sort_values(by = dr_pattern, ascending = False)

      # ## Only filter files if WT or 7KO
      # if (
      #       re.search(fr"WT.*-NBS$", str(folder_name))
      #       or
      #       re.search(fr"7KO.*-BS$", str(folder_name))
      # ):
      #    df_draft = df_draft[df_draft[dr_pattern].le(0.3)]
      
      ## Output final TSV
      df_final = (
         pd.merge(df_draft, df_transcripts,
                  how = "left",
                  on = ["Chrom", "GenomicModBase"])
      )
      df_final.insert(17, "AllAssocTranscripts", df_final.pop("AllAssocTranscripts"))
      
      df_final.to_csv(output_tsv_name, sep = "\t", index = False)

   def calc_rate(self, df_original: pd.DataFrame, df_count: pd.DataFrame, 
                 key: dict) -> pd.DataFrame:
      ## Calculate observed deletion rates
      coverage_list = [col for col in df_count.columns
                       if re.match("(A|C|G|T|Deletions)$", col)]
      df_count["TotalCoverage"] = df_count[coverage_list].sum(axis = 1)
      df_count["DeletionRate"] = df_count["Deletions"] / df_count["TotalCoverage"]
      
      ## Calculate real deletion rates
      count_cols = (df_count.columns.tolist())[3:]
      df_calc = (
         pd.merge(df_original, df_count, how = "left", 
                  on = ["Chrom", "GenomicModBase"])
         .dropna(subset = count_cols, how = "all")
      )
      df_calc["Deletions"] = df_calc["Deletions"].fillna(0) 
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

   def count_single_dels(self, chrom: str, 
                         pos: int, 
                         bamfile: pysam.AlignmentFile) -> int:
      """
      PURPOSE:
      1. Obtain all reads at specific genomic coordinate.
      2. If read has 1D (single-nt deletion) in its CIGAR 
         string, append deletion count.
         -> We only consider deletions of the central U in
            the UNUAR motif.
      3. Skip all reads that (i) don't have deletions or
         (ii) have multi-nt deletions.
      ---
      NOTES:
      * pos is 1-based, but base pysam is 0-based.
      * bamfile.fetch() obtains ALL reads overlapping the 
        specified position.
      * We can extract deletion info from CIGAR strings by 
        searching for digits followed by a D.
        -> EXAMPLE: Given CIGAR string 4M1I3M4D9M, we can
           extract 4D which indicates there are 4 deletions.
      * If deletions == 0, return None so null rows can be
        properly dropped later (outside this func).
      """
      deletions = 0

      for read in bamfile.fetch(chrom, pos - 1, pos):
         if read.is_secondary or read.is_duplicate:
            continue
         cigar = read.cigarstring
         del_info = re.search(r"(\d+)D", str(cigar))
         if del_info and int(del_info.group(1)) == 1:
            deletions += 1
         
      if (deletions == 0):
         return None
      
      return deletions

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
                     # "Deletions": self.count_single_dels(chrom, pos, bamfile)
                  }
               )
            
         bamfile.close()
      except Exception as e:
         print(f"Failed to count bases/deletions in UNUAR sites"
               f" within chromosome {chrom}: {e}")
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

   def create_agg_dict(self, df: pd.DataFrame, 
                       group1_col: str|list, group2_col: str) -> dict:
      if group2_col == "GenomicModBase":
         df = (df[[group1_col, group2_col]])
      genome_coord = df.groupby(group1_col)
      unuar_sites = genome_coord[group2_col].agg(set).to_dict()
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
         match = re.match(r"(7KO|7LKO|WT)(?:.*)(Cyto|Nuc)", folder_name)
      except Exception as e:
         print(f"Failed to RegEx match input folder to group: {e}")
         traceback.print_exc()
         raise
      return match.group(1) + "-" + match.group(2)

def main(folder_name: str, fasta: str, rep_index: int):
   """
   PURPOSE: 
   Opens .bam in folder and runs calculations
   """
   current_path = Path.cwd()
   input_dir = current_path/"realignments"/folder_name
   fasta_dir = Path(fasta).expanduser()
   
   fwd = (
      pd.read_csv(Path("/nfs/turbo/umms-RNAlabDATA/Software/B-PsiD_tools"
                  "/B-PsiD_UNUAR_motif_sites_mRNA_hg38_fwd.tsv"), 
                  sep = "\t")
   )
   rev = (
      pd.read_csv(Path("/nfs/turbo/umms-RNAlabDATA/Software/B-PsiD_tools"
                  "/B-PsiD_UNUAR_motif_sites_mRNA_hg38_rev.tsv"), 
                  sep = "\t")
   )
   
   ## Initialize classes
   prep = PrepData()
   counter = BaseDelCounter()
   
   ## Define full dataframe containing sites on both +/- strand
   df_original = pd.concat([fwd, rev], ignore_index = True)

   ## Group sites by chromosome, then convert to dict
   fwd_unuar = prep.create_agg_dict(fwd, "Chrom", "GenomicModBase")
   rev_unuar = prep.create_agg_dict(rev, "Chrom", "GenomicModBase")
   
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
            bam = input_dir/f"{type}.sorted.bam"
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
   parser.add_argument("--folder_name", help = "Name of realignments folder", 
                       required = True)
   parser.add_argument("--fasta", help = "Directory to FASTA file", 
                       required = True)
   parser.add_argument("--rep_index", help = "Index number of replicate in file naming structure", 
                       required = True)
   args = parser.parse_args()

   print("Calculating deletion rates...")
   main(args.folder_name, args.fasta, args.rep_index)
   print("Process finished.")
