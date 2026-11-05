"""
@author: ling-sn
Last updated: November 2026
Use B-PSID conda environment
MAIN CHANGE: Use rate ratio distribution instead of Fisher's Exact Test
"""
from pathlib import Path
import traceback
import pandas as pd
import numpy as np
import re
import requests
import concurrent.futures
from bs4 import BeautifulSoup
pd.options.mode.chained_assignment = None

class MergeTSV:
   def obtain_gene(self, url: str):
      page = requests.get(url)
      soup = BeautifulSoup(page.content, "html.parser")
      gene = soup.find("span", class_= "gn")
      return gene.text
   
   def merge_all(self, matching_name, pvals_tsv, final_dir):
      matches = [tsv for tsv in pvals_tsv if re.search(matching_name, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
      df_merged = self.iteratively_merge(df_list, "inner")
      
      ## Add NCBI links
      df_merged["TranscriptID"] = (
         df_merged["TranscriptID"]
         .str.replace("rna-", "", regex = False)
      )
      df_merged["NCBILink"] = "https://www.ncbi.nlm.nih.gov/gene/?term=" + df_merged["TranscriptID"]
      df_merged.insert(18, "NCBILink", df_merged.pop("NCBILink"))
      
      ## Obtain gene name
      MAX_THREADS = 20
      threads = min(MAX_THREADS, len(df_merged["NCBILink"]))
      with concurrent.futures.ThreadPoolExecutor(max_workers = threads) as executor:
         df_merged["Gene"] = list(executor.map(self.obtain_gene, df_merged["NCBILink"]))
      df_merged.insert(0, "Gene", df_merged.pop("Gene"))
      
      ## 1. Sort rows by TotalAvgDeletionRate for WT-BS sample (descending order)
      ## 2. Move all important rows (TotalAvg, Std, PvaluePass) to end of dataframe
      df_merged.sort_values(
         by = "WT_TotalAvgDeletionRate_BS",
         ascending = False, 
         inplace = True
      )
      important_cols = [
         col for col in df_merged.columns 
         if re.search(
            "(TotalAvgDeletionRate|StdDeletionRate|Pvalue_Pass)", col
         )
      ]
      diff_cols = list(df_merged.columns.difference(important_cols, sort = False))
      new_col_order = diff_cols + important_cols
      df_merged = df_merged[new_col_order]
      
      ## 1. Create output name
      ##    e.g., 7KO-Cyto-Pvals + 7LKO-Cyto-Pvals + WT-Cyto-Pvals -> Cyto
      ## 2. Save merged dataframe as TSV
      output_name = (matches[0].stem).split("-")[1]
      merged_dir = final_dir/f"{output_name}.tsv"
      df_merged.to_csv(merged_dir, sep = "\t", index = False)
      
   def iteratively_merge(self, list_of_dfs: list, merge_type: str):
      df1_colnames = list_of_dfs[0].columns.tolist()
      selected_colnames = df1_colnames[0:18]
      merged = list_of_dfs[0]

      for df in list_of_dfs[1:]:
         if not df.empty:
            merged = pd.merge(
               merged, df,
               on = selected_colnames,
               how = merge_type
            )
   
      return merged

   def merge_reps(
      self, suffix: str, 
      tsv_list: list
   ) -> pd.DataFrame:
      """
      1. Search TSVs for matching suffix in filename
      2. Put them in list
      3. Read in as pandas dataframes
      4. Iteratively merge dataframes with outer join (do not drop rows yet)
      5. Drop rows that are completely identical
      """
      matches = [tsv for tsv in tsv_list if re.search(suffix, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
      merged = self.iteratively_merge(df_list, "outer").drop_duplicates()
      return merged
      
class FilterTSV:
   def calc_rr(self, df, sample_type):
      """
      PURPOSE:
      Calculate rate ratios for each UNUAR site, 
      and return the updated dataframe.
      """
      ## Collect lists of corresponding counts (Deletions/TotalCoverage
      ## for BS/NBS samples) and store them in a dictionary
      all_cols = {}
      
      for count_type, treatment in zip(
         ["Deletions", "TotalCoverage"], ["BS", "NBS"]
      ):
         list_name = (f"{treatment}_{count_type}_cols").lower()
         selected_cols = [
            col for col in 
            re.search(f"^{sample_type}_Rep(\d+)_{count_type}_{treatment}$")
         ] 
         all_cols[list_name] = selected_cols
      
      ## Create required columns
      # (i) BS
      df["num_bs_dels"] = df[all_cols["bs_deletions_cols"]].sum(axis = 1)
      df["total_bs"] = df[all_cols["bs_totalcoverage_cols"]].sum(axis = 1)
      df["bs_ratio"] = df["num_bs_dels"] / df["total_bs"]
      
      # (ii) NBS
      df["num_nbs_dels"] = df[all_cols["nbs_deletions_cols"]].sum(axis = 1)
      df["total_nbs"] = df[all_cols["nbs_totalcoverage_cols"]].sum(axis = 1)
      df["nbs_ratio"] = df["num_nbs_dels"] / df["total_nbs"]
      
      ## Calculate rate ratios
      df["RateRatio"] = df["bs_ratio"] / df["nbs_ratio"]
      
      if (sample_type == "7KO"):
         df.sort_values(
            by = "RateRatio",
            ascending = True, 
            inplace = True
         )
         return df
      
      df.sort_values(
         by = "RateRatio",
         ascending = False, 
         inplace = True
      )
      
      return df
      
def main():
   """
   PURPOSE:
   Filters .tsv files in grouped folders
   """
   current_path = Path.cwd()
   start_dir = current_path/"calculations"

   ## Initialize classes
   mergetsv = MergeTSV()
   filtertsv = FilterTSV()

   try: 
      ## Create output folder directories
      processed_folder = current_path/"merged"
      final_dir = processed_folder/"final_outputs"
      for dirname in [processed_folder, final_dir]:
         dirname.mkdir(exist_ok = True, parents = True)

      all_merged = {}
      for subfolder in start_dir.iterdir():
         if subfolder.is_dir():
            ## Collect paths of .tsv files and put in list
            tsv_list = sorted(
               subfolder.glob("*.tsv"),
               # Order by replicate integer
               key = lambda x: int(re.search(r"Rep(\d+)", x.name).group(1))
            )

            ## Merge by replicate (Rep 1-3)
            merged_bs_nbs = []
            for suffix in ["-BS", "-NBS"]:
               merged = mergetsv.merge_reps(suffix, tsv_list, subfolder)
               merged_bs_nbs.append(merged)
            
            ## Merge by BS-NBS
            """
            EXAMPLE: 
            WT-Cyto-BS + WT-Cyto-NBS -> WT-Cyto
            """
            df1 = merged_bs_nbs[0]
            df2 = merged_bs_nbs[1]
            selected_colnames = (df1.columns.tolist())[0:18]
            
            if not (df1.empty and df2.empty):
               df_merged = pd.merge(df1, df2, on = selected_colnames, how = "outer")
            elif df1.empty:
               df_merged = df2
            else:
               df_merged = df1
            
            """
            EXAMPLE:
            subfolder.name = 7KO-Cyto
            """
            all_merged[subfolder.name] = df_merged
            
      ## Calculate p-value for each merged subfolder dataframe
      for location in ["Cyto", "Nuc"]:
         # Open up 7KO dataframe and calculate rate ratios
         knockout_df = filtertsv.calc_rr(
            all_merged["7KO-" + location], "7KO"
         )
         knockout_df 
         
         # Open up WT datafarme
         wt_df = all_merged["WT-" + location]

      ## After p-value calculations, create final merged ouputs
      pvals_tsv = list(pvals_dir.glob("*.tsv"))

      for matching_name in ["-Cyto-Pvals", "-Nuc-Pvals"]:
         mergetsv.merge_all(matching_name, pvals_tsv, final_dir)

   except Exception as e:
      print(f"Failed to create merged .tsv file: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")