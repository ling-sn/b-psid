## use RNA-STAR conda environment
from pathlib import Path
import traceback
import argparse
import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact

class FilterTSV:
   def merge_WT_7KO(self, matching_name, pvals_tsv, final_dir):
      matches = [tsv for tsv in pvals_tsv if re.search(matching_name, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
      
      """
      1. Ensure 7KO is merged with WT, so WT columns appear first
      2. If either dataframe is not empty, then merge w/ outer join
      3. No need to iteratively merge because there are only 2 files
      """  
      first_cols = df_list[0].columns.tolist()
      
      if any(re.search("WT", col) for col in first_cols):
         df1 = df_list[0]
         df2 = df_list[1]
      else:
         df1 = df_list[1]
         df2 = df_list[0]
            
      """
      Merge df1 & df2
      """
      selected_colnames = (df1.columns.tolist())[0:18]
      if not df1.empty and not df2.empty:
         df_merged = pd.merge(df1, df2, on = selected_colnames, how = "outer")
      elif df1.empty:
         df_merged = df2
      else:
         df_merged = df1
      
      """
      Add NCBI links
      """
      df_merged["TranscriptID"] = df_merged["TranscriptID"].str.replace("rna-", "", regex = False)
      df_merged["NCBILink"] = "https://www.ncbi.nlm.nih.gov/gene/?term=" + df_merged["TranscriptID"]
      df_merged.insert(18, "NCBILink", df_merged.pop("NCBILink"))
      
      """
      1. Create output name
         e.g., 7KO-Cyto-Pvals + WT-Cyto-Pvals -> Cyto
      2. Save merged dataframe as TSV
      """
      output_name = (matches[0].stem).split("-")[1]
      merged_dir = final_dir/f"{output_name}.tsv"
      df_merged.to_csv(merged_dir, sep = "\t", index = False)
  
   def match_cols(self, merged_colnames: list, rep: str, basedel: str) -> list:
      col_list = []
      for sample in ["BS", "NBS"]:
         match = next(col for col in merged_colnames
                      if re.search(f"{rep}_{basedel}_{sample}", col))
         col_list.append(match)
      return col_list
  
   def fisher_test(self, df_merged: pd.DataFrame,
                   wt_7ko_7lko: str) -> pd.DataFrame:
      try:
         ## Rename columns with WT/7KO prefix (excluding last 2, which already contain it)
         selected_cols = (df_merged.columns.tolist())[18:-2]   
         for old_name in selected_cols:
            new_name = wt_7ko_7lko + "_" + old_name
            df_merged.rename(columns = {old_name: new_name}, inplace = True)
         
         ## Collect column names and replicates in dataframe
         merged_colnames = df_merged.columns.tolist()
         rep_list = sorted(
            set([re.search(r"(Rep\d+)", col).group(1) 
                  for col in merged_colnames 
                  if re.search(r"(Rep\d+)", col)]), 
            key = lambda x: int(re.search(r"Rep(\d+)", x).group(1))
         )
         
         ## Split dataframe by strand (+/-)
         fwd = df_merged[df_merged["Strand"] == "+"]
         rev = df_merged[df_merged["Strand"] == "-"]
         
         ## Adjust p-value calc. based on which strand UNUAR site is on
         results = []
         for df, strand in zip([fwd, rev], ["+", "-"]):
            for rep in rep_list:
               del_cols = self.match_cols(merged_colnames, rep, "Deletions")
               t_cols = self.match_cols(merged_colnames, rep, "T")
               a_cols = self.match_cols(merged_colnames, rep, "A")

               fisher_cols = [t_cols[0],     ## BS
                              del_cols[0],   ## BS
                              t_cols[1],     ## NBS
                              del_cols[1]]   ## NBS
               
               if strand == "-":
                  ## If UNUAR site on - strand, use A counts instead of T
                  fisher_cols[0] = a_cols[0] ## BS
                  fisher_cols[2] = a_cols[1] ## NBS
               
               ## Create copy to disable SettingWithCopyWarning
               df = df.copy()

               ## Calculate p-values
               if set(fisher_cols).issubset(df.columns):
                  df = df.dropna(subset = fisher_cols)
                  arr = df[fisher_cols].values.reshape(-1, 2, 2) 
                  pvals = [fisher_exact(table, alternative = "less")[1] 
                           for table in arr]
                  df[f"{rep}_Pvalue"] = pvals
            results.append(df)
         
         df_pval = self.iteratively_merge(results)
                  
         return df_pval
      except Exception as e:
         print(f"Failed to calculate p-value for {rep}: {e}")
         traceback.print_exc()
         raise

   def calc_pval(self, df_merged: pd.DataFrame, sample: str, pvals_dir: Path):
      ## Calculate p-values for each replicate
      wt_7ko_7lko = sample.split("-")[0] # Determine if sample is WT, 7KO, 7LKO
      df_pval = self.fisher_test(df_merged, wt_7ko_7lko)

      ## Filter by p-value (at least 2/3 replicates pass cutoff)
      pval_cutoff_name = f"{wt_7ko_7lko}_Pvalue_Pass"
      df_pval[pval_cutoff_name] = 0

      pval_list = [col for col in df_pval.columns 
                   if re.search("_Pvalue$", col)]

      ## Count p-vals that pass cutoff
      if re.match(fr"WT.*", str(sample)):
         df_pval[pval_cutoff_name] = (df_pval[pval_list] <= 0.05).sum(axis = 1)               
      else:
         df_pval[pval_cutoff_name] = (df_pval[pval_list] > 0.05).sum(axis = 1)

      ## Select only UNUAR sites that have >=1 p-vals that pass cutoff
      count_cutoff = df_pval[pval_cutoff_name].ge(1)
      df_final = df_pval.loc[count_cutoff]

      ## Save as output
      output_dir = pvals_dir/f"{sample}-Pvals.tsv"
      df_final.to_csv(output_dir, sep = "\t", index = False)

   def calc_avg_std(self, df: pd.DataFrame, 
                    avg_col: str, 
                    std_col: str) -> pd.DataFrame:
      dr_col = [col for col in df.columns 
                if re.search("_DeletionRate_", col)]
      df[avg_col] = df[dr_col].mean(axis = 1)
      df[std_col] = df[dr_col].std(axis = 1)
      return df

   def iteratively_merge(self, list_of_dfs: list):
      df1_colnames = list_of_dfs[0].columns.tolist()
      selected_colnames = df1_colnames[0:18]
      merged = list_of_dfs[0]

      for df in list_of_dfs[1:]:
         if not df.empty:
            merged = pd.merge(merged, df,
                              on = selected_colnames,
                              how = "outer")
   
      return merged

   def merge_reps(self, suffix: str, 
                  tsv_list: list, 
                  subfolder: Path) -> pd.DataFrame:
      """
      1. Search TSVs for matching suffix in filename
      2. Put them in list
      3. Read in as pandas dataframes
      4. Iteratively merge dataframes
      """
      matches = [tsv for tsv in tsv_list if re.search(suffix, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
      merged = self.iteratively_merge(df_list)
      
      """
      1. Define col_start and col_end so that concatenation
         results in examples like:
         a. 7KO_AvgDeletionRate_BS
         b. 7KO_StdDeletionRate_BS
      2. Create AvgDeletionRate and StdDeletionRate columns
         in merged df
      """
      col_start = subfolder.name.split("-")[0]
      col_end = suffix.split("-")[1]
      
      avg_std_colnames = []
      for type in ["Avg", "Std"]:
         name = col_start + f"_{type}DeletionRate_" + col_end
         avg_std_colnames.append(name)

      calc_merged = self.calc_avg_std(merged, 
                                      avg_std_colnames[0], 
                                      avg_std_colnames[1])
      calc_merged = (
         calc_merged.drop_duplicates()
         .sort_values(by = avg_std_colnames[0], ascending = False)
      )
      return calc_merged

def main():
   """
   PURPOSE:
   Filters .tsv files in grouped folders
   """
   current_path = Path.cwd()
   start_dir = current_path/"calculations"

   ## Initialize class
   filtertsv = FilterTSV()

   try: 
      ## Create output folder directories
      processed_folder = current_path/"merged"
      pvals_dir = processed_folder/"pvals"
      final_dir = processed_folder/"final_outputs"
      for dirname in [processed_folder, pvals_dir, final_dir]:
         dirname.mkdir(exist_ok = True, parents = True)

      all_merged = []
      for subfolder in start_dir.iterdir():
         if subfolder.is_dir():
            ## Collect paths of .tsv files and put in list
            tsv_list = sorted(
               subfolder.glob("*.tsv"),
               ## Order by replicate integer
               key = lambda x: int(re.search(r"Rep(\d+)", x.name).group(1))
            )

            ## Merge by replicate, and calc. D.R. avg/std
            merged_bs_nbs = []
            for suffix in ["-BS", "-NBS"]:
               merged = filtertsv.merge_reps(suffix, tsv_list, subfolder)
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
            
            all_merged.append(df_merged)
            
      ## Calculate p-value
      for sample in ["WT-Cyto", "WT-Nuc", "7KO-Cyto", "7KO-Nuc"]:
         filtertsv.calc_pval(df_merged, sample, pvals_dir)

      ## After p-value calculations, create final merged ouputs
      pvals_tsv = list(pvals_dir.glob("*.tsv"))

      for matching_name in ["-Cyto-Pvals", "-Nuc-Pvals"]:
         filtertsv.merge_WT_7KO(matching_name, pvals_tsv, final_dir)

   except Exception as e:
      print(f"Failed to create merged .tsv file: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")