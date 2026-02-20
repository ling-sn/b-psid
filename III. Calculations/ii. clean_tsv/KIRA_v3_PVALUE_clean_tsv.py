## use RNA-STAR conda environment
from pathlib import Path
import traceback
import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact

class FilterTSV:
   def merge_reps(self, suffix, tsv_list, subfolder, reps_dir):
      """
      1. Search TSVs for matching suffix in filename
      2. Put them in list
      3. Read in as pandas dataframes
      """
      matches = [tsv for tsv in tsv_list if re.search(suffix, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]

      """
      Copy + paste iterative merging code from original clean_tsv
      because there are 3 replicates
      """
      df1_colnames = df_list[0].columns.tolist()
      selected_colnames = df1_colnames[0:17]
      merged = df_list[0]

      for df in df_list[1:]:
         if not df.empty:
            merged = pd.merge(merged, df,
                              on = selected_colnames,
                              how = "outer")
      
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
      avg_col = col_start + "_AvgDeletionRate_" + col_end
      std_col = col_start + "_StdDeletionRate_" + col_end

      calc_merged = self.calc_avg_std(merged, avg_col, std_col)
      calc_merged = calc_merged.drop_duplicates().sort_values(by = avg_col, ascending = False)

      """
      Save merged dataframe as TSV
      """
      merged_dir = reps_dir/f"{subfolder.name}{suffix}.tsv"
      calc_merged.to_csv(merged_dir, sep = "\t", index = False)

   def calc_avg_std(self, df, avg_col, std_col):
      dr_col = [col for col in df.columns if re.search("_DeletionRate_", col)]
      df[avg_col] = df[dr_col].mean(axis = 1)
      df[std_col] = df[dr_col].std(axis = 1)
      return df

   def merge_BS_NBS(self, sample, merged_reps_tsv, pvals_dir):
      matches = [tsv for tsv in merged_reps_tsv if re.search(sample, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
      
      """
      1. Ensure NBS is merged with BS, so BS columns appear first
      2. If either dataframe is not empty, then merge w/ inner join
      3. No need to iteratively merge because there are only 2 files
      """
      first_cols = df_list[0].columns.tolist()

      if any(re.search("-BS", col) for col in first_cols):
         df1 = df_list[0]
         df2 = df_list[1]
      else:
         df1 = df_list[1]
         df2 = df_list[0]

      """
      Merge df1 & df2
      """
      selected_colnames = (df1.columns.tolist())[0:17]
      if not df1.empty and not df2.empty:
         df_merged = pd.merge(df1, df2, on = selected_colnames, how = "outer")
      elif df1.empty:
         df_merged = df2
      else:
         df_merged = df1

      ## Calculate p-values in each TSV
      merged_colnames = df_merged.columns.tolist()
      rep_list = sorted(
                  set([re.search(r"(Rep\d+)", col).group(1) for col in merged_colnames 
                  if re.search(r"(Rep\d+)", col)]), 
                  key = lambda x: int(re.search(r"Rep(\d+)", x).group(1))
                 )
      df_pval = self.calc_pval(df_merged, merged_colnames, rep_list)

      ## Filter by p-value (at least 2/3 replicates pass cutoff)
      wt_7ko = sample.split("-")[0]
      pval_cutoff_name = f"{wt_7ko}_Pvalue_Pass"
      df_pval[pval_cutoff_name] = 0

      pval_list = [col for col in df_pval.columns 
                   if re.search("_Pvalue$", col)]

      for col in pval_list:
         if re.match(fr"WT.*", str(sample)):
            pval_condition = df_pval[col] <= 0.05                  
         elif re.match(fr"7KO.*", str(sample)):
            pval_condition = df_pval[col] > 0.05

         df_pval.loc[pval_condition, pval_cutoff_name] += 1

      count_cutoff = df_pval[pval_cutoff_name].ge(2)
      df_final = df_pval.loc[count_cutoff]

      ## Save as output
      output_dir = pvals_dir/f"{sample}-Pvals.tsv"
      df_final.to_csv(output_dir, sep = "\t", index = False)

   def calc_pval(self, df_merged, merged_colnames, rep_list):
      try:
         for rep in rep_list:
            bs_del_col = [col for col in merged_colnames 
                          if re.search(f"{rep}_Deletions_BS", col)]
            nbs_del_col = [col for col in merged_colnames 
                           if re.search(f"{rep}_Deletions_NBS", col)]

            bs_t_col = [col for col in merged_colnames
                        if re.search(f"{rep}_T_BS", col)]
            nbs_t_col = [col for col in merged_colnames
                         if re.search(f"{rep}_T_NBS", col)]

            t_cols = [bs_t_col[0], nbs_t_col[0]]
            del_cols = [bs_del_col[0], nbs_del_col[0]]

            fisher_cols = [t_cols[0], 
                           del_cols[0], 
                           t_cols[1], 
                           del_cols[1]]
            
            ## Create copy to disable SettingWithCopyWarning
            df_merged = df_merged.copy()

            ## Calculate p-values
            if set(fisher_cols).issubset(df_merged.columns):
               df_merged = df_merged.dropna(subset = fisher_cols)
               arr = df_merged[fisher_cols].values.reshape(-1, 2, 2) 
               pvals = [fisher_exact(table, alternative = "less")[1] for table in arr]
               df_merged[f"{rep}_Pvalue"] = pvals
                  
         return df_merged
      except Exception as e:
         print(f"Failed to calculate p-value for {rep}: {e}")
         traceback.print_exc()
         raise
   
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
      Rename differing columns (except {WT|7KO}_Pvalue_Pass column) with prefix
      """
      selected_colnames = (df1.columns.tolist())[0:17]
      diff_cols = (df1.columns.difference(selected_colnames, sort = False))[:-1]
      
      for df, prefix in zip([df1, df2], ["WT_", "7KO_"]):
         for old_name in diff_cols:
            new_name = prefix + old_name
            df.rename(columns = {old_name: new_name}, inplace = True)
            
      """
      Merge df1 & df2
      """
      if not df1.empty and not df2.empty:
         df_merged = pd.merge(df1, df2, on = selected_colnames, how = "inner")
      elif df1.empty:
         df_merged = df2
      else:
         df_merged = df1
         
      """
      1. Create output name
         e.g., 7KO-Cyto-Pvals + WT-Cyto-Pvals -> Cyto
      2. Save merged dataframe as TSV
      """
      output_name = (matches[0].stem).split("-")[1]
      merged_dir = final_dir/f"{output_name}.tsv"
      df_merged.to_csv(merged_dir, sep = "\t", index = False)

def main():
   """
   PURPOSE:
   Filters .tsv files in grouped folders
   """
   current_path = Path.cwd()
   input_dir = current_path/"calculations"

   ## Initialize class
   filtertsv = FilterTSV()

   try: 
      ## Main folder that stores subfolders of merged outputs
      processed_folder = current_path/"merged"
      processed_folder.mkdir(exist_ok = True, parents = True)

      reps_dir = processed_folder/"merged_reps"
      reps_dir.mkdir(exist_ok = True, parents = True)

      pvals_dir = processed_folder/"pvals"
      pvals_dir.mkdir(exist_ok = True, parents = True)

      final_dir = processed_folder/"final_outputs"
      final_dir.mkdir(exist_ok = True, parents = True)

      for subfolder in input_dir.iterdir():
         tsv_folder = input_dir/subfolder/"individual_tsv"

         if subfolder.is_dir():
            ## Collect paths of .tsv files and put in list
            tsv_list = sorted(
               tsv_folder.glob("*.tsv"),
               ## Order by replicate integer
               key = lambda x: int(re.search(r"Rep(\d+)", x.name).group(1))
            ) 

            ## Merge replicates for each sample type, and calc. avg/std D.R.
            for suffix in ["-BS", "-NBS"]:
               filtertsv.merge_reps(suffix, tsv_list, subfolder, reps_dir)
      
      ## Collect all TSVs in reps_dir
      merged_reps_tsv = list(reps_dir.glob("*.tsv"))

      ## Merge TSV pairs by BS/NBS
      """
      Example: 
      * WT-Cyto-BS + WT-Cyto-NBS <-> WT-Cyto
      """
      for sample in ["WT-Cyto", "WT-Nuc", "7KO-Cyto", "7KO-Nuc"]:
         ## Merge matching files, then calculate and filter by p-value
         filtertsv.merge_BS_NBS(sample, merged_reps_tsv, pvals_dir)

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