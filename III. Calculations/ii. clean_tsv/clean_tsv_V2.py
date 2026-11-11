"""
@author: ling-sn
Last updated: November 2026
Use B-PSID conda environment
MAIN CHANGE: Use rate ratio distribution instead of Fisher's Exact Test
"""
from pathlib import Path
import traceback
import argparse
import pandas as pd
import re
import concurrent.futures
import Bio
from Bio import Entrez
pd.options.mode.chained_assignment = None

class MergeData:
   def obtain_mane(self, email: str, genes_str: str):
      """
      PURPOSE:
      Obtain all MANE transcripts given a list of genes.
      * MANE = Transcript that is most representative of a certain gene
      """
      ## Provide identity to access Entrez
      Entrez.email = email
      
      ## Iterate through each gene to obtain all associated MANE transcripts
      
   
   def obtain_all_genes(self, email: str, transcripts_str: str):
      """
      PURPOSE:
      Use list of transcripts associated with the UNUAR site to find
      all associated genes.
      """
      ## Provide identity to access Entrez
      Entrez.email = email
      
      ## Iterate through each transcript to obtain all associated genes
      all_transcripts = transcripts_str.split(", ")
      all_genes = []
      for transcript in all_transcripts:
         """
         I. Obtain nucleotide ID
         """
         with Entrez.esearch(
            db = "nucleotide", term = transcript
         ) as handle:
            nucleotide_record = Entrez.read(handle)
            nucleotide_id = nucleotide_record["IdList"][0]

         """
         II. Obtain gene ID
            * gene_record[0]: gene_record is a list with one item
            * LinkSetDb[0]: LinkSetDb is a list with one item
            * ["Link"][0]["Id"]: Within the Link dict, select the value for Id
         """
         with Entrez.elink(
            dbfrom = "nucleotide", db = "gene", id = nucleotide_id
         ) as handle:
            gene_record = Entrez.read(handle)
            gene_id = gene_record[0]["LinkSetDb"][0]["Link"][0]["Id"]
            
         """
         III. Query the database
            * gene_summary["DocumentSummarySet"]: Select DocumentSummarySet key from dict
            * ["DocumentSummary"][0]: DocumentSummary key returns a list with one item
            * ["Name"]: Obtain gene name
         """
         with Entrez.esummary(
            db = "gene", id = gene_id
         ) as handle:
            gene_summary = Entrez.read(handle)
            gene_name = gene_summary["DocumentSummarySet"]["DocumentSummary"][0]["Name"]

         if not gene_name in all_genes:
            all_genes.append(gene_name)
         
      combined_genes = ", ".join(all_genes)
      return combined_genes
   
   def clean_tsv(
      self, df: pd.DataFrame, 
      location: str, 
      email: str,
      final_dir: Path
   ):
      """
      PURPOSE:
      Explicitly list all genes associated with transcripts at UNUAR site.
      """
      ## Obtain all associated genes
      MAX_THREADS = 20
      threads = min(MAX_THREADS, len(df["AllAssocTranscripts"]))
      with concurrent.futures.ThreadPoolExecutor(max_workers = threads) as executor:
         df["AllAssocGenes"] = list(
            executor.map(
               self.obtain_all_genes, email, df["AllAssocTranscripts"]
            )
         )
         df["MANETranscripts"] = list(
            executor.map(
               self.obtain_mane, email, df["AllAssocGenes"]
            )
         )
      df.insert(0, "AllAssocGenes", df.pop("AllAssocGenes"))
      df.insert(2, "MANETranscripts", df.pop("MANETranscripts"))
      
      ## Save merged dataframe as TSV
      merged_dir = final_dir/f"{location}.tsv"
      df.to_csv(merged_dir, sep = "\t", index = False)
      
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
      
def main(email: str):
   """
   PURPOSE:
   Filters .tsv files in grouped folders
   """
   current_path = Path.cwd()
   start_dir = current_path/"calculations"

   ## Initialize classes
   mergedata = MergeData()
   filtertsv = FilterTSV()

   try: 
      ## Create output folder directories
      final_dir = current_path/"final_outputs"
      final_dir.mkdir(exist_ok = True, parents = True)

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
               merged = mergedata.merge_reps(suffix, tsv_list, subfolder)
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
            
      ## Collect WT dataframes with filtered rate ratios
      for location in ["Cyto", "Nuc"]:
         ## Implicitly ignore 7LKO and only handle 7KO and WT data
         # (i) Open up 7KO dataframe 
         # (ii) Calculate rate ratios
         knockout_df = filtertsv.calc_rr(
            all_merged["7KO-" + location], "7KO"
         )
         top5_perc = knockout_df["RateRatio"].quantile(0.95)
         
         # (i) Open up WT dataframe
         # (ii) Calculate rate ratios
         # (iii) Only keep rate ratios >= 95th percentile from 7KO distribution
         wt_df = filtertsv.calc_rr(
            all_merged["WT-" + location], "WT"
         )
         filtered_wt_df = wt_df[wt_df["RateRatio"] >= top5_perc]
         
         ## After rate ratio filtering, create final outputs
         mergedata.clean_tsv(filtered_wt_df, location, email, final_dir)

   except Exception as e:
      print(f"Failed to create merged .tsv file: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   parser = argparse.ArgumentParser(
      description = "Merges and cleans dataframes from deletion rate calculation."
   )
   parser.add_argument(
      "--email", 
      default = "<uniqname>@umich.edu"
   )
   args = parser.parse_args()
   
   print("Cleaning .tsv files...")
   main(args.email)
   print("Process finished.")