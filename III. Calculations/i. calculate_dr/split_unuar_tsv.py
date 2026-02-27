## Use RNA-STAR conda environment
from pathlib import Path
import traceback
import argparse
import pandas as pd

## Disable .loc indexing warning
pd.options.mode.chained_assignment = None

def main():
   try:
      """
      PURPOSE: 
      1. Merges dataframes (UNUAR sites & fitted values from BID-Seq)
      2. Keeps only 3UTR sites and splits merged dataframe by strand [+/-]
      """
      current_path = Path.cwd()
      left = pd.read_csv(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools"
                              "/UNUAR_motif_sites_mRNA_hg38p14.tsv").expanduser(), sep = "\t")
      right = pd.read_excel(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools"
                                 "/Zhang_HE_NatureProtocols_2023_SupplementaryTable1.xlsx").expanduser())

      df = pd.merge(left, right, how = "left", on = "Motif")
      df = df[df["Region"] == "3UTR"]

      ## Split into two dataframes based on strand (+/-)
      fwd_condition = df["Strand"] == "+"
      fwd = df[fwd_condition]
      rev = df[~fwd_condition]

      ## Output as TSV files
      base = "B-PsiD_UNUAR_motif_sites_mRNA_hg38"
      types = ["fwd", "rev"]
      unuar_sites = [fwd, rev]

      for type, df in zip(types, unuar_sites):
         output_name = current_path/f"{base}_{type}.tsv"
         df.to_csv(output_name, sep = "\t", index = False)
   
   except Exception as e:
      print(f"Failed to merge dataframes and output separate .TSVs based on strand: {e}")
      traceback.print_exc()
      raise

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description = ("Merges dataframes for B-PsiD and splits "
                                                   "based on strand [+/-] for UNUAR site"))
   print("Splitting UNUAR .tsv file...")
   main()
   print("Process finished.")