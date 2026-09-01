"""
@author: ling-sn
Last updated: August 2026
Use B-PSID conda environment
"""
from pathlib import Path
import traceback
import argparse
import subprocess

class SplitBAM:
   def merge_bams(self, bam1: Path, bam2: Path, output: Path):
      try:
         sorted_output = output.with_name(f"{output.stem}_sorted.bam")
         
         subprocess.run(
            [
               "samtools", "merge",
               "-f", str(output),
               str(bam1), str(bam2)
            ],
            check = True
         )
         subprocess.run(["samtools", "sort", str(output), "-o", 
                         str(sorted_output)], check = True)
         subprocess.run(["samtools", "index", str(sorted_output)], 
                        check = True)

      except Exception as e:
         print(f"Failed to merge {bam1.name} and {bam2.name}: {e}")
         traceback.print_exc()
         raise

   def map_reverse(self, output: Path, small_f: int, bam: str):
      try:
         if not output.exists():
            with open(output, "wb") as f:
               subprocess.run(
                  [
                     "samtools", "view", "-b",
                     "-f", str(small_f),
                     bam
                  ],
                  stdout = f,
                  check = True
               )

      except Exception as e:
         print(f"Failed to collect alignments for {output.name}: {e}")
         traceback.print_exc()
         raise

   def map_forward(self, output: Path, small_f: int, big_F: int, bam: str):
      try:
         if not output.exists():
            with open(output, "wb") as f:
               subprocess.run(
                  [
                     "samtools", "view", "-b",
                     "-f", str(small_f),
                     "-F", str(big_F),
                     bam
                  ],
                  stdout = f,
                  check = True
               )

      except Exception as e:
         print(f"Failed to collect alignments for {output.name}: {e}")
         traceback.print_exc()
         raise

   def collect_rev(self, rev1: Path, rev2: Path, rev: Path, bam: str):
      """
      PURPOSE:
      * Collect alignments for rev.bam
      """
      ## 1) Second-in-pair maps to the reverse strand
      self.map_reverse(rev1, 144, bam)

      ## 2) First-in-pair maps to the forward strand
      self.map_forward(rev2, 64, 16, bam)

      ## 3) Combine alignments that originate on the reverse strand
      self.merge_bams(rev1, rev2, rev)

   def collect_fwd(self, fwd1: Path, fwd2: Path, fwd: Path, bam: str):
      """ 
      PURPOSE:
      * Collect alignments for fwd.bam
      """
      ## 1) Second-in-pair maps to the forward strand
      self.map_forward(fwd1, 128, 16, bam)

      ## 2) First-in-pair maps to the reverse strand
      self.map_reverse(fwd2, 80, bam)

      ## 3) Combine alignments that originate on the forward strand
      self.merge_bams(fwd1, fwd2, fwd)

def main(bam_folder: str, library_type: str):
   """
   PURPOSE:
   * Obtain specific alignments from BAM files, and then combine 
     them into separate forward and reverse strand files.
   ---
   LIBRARY TYPE:
   * Default: Assumes paired end reads with RF library type (NEBNext).
   * If FR specified, then whatever is in fwd -> (goes to) rev, 
     and whatever is in rev -> (goes to) fwd.
   """
   try:
      ## Obtain BAM file in folder
      folder = Path("realignments")/bam_folder
      bam = str(next(file for file in folder.glob("*.bam")))

      ## Initialize class
      splitbam = SplitBAM()

      ## Define directories
      """
      Explanation of directories list
      [0] = fwd1.bam
      [1] = fwd2.bam
      [2] = fwd.bam
      [3] = rev1.bam
      [4] = rev2.bam
      [5] = rev.bam
      """
      directories = []
      for type in ["fwd", "rev"]:
         r1 = folder/f"{type}1.bam"
         r2 = folder/f"{type}2.bam"
         output = folder/f"{type}.bam"
         all_files = [r1, r2, output]
         directories.extend(all_files)

      if library_type == "RF":
         splitbam.collect_fwd(directories[0], directories[1], directories[2], bam)
         splitbam.collect_rev(directories[3], directories[4], directories[5], bam)
      else:
         ## Swap around output file names
         splitbam.collect_fwd(directories[3], directories[4], directories[5], bam)
         splitbam.collect_rev(directories[0], directories[1], directories[2], bam)

      """
      PART III:
      * Clean up intermediate files
      """
      subprocess.run(
         [
            "rm",
            str(directories[0]), str(directories[1]),
            str(directories[2]), str(directories[3]),
            str(directories[4]), str(directories[5])
         ]
      )

   except Exception as e:
      print(f"Failed to split BAM in {folder.stem}: {e}")
      traceback.print_exc()
      raise

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description = ("Splits alignment BAM into separate files"
                                                   " (forward and reverse) by detecting strandedness"))
   parser.add_argument("--bam_folder", help = "Path to folder with BAM file", required = True)
   parser.add_argument("--library_type", choices = ["RF", "FR"], default = "RF")
   args = parser.parse_args()

   print("Splitting BAM files...")
   main(args.bam_folder, args.library_type)
   print("Process finished.")