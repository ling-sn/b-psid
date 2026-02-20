## use RNA-SEQ conda environment
from pathlib import Path
import traceback
import argparse
import subprocess

def merge_bams(bam1: Path, bam2: Path, output: Path):
   try:
      subprocess.run(
         [
            "samtools", "merge",
            "-f", str(output),
            str(bam1), str(bam2)
         ],
         check = True
      )
      subprocess.run(["samtools", "index", str(output)], check = True)

   except Exception as e:
      print(f"Failed to merge {bam1.name} and {bam2.name}: {e}")
      traceback.print_exc()
      raise

def map_reverse(output: Path, small_f: int, bam: str):
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

def map_forward(output: Path, small_f: int, big_F: int, bam: str):
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

def main(bam_folder: str):
   """
   PURPOSE:
   * Obtain specific alignments from BAM files, and then combine 
     them into separate forward and reverse strand files.
   ---
   REQUIRES:
   * Paired end reads with RF library type (NEBNext).
   """
   try:
      ## Obtain BAM file in folder
      folder = Path(bam_folder)
      bam = str(next(file for file in folder.glob("*.bam")))

      """ 
      PART I:
      * Collect alignments for fwd.bam
      """
      ## 1) Second-in-pair maps to the forward strand
      fwd1 = folder/"fwd1.bam"
      map_forward(fwd1, 128, 16, bam)
      
      ## 2) First-in-pair maps to the reverse strand
      fwd2 = folder/"fwd2.bam"
      map_reverse(fwd2, 80, bam)

      ## 3) Combine alignments that originate on the forward strand
      fwd = folder/"fwd.bam"
      merge_bams(fwd1, fwd2, fwd)

      """
      PART II:
      * Collect alignments for rev.bam
      """
      ## 1) Second-in-pair maps to the reverse strand
      rev1 = folder/"rev1.bam"
      map_reverse(rev1, 144, bam)

      ## 2) First-in-pair maps to the forward strand
      rev2 = folder/"rev2.bam"
      map_forward(rev2, 64, 16, bam)

      ## 3) Combine alignments that originate on the reverse strand
      rev = folder/"rev.bam"
      merge_bams(rev1, rev2, rev)

      """
      PART III:
      * Clean up intermediate files
      """
      subprocess.run(
         [
            "rm",
            str(fwd1), str(fwd2),
            str(rev1), str(rev2)
         ]
      )

   except Exception as e:
      print(f"Failed to split BAM in {folder.stem}: {e}")
      traceback.print_exc()
      raise

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description = ("Splits alignment BAM into separate files (forward and reverse)"
                                                   " by detecting strandedness"))
   parser.add_argument("--bam_folder", help = "Path to folder with BAM file", required = True)
   # parser.add_argument("--library", help = "Library type. Choices: RF or FR")
   args = parser.parse_args()

   print("Splitting BAM files...")
   main(args.bam_folder)
   print("Process finished.")