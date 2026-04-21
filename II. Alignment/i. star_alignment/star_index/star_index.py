from pathlib import Path
import subprocess
import traceback
import argparse

def star_index(overhang):
    """
    Before proceeding, the `star_index_hg38` folder must be
    located in the directory containing the other data 
    folders you want to process.
    """
    current_path = Path.cwd()
    genome_dir = f"{current_path}/star_index_hg38"
    genome_fasta = f"{genome_dir}/GCF_000001405.40_GRCh38.p14_genomic.fa"
    gtf_file = f"{genome_dir}/GCF_000001405.40_GRCh38.p14_genomic.gtf"
    threads = 2

    if len(list(Path(genome_dir).glob("*"))) == 4: ## ensures this only runs if index hasn't been created yet
        try:
            cmd = ["STAR", "--runThreadN", str(threads),
                   "--runMode", "genomeGenerate",
                   "--genomeDir", str(genome_dir),
                   "--genomeFastaFiles", str(genome_fasta),
                   "--sjdbGTFfile", str(gtf_file),
                   "--sjdbOverhang", str(overhang)]
            subprocess.run(cmd, 
                           check = True,
                           capture_output = True, 
                           text = True)
            
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to build STAR index: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
    else:
        print("Index already exists.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overhang", type = int, default = 100, help = "Overhang value (default: 100)")
    args = parser.parse_args()
    print("Creating STAR index...")
    star_index(args.overhang)
    print("Index created.")