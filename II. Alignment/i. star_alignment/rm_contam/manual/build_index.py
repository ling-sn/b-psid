import subprocess
from pathlib import Path
import traceback

def build_bowtie2_index():
    """
    Builds bowtie2 index once.
        1. Requires contaminants.fa in current working directory
        3. Call build_index.sbatch after activating B-PSID conda environment
        3. Use before rm_contam.sbatch
    """
    current_path = Path.cwd()
    contaminants_dir = current_path/"contaminants.fa"
    bowtie2_index = current_path/"contaminants_index"
    bt2_files = list(current_path.glob("*.bt2")) ## produces list of files
    
    if not bt2_files: ## checks if list is empty; if so, proceed
        try:
            cmd = ["bowtie2-build",
                    str(contaminants_dir),
                    str(bowtie2_index)]
            result = subprocess.run(cmd, 
                                    check = True, ## if command returns non-zero exit status, raise error
                                    capture_output = True, 
                                    text = True)
            return result
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to build bowtie2 index: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise

if __name__ == "__main__":
    print("Creating bowtie2 index...")
    build_bowtie2_index()
    print("Index created.")