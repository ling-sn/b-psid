# PART I: Adapter Trimming
### Necessary files
<img src="https://github.com/user-attachments/assets/b281ca65-3535-4770-b361-f69a611ed3e5" width="400"/>

### Overview
* Performs quality control (adapter trimming and merging) on raw fastqs.
* The raw fastq files contain information for each **Read 1** and **Read 2**, such as their sequences and base quality scores. If one of the reads has poor quality scores, then it will be discarded and the read will be unpaired. Otherwise, depending on if the reads overlap, they will be merged or unmerged.

  <img src="https://github.com/user-attachments/assets/96b31246-c87b-482b-8649-53b5c5ded5bb" width="400"/>


### Instructions
1. To bypass manually writing SLURM tasks for each sample in the 📁 `raw_fastqs` directory, run `write_slurm_cutadapt.py` directly in Bash with the following input commands:
   * **--input_folder:** Folder for input raw fastqs.
   * **--output_folder:** Folder for trimmed fastqs.
   * **--email:** Email that will be notified when SLURM task begins/ends.
   * **--slurm_acct:** SLURM account.
   * **--walltime:** Amount of time allocated for job.
   * **--mem:** Amount of memory allocated for job.
   
   See example:
   ```
   python3 write_slurm.py --input_folder raw_fastqs --output_folder trimmed_reads --email uniqname@umich.edu --slurm_acct cweidman99 --walltime 1:00:00 --mem 10000
   ```
   **Output:** 📄 `SBATCHSubArr-CUT_FASTP.sbatch`.
2. In Bash, run the following commands:
   ```
   conda activate B-PSID
   sbatch SBATCHSubArr-CUT_FASTP.sbatch
   ```
   **Output:** 📁 `trimmed_reads`

### When do I use this script?
* Run after creating B-PSID conda environment and moving all sequencing files to a 📁 `raw_fastqs` folder.

### Understanding the SBATCH
```
python3 run_cutadapt_fastp.py --input raw_fastqs --output trimmed_reads -C 2 -U 12 -S KEH-Rep1-WT-HEK293T-Nuc-BS
```
* **--input:** Name of folder containing raw fastq.gz sequencing files.
* **--output:** Name of output folder for trimmed reads.
* **-C:** Number of CPUs. (Default = 2)
* **-U:** Length of UMI at 5' end of reads.
* **-S:** Only process files with sample prefix within the input directory.

### Citations
* `run_cutadapt_fastp.py` by Chase Weidmann. If you have any questions, please reach out to [chaseaw](https://github.com/chaseaw).
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5
