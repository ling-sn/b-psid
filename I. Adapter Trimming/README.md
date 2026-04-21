# PART I: Adapter Trimming
### Necessary files
<img src="https://github.com/user-attachments/assets/b281ca65-3535-4770-b361-f69a611ed3e5" width="400"/>

### Overview
* Performs quality control (adapter trimming and merging) on raw fastqs.
* The raw fastq files contain information for each Read 1 and Read 2, such as their sequences and base quality scores. If one of the reads has poor quality scores, then it will be discarded and the read will be unpaired. Otherwise, depending on if the reads overlap, the reads will be merged or unmerged.

<img src="https://github.com/user-attachments/assets/fbc286ad-1a8c-473e-82e7-5c672a42119e" width="400"/>

### Instructions
1. Run `write_slurm_cutadapt.py` directly in Bash with the following input commands:
   * --input_folder
   * s
   
   See example:
   ```
   python3 write_slurm.py --input_folder raw_fastqs --output_folder trimmed_reads --email uniqname@umich.edu --slurm_acct cweidman99 --walltime 1:00:00 --mem 10000
   ```

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
