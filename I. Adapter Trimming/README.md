# PART I: Adapter Trimming
### Starter files
<img src="https://github.com/user-attachments/assets/b281ca65-3535-4770-b361-f69a611ed3e5" width="400"/>

### Overview
* Performs quality control, or adapter trimming and fastq pre-processing (fastp), on raw fastqs.
  
  i. **Adapter trimming**
    * STAR two-pass adapter removal from 3' ends.
      * 1st pass = Targets barcode primer sequences
      * 2nd pass = Targets remaining tn5 chimeric sequences
    * Reads <30 nt after trimming are discarded.
    * Fastq outputs are in the R1 and R2 format (**suffix:** `cutadapt_R#.fastq.gz`), and are retained for debugging purposes.
  
  ii. **Fastp**
    * Extracts UMI from 5' end of both reads.
    * Performs quality trimming from 5' → 3' ends, and corrects mismatches.
    * Orphaned/unpaired reads and reads with >30% remaining bases below Q15 are discarded.
    * Fastq outputs remain in the R1 and R2 format (**suffix:** `cutadapt_fastp_R#.fastq.gz`) because no merging occurs.
      * These files contain sequences and base quality scores.

### Instructions
1. Create `raw_fastqs` folder if it doesn't already exist, and upload the remaining starter files to your GLC directory.
2. To bypass manually writing SLURM tasks for each sample in the 📁 `raw_fastqs` directory, create an SBATCH by running `write_slurm_cutadapt.py` in Bash with the following input commands:
   * **--input_folder:** Name of folder containing input raw fastqs.
   * **--output_folder:** Name of folder for trimmed fastqs.
   * **--email:** Email that will be notified when SLURM task begins/ends.
   * **--slurm_acct:** SLURM account.
   * **--walltime:** Amount of time allocated for job.
   * **--mem:** Amount of memory allocated for job.
   
   See example:
   ```
   python3 write_slurm_cutadapt.py --input_folder raw_fastqs --output_folder trimmed_reads --email uniqname@umich.edu --slurm_acct cweidman99 --walltime 1:00:00 --mem 10000
   ```
   **Output:** 📄 `SBATCHSubArr-CUT_FASTP.sbatch`
3. In Bash, run the following commands:
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
---
### Citations
* `run_cutadapt_fastp.py` by Chase Weidmann. If you have any questions, please reach out to [chaseaw](https://github.com/chaseaw).
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5
