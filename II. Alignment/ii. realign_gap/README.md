# PART II: Alignment
> ### Table of contents
> 1. [STAR alignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/i.%20star_alignment#table-of-contents)
> 2. [Realignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/ii.%20realign_gap#table-of-contents) :star: **YOU ARE HERE!**
> 3. [BAM splitting](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/iii.%20bam_splitting#table-of-contents)
## Realignment
### Starter files
<img src="https://github.com/user-attachments/assets/a35085c4-a14d-40cf-a4e9-832868e20e7a" width="400"/>

### Overview
* STAR alignment optimizes both speed and efficiency, which can lead to some misalignments. Therefore, we use the Smith-Waterman algorithm for slower, but more accurate mapping around areas we expect to see deletions.
  
  > "Alignment software usually uses a seed alignment algorithm to increase alignment speed; however, this also affects pairwise alignment accuracy, especially for bases near deletion signatures. To solve this, we integrated the Smith-Waterman local alignment algorithm into the pipeline for realignment. Reads that contained any mismatch, deletion, insertion, soft-clip or splicing were further processed by the realignment tool in the BID-pipe package. By setting the penalty of gap open and gap extension as −3 and −2, respectively, deletion signatures can have a higher priority in the alignment" (_Zhang et al., 520_).

  <img src="https://github.com/user-attachments/assets/a64d90a1-1cf7-4e93-8205-058c1c15279b" width="400"/>
### Instructions
1. Excluding the 📁 `star_aligned` folder (which should already exist), upload the remaining starter files to your GLC directory.
2. Create an SBATCH by running `write_slurm_realignment.py` in Bash with the following input commands:
   * **--star_folder:** Name of folder containing deduplicated BAMs.
   * **--fasta:** Path to reference FASTA.
   * **--email:** Email that will be notified when SLURM task begins/ends.
   * **--slurm_acct:** SLURM account.
   * **--walltime:** Amount of time allocated for job.
   * **--mem:** Amount of memory allocated for job.

   See example:
   ```
   python3 write_slurm_realignment.py --star_folder star_aligned --fasta ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --email uniqname@umich.edu --slurm_acct cweidman99 --walltime 30:00 --mem 4000
   ```
   **Output:** 📄 `realignGap.sbatch`
3. In Bash, run the following commands to obtain the realigned files:
   ```
   conda activate B-PSID
   sbatch realignGap.sbatch
   ```
   **Output:** 📁 `realignments`

### When do I use this script?
* After obtaining deduplicated BAMs, run this in the same directory.

### Understanding the SBATCH
```
python3 realignGap.py --star_folder star_aligned --subf KEH-Rep1-7LKO-HEK293T-Cyto-NBS --fasta_dir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa
```
* **--star_folder:** Name of folder containing deduplicated BAMs.
* **--subf:** Name of sample subfolder.
* **--fasta_dir:** Path to reference FASTA that was used to build the STAR genome index. This will always be `~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa`.

---
### Citations
* Realignment code adapted from [realignGap](https://github.com/y9c/pseudoU-BIDseq/blob/main/bin/realignGap) by Ye Chang.
