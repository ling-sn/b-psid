# PART II: Alignment
> ### Table of contents
> 1. [STAR alignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/i.%20star_alignment#star-alignment) :star: **YOU ARE HERE!**
> 2. Realignment
> 3. BAM splitting
## STAR alignment
### Starter files
<img src="https://github.com/user-attachments/assets/a952f05f-8174-4efc-bf71-19d750776d53" width="400"/>

### Overview
* Aligns trimmed reads with [STAR](https://github.com/alexdobin/STAR), using the hg38 genome as a reference.
* During alignment, we map contaminant RNAs using a FASTA and proceed downstream with the unmapped reads (mRNA) only.
  
  <img src="https://github.com/user-attachments/assets/8de821ce-dce6-4050-958b-eca6f8a94a0a" width="200"/>



### Instructions
1. Excluding the 📁 `trimmed_reads` folder (which should already exist), upload the remaining starter files to your GLC directory.
2. Create an SBATCH by running `write_slurm_alignment.py` in Bash with the following input commands:
   * **--input_folder:**
   * **--output_folder:**
   * **--aligner_type:**
   * **--genome_idx:**
   * **--filter_idx:**
   * **--library:**
   * **--two_pass:**
   * **--emit_dedup:**
   * **--email:**
   * **--slurm_acct:**
   * **--walltime:**
   * **--mem:**

   See example:
   ```
   python3 write_slurm_alignment.py --input_folder trimmed_reads --output_folder star_aligned --aligner_type star --genome_idx ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --filter_idx ~/umms-RNAlabDATA/Software/genome_indices/contaminants/bowtie2_contam_index/SL_contaminants/contaminants_index --library RF --two_pass --emit_dedup SBATCHSubArr-DEDUP.sbatch --email uniqname@umich.edu --slurm_acct cweidman99 --walltime 10:00:00 --mem 50000
   ```
   **Output:** 📄 `SBATCHSubArr-Align-STAR.sbatch`
3. In Bash, run the following commands to obtain STAR aligned files:
   ```
   
   ```
   **Output:** 📁 `star_aligned`
4. An SBATCH file called 

### When do I use this script?
* After obtaining trimmed reads, run this in the same directory.

### Understanding the SBATCH
```
python3 -u run_align.py --input trimmed_reads --output star_aligned --aligner star --index ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 -C 2 -L RF -S KEH-Rep1-7LKO-HEK293T-Cyto-NBS --filter_index ~/umms-RNAlabDATA/Software/genome_indices/contaminants/bowtie2_contam_index/SL_contaminants/contaminants_index -T --emit_dedup_slurm SBATCHSubArr-DEDUP.sbatch
```
* Text

### Creating STAR and contaminants index (optional)
