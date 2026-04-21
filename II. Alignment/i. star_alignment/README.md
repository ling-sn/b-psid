# PART II: Alignment
> ### Table of contents
> 1. [STAR alignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/i.%20star_alignment#table-of-contents) :star: **YOU ARE HERE!**
> 2. Realignment
> 3. BAM splitting
## STAR alignment
### Starter files
<img src="https://github.com/user-attachments/assets/a952f05f-8174-4efc-bf71-19d750776d53" width="400"/>

### Overview
* Aligns trimmed reads with [STAR](https://github.com/alexdobin/STAR) using hg38 as the reference genome.
* During alignment, mapped contaminant RNAs are obtained with a FASTA and only the unmapped reads (mRNA) are used for downstream analysis.
  
  <img src="https://github.com/user-attachments/assets/8de821ce-dce6-4050-958b-eca6f8a94a0a" width="200"/>

### Instructions
1. Excluding the 📁 `trimmed_reads` folder (which should already exist), upload the remaining starter files to your GLC directory.
2. Create an SBATCH by running `write_slurm_alignment.py` in Bash with the following input commands:
   * **--input_folder:** Name of folder containing trimmed reads.
   * **--output_folder:** Name of folder for STAR aligned outputs.
   * **--aligner_type:** Specify aligner type. (Choices = hisat2, star)
   * **--genome_idx:** Directory to STAR genome index.
   * **--filter_idx:** Directory to contaminants index, including prefix of index (e.g., everything up to .1.bt2).
   * **--library:** Specify which strand is on Read 1. Default is unstranded, but NEBNext uses RF. (Choices = RF, FR, unstranded)
   * **--two_pass:** Enable two-pass STAR alignment.
   * **--emit_dedup:** Name of output SBATCH file for deduplication.
   * **--email:** Email that will be notified when SLURM task begins/ends.
   * **--slurm_acct:** SLURM account.
   * **--walltime:** Amount of time allocated for job.
   * **--mem:** Amount of memory allocated for job.

   See example:
   ```
   python3 write_slurm_alignment.py --input_folder trimmed_reads --output_folder star_aligned --aligner_type star --genome_idx ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --filter_idx ~/umms-RNAlabDATA/Software/genome_indices/contaminants/bowtie2_contam_index/SL_contaminants/contaminants_index --library RF --two_pass --emit_dedup SBATCHSubArr-DEDUP.sbatch --email uniqname@umich.edu --slurm_acct cweidman99 --walltime 10:00:00 --mem 50000
   ```
   **Output:** 📄 `SBATCHSubArr-Align-STAR.sbatch`
3. In Bash, run the following commands to obtain STAR aligned files:
   ```
   conda activate B-PSID
   sbatch SBATCHSubArr-Align-STAR.sbatch
   ```
   **Output:** 📁 `star_aligned`
4. After alignment concludes, deduplicate reads with the SBATCH that was automatically generated in the same directory:
   ```
   sbatch SBATCHSubArr-DEDUP.sbatch
   ```
   **Output:** Deduplicated BAM files with the "dedup" suffix (excluding file extension). These are saved in 📁 `star_aligned`.

### When do I use this script?
* After obtaining trimmed reads, run this in the same directory.

### Understanding the SBATCH
```
python3 -u run_align.py --input trimmed_reads --output star_aligned --aligner star --index ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 -C 2 -L RF -S KEH-Rep1-7LKO-HEK293T-Cyto-NBS --filter_index ~/umms-RNAlabDATA/Software/genome_indices/contaminants/bowtie2_contam_index/SL_contaminants/contaminants_index -T --emit_dedup_slurm SBATCHSubArr-DEDUP.sbatch
```
* **--input:** Name of folder containing trimmed reads.
* **--output:** Name of folder for STAR aligned outputs.
* **--aligner:** Specify HISAT2 or STAR aligner.
* **--index:** HISAT2 index prefix or directory to STAR genome index.
* **-C:** Number of CPUs
* **-L:** Library type.
* **-S:** Sample prefix.
* **--filter_index:** If set, align reads to given bowtie2 index and retain unaligned reads for downstream processing.
* **-T:** Enable two-pass STAR alignment, which is recommended for transcript discovery (see [manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)).
* **--emit_dedup_slurm:** Name of output SBATCH file for deduplication. Use a custom name or simply just `SBATCHSubArr-DEDUP.sbatch`.

### Contaminant removal
* `contaminants.fa` contains human rRNA, tRNA, snoRNA, and snRNA sequences sourced from anRNAlab and various online sources. You can view the individual FASTA files used in `fasta.zip`.
  
  * rRNA online source: [fallerlab](https://github.com/fallerlab/ARF/blob/main/rRNAs/4V6X_human_rRNAs.fa)
  * tRNA online source: [GtRNAdb](https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html)
  * snoRNA online source: [snoRNABase](https://www-snorna.biotoul.fr/browse.php)
 
* The contaminants index is pre-built and located in its permanent directory `~/umms-RNAlabDATA/Software/genome_indices/contaminants/bowtie2_contam_index/SL_contaminants`.
  * If these files are lost, they can be accessed via 📁 `rm_contam` > `contaminants_index` in this repository.
 
* Alternatively, the contaminants index can be manually created.
  * Go into 📁 `rm_contam` > `manual`, then copy `contaminants.fa`, `build_index.py`, and `build_index.sbatch` into your directory. Finally, run the SBATCH file.

### STAR index
* The STAR hg38 genome index is pre-built and located in its permanent directory `~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38`.
* Alternatively, the genome index can be manually created.
  1. Create an empty folder called 📁 `star_index_hg38` in your directory.
  2. In GLC, navigate to `~/umms-RNAlabDATA/Software/genome_indices/hisat2_hg38/hg38p14_tran/` and copy the following 4 files into 📁 `star_index_hg38`:

     * `GCF_000001405.40_GRCh38.p14_exons.txt`
     * `GCF_000001405.40_GRCh38.p14_genomic.fa`
     * `GCF_000001405.40_GRCh38.p14_genomic.gtf`
     * `GCF_000001405.40_GRCh38.p14_splice_sites.txt`
       
  3. Go to 📁 `star_index` in this repository, then copy `star_index.py` and `star_index.sbatch` into your directory. Finally, run the SBATCH file. 

---
### Citations
* `run_align.py` by Chase Weidmann. If you have any questions, please reach out to [chaseaw](https://github.com/chaseaw).
