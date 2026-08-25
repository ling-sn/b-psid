# PART II: Alignment
> ### Table of contents
> 1. [STAR alignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/i.%20star_alignment#table-of-contents)
> 2. [Realignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/ii.%20realign_gap#table-of-contents)
> 3. [BAM splitting](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/iii.%20bam_splitting#table-of-contents) :star: **YOU ARE HERE!**
## BAM splitting
### Starter files
<img src="https://github.com/user-attachments/assets/82916a05-55eb-431b-a59d-f2895005d409" width="400"/>

### Overview
* Since UNUAR sites are located on either the forward or reverse strand, we employ BAM splitting to <ins>**only**</ins> count base/deletions from one strand for each site.

* We obtain `fwd.bam` and `rev.bam` for each BAM that we are interested in investigating downstream.

### Instructions
1. Excluding the 📁 `realignments` folder (which should already exist), upload the remaining starter files to your GLC directory.
2. Create an SBATCH by running `write_slurm_bam_split.py` in Bash with the following input commands:
   * **--email:** Email that will be notified when SLURM task begins/ends.
   * **--slurm_acct:** SLURM account.
   * **--walltime:** Amount of time allocated for job.
   * **--mem:** Amount of memory allocated for job.
   * **--library:** Specify library type. (Default: RF)
  
   See example:
   ```
   python3 write_slurm_bam_split.py --email uniqname@umich.edu --slurm_acct cweidman99 --walltime 3:00 --mem 3 --library RF
   ```
   **Output:** 📄 `split.sbatch`
3. In Bash, run the following commands to obtain the split BAMs:
   ```
   conda activate B-PSID
   sbatch split.sbatch
   ```
   **Output:** Split BAM files (`fwd.bam` and `rev.bam`). These will be saved in the same locations as the original BAMs, _i.e.,_ in the subfolders of 📁 `realignments`.

### When do I use this script?
* After obtaining realignments, run this in the same directory.

### Understanding the SBATCH
```
python3 split.py --bam_folder KEH-Rep1-7LKO-HEK293T-Cyto-NBS --library_type RF
```
* **--bam_folder:** Name of folder containing realigned reads.
* **--library_type:** Specify type based on library preparation method. (Choices = RF, FR)

### Explanation of fwd/rev outputs
* Our libraries are prepared with the NEBNext strand-specific protocol. Because this uses the dUTP method, our reads are in the form "RF".
  * In other words, the first-in-pair is recognized as reverse, while the second-in-pair is recognized as forward or in the direction of the transcript (reference/FASTA).

* We can use samtools flags to obtain specific alignments from each BAM file, and then combine them into separate "forward strand" and "reverse strand" files.
  * <ins>Forward strand</ins>
    * Second-in-pair maps to the forward strand
      ```
      samtools view -b -f 128 -F 16 $DATA > fwd1.bam
      ```
      * `-f 128` = Second-in pair
      * `-F 16` = Do not include read reverse strand flag
        
    * First-in-pair maps to the reverse  strand
      ```
      samtools view -b -f 80 $DATA > fwd2.bam
      ```
      * `-f 80` = Read reverse strand, first-in-pair
        
  * <ins>Reverse strand</ins>
    * Second-in-pair maps to the reverse strand
      ```
      samtools view -b -f 144 $DATA > rev1.bam
      ```
      * `-f 144` = Read reverse strand, second-in-pair
        
    * First-in-pair maps to the forward strand
      ```
      samtools view -b -f 64 -F 16 $DATA > rev2.bam
      ```
      * `-f 64` = First-in-pair
      * `-F 16` = Do not include read reverse strand flag
        
* Altogether, we obtain these two outputs after BAM splitting:
  * `fwd.bam` = All reads have "Pair orientation = F2R1" in IGV
  * `rev.bam` = All reads have "Pair orientation = F1R2" in IGV

---
### Citations
* BAM splitting code adapted from a [Biostars](https://www.biostars.org/p/92935) forum post by Istvan Albert.
* UC Davis DNA Technologies Core. _Which strand is sequenced for my strand-specific RNA-seq data?_ Retrieved April 27, 2026, from https://dnatech.ucdavis.edu/faqs/which-strand-is-sequenced-for-my-strand-specific-rna-seq-data
