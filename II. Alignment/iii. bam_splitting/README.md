# PART II: Alignment
> ### Table of contents
> 1. [STAR alignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/i.%20star_alignment#table-of-contents)
> 2. [Realignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/ii.%20realign_gap#table-of-contents)
> 3. BAM splitting :star: **YOU ARE HERE!**
## BAM splitting
### Starter files
<img src="https://github.com/user-attachments/assets/82916a05-55eb-431b-a59d-f2895005d409" width="400"/>

### Overview
* Text

### Instructions
1. Create SBATCH
   ```
   python3 write_slurm_bam_split.py --email lingsn@umich.edu --slurm_acct cweidman99 --walltime 5:00 --mem 2000 --library RF
   ```

### When do I use this script?
* Text

### Understanding the SBATCH
```
Text
```
* Text

---
### Citations
* https://dnatech.ucdavis.edu/faqs/which-strand-is-sequenced-for-my-strand-specific-rna-seq-data
* https://www.biostars.org/p/92935/
---
### DELETE LATER
* Objective: If UNUAR site is on reverse strand, only consider base/del counts from the reverse strand?
* We use the NEBNext workflow for library preparation. The description for NEBNext Ultra II RNA states that directional (strand-specific) library prep uses the “dUTP method” (p. 12).
* When we use the “dUTP method” during second-strand synthesis, the strand with dUTP will not be amplified and thus preserve the strand information for RNA-seq.
* For paired-end sequencing, the forward read of the resulting sequencing data represents the “anti-sense strand” and the reverse read the “sense strand” of the gene.

### Overview
* Due to NEBNext library preparation, our reads are -RF. Meaning that first-in-pair → “reverse” and second-in-pair is in the direction of the transcript (i.e., reference/FASTA) → “forward”
* Using samtools flags allows us to obtain specific alignments from BAM files, and then combine them into separate “forward strand” and “reverse strand” files.
  * Forward strand:
    * Second-in-pair maps to the forward strand
      ```
      samtools view -b -f 128 -F 16 $DATA > fwd1.bam
      ```
      * -f 128 $=$ Second-in pair
      * -F 16 $=$ Do not include read reverse strand flag
    * First-in-pair maps to the reverse  strand
      ```
      samtools view -b -f 80 $DATA > fwd2.bam
      ```
      * -f 80 $=$ Read reverse strand, first-in-pair
  * Reverse strand:
    * Second-in-pair maps to the reverse strand
      ```
      samtools view -b -f 144 $DATA > rev1.bam
      ```
      * -f 144 $=$ Read reverse strand, second-in-pair
    * First-in-pair maps to the forward strand
      ```
      samtools view -b -f 64 -F 16 $DATA > rev2.bam
      ```
      * -f 64 $=$ First-in-pair
      * -F 16 $=$ Do not include read reverse strand flag
* Outputs
  * fwd.bam = Reads all have Pair orientation = F2R1
  * rev.bam = Reads all have Pair orientation = F1R2
