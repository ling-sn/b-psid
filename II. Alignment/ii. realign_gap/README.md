# PART II: Alignment
> ### Table of contents
> 1. [STAR alignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/i.%20star_alignment#table-of-contents)
> 2. [Realignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/ii.%20realign_gap#table-of-contents) :star: **YOU ARE HERE!**
> 3. BAM splitting
## Realignment
### Starter files
<img src="https://github.com/user-attachments/assets/a35085c4-a14d-40cf-a4e9-832868e20e7a" width="400"/>

### Overview
* STAR alignment optimizes both speed and efficiency, which can lead to some misalignments. Therefore, we use the Smith-Waterman algorithm for slower but more accurate mapping around areas we expect to see deletions.
> "Alignment software usually uses a seed alignment algorithm to increase alignment speed; however, this also affects pairwise alignment accuracy, especially for bases near deletion signatures. To solve this, we integrated the Smith-Waterman local alignment algorithm into the pipeline for realignment. Reads that contained any mismatch, deletion, insertion, soft-clip or splicing were further processed by the realignment tool in the BID-pipe package. By setting the penalty of gap open and gap extension as −3 and −2, respectively, deletion signatures can have a higher priority in the alignment" (_Zhang et al., 520_). 

### Instructions
1. Excluding the 📁 `star_aligned` folder (which should already exist), upload the remaining starter files to your GLC directory.
2. Create an SBATCH by running `write_slurm_realignment.py` in Bash with the following input commands:
   * **--input_folder:** 

   See example:
   ```
   Text
   ```
   **Output:** 📄 Text
3. Text

### When do I use this script?
* After obtaining deduplicated BAMs, run this in the same directory.

### Understanding the SBATCH
```
text
```
* **--input:**

---
### Citations
* Realignment code adapted from [realignGap](https://github.com/y9c/pseudoU-BIDseq/blob/main/bin/realignGap) by Ye Chang
