# PART II: Alignment
> ### Table of contents
> 1. [STAR alignment](https://github.com/ling-sn/b-psid/tree/main/II.%20Alignment/i.%20star_alignment#table-of-contents)
> 2. Realignment :star: **YOU ARE HERE!**
> 3. BAM splitting
## Realignment
### Starter files
<img src="https://github.com/user-attachments/assets/a952f05f-8174-4efc-bf71-19d750776d53" width="400"/>

### Overview
* STAR alignment optimizes both speed and efficiency, but can lead to some misalignments. Therefore, we use the Smith-Waterman algorithm for slower but more accurate mapping around areas we expect to see deletions.

### Instructions
1. Excluding the 📁 `star_aligned` folder (which should already exist), upload the remaining starter files to your GLC directory.
2. Create an SBATCH by running `write_slurm_alignment.py` in Bash with the following input commands:
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
