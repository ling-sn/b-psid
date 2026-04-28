# PART III: Calculations
> ### Table of contents
> 1. [Calculating deletion rates](https://github.com/ling-sn/b-psid/tree/main/III.%20Calculations/i.%20calculate_dr#table-of-contents)
> 2. [Cleaning TSVs](https://github.com/ling-sn/b-psid/blob/main/III.%20Calculations/ii.%20clean_tsv/README.md#table-of-contents) :star: **YOU ARE HERE!**
## Cleaning TSVs
### Starter files
<img src="https://github.com/user-attachments/assets/9d3f6eef-65df-42f3-9490-203167539b1d" width="400"/>

### Overview
* Cleans TSV outputs from D.R. calculations via merging and calculation of summary statistics / $p$-values.
  * Merges replicates for each sample and calculates weighted averages/standard deviations.
  * Merges BS and NBS samples, calculates $p$-values with Fisher's Exact Test, and filters UNUAR sites by WT p-value.
  * Merges sample groups (WT/7KO/7LKO) together to obtain all UNUAR sites in cytoplasm (`Cyto.tsv`) or nucleus (`Nuc.tsv`).

### Instructions
1. Excluding the 📁 `calculations` folder (which should already exist), upload the remaining starter files to your GLC directory.
   
3. Edit the following lines in 📄 `run_clean_tsv.sbatch` to match your experiment:
   * `#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu`
   * `#SBATCH --time=30:00`
   * `#SBATCH --mem=30000m`
     
4. In Bash, run the following commands to obtain the output TSVs:
   ```
   conda activate B-PSID
   sbatch run_clean_tsv.sbatch
   ```
   **Output:** 📁 `merged` > 📁 `final_outputs` > 📄 `Cyto.tsv` and `Nuc.tsv`

### When do I use this script?
* After obtaining the deletion rate calculations, run this in the same directory.

### Explanation of merging

### Fisher's Exact Test
* Fisher’s Exact Test assumes that the marginal sums are fixed, but these assumptions are not met in our experiment because sequencing counts from the genome are random.
  * Therefore, an explanation of this test will be omitted as it will be replaced with rate ratio calculations (July 2026). 

### Citations
* `clean_tsv.py` by Sonia Ling. If you have any questions, please reach out to [ling-sn](https://github.com/ling-sn).
* Heger, A., Marshall, J., & Jacobs, K. _Pysam manual_. Retrieved April 28, 2026, from https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.pileup

---
### Overview
* For Fisher’s Exact Test, detect strand type and use correct T-bases. Check if UNUAR site is on negative or positive strand. If on negative strand, use A count for test; if on positive strand, use T count for test
* I will omit explanation about Fisher's Exact Test since it will be replaced when I come back in July.
* Secondary alignments are automaticalyl filtered out. Secondary alignments are reads that “align reasonably well to more than one place.” they are  ignored by downstream tools in order to avoid overcounting. In pysam (base tool of pysamstats), secondary alignments are ignored by DEFAULT because they could artificially inflate base/del counts and we only want high accuracy in our counts.
  * During STAR alignment, the primary alignment is the first top-scoring alignment that STAR finds. It is deterministic, meaning that “everytime it
    maps to the same read, the same primary alignment will be marked.”
* TODO: Only count single-nt deletions, use rate ratios instead of Fisher's Exact Test
* The problem with dropping all duplicate Chrom and GenomicModBase pairs and only selecting the first one is that we lose all of the other TranscriptID information from the “duplicates.”
  * As a solution, add a column that stores the remaining TranscriptIDs in a set/list, separated by commas, and name this column “OtherAssocTranscripts.”
  * We drop these pairs when we read in the UNUAR TSV
  * Go back to my presentations to see justifications for hwy
* Only kept WT p-value filtering (p ≤ 0.05)
* Ordered UNUAR sites in descending order by TotalAvgDeletionRate, a weighted average calculation (sum of deletions / sum of coverage)

### Explanation of joining
* Reference 4/17/26
