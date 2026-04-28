# PART III: Calculations
> ### Table of contents
> 1. [Calculating deletion rates](https://github.com/ling-sn/b-psid/tree/main/III.%20Calculations/i.%20calculate_dr#table-of-contents)
> 2. [Cleaning TSVs](https://github.com/ling-sn/b-psid/blob/main/III.%20Calculations/ii.%20clean_tsv/README.md#table-of-contents) :star: **YOU ARE HERE!**
## Cleaning TSVs
### Starter files
<img src="https://github.com/user-attachments/assets/0010b96f-b977-4f74-b812-79632984dc1d" width="400"/>

### Overview
* Text

### Instructions
1. Text

### When do I use this script?
* After obtaining calculated deletion rates, run this in the same directory.

### Understanding the SBATCH
```
Text
```
* **Text:** Text
---
### Citations
* `clean_tsv.py` by Sonia Ling. If you have any questions, please reach out to [ling-sn](https://github.com/ling-sn).



---

### Summary
- Cleans TSV outputs from D.R. calculations
- Formats data for graphing
- Outputs graphs

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

### Citations
* https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.pileup
* https://github.com/alexdobin/STAR/issues/1475
