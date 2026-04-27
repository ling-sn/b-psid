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
