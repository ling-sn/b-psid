### Overview
* Use pysamstats to obtain base/del counts for each bam file, and store in dictionary
* stat_pileup() returns the iterator recs, where “each dict holds data for a single genome position”
* Call pileup traversal once per chromosome by leveraging min/max values of genomic coordinates for that chrom 
  * Group sites by chromosome like {Chr1: (Site1, Site2, Site3)}
  * Store base counts only at sites of interest in a dictionary
  * Append dictionary to list
* Repeat for all chromosomes, and at end, convert list to dataframe

### Citations
https://github.com/alimanfoo/pysamstats/blob/master/pysamstats/pileup.py
