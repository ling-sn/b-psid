### Overview
* Use pysamstats to obtain base/del counts for each bam file, and store in dictionary
* stat_pileup() returns the iterator recs, where “each dict holds data for a single genome position”
* Call pileup traversal once per chromosome by leveraging min/max values of genomic coordinates for that chrom 
  * Group sites by chromosome like {Chr1: (Site1, Site2, Site3)}
  * Store base counts only at sites of interest in a dictionary
  * Append dictionary to list
* Repeat for all chromosomes, and at end, convert list to dataframe
* Permanent directory for UNUAR tsvs
* How did we obtain the fwd and rev files? -> split_unuar_tsv
  * fwd_unuar = Dataframe with Ψ sites on forward strand (+)
  * rev_unuar = Dataframe with Ψ sites on reverse strand (-)
* Open up fwd and rev BAMs separately and obtain base/del counts for each
* Afterwards, join fwd and rev counts into one dataframe

### Citations
https://github.com/alimanfoo/pysamstats/blob/master/pysamstats/pileup.py
