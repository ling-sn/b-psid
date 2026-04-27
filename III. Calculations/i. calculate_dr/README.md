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

<img width="879" height="573" alt="image" src="https://github.com/user-attachments/assets/75a61bd6-c499-4db0-88f0-1cdb91848aca" />

* Use stat_variation_strand() to return an iterator. Only append pileup statistics to results dictionary if pos matches

* Many rows with duplicate Chrom and GenomicModBase pairs, leading to rows with duplicate base/del counts. Hence, we removed duplicate GenomicModBase + Chrom pairs when we read in the UNUAR tsvs
  * Start with original left-joined TSV: 3425693 rows
  * Obtain TSV of rows in region "3UTR" with strand “+” (fwd): 1057146 rows
  * Obtain TSV of rows in region "3UTR" with strand “-” (rev): 1014095 rows
* After dropping duplicate rows in fwd and rev:
  * Updated fwd: 187048 rows
  * Updated rev: 181396 rows
* Every position is associated with multiple TranscriptIDs, so we group them in a separate column "OtherAssocTranscripts" before dropping the pairs and keeping only the first one



### Citations
https://github.com/alimanfoo/pysamstats/blob/master/pysamstats/pileup.py
