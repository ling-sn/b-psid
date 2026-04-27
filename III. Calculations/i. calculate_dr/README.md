# PART III: Calculations
> ### Table of contents
> 1. [Calculating deletion rates](https://github.com/ling-sn/b-psid/tree/main/III.%20Calculations/i.%20calculate_dr#table-of-contents) :star: **YOU ARE HERE!**
> 2. Cleaning TSVs
## Calculating deletion rates
### Starter files
<img src="https://github.com/user-attachments/assets/0010b96f-b977-4f74-b812-79632984dc1d" width="400"/>

### Overview
* Opens up `fwd.bam` and `rev.bam` separately, and then obtains base and deletion counts for each.
* After each process concludes, count results are concatenated and deletion rates are calculated.

### Instructions
1. Excluding the 📁 `realignments` folder (which should already exist), upload the remaining starter files to your GLC directory.
2. Create an SBATCH by running `write_slurm_dr.py` in Bash with the following input commands:
   * **--email:** Email that will be notified when SLURM task begins/ends.
   * **--slurm_acct:** SLURM account.
   * **--walltime:** Amount of time allocated for job.
   * **--mem:** Amount of memory allocated for job.
   * **--fa:** Path to reference FASTA.
  
   See example:
   ```
   python3 write_slurm_dr.py --email uniqname@umich.edu --slurm_acct cweidman99 --walltime 1:00:00 --mem 5000 --fa ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa
   ```
   **Output:** 📄 `run_calculate_dr.sbatch`
3. In Bash, run the following commands to obtain the output TSVs:
   ```
   conda activate B-PSID
   sbatch run_calculate_dr.sbatch
   ```
   **Output:** 📁 `calculations`

### When do I use this script?
* After obtaining split BAMs, run this in the same directory.

### Understanding the SBATCH
```
Text
```
* **Text:** Path to reference FASTA that was used to build the STAR genome index. This will always be `~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa`.

### Datasets & calculations
* The following two datasets were merged with a left-join:
  * `UNUAR_motif_sites_mRNA_hg38p14.tsv` contains the GenBank accession number (_Chrom_) and genomic coordinate of the modified base (_GenomicModBase_) for all UNUAR sites in the human genome.
  * `SupplementaryTable1.xlsx` contains the best-fit parameters for the calibration curves of 256 UNUAR motifs (_Zhang et al., 533_). They are plugged into the equation below to estimate the fraction of actual Ψ modification, which is also referred to as "RealRate" in the script.
  
  $$f = \large{\frac{b-r}{c \cdot (b+s-s \cdot r -1)}}$$

  in which:

  <img src="https://github.com/user-attachments/assets/0f207398-ca2f-4963-9f65-c554dccd4a46" width="300"/>

  | Variable | Name | Meaning
  | --- | --- | --- |
  | $$f$$ | Real deletion rate $$(f > 0)$$ | Ψ modification stoichiometry |
  | $$r$$ | Observed deletion rate $$(0 < r < 1)$$ | Percentage of deletions at given `GenomicModBase` |
  | $$b$$ | Background deletion rate | Baseline deletion rates due to experimental conditions; `fit_b` |
  | $$s$$ | RT dropoff ratio | Ratio of times a site “falls out” when it gets hit by bisulfite; `fit_s` |
  | $$c$$ | Conversion ratio | Maximum amount of times a site can be turned into a deletion; `fit_c` |
  
* The observed deletion rate ($$r$$) at each UNUAR site is calculated with the following formula:

  $$\large{\frac{\text{Number of deletions}}{\text{Total amount of A, C, T, G, and deletions}}}$$
---
### Citations
* `calculate_dr.py` by Sonia Ling. If you have any questions, please reach out to [ling-sn](https://github.com/ling-sn).


---
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
