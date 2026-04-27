# PART III: Calculations
> ### Table of contents
> 1. [Calculating deletion rates](https://github.com/ling-sn/b-psid/tree/main/III.%20Calculations/i.%20calculate_dr#table-of-contents) :star: **YOU ARE HERE!**
> 2. Cleaning TSVs
## Calculating deletion rates
### Starter files
<img src="https://github.com/user-attachments/assets/0010b96f-b977-4f74-b812-79632984dc1d" width="400"/>

### Overview
* Opens up `fwd.bam` and `rev.bam` separately, and then obtains base and deletion counts for each.
* After each process concludes, count results are concatenated into one dataframe and deletion rates are calculated.
* No TSV filtering occurs at this step of the pipeline.

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
python3 calculate_dr.py --folder_name KEH-Rep1-7LKO-HEK293T-Cyto-NBS --fasta ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa
```
* **--folder_name:** Name of folder containing split BAMs.
* **--fasta:** Path to reference FASTA that was used to build the STAR genome index. This will always be `~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa`.

### UNUAR datasets
* The following datasets were first merged with a left-join:
  * `UNUAR_motif_sites_mRNA_hg38p14.tsv` contains the GenBank accession number (_Chrom_) and genomic coordinate of the modified base (_GenomicModBase_) for all UNUAR sites in the human genome.
  * `SupplementaryTable1.xlsx` contains the best-fit parameters for the calibration curves of 256 UNUAR motifs (_Zhang et al., 533_).
* To accommodate the split BAM structure, this merged dataset was filtered to only contain sites in the 3UTR region, and then split into two datasets by strand type (+/-). 
  * <ins>**Permanent directories**</ins>: These paths are already included in the code by default, so nothing additional needs to be done.
    * **fwd_unuar:** `~/umms-RNAlabDATA/Software/B-PsiD_tools/B-PsiD_UNUAR_motif_sites_mRNA_hg38_fwd.tsv`
    * **rev_unuar:** `~/umms-RNAlabDATA/Software/B-PsiD_tools/B-PsiD_UNUAR_motif_sites_mRNA_hg38_fwd.tsv`
  * <ins>**Manual**</ins>: If the fwd_unuar and rev_unuar files are missing but you have the original files, run the following code in Bash to create them in your current directory.
    
    ```
    python3 split_unuar_tsv.py
    ```

* However, there were too many rows with duplicate Chrom and GenomicModbase pairs, so this led to duplicate counts in the output TSVs based on how the counting algorithm was designed.
  * As a solution, these duplicate pairs were dropped when **fwd_unuar** and **rev_unuar** were read in as dataframes in `calculate_dr.py`.
  * Because every Chrom and GenomicModBase pair was associated with multiple TranscriptIDs, the transcripts were grouped in a separate column "OtherAssocTranscripts" prior to the pairs (excluding the first occurrence) being dropped.
* Filtering statistics:
  * **Original left-joined TSV:** 3425693 rows
  * **fwd_unuar:** 1057146 rows
  * **rev_unuar:** 1014095 rows
  * **Updated fwd_unuar (no dup):** 187048 rows
  * **Updated rev_unuar (no dup):** 181396 rows

### Real rate calculation
* The best-fit parameters from `SupplementaryTable1.xlsx` were plugged into the following equation to estimate the fraction of actual Ψ modification, which is also referred to as "RealRate" in the script.
  
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
  
* The observed deletion rate ($$r$$) at each UNUAR site was calculated with the following formula:

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
* Afterwards, join fwd and rev counts into one dataframe

<img width="879" height="573" alt="image" src="https://github.com/user-attachments/assets/75a61bd6-c499-4db0-88f0-1cdb91848aca" />

* Use stat_variation_strand() to return an iterator. Only append pileup statistics to results dictionary if pos matches



### Citations
https://github.com/alimanfoo/pysamstats/blob/master/pysamstats/pileup.py
