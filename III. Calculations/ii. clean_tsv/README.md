# PART III: Calculations
> ### Table of contents
> 1. [Calculating deletion rates](https://github.com/ling-sn/b-psid/tree/main/III.%20Calculations/i.%20calculate_dr#table-of-contents)
> 2. [Cleaning TSVs](https://github.com/ling-sn/b-psid/blob/main/III.%20Calculations/ii.%20clean_tsv/README.md#table-of-contents) :star: **YOU ARE HERE!**
## Cleaning TSVs
### Starter files
<img src="https://github.com/user-attachments/assets/9d3f6eef-65df-42f3-9490-203167539b1d" width="400"/>

### Overview
* Cleans TSV outputs from D.R. calculations via merging and calculation of summary statistics / $p$-values.
  1. Merges replicates for each sample and calculates weighted averages/standard deviations.
  2. Merges BS and NBS samples, calculates $p$-values with Fisher's Exact Test, and filters UNUAR sites by WT $p$-value ($p \le 0.05$).
  3. Merges sample groups (WT/7KO/7LKO) together to obtain all UNUAR sites in cytoplasm (`Cyto.tsv`) or nucleus (`Nuc.tsv`).

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
**Example:** From 📁 `calculations`, start from the 📁 `7KO-Cyto` subfolder.

<img src="https://github.com/user-attachments/assets/d5abe6a1-9b5e-4134-ba77-47c581f9006d" width="400"/>

1. Merge replicates with a full join to obtain separate files for the $\color{darkseagreen}BS$ and $\color{lightcoral}NBS$ samples.
   
   <img src="https://github.com/user-attachments/assets/57eea18d-98c5-4bc8-ab32-c4a89fb27139" width="400"/>

   * For each file (7KO-Cyto-BS & 7KO-Cyto-NBS), calculate StdDeletionRate and TotalAvgDeletionRate (weighted average) columns.
     
     $$\text{Weighted avg.} = \frac{\text{Sum of replicate deletions}}{\text{Sum of replicate total coverage}}$$
2. Merge BS and NBS files with a full join to obtain one file per sample group (_e.g.,_7KO-Cyto).
   * For each sample group, calculate and filter by WT $p$-values ($p \le 0.05$).
3. Merge WT, 7KO, and 7LKO files with an inner join to obtain separate outputs `Cyto.tsv` and `Nuc.tsv`.
   <img src="https://github.com/user-attachments/assets/0e191601-86e1-4cb0-beb5-946267bc97e1" width="400"/>

### Fisher's Exact Test
<img src="https://github.com/user-attachments/assets/d0150acf-56d2-491e-ae25-692627d1c275" width="300"/>

* For each UNUAR site, we used this test to determine if there were significant differences between $\color{darkseagreen} \verb|%| \textbf{ deletions in BS}$ vs. $\color{lightcoral}\verb|%|\textbf{ deletions in NBS}$.
  * The test was in the "less-than" direction because we wanted to minimize T bases (and implicitly maximize deletions) in our BS-treated samples.
* For sites on the reverse strand (-), T bases were counted as A bases.
  * If a site was on the reverse strand, the A count was used for the test.
  * If a site was on the forward strand, the T count was used for the test.
* However, while Fisher’s Exact Test assumes that marginal sums are fixed, these assumptions were not met in our experiment because sequencing counts from the genome were random.
  * Therefore, this test will be replaced with rate ratio calculations in the future (sometime later in 2026). 

---

### Citations
* `clean_tsv.py` by Sonia Ling. If you have any questions, please reach out to [ling-sn](https://github.com/ling-sn).
