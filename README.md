# B-PsiD: A pipeline for bisulfite pseudouridine detection
### Abstract
<img src="https://github.com/user-attachments/assets/caba8228-ab7d-4fe8-bbeb-5ce738c4f23e" width="600"/>

> Pseudouridine (Ψ), an isomer of uridine, is a prevalent RNA modification driven by pseudouridine synthase (PUS) enzymes. Dysregulation of one PUS enzyme, PUS7, is implicated in various disorders such as autism and microcephaly. PUS7 modifies the central U within UNUAR motifs of RNA, but only a small fraction is modified within cells, and little is known about how PUS7 chooses which UNUAR substrates to modify. In turn, our laboratory is interested in defining PUS7 substrate specificity for mRNAs, which requires quantification of the number of Ψ sites in cells. The Bisulfite Pseudouridine Detection (B-PsiD) pipeline is therefore being developed for integration into our laboratory’s workflows. The core of this pipeline extends upon the bisulfite-induced deletion sequencing (BID-Seq) method from Zhang et al., which relies on bisulfite-based conversion of pseudouridines into ribose ring-openings that are read as deletions during reverse transcription [*Nature Protocols, 19, 517–538 (2024)*]. B-PsiD adapts BID-Seq using alternative library preparation methods, enabling orthogonal validation and quantification of the ratio of Ψ:U at every UNUAR site. Logistically, B-PsiD mirrors the five main parts of BID-Seq: (i) quality control and merging of the Illumina sequencing reads; (ii) contaminant RNA removal to isolate the mRNA; (iii) alignment to the human genome; (iv) realignment for more accurate deletion-only mapping; and finally, (v) the calculation of real Ψ modification stoichiometry at every UNUAR site. Once all high confidence modified Ψ sites are identified within mRNAs, our laboratory will compare the RNA structure around modified versus unmodified sites via chemical probing and mutational profiling.

### Instructions
1. Set-up the B-PSID conda environment in GLC.
   
   * Download `create_env.sbatch` from this repository and upload it to your GLC directory.
   * In Bash, run the command `sbatch create_env.sbatch`. This should take a maximum of 5 minutes to run.
     
3. In your GLC directory, create a 📁 `raw_fastqs` folder containing **all** the raw fastq.gz files from Illumina sequencing.
   
   * Do not add any subfolders.
     
5. Come back to this repository and go to each folder in the order listed below. Then, follow the instructions in their READMEs.
   
   * 📁 `I. Adapter Trimming`
   * 📁 `II. Alignment`
   * 📁 `III. Calculations`
     
### Citations
* Kim, K. Q., Burgute, B. D., Tzeng, SC, et al. (2022). N1-methylpseudouridine found within COVID-19 mRNA vaccines produces faithful protein products. _Cell Reports_, 40(9), 1–14. https://doi.org/10.1016/j.celrep.2022.111300
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5
---
### To-do list
* **bam_splitting**
  - [ ] Add visual in README
* **calculate_dr**
  - [x] Only count single-nt deletions instead of multi-nt ones (default behavior)
* **clean_tsv**
  - [ ] Replace Fisher's Exact Test with rate ratio calculations (Poisson) and LASSO regression
* **Other**
  - [ ] Create Snakemake pipeline for B-PsiD
