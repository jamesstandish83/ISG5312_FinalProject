# ISG5312 Final Project: Fire Blight QTL Mapping in Apple

## Overview
This project recreates and extends the fire blight resistance mapping workflow from:

Peil A, Hübert C, Wensing A, Horner M, Emeriewen OF, Richter K, Wöhner T, Chagné D, Orellana-Torrejon C, Saeed M, Troggio M, Stefani E, Gardiner SE, Hanke MV, Flachowsky H, Bus VGM. *Mapping of fire blight resistance in Malus ×robusta 5 flowers following artificial inoculation.* BMC Plant Biology. 2019;19:532. https://doi.org/10.1186/s12870-019-2154-7

The original study identified a major fire blight resistance QTL in *Malus ×robusta 5* and reported that the resistance signal was consistent across inoculation modes, including flower and shoot inoculation. This project reproduces the analysis goal in R using **R/qtl** and **ASMap** instead of the JoinMap/MapQTL workflow used in the paper. citeturn730037view0turn730037view3turn730037view2

The biological question is whether fire blight resistance in this apple population is linked to one or more genomic regions and whether the strongest signal agrees with the published result. citeturn730037view0

## Workflow summary
The project is organized into two scripts:

1. `GenomeStudio_OutPut_QC.R`
   - reads the raw GenomeStudio genotype file and phenotype mapping file
   - standardizes sample IDs
   - matches genotype and phenotype records
   - converts missing genotype calls to `NA`
   - calculates sample-level and SNP-level call rates
   - writes cleaned mapping and QC summary files

2. `FireBlight_Rqtl_ASMap_Cleaned_Workflow.R`
   - reads the cleaned/filtered input files
   - filters low-quality SNPs by call rate
   - builds a physical-position map and ASMap genetic map
   - checks for allele flips
   - runs marker regression, manual Kruskal-Wallis, Haley-Knott, and EM scans
   - runs permutation tests for marker regression and EM interval mapping
   - extracts the primary QTL peak and genotype-effect summaries
   - attempts MQM analysis when the cross type supports it

R/qtl provides the core QTL scanning functions, including `scanone`, and ASMap provides linkage-map construction and diagnosis for R/qtl cross objects using the MSTmap algorithm. citeturn730037view2turn422806search10turn730037view3

## Repository structure
A clean repo layout should look like this:

```text
ISG5312_FinalProject/
├── README.md
├── scripts/
│   ├── GenomeStudio_OutPut_QC.R
│   └── FireBlight_Rqtl_ASMap_Cleaned_Workflow.R
├── data/
│   ├── raw/
│   └── processed/
├── results/
│   ├── qc/
│   └── qtl/
├── docs/
└── logs/
```

## Software requirements
The analysis was written in R and uses these packages:

- `tidyverse`
- `ggplot2`
- `qtl`
- `ASMap`
- `cli`

## Input files
The scripts expect these input files:

- `R5xRG_Genome_Studio_Carolina_08_2016.txt` — raw GenomeStudio genotype matrix
- `Phenotype_Genotype_Mapping.txt` — phenotype file with sample IDs and disease scores
- `SNP_CallRates.txt` — SNP call-rate summary used to filter low-quality markers
- `All_Markers.tsv` — marker annotation table
- `9K_SNP_positions_DH_genome.txt` — SNP physical-position reference table
- `Apple-IM-F1-Mr5_cM-Map.tsv` — published reference map file used for annotation support

## Reproducibility notes
This is not a byte-for-byte reproduction of the original JoinMap/MapQTL workflow. It is a transparent R-based recreation of the same biological analysis goal, using `R/qtl` and `ASMap` to map fire blight resistance in the Mr5-derived population. The scripts are written so another person can follow the logic, inspect intermediate files, and reconstruct the analysis.

## How to run
Run the QC script first, then the QTL script.

```r
source("scripts/GenomeStudio_OutPut_QC.R")
source("scripts/FireBlight_Rqtl_ASMap_Cleaned_Workflow.R")
```

If you keep the current absolute paths, the scripts will read and write in the Windows folders shown inside the scripts. For a cleaner repository, update the paths so they point to `data/raw/`, `data/processed/`, and `results/`.

## What the QC outputs mean
- `Sample_CallRate_Histogram.png` / `Sample_CallRate_Histogram.pdf` — distribution of call rates across samples; used to spot poor-quality samples.
- `SNP_CallRate_Histogram.png` / `SNP_CallRate_Histogram.pdf` — distribution of call rates across SNP markers; used to spot low-quality markers.
- `Sample_CallRates.txt` — per-sample call rate table.
- `SNP_CallRates.txt` — per-marker call rate table.
- `Genotype_Distribution.txt` — counts and proportions of AA, AB, BB, and missing calls.
- `Phenotype_Genotype_Mapping_Cleaned.txt` — cleaned sample-to-phenotype mapping used for downstream analysis.
- `QC_Summary.txt` — overall summary including sample count, SNP count, and call-rate statistics.

The QC stage shows whether the data are strong enough for downstream mapping and whether any samples or markers should be removed before analysis.

## What the QTL outputs mean
### Cross objects and map objects
- `cross_phys_final.rds` — R/qtl cross object using the physical-position map.
- `cross_qtl_f2_final.rds` — preserved F2 cross object used to keep the QTL analysis compatible with F2 methods.
- `cross_qtl_for_analysis_final.rds` — final analysis object after map alignment and genotype-probability calculation.
- `cross_phys_summary.txt` — summary of the physical-position cross.
- `cross_asmap_summary.txt` — summary of the ASMap cross, if ASMap completed successfully.
- `cross_qtl_summary.txt` — summary of the final QTL analysis object.
- `physical_map_tbl.rds` — physical map table as an R object.
- `genetic_map_tbl.rds` — final ASMap genetic map table as an R object.
- `ASMap_Genetic_Map_Table.csv` — tabular genetic map output.
- `ASMap_Genetic_Map_Summary.txt` — text summary of the ASMap map.
- `ASMap_Chromosome_Summary.csv` — marker count by chromosome.
- `ASMap_Chromosome_Subgroup_Summary.csv` — marker count by chromosome subgroup.
- `Physical_Map_Table.csv` — marker coordinates in base pairs and megabases.

### Mapping and scan results
- `ScanOne_MarkerRegression_Results.csv` — full marker-regression genome scan.
- `ScanOne_KruskalWallis_Results.csv` / `Manual_KruskalWallis_Results.csv` — manual Kruskal-Wallis scan results.
- `ScanOne_HaleyKnott_Results.csv` — Haley-Knott interval-mapping results.
- `ScanOne_EM_Results.csv` — EM interval-mapping results.
- `perm_scan_mr.rds` — marker-regression permutation test object.
- `perm_scan_em.rds` — EM permutation test object.
- `Permutation_Threshold_MR.txt` — 5% genome-wide threshold for marker regression.
- `Permutation_Threshold_EM.txt` — 5% genome-wide threshold for EM interval mapping.
- `SingleQTL_MarkerRegression.png` — marker-regression scan plot.
- `SingleQTL_EM.png` — EM scan plot.
- `SingleQTL_HK.png` — Haley-Knott scan plot.
- `SingleQTL_KruskalWallis.png` / `SingleQTL_KruskalWallis.pdf` — Kruskal-Wallis scan plot.
- `R_qtl_QTL_Scans.pdf` — combined multi-panel QTL scan report.

### Peak and effect summaries
- `Primary_MarkerRegression_Peak.csv` — strongest marker-regression peak.
- `Significant_QTL_Peaks_EM.csv` — significant EM peaks based on permutations.
- `Top_LOD_Sites_EM.csv` — top EM LOD sites.
- `Chromosome_MaxLOD_Summary_EM.csv` — maximum EM LOD per chromosome.
- `QTL_Peak_Genotype_Effects.csv` — genotype-group means and variation at the primary MR peak.
- `QTL_Peak_Quartile_Association.csv` — genotype counts split by phenotype quartiles.
- `QTL_Peak_Genotype_Association.pdf` — boxplot/jitter plot of phenotype by genotype at the primary QTL peak.
- `qtl_effects.rds` — serialized object containing all peak-effect summaries and plotting data.

### Marker annotation and candidate marker outputs
- `Marker_ID_Lookup.csv` — maps SNP index values to the generated marker IDs.
- `Marker_Annotation_Merged.csv` — merged marker annotation table combining the annotation files and the QTL map.
- `Candidate_Markers_For_Figures.csv` — markers near the strongest QTL peak used for reporting and figure support.
- `Physical_Map_Table.csv` — marker positions in physical coordinates.
- `Top_ErrorLOD_Genotypes.csv` — markers with high genotype error LOD values.
- `CheckAlleles_Suspects.csv` — markers flagged as possible allele flips.
- `Dropped_Markers_Invalid_Chr_or_Position.csv` — markers removed because of invalid chromosome or position values.

### Phenotype summaries and plots
- `pheno.csv` — cleaned phenotype file used by R/qtl.
- `Phenotype_Class_Assignment.csv` — individual-level resistant/intermediate/susceptible assignment.
- `Phenotype_Class_Summary.csv` — count table for the phenotype classes.
- `Phenotype_Class_Boxplot.png` — boxplot of disease incidence by phenotype class.
- `Phenotype_Ranked_Barplot.png` — ranked disease incidence plot for the progeny.

### MQM outputs
- `mqm_cofactors.rds` — cofactor list for MQM, if the model stabilized.
- `mqm_result.rds` — MQM scan result, if the model stabilized.
- `MQM_Scan_Results.csv` — MQM scan table, if produced.
- `MQM_Scan.png` — MQM scan plot, if produced.
- `MQM_Status.txt` — plain-text note describing whether MQM completed or was skipped.

### Why these files were made
These outputs separate the project into clear stages. The QC files justify which data were kept. The cross and map files preserve the analysis state. The scan and peak files identify candidate QTL. The annotation and candidate-marker tables connect the statistical results back to the biology and the published map.

## Notes on the analysis strategy
- Marker regression was kept because it gives a direct genotype-by-phenotype association scan.
- Manual Kruskal-Wallis was added as a nonparametric complement for the F2 cross.
- Haley-Knott and EM interval mapping provide interval-based support for QTL location.
- ASMap is used for genetic map construction.
- MQM was attempted as a more advanced mapping step, but it may fail when the model becomes unstable; the script records that outcome instead of stopping the workflow.

## Troubleshooting
If the scripts fail, check the following first:

1. The file paths at the top of each script.
2. Whether the input files are present and named exactly as expected.
3. Whether the phenotype file contains the `CleanID` and `MeanDiseaseIncidence` columns.
4. Whether the genotype file has the expected GenomeStudio structure.
5. Whether all required R packages are installed.

## GitHub update workflow
To update the repository after editing the scripts and README locally:

```bash
cd /path/to/ISG5312_FinalProject
git status
git add README.md scripts/GenomeStudio_OutPut_QC.R scripts/FireBlight_Rqtl_ASMap_Cleaned_Workflow.R
git add data/raw data/processed results docs logs
git commit -m "Update fire blight QC and QTL workflow"
git push origin main
```

GitHub documents `git push` as the command used to send local commits to the remote repository, commonly with `git push origin main`, and the repository README is intended to explain what the project does, how to get started, and how to use it. citeturn730037view4turn730037view5

If you prefer the GitHub website instead of the command line, open the repository, use **Add file** or the pencil icon to edit `README.md`, and then commit the changes through the web interface. citeturn546190search3turn546190search6

## Final project reflection
If I were starting over, I would:

- move all scripts to relative paths so the project runs from the repository root
- place raw data, cleaned data, and final results into separate folders
- document every package dependency in a single setup section
- save a small workflow diagram alongside the README
- keep a short changelog of each analysis decision, especially around F2 handling, ASMap map replacement, and MQM fallback behavior

That would make the project easier to reproduce, easier to review, and easier to extend.
