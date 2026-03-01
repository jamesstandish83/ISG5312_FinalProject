# ISG5312_FinalProject

This project reproduces and analyzes the computational workflow from the published paper:

Shen, F., Huang, Z., Zhang, B., Wang, Y., Zhang, X., Wu, T., Xu, X., Zhang, X., & Han, Z. (2019). 
Mapping Gene Markers for Apple Fruit Ring Rot Disease Resistance Using a Multi-omics Approach. 
G3: Genes|Genomes|Genetics, 9(5), 1663–1674. https://doi.org/10.1534/g3.119.400167

The study identifies candidate genes and QTL associated with resistance to apple fruit ring rot (FRR) 
using a combination of whole-genome re-sequencing, bulk segregant analysis (BSA-Seq), and RNA-Seq 
meta-analysis.

## Data Description

- Organism: Apple (Malus domestica)
- Phenotypes: Resistance/susceptibility to Botryosphaeria dothidea isolates
- Samples: 1,773 F1 hybrids from 'Jonathan' × 'Golden Delicious', plus 60 Malus germplasm accessions
- Sequencing type: Illumina paired-end sequencing (HiSeq 2500 and HiSeq X Ten)
- Data size: ~1.4 billion reads total (parental + bulks)
- Publicly available at: NCBI SRA accession PRJNA392908
- Supplemental material: https://doi.org/10.25387/g3.7819688

## Project Structure

- `data/` – Processed data files
- `raw_data/` – Original raw sequencing data (optional download)
- `scripts/` – Analysis scripts
- `results/` – Output from analysis (plots, tables)
- `docs/` – Notes, figures, and supplementary info
- `logs/` – Logs from workflow runs
- `README.md` – This file


