# FUR Cat Driver analysis 
## Summary
Significantly mutated genes within each FUR-CAT cohort were identified using dNdScv. 

To perform dNdScv analysis:
- 1) The Ensembl v104 canonical transcript for each gene in the FUR cat baitset was obtained (`scripts/01_fetch104_canonical_transcripts.sh`)
- 2) A dNdScv reference database was created using the FelCat9 genome and the FelCat 104 canonical transcripts (`scripts/01_fetch104_canonical_transcripts.sh`)
- 3) The cohort MAF files were filtered and reformatted to input into dNdScv (`scripts/03_prepare_mutation_sets.R`)

Each cohort's "Keep" variants - those passing the QC filtering steps for the FUR Cat project – were subject to two filters.
First, indels with VAF < 0.1 were removed from the cohort MAF files to reduce the impact of possible FFPE artefacts. 
Then, any variants detected in the "99 Lives" set (likely germline variants) were removed from the analysis.


## Project organistion
```
.
├── metadata
│   ├── felis_catus_104_biomart.txt
│   └── FUR_BAITSET_GENES.xlsx
├── README.md
├── renv
│   ├── activate.R
│   ├── library
│   ├── settings.json
│   └── staging
├── renv.lock
├── results
│   ├── dndscv
│   └── inputs
└── scripts
    ├── 01_fetch104_canonical_transcripts.sh
    ├── 02_generate_refcds_object.R
    ├── 03_prepare_mutation_sets.R
    └── get_canonical_transcripts.pl
```
