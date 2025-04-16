# FUR Cat Driver analysis 
## Summary
Significantly mutated genes within each FUR Felis Catus cohort were identified using dNdScv. 

To perform dNdScv analysis:
- 1) The Ensembl v104 canonical transcript for each gene in the FUR cat baitset was obtained (`scripts/01_fetch104_canonical_transcripts.sh`)
- 2) A dNdScv reference database was created using the FelCat9 genome and the FelCat 104 canonical transcripts (`scripts/01_fetch104_canonical_transcripts.sh`)
- 3) Cohort MAF files generated in upstream parts of the study were filtered and reformatted to input into dNdScv and genes were assessed for significance (`scripts/03_prepare_mutation_sets.R`)

Each cohort's "Keep" variants - those passing the QC filtering steps for the FUR Cat project – were subject to two filters.
First, indels with VAF < 0.1 were removed from the cohort MAF files to reduce the impact of possible FFPE artefacts. 
Then, any variants detected in the "99 Lives" set (likely germline variants) were removed from the analysis.

The max genes per sample for dNdScv was set to 2 in order to reduce the contribution of multi-hit genes in a single sample on the dNdScv analysis.


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

## Dependencies

This analysis was conducted within a development container with docker and R (4.4.0) installed (`.devcontainer/devcontainer.json`). All scripts can therefore be run in a Github codespace or Vscode session with docker available. All R dependencies for this project are detailed within the project renv.lock file and can be installed by running `renv::restore()` in a R terminal.



## Contact 
- jb63@sanger.ac.uk