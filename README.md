# FUR Cat Driver analysis 

Significantly mutated genes within each FUR cohort were identified using dNdscv. To do so, the studies "Keep" variants - those passing QC filtering steps in the FUR project were subject to two filtering steps.First indels with VAF < 0.1 were removed from the cohort MAF files to reduce the impact of possible  FFPE artefacts in dNdSCv analysis. Then, any variants with matches to the 99 Lives set were removed to remove likely Germline variants from the cohort MAFs. 

In brief, this process involved:
- 1) Obtaining the canonical transcript for each of the genes included within the FUR cat baiset.
- 2) Generating a reference db from the FelCat9 genome and the canonical transcripts.
- 3) Reformating + filtering the input MAF files and then running the variants within them through dNdScv
```
.
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
    ├── 01_generate_refcds_object.R
    ├── 02_prepare_mutation_sets.R
    └── get_canonical_transcripts.pl
```
