library(dndscv)
library(biomaRt)
options(timeout = max(600, getOption("timeout")))
library(tibble)
library(dplyr)
library(here)
library(readxl)
library(readr)
library(stringr)

# List available datasets to find your species of interest


path_genome_fasta <- here("felCat9/ncbi_dataset/data/GCF_000181335.3/GCF_000181335.3_Felis_catus_9.0_genomic.fna")

prim_contigs <- c(
  "A1", "A2", "A3", "B1", "B2", "B3", "B4",
  "C1", "C2", "D1", "D2", "D3", "D4", "E1",
  "E2", "E3", "F1", "F2", "X"
)

all_genes <- read_tsv(here("metadata/felis_catus_104_ biomart.txt"))
canonical_transcripts <- read_tsv(
  here("feline_104_canonical_transcripts.tsv"),
  col_names = c("Symbol", "Gene", "Transcript", "Transcript Stable ID")
)

canonical_genes <- all_genes |>
  filter(`Transcript stable ID version` %in% canonical_transcripts[["Transcript Stable ID"]]) |>
  select(-`Transcript stable ID version`)


write_tsv(canonical_genes, here("results/inputs/felis_catus_104_canonical_refcds_inputs.txt"))


buildref(
  cdsfile = here("results/inputs/felis_catus_104_canonical_refcds_inputs.txt"),
  genomefile = path_genome_fasta,
  onlychrs = prim_contigs,
  useids = TRUE,
  outfile = here("results/inputs/feline_transcript_104_canon_dset.rda")
)
