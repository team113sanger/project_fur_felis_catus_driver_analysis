library(dndscv)
library(biomaRt)
options(timeout = max(600, getOption("timeout")))
library(tibble)
library(dplyr)
library(here)
library(readxl)
library(readr)
library(stringr)

path_genome_fasta <- "/lustre/scratch124/casm/team78pipelines/reference/Cat/Felis_catus_9.0/genome.fa"

prim_contigs <- c("A1","A2","A3","B1","B2","B3","B4",
                  "C1","C2","D1","D2","D3","D4","E1",
                  "E2","E3","F1","F2","X")

buildref(cdsfile= here("felis_catus_104.txt"), 
         genomefile=path_genome_fasta, 
         onlychrs = prim_contigs,
         useids = TRUE,
         outfile = "feline_transcript_104_dset.rda")



