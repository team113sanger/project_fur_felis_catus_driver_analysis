library(dndscv)
library(biomaRt)
options(timeout = max(600, getOption("timeout")))
library(tibble)
library(dplyr)
library(here)
library(readxl)
library(readr)
library(stringr)

ensembl <- useMart("ensembl")

# List available datasets to find your species of interest
datasets <- listDatasets(ensembl)
datasets |> 
dplyr::filter(grepl(version, pattern = "Felis_catus_9.0"))

ensembl <- useMart("ensembl", dataset = "fcatus_gene_ensembl")

# Define the 10 attributes needed by buildref in the correct order
attributes <- c(
    "ensembl_gene_id",       # Gene stable ID
    "external_gene_name",    # Gene name
    "ensembl_peptide_id",    # Protein stable ID
    "chromosome_name",       # Chromosome/scaffold name
    "exon_chrom_start",      # Exon chromosome start
    "exon_chrom_end",        # Exon chromosome end
    "cds_start",             # Genomic coding start
    "cds_end",               # Genomic coding end
    "cds_length",            # Length of CDS 
    "strand"                 # Strand
)

# Retrieve the data from Ensembl BioMart
transcript_data <- getBM(attributes = attributes, mart = ensembl)



cat_baitset <- read_xlsx(here("FUR_BAITSET_GENES.xlsx"),sheet = 3)
head(transcript_data)


# These genes are causing problems with RefCDS construction.. Not clear why 


feline_transcripts <- transcript_data |> 
filter(ensembl_gene_id %in% cat_baitset[["Feline_Ensembl_ID"]]) |>
filter(!ensembl_gene_id %in% c("ENSFCAG00000001254","ENSFCAG00000007992", "ENSFCAG00000011205"))

write_tsv(feline_transcripts,"felis_catus_baitset_transcripts.tsv")

path_genome_fasta <- "/lustre/scratch124/casm/team78pipelines/reference/Cat/Felis_catus_9.0/genome.fa"


prim_contigs <- c("A1","A2","A3","B1","B2","B3","B4",
                  "C1","C2","D1","D2","D3","D4","E1",
                  "E2","E3","F1","F2","X")

buildref(cdsfile= here("felis_catus_104.txt"), 
         genomefile=path_genome_fasta, 
         onlychrs = prim_contigs,
         useids = TRUE,
         outfile = "feline_transcript_104_dset.rda")



