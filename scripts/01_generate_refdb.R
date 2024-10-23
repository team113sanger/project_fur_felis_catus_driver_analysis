library(dndscv)
library(biomaRt)
options(timeout = max(600, getOption("timeout")))
library(tibble)
library(dplyr)
library(here)
library(readxl)
library(readr)

ensembl <- useMart("ensembl")

# List available datasets to find your species of interest
datasets <- listDatasets(ensembl)
datasets |> dplyr::filter(grepl(version, pattern = "Felis_catus_9.0"))

# Replace 'hsapiens_gene_ensembl' with the dataset corresponding to your species
ensembl <- useMart("ensembl", dataset = "fcatus_gene_ensembl")

# Define the 10 attributes needed by buildref in the correct order
attributes <- c(
    "ensembl_gene_id",       # Gene stable ID
    'external_gene_name',    # Gene name
    'ensembl_transcript_id', # Transcript stable ID
    'rank',             # Exon rank in transcript
    'exon_chrom_start',      # Exon chromosome start
    'exon_chrom_end',        # Exon chromosome end
    'cds_start',             # Genomic coding start
    'cds_end',               # Genomic coding end
    'strand',                # Strand
    'chromosome_name'        # Chromosome/scaffold name
)

# Retrieve the data from Ensembl BioMart
transcript_data <- getBM(attributes = attributes, mart = ensembl)

# # Check for missing gene names and replace them with gene stable IDs if necessary
# transcript_data$external_gene_name[transcript_data$external_gene_name == ""] = 
#     transcript_data$ensembl_gene_id[transcript_data$external_gene_name == ""]



cat_baitset <- read_xlsx(here("FUR_BAITSET_GENES.xlsx"),sheet = 3)
head(transcript_data)

feline_transcripts <-transcript_data |> 
filter(ensembl_gene_id %in% cat_baitset[["Feline_Ensembl_ID"]])
write_tsv(feline_transcripts,"felis_catus_baitset_transcripts.tsv")

path_genome_fasta <- "/lustre/scratch124/casm/team78pipelines/reference/Cat/Felis_catus_9.0/genome.fa"


buildref(cdsfile=here("felis_catus_baitset_transcripts.tsv"), 
         genomefile=path_genome_fasta, 
         outfile = "feline_transcripts.rda")

# Export the data as a tab-delimited file
write.table(transcript_data, file = "felis_catus_9.0_transcripts.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)