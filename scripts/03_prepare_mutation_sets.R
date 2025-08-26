################################################################################
# Setup
################################################################################
library(readr)
library(dndscv)
library(fs)
library(tibble)
library(dplyr)
library(here)
library(readxl)
library(stringr)
library(purrr)
library(logger)



################################################################################
# Dndscv helper
################################################################################

run_dndscv <- function(cohort_file, targeted = TRUE, target_genes, suffix, tagged) {
  test <- read_tsv(cohort_file)

  mutations <- test |>
    filter(
      !(Variant_Type %in% c("INS", "DEL")) |
        (VAF_tum >= 0.1 & Variant_Type %in% c("INS", "DEL")) # Filter out likely indel artefacts
    )
    
    if (tagged){
    mutations <- mutations |> 
    filter(`99_Lives` == "-")
    }
    mutations <- mutations |>
    select(Tumor_Sample_Barcode, Chromosome, POS_VCF, REF_VEP, ALT_VEP) |>
    rename(sampleID = Tumor_Sample_Barcode, chr = Chromosome, pos = POS_VCF, ref = REF_VEP, mut = ALT_VEP) |>
    mutate(chr = str_replace(chr, "chr", "")) |>
    distinct()
  if (targeted){
    logger::log_info("Running dndscv on targeted mutations")
  dndsout <- dndscv(mutations,
    gene_list = target_genes,
    refdb = here("results/inputs/feline_transcript_104_canon_dset.rda"), 
    max_muts_per_gene_per_sample = 2,
    cv = NULL
  )}
  else {
    logger::log_info("Running dndscv on WES/WGS mutations")
    dndsout <- dndscv(mutations,
    refdb = here("results/inputs/feline_transcript_104_canon_dset.rda"), 
    max_muts_per_gene_per_sample = 2,
    cv = NULL
  )
  }
  dndsout

  sel_cv <- dndsout$sel_cv
  print("Sig. genes")
  print(head(sel_cv, 20), digits = 3)
  tibble(sel_cv) |>
    select(gene_name, contains("q"))
  write.table(sel_cv, file = paste0("results/dndscv/dndscv_genes_", suffix, ".tsv"), quote = F, sep = "\t", row.names = F)



  print("Global dnds")
  print(dndsout$globaldnds)

  print("Annot muts")
  head(dndsout$annotmuts)

  print("Theta")
  print(dndsout$nbreg$theta)

  sel_out <- sel_cv |> mutate(theta = dndsout$nbreg$theta)
  return(sel_out)
}

################################################################################
# Read in cohort MAFs
################################################################################
cohort_mafs <- fs::dir_ls("inputs",
  recurse = TRUE, regexp = "keep_vaf_size_filt_matched.*[0-9]+.*.maf$"
)
names(cohort_mafs) <- str_extract(cohort_mafs, pattern = "([0-9]+_[0-9]+)")


################################################################################
# Prep targeted gene list
################################################################################
cat_baitset <- read_xlsx(here("metadata/FUR_BAITSET_GENES.xlsx"), sheet = 3)
gene_list <- cat_baitset |>
  filter(Feline_Gene_Symbol != "-") |>
  mutate(id = paste0(Feline_Ensembl_ID, ":", Feline_Gene_Symbol)) |>
  pull(id)


missing_genes <- c(
  "ENSFCAG00000011854:ALK", "ENSFCAG00000006946:CARS",
  "ENSFCAG00000005106:FGFR1OP", "ENSFCAG00000034898:GNAS",
  "ENSFCAG00000030424:H3F3B", "ENSFCAG00000006079:HIST1H1B",
  "ENSFCAG00000003805:HIST1H1C", "ENSFCAG00000011368:YWHAE",
  "ENSFCAG00000006768:HIST1H1D", "ENSFCAG00000051764:HIST1H1E",
  "ENSFCAG00000053015:HIST1H2AC", "ENSFCAG00000052405:HIST1H2BD",
  "ENSFCAG00000047792:HIST2H3C", "ENSFCAG00000044858:HOXC13",
  "ENSFCAG00000009764:HSP90AB1", "ENSFCAG00000005238:ICK",
  "ENSFCAG00000015604:MAGED1", "ENSFCAG00000024092:PIK3R3",
  "ENSFCAG00000044860:RIT1"
)




filtered_gene_list <- setdiff(gene_list, missing_genes)


################################################################################
# Check genes for exome data
################################################################################

all_genes <- read_tsv(here("metadata/felis_catus_104_biomart.txt"))
observed_genes <- read_tsv("inputs/keep_vaf_size_filt_matched_7834_feline_mammary_wes.filtered.tissue.maf") |>
                    filter(Hugo_Symbol != "-") |> 
                    pull(Hugo_Symbol)



observed_exome <- all_genes |> 
filter(`Gene name` %in% observed_genes) |>
transmute(Gene =  paste0(`Gene stable ID`, ":", `Gene name`)) |> 
distinct() |> 
pull()




################################################################################
# Run dndsvc on targeted cohorts
################################################################################

sig_genes_per_cohort <- imap_dfr(cohort_mafs[1:15], 
                        ~ run_dndscv(.x, targeted = TRUE, 
                                     filtered_gene_list, .y, TRUE), 
                                     .id = "cohort")

sig_genes_per_cohort |>
  tibble() |>
  filter(qglobal_cv < 0.05) |>
  group_by(cohort) |>
  select(cohort, gene_name, qglobal_cv, theta) |>
  group_split() |>
  map(knitr::kable)



knitr::kable(sig_genes_per_cohort |> tibble() |>
  filter(qglobal_cv < 0.05) |>
  select(cohort, gene_name, qglobal_cv, theta) |>
  arrange(qglobal_cv) |>
  print(n = 100))

################################################################################
# Run dndsvc on WES cohorts
################################################################################
names(cohort_mafs)[16] <- "7834_3518_organoids"
names(cohort_mafs)[17] <- "7834_3518_tumours"

sig_genes_wes <- imap_dfr(cohort_mafs[14:15],
~ run_dndscv(.x, targeted = FALSE, 
             observed_exome, suffix = .y , FALSE), 
             .id = "cohort")

sig_genes_wes |>
  tibble() |>
  filter(qglobal_cv < 0.05) |>
  group_by(cohort) |>
  select(cohort, gene_name, qglobal_cv, theta) 


knitr::kable(sig_genes_wes |> tibble() |>
  filter(qglobal_cv < 0.05) |>
  select(cohort, gene_name, qglobal_cv, theta) |>
  arrange(cohort,qglobal_cv) |>
  print(n = 100))
