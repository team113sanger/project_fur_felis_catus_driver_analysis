library(readr)
library(dndscv)
library(fs)
library(tibble)
library(dplyr)
library(here)
library(readxl)
library(stringr)
library(purrr)

cohort_mafs <- fs::dir_ls("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/phased_mafs/outputs", glob = "**baitset_only.whole_cohort.keep.maf")
names(cohort_mafs) <- str_extract(cohort_mafs, pattern = "([0-9]+_[0-9]+)")


cat_baitset <- read_xlsx(here("FUR_BAITSET_GENES.xlsx"), sheet = 3)
gene_list <- cat_baitset |>
  filter(Feline_Gene_Symbol != "-") |>
  mutate(id = paste0(Feline_Ensembl_ID, ":", Feline_Gene_Symbol)) |>
  pull(id)

cat_baitset
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


run_dndscv <- function(cohort_file, target_genes, suffix) {
  test <- read_tsv(cohort_file)

  mutations <- test |>
    filter(!(Variant_Type %in% c("DNP", "TNP", "ONP"))) %>%
    select(Tumor_Sample_Barcode, Chromosome, POS_VCF, REF_VEP, ALT_VEP) %>%
    rename(sampleID = Tumor_Sample_Barcode, chr = Chromosome, pos = POS_VCF, ref = REF_VEP, mut = ALT_VEP) %>%
    mutate(chr = str_replace(chr, "chr", "")) %>%
    distinct()

  dndsout <- dndscv(mutations,
    gene_list = target_genes,
    refdb = here("feline_transcript_104_dset.rda"), cv = NULL
  )
  dndsout

  sel_cv <- dndsout$sel_cv
  print("Sig. genes")
  print(head(sel_cv, 20), digits = 3)
  tibble(sel_cv) |>
    select(gene_name, contains("q"))
  write.table(sel_cv, file = paste0("results/dndscv_genes_", suffix, ".tsv"), quote = F, sep = "\t", row.names = F)



  print("Global dnds")
  print(dndsout$globaldnds)

  print("Annot muts")
  head(dndsout$annotmuts)

  print("Theta")
  print(dndsout$nbreg$theta)
  return(sel_cv)
}


filtered_gene_list <- setdiff(gene_list, missing_genes)
sig_genes_per_cohort <- imap_dfr(cohort_mafs, ~ run_dndscv(.x, filtered_gene_list, .y), .id = "cohort")

sig_genes_per_cohort |>
  tibble() |>
  filter(qallsubs_cv < 0.05) |>
  group_by(cohort) |>
  group_split()


knitr::kable(sig_genes_per_cohort |> tibble() |>
  filter(qallsubs_cv < 0.05) |>
  select(cohort, gene_name, qallsubs_cv) |>
  arrange(qallsubs_cv) |>
  print(n = 100))
