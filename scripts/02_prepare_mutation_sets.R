library(readr)
library(dndscv)
library(fs)
library(purrr)

cohort_mafs <- fs::dir_ls("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/phased_mafs/outputs", glob = "**baitset_only.whole_cohort.keep.maf")
names(cohort_mafs) <- str_extract(cohort_mafs, pattern = "([0-9]+_[0-9]+)")

mutations <- read_tsv("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/phased_mafs/outputs/7040_3064.baitset_only.whole_cohort.keep.maf") |> 
filter(Hugo_Symbol == "TP53") |>
filter(!(Variant_Type %in% c("DNP", "TNP", "ONP"))) %>%
  select(Tumor_Sample_Barcode, Chromosome, POS_VCF, REF_VEP, ALT_VEP) %>%
  rename(sampleID = Tumor_Sample_Barcode, chr = Chromosome, pos = POS_VCF, ref = REF_VEP, mut = ALT_VEP) %>%
  mutate(chr = str_replace(chr, "chr", "")) %>%
  distinct()

run_dndscv <- function(cohort_file, suffix){

test <- read_tsv(cohort_file)

mutations <- test |>
  filter(!(Variant_Type %in% c("DNP", "TNP", "ONP"))) %>%
  select(Tumor_Sample_Barcode, Chromosome, POS_VCF, REF_VEP, ALT_VEP) %>%
  rename(sampleID = Tumor_Sample_Barcode, chr = Chromosome, pos = POS_VCF, ref = REF_VEP, mut = ALT_VEP) %>%
  mutate(chr = str_replace(chr, "chr", "")) %>%
  distinct()

dndsout <- dndscv(mutations, refdb = here("feline_transcript_dset.rda"), cv = NULL)
dndsout

sel_cv <- dndsout$sel_cv
print("Sig. genes")
print(head(sel_cv, 20), digits = 3)
tibble(sel_cv)  |>
select(gene_name, contains("q"))
write.table(sel_cv, file = paste0("results/dndscv_genes_", suffix,".tsv"), quote = F, sep = "\t", row.names = F)



print("Global dnds")
print(dndsout$globaldnds)

print("Annot muts")
head(dndsout$annotmuts)

print("Theta")
print(dndsout$nbreg$theta)
return(sel_cv)
}

sig_genes_per_cohort <- imap_dfr(cohort_mafs,~run_dndscv(.x,.y), .id = "cohort")

sig_genes_per_cohort |> tibble() |> filter(qallsubs_cv < 0.05)


knitr::kable(sig_genes_per_cohort |> tibble() |> 
filter(pallsubs_cv < 0.05) |> 
select(cohort, gene_name,qallsubs_cv) |> 
arrange(qallsubs_cv) |> 
print(n = 100))