module load perl-5.38.0
module load ensembl-vep/104


ensembl-vep/104-exec \
perl get_canonical_transcripts.pl > feline_104_canonical_transcripts.tsv