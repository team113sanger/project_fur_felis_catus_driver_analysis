
docker run -v $PWD:/$PWD \
ensemblorg/ensembl-vep:release_104.1 \
perl scripts/get_canonical_transcripts.pl > results/inputs/feline_104_canonical_transcripts.tsv