
# Fetch canonical transcripts for feline genome assembly GCF_000181335.3
docker run -v $PWD:$PWD \
--workdir $PWD \
ensemblorg/ensembl-vep:release_104.1 \
perl scripts/get_canonical_transcripts.pl > results/inputs/feline_104_canonical_transcripts.tsv

docker run -v $PWD:$PWD \
--workdir $PWD \
quay.io/staphb/ncbi-datasets:16.41.0 \
datasets download genome accession GCF_000181335.3 \
--include genome --filename GCF_000181335.3.fna.gz

unzip GCF_000181335.3.fna.gz -d felCat9


