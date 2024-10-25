#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";
 
$reg->load_registry_from_db(
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
   -dbname  => 'felis_catus_core_104_9'
);
 
#open(OUTFILE, ">human_canonical_transcripts_v103.fa");
 
my $gene_adaptor = $reg->get_adaptor('cat', 'core', 'gene');
my @genes = @{$gene_adaptor->fetch_all};
 
while(my $gene = shift @genes) {
#	print STDERR "$gene\n";
	my $symbol = $gene->external_name ? $gene->external_name : "-";
	print $symbol . "\t" . $gene->stable_id . "\t" . $gene->canonical_transcript->stable_id . "\t" . $gene->canonical_transcript->stable_id . "." . $gene->canonical_transcript->version . "\n";
}

