#!/usr/bin/perl
#GetSeq: iPlant tool to retrieve sequences from NCBI given a list on IDs
#Author: Naim Matasci <nmatasci@iplantcollaborative.org>
#
# The contents of this file are subject to the terms listed in the LICENSE file you received with this code.
# Copyright (c) 2011, The Arizona Board of Regents on behalf of
# The University of Arizona
#
###############################################################################


use 5.008000;

use Bio::DB::SoapEUtilities;
use Bio::Seq;
use Bio::SeqIO;
use XML::Simple;
use Data::Dumper;
use Bio::Species;
use Getopt::Long;
use JSON;

our $VERSION = '0.0.1';

my ( $mode, $inputf, $outputf, $use_json, $version, $help, $ids );
my $EMAIL = 'nmatasci@iplantcollaborative.org';

$outputf = 'sequences.fa';

GetOptions(
	'i=s' => \$inputf,
	'o=s' => \$outputf,
	'm=s' => \$mode,
	'j'   => \$use_json,
	'v'   => \$version,
	'h'   => \$help
);

if ($version) {
	print "	
	GetSeq 0.0.1
	Copyright (c) 2011, The Arizona Board of Regents on behalf of The University of Arizona
	This software is released under a BSD 3-clause license. For details see http://iplantcollaborative.org/opensource\n\n";
	exit;
}

if ($help) {
	print "
	GetSeq 0.0.1
	Retrieves sequences from GenBank given a list of ids or accessions.
	Usage: ./GetSeq.pl  [options] -i input_file | id_list
	-h: print this help
	-v: prints the current version number
	-j: send output to STDOUT in JSON format 
	-m mode: Controls how much info is stored into the fasta header. One of 'full','species','short'. Default uses ids only
	-o output_file: Name of the output file. Defaults to \"sequences.fa\"
	-i input_file: name of a file containing a list of ids to retrieve
	id_list: list of ids to retrieve\n\n";
	exit;

}

#Obtains the IDs
if ($inputf) {
	$ids = read_list($inputf);    #From a file
}
else {
	$ids = \@ARGV;                #Directly from the command line
}
if ( !$ids->[0] ) {    #If no ids are specified, there is nothing to do
	warn "Error: Please provide a list of sequence identifiers.\n";
	exit;
}
print $ids;
for(@{$ids}){
	print "$_\n";	
}
#Pulls the sequence records from NCBI
my $records = retrive_seqs( $ids, $EMAIL );

#Transform the records into Bio::Seq objects
my $seqs = make_seq($records);

#Generate the output
if ($use_json) {
	as_json($seqs);
}
else {
	write_file( $seqs, $mode, $outputf );
}
exit;

#############
##FUNCTIONS##
#############

#Read the ID list
sub read_list {
	my $input = shift;
	open my $IDLIST, "<$inputf"
	  or die "Unable to load the sequence IDs list $inputf: $!\n";
	my @ids = <$IDLIST>;
	close $IDLIST;

	if ( @ids == 1 ) {
		@ids =
		  split( /\s+|\t+|,\s*/, $ids[0] )
		  ; #In case the ids are not on different lines but separated by whitespace or commas.
	}
	%seen = (); 
	foreach $item (@ids) { 
		$seen{$item}++; 
	}
	return [keys %seen];
}

#Retrieve the sequence records from NCBI
sub retrive_seqs {
	my $ids   = shift;
	my $email = shift;

	my $fac = Bio::DB::SoapEUtilities->new();

	my $xs_parser =
	  XML::Simple->new( ForceArray => ['TSeq'] )
	  ; #The fetch adaptor for fasta files doesn't export all the available info

	my $records = $xs_parser->XMLin(
		$fac->efetch(
			-db      => 'nucleotide',
			-id      => $ids,
			-rettype => 'fasta',
			-email   => $email,
			-tool    => 'GetSeq_dev'
		  )->run( -raw_xml => 1 )
	);

	if ( $records =~ m/Error/g ) {
		warn "Nothing found\n";
		print Dumper($records);
		exit;
	}
	return $records;
}

#write sequence in a multiple fasta file
sub write_file {
	my ( $seqs, $mode, $out ) = @_;
	my $outfh = Bio::SeqIO->new(
		'-file'   => ">$out",
		'-format' => 'fasta'
	);
	for ( @{$seqs} ) {
		$outfh->write_seq( _compile_id( $_, $mode ) )
		  or die "Cannot write to output file $out: $!\n";
	}
}

#make a json string
sub as_json {
	my $seqs = shift;

	my %rep;
	for ( @{$seqs} ) {
		$rep{ $_->id } = {
			id               => $_->id,
			sequence         => $_->seq,
			accession_number => $_->accession_number,
			species          => $_->species->binomial
			  . (
				$_->species->sub_species ? " " . $_->species->sub_species : ''
			  ),
			ncbi_taxid => $_->species->ncbi_taxid,
		};
	}
	print JSON::to_json( \%rep, { pretty => 1, allow_nonref => 1 } );
}

#Create new sequence
#A parser for the XML returned by NCBI
sub make_seq {
	my $records  = shift;
	my @returned =
	  @{ $records->{'SOAP-ENV:Body'}->{'eFetchResult'}->{'TSeqSet'}->{'TSeq'} };
	my @seqs;

	for (@returned) {
		my %seq_data = (
			-id               => $_->{'TSeq_gi'},
			-accession_number => $_->{'TSeq_accver'},
			-species          =>
			  _make_species( $_->{'TSeq_orgname'}, $_->{'TSeq_taxid'} )
			,    #the species is a Bio::Species object
			-seq => $_->{'TSeq_sequence'}
		);

		push @seqs, Bio::Seq->new(%seq_data);
	}
	return \@seqs;
}

#Generates the Bio::Species object
sub _make_species {
	my ( $name_string, $taxon_id ) = @_;

	#In case the species is unknown
	if ( !$name_string ) {
		my $speciesobj = Bio::Species->new();
		$speciesobj->classification( 'unknown', 'Unknown' );
		$speciesobj->ncbi_taxid($taxon_id);
		return $speciesobj;
	}

#Splits the name string into individual component: the species name, the genus name and the rest (subspecies, varieties, etc).
	my @bin_comp = split " ", $name_string;
	my $genus    = shift @bin_comp;
	my $species  = shift @bin_comp;
	my $subepi = join " ", @bin_comp;    #The rest is turned back into a string

	#Create the Species object
	my $speciesobj = Bio::Species->new();
	$speciesobj->classification( $species, $genus )
	  ; #The species name needs to be passed according to a classification: SPECIES GENUS FAMILY---->KINGDON
	if ($subepi) {
		$speciesobj->sub_species($subepi);
	}
	$speciesobj->ncbi_taxid($taxon_id);

	return $speciesobj;
}

#The default output of Bio::SeqIO for fasta file puts only the sequence id in the sequence header.
#This adds additional information
sub _compile_id {
	my ( $seq, $mode ) = @_;

	if ( !$mode ) {
		return $seq;    #In case id only is sufficient
	}

	#Extract the information from the Species object
	#	my$species_string=$seq->species;
	my $nid;
	my $subsp =
	    $seq->species->sub_species()
	  ? $seq->species->sub_species
	  : '';             #If a subspecies is present is added
	my $species = $seq->species->species;
	my $genus   = $seq->species->genus;
	$subsp =~ s/\s+/_/g; #Spaces in the subspecies string are transformed into _

	#Puts the header text together
	if ( $mode eq 'full' ) {
		$nid = join '_',
		  (
			$seq->id, $seq->accession_number, $genus, $species, $subsp,
			$seq->species->ncbi_taxid
		  );             #All available info
	}
	if ( $mode eq 'short' ) {

		$nid = join '_',
		  (
			$seq->id,
			substr( $genus, 0, 3 )
			  . uc( substr( $species, 0, 1 ) )
			  . substr( $species, 1, 2 ),
			substr( $subsp,       0, 3 )
		  )
		  ; #Id and taxon name abbreviated to 3 character for genus and species (and subspecies)
	}
	if ( $mode eq 'species' ) {
		$nid = join '_', ( $seq->id, $genus, $species, $subsp )
		  ;    #Only id and full length species info
	}

	$nid =~ s/_+/_/mg;     #Removes double underscores
	$nid =~ s/\./__/mg;    #Transforms the . character into a double underscore
	$nid =~ s/_$//mg;      #Removes trailing underscores

	$seq->id($nid);

	return $seq;
}

