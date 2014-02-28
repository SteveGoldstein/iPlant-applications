#!/usr/local2/perl-5.16.2/bin/perl
#blast_ortho.pl: Blasts nucleotide or protein sequences against a database of reference Arabidopsis sequences and retrieves the matching files from iRODS
#Author: Naim Matasci <nmatasci@iplantcollaborative.org>

use strict;
use version;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile tmpnam tempdir);
use Config::Simple;
use FindBin qw($Bin);

#use HTTP::Tiny;
#use Data::Dumper;
#use JSON qw(from_json to_json);

our $VERSION = qv(0.1.4);

main();

sub main {
	my $params = init();
	my $accs = blast( $params->{'in'}, $params->{'type'}, $params );
	if ( !$accs ) {
		print "No matches\n";
		exit 0;
	}
	my $acc2hmmid = retrieve_hmmid( $params->{'hmmid_path'} );

	for ( @{$accs} ) {
		my ( $oid, $acc ) = split /\t/, $_;
		if ( !$acc2hmmid->{$acc} ) {
			next;
		}
		my $files = query_ids( $acc2hmmid->{$acc}, $params );
		if ( !$files ) {
			next;
		}
		mkdir $oid;
		get_files( $files, $oid, $params );
	}

}

sub init {
	my $sdir = $Bin;

	my $cfg    = new Config::Simple("$sdir/blast_ortho.cfg");
	my %params = $cfg->vars();

	my ( $file, $type, $query, @outf );
	GetOptions(
		"in=s"          => \$file,
		"type=s"        => \$type,
		"query=s"       => \$query,
		"format=s{1,5}" => \@outf
	);

	if ( !$file ) {
		if ( !$query ) {
			die "Please provide either an input file or an query string";
		} 
		#otherwise write input file as temp
		my ( $TF, $filename ) = tempfile( UNLINK => 0 );    #
		print $TF $query;
		$file = $filename;
	}
	
	if ( !$type || lc($type) eq 'auto' ) {
		$type = guess_type($file);
	}
	
	$params{'blastdb_path'} = "$sdir/$params{'blastdb_path'}";
	$params{'hmmid_path'}   = "$sdir/$params{'hmmid_path'}";
	$params{'in'}           = $file;
	$params{'type'}         = $type;
	$params{'format'}       = \@outf;

	return \%params;
}

sub get_files {
	my ( $file_loc, $oid, $params ) = @_;
	for my $format ( @{$params->{'format'}} ) {
		if ( !$file_loc->{$format} ) {
			warn "Unable to retrieve files for format $format.\n";
			next;
		}
		my @tn = split /\//, $file_loc->{$format};
		my $fn = pop @tn;
		mkdir "$oid/$format";
		my$exec=$params->{'icommands'}.'/iget';
		my$vare= 'export irodsEnvFile='.$params->{'irods_env'};
		my $cmd = "$vare; $exec -f " . $file_loc->{$format} . " $oid/$format/$fn";
		my $k   = `$cmd`;
		if ( $k =~ /ERROR/ ) {
			warn "Unable to retrieve files for "
			  . $file_loc->{$format} . "\n$k";
		}
	}
}

sub query_ids {
	my ( $hmmid, $params ) = @_;
	my$exec=$params->{'icommands'}.'/iquest';
	my$vare= 'export irodsEnvFile='.$params->{'irods_env'};
	my $IQUERY_L =
	   $exec.' --no-page "select COLL_NAME, DATA_NAME where DATA_NAME like \'';
	my $IQUERY_R =
	  '.%\' AND COLL_NAME like \'' . $params->{'irods_path'} . '%\'"';
	my %file_loc;
	$file_loc{'id'} = $hmmid;
	my $k = "$vare; ".$IQUERY_L . $hmmid . $IQUERY_R;
	my $r = `$k`;
	if ( !$r || $r =~ /CAT_NO_ROWS_FOUND/ ) {
		warn "Files for HMMID $hmmid could not be found:\n$r\n";
		return 0;
	}
	my @res = split /^-+$/m, $r;
	for (@res) {
		if ( $_ =~ /COLL_NAME = (.*)\nDATA_NAME = (.*)/ ) {
			my $file = "$1/$2";
			if ( $file =~ /alignments\/FAA/ ) {
				$file_loc{'aa_aln'} = $file;
			}
			elsif ( $file =~ /alignments\/FNA/ ) {
				$file_loc{'nt_aln'} = $file;
			}
			elsif ( $file =~ /gene_trees\/FAA/ ) {
				$file_loc{'aa_tre'} = $file;
			}
			elsif ( $file =~ /12_codon/ ) {
				$file_loc{'12_tre'} = $file;
			}
			elsif ( $file =~ /all_codon/ ) {
				$file_loc{'nt_tre'} = $file;
			}

		}
	}
	return \%file_loc;
}

sub retrieve_hmmid {
	my $hmmid_f = shift;
	my %acc2hmmid;
	open my $HMMID, "<$hmmid_f" or die "Cannot open $hmmid_f: $!\n";
	while (<$HMMID>) {
		chomp;
		my ( $acc, $id ) = split /\t/, $_;
		$acc2hmmid{$acc} = $id
		  ; #NOTE: there are 15 case where a gene id belongs to multiple hmms. THESE ARE IGNORED
	}
	return \%acc2hmmid;

}

#sub query_db {
#	my $acc = shift;
#
#	my $call     = "http://$HOST/$EP\?accession=$acc";
#	my $response = HTTP::Tiny->new->get($call);
#	if ( !$response->{'success'} ) {
#		die "Unable to continue - orthoDB Error: ", $response->{'status'}, " ",
#		  $response->{'reason'}, "\n";
#	}
#	my $result = from_json( $response->{'content'} );
#	my $hmmid  = @{ $result->{'results'} }[0]->{'hmmId'};
#	return $hmmid;
#}

sub blast {
	my ( $file, $type, $params ) = @_;
	my $blastdb  = $params->{'blastdb_path'};
	my $nthreads = $params->{'nthreads'};
	my $blast_exe = $params->{'blast_exe'};
	my $exec  = "$blast_exe/blastp";
	if ( $params->{'type'} eq 'nuc' ) {
		$exec = "$blast_exe/blastx";
	}
	my $cmd = eval {
"$exec -db $blastdb -query $file -num_threads $nthreads -outfmt \"6 qseqid sseqid \" -max_target_seqs 1";
	};
	my $blout = `$cmd`;    #System call to blast. The result is stored in $blout
	if ( !$blout ) {
		return 0;
	}
	my @res = split /\n/, $blout;
	my @acc_list;
	for (@res) {
		my ( $qname, $def ) = split /\t/, $_;

		my $acc = ( split /\|/, $def )[2];
		chomp $acc;
		$acc =~ s/\.\d+$//;
		if ( !$qname || $qname eq '' || $qname eq 'unnamed' ) {
			$qname = $acc;
		}
		push @acc_list, "$qname\t$acc";
	}
	return \@acc_list;
}

sub guess_type {
	my $file = shift;
	my $type = 'nuc';
	open my $IF, "<$file" or die "Cannot load input query $file: $!\n";
	while (<$IF>) {
		if ( $_ =~ m/[efilpq]/i ) {
			$type = 'prot';
			close $IF;
			return $type;
		}
	}
	close $IF;
	return $type;
}
