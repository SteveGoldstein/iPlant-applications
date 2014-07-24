#!/usr/bin/perl
#readseq_wrapper.pl: wrapper for iubio's Readseq: http://iubio.bio.indiana.edu/soft/molbio/readseq/java/
#usage: readseq_wrapper.pl -f[ormat] <format_id> [-noclean -nodedup -l[imit] <name_length>] [-o[ut] <output_name>] <input_file> [<input_files>...]
#
#Author: Naim Matasci <nmatasci@iplantcollaborative.org>

use strict;
use POSIX qw(ceil);
use Digest::MD5 qw(md5_hex);
use Getopt::Long qw(GetOptions);

our $VERSION = 1.2;

my $MAX_LEN = 10;

my %SUFFIX = (
	1  => '.ig',
	2  => '.gb',
	3  => '.nbrf',
	4  => '.embl',
	5  => '.gcg',
	6  => '.strider',
	8  => '.fas',
	11 => '.phylip2',
	12 => '.phylip',
	13 => '.seq',
	14 => '.pir',
	15 => '.msf',
	17 => '.nexus',
	18 => '.pretty',
	19 => '.xml',
	22 => '.aln',
	23 => '.fff',
	24 => '.gff',
	25 => '.ace'
);

my $FORMATS = "
ID\tFormat
1	IG|Stanford
2	GenBank|gb
3	NBRF
4	EMBL|em
5	GCG
6	DNAStrider
8	Pearson|Fasta|fa
11	Phylip3.2 (sequential)
12	Phylip|Phylip4 (interleaved)
13	Plain|Raw
14	PIR|CODATA
15	MSF
17	PAUP|NEXUS
18	Pretty
19	XML
22	Clustal
23	FlatFeat|FFF
24	GFF
25	ACEDB
";

my $DIR = '.';

main(@ARGV);

sub main {
	my $inf;
	my $out;
	my $format;
	my $clean = 1;
	my $dedup = 1;
	my $limit = 0;

	GetOptions(
		"format=i" => \$format,
		"out=s"    => \$out,
		"clean!"   => \$clean,
		"dedup!"   => \$dedup,
		"limit:i"  => \$limit
	);
	print "Running on ", join "\t", @ARGV;
	if ( $limit == 1 ) {
		$limit = 10;
	}
	my $c   = 'cd "$( dirname "$0" )" && pwd';
	my $dir = `$c`;
	if ($dir) {
		chomp $dir;
		$DIR = $dir;
	}
	if ( !$format || !@ARGV || !$SUFFIX{$format} ) {
		die
"Please provide an input file and a valid output format ID. Accepted formats are $FORMATS\n";
	}
	elsif ( @ARGV == 1 ) {
		if ( !$out ) {
			$out = 'sequences';
		}
		my $text = run( $ARGV[0], $format, $clean, $dedup, $limit );
		write_out( $text, $out );
	}
	else {
		if ($out) {
			warn
"Option -o not valid with multiple files. Output files will be named after input name.";
		}
		for my $inf (@ARGV) {
			my $text = run( $inf, $format, $clean, $dedup, $limit );
			write_out( $text, $inf . $SUFFIX{$format} );
		}
	}
}

sub run {

	my ( $inf, $format, $clean, $dedup, $limit ) = @_;
	open my $INF, "<$inf" or die "Cannot load input file $inf:$!\n";
	my $text_file = join '', <$INF>;
	close $INF;

	my $names = get_names($inf);

	my ( $clean_text, $subst_table ) = names2nums( $text_file, $names );

	my $c =
	  "echo \"$clean_text\"|java -cp \"$DIR/readseq.jar\" run -f $format -p";
	my $new_form = `$c`;

	my $out_text =
	  nums2names( $new_form, $subst_table, $clean, $dedup, $limit );
	return $out_text;
}

sub write_out {
	my ( $text, $ofn ) = @_;
	open my $OFN, ">$ofn" or die "Cannot write output file $ofn:$!\n";
	print $OFN $text;
	close $OFN;
}

sub get_names {
	my $fn       = shift;
	my $c        = "java -cp \"$DIR/readseq.jar\" run -p -l \"$fn\"";
	my @id_names = `$c`;
	print join "", @id_names;
	my @names;
	for my $name (@id_names) {
		chomp $name;
		$name =~ s/^\s*\d+\)\s+//g;
		push @names, $name;
	}
	return \@names;
}

sub clean_names {
	my $dirty_names = shift;
	my @clean_names;
	for ( @{$dirty_names} ) {
		push @clean_names, _clean_name($_);
	}
	return \@clean_names;
}

sub _clean_name {
	my $dirty_name = shift;
	$dirty_name =~ s/[<> _=',;:\]\[\\\/\(\)]+/_/g;
	$dirty_name =~ s/_*(\W)_*/$1/g;
	$dirty_name =~ s/[\s\-_]+$//g;
	return $dirty_name;
}

sub dedup_names {
	my ( $names, $limit ) = @_;
	my @new_names;
	my %seen;
	for my $name ( @{$names} ) {
		if ( $limit && length($name) >= $limit ) {
			$name = substr $name, 0, $limit;
		}
		if ( $name =~ m/\#\d+$/ ) {
			$name .= "#o";
		}
		if ( !$seen{$name} ) {
			push @new_names, $name;
			if ( $limit && length($name) >= $limit - 2 ) {
				$name = substr $name, 0, $limit - 2;
			}
			$seen{$name} = 1;
		}

		else {
			$name =~ s/\#\d+$//;
			if ( $limit && length($name) >= $limit - 2 ) {
				$name = substr $name, 0, $limit - 2;
			}
			my $nname = $name . "#" . ( $seen{$name} + 1 );

			my $i = 1;
			while ( $seen{$nname} ) {
				$nname = $name . "#" . ( $seen{$name} + $i );
				$i++;
			}
			$seen{$name}++;
			push @new_names, $nname;

		}
	}
	return \@new_names;
}

sub scrub_names {
	my ( $text_file, $names, $clean, $dedup, $limit ) = @_;
	my %name_map;
	my @onames = @{$names};

	if ($clean) {
		$names = clean_names($names);

		#		print join "\n",@{$names},"\n";
	}
	if ($dedup) {
		$names = dedup_names( $names, $limit );

		#		print join "\n",@{$names},"\n";
	}
	for ( my $i = 0 ; $i < @onames ; $i++ ) {
		push @{ $name_map{ $onames[$i] } }, $names->[$i];

		#		print "$onames[$i] -> ",$names->[$i],"\n";
	}
	my $k = 0;
	for ( my $i = 0 ; $i < @onames ; $i++ ) {
		my $n_name = shift @{ $name_map{ $onames[$i] } };

		#		print "$num -> $name -> $n_name\n";
		$k += $text_file =~ s/\Q$onames[$i]\E/$n_name/;
	}
	$text_file =~ s/\s+$//ig;
}

sub names2nums {
	my ( $text_file, $names, $max_len ) = @_;
	my $n_names = scalar( @{$names} );
	my $post    = sig_n($n_names);
	if ( !$max_len ) {
		$max_len = $MAX_LEN;
	}
	my $pre = $max_len - ( 1 + $post );

	my $pattern = 0;

  OUTER:
	for ( my $i = $pre ; $i > 0 ; $i-- ) {
		for my $k ( '_', '-', '=', 'a' .. 'z', 'A' .. 'Z' ) {

			if ( $text_file =~ m/\d{$i}$k\d{1,$post}/ ) {
				next;
			}
			else {
				$pattern = [ $i, $k, $post ];
				last OUTER;
			}
		}
	}

	if ( !$pattern ) {
		$pattern = [$max_len];
	}

	my $subst_table = _make_nums( $names, $pattern );
	my $k = 0;
	while ( my ( $num, $name ) = each %{$subst_table} ) {
		my $w = $text_file =~ s/\s*\Q$name\E(\s+|\r?\n|\r)/$num$1/;
		if ( !$w ) {
			$w = $text_file =~ s/\s*\Q$name\E/$num/
			  ; #strict phylip allows the first 10 chars to be label and the 11th to be the start of the sequence
		}
		$k += $w;
	}
	warn "Substituted $k names for $n_names\n";
	return ( $text_file, $subst_table );
}

sub nums2names {
	my ( $text_file, $subst_table, $clean, $dedup, $strict ) = @_;
	my %name_map;
	my @onames = values %{$subst_table};
	my $names  = \@onames;
	if ($clean) {
		$names = clean_names($names);

		#		print join "\n",@{$names},"\n";
	}
	if ($dedup) {
		$names = dedup_names( $names, $strict );

		#		print join "\n",@{$names},"\n";
	}
	for ( my $i = 0 ; $i < @onames ; $i++ ) {
		push @{ $name_map{ $onames[$i] } }, $names->[$i];

		#		print "$onames[$i] -> ",$names->[$i],"\n";
	}
	my $k = 0;
	while ( my ( $num, $name ) = each %{$subst_table} ) {
		my $n_name = shift @{ $name_map{$name} };

		#		print "$num -> $name -> $n_name\n";
		$k += $text_file =~ s/\Q$num\E/$n_name/;
	}
	$text_file =~ s/\s+$//ig;
	return $text_file;
}

sub _make_nums {
	my ( $names, $pattern ) = @_;

	my $subst_table;
	if ( @{$pattern} > 1 ) {
		my $pre = sprintf( "%0" . $pattern->[0] . "d", 0 );
		my $pattern = $pre . $pattern->[1];
		for ( my $i = 0 ; $i < @{$names} ; $i++ ) {
			$subst_table->{ $pattern . $i } = $names->[$i];
		}
	}
	else {
		my @alphabet = ( 0 .. 9, "a" .. "f" );
		for my $name ( @{$names} ) {
			my $chk   = md5_hex($name);
			my $nname = $chk;
			if ( $pattern->[0] < 32 ) {
				$nname = substr $chk, 0, $pattern->[0];
			}
			while ( $subst_table->{$nname} ) {
				my $schar = $alphabet[ int( rand() * 16 ) ];
				substr $nname, int( rand() * length $nname ), 1, $schar;
			}
			$subst_table->{$nname} = $name;
		}
	}
	return $subst_table;
}

sub sig_n {
	my $x = shift;
	if ( $x <= 0 ) {
		die "Log($x) is undefined:$!\n";
	}
	return ceil( log($x) / log(10) );
}
