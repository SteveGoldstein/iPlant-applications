#!/usr/bin/perl -w
#paml_run.pl v. 1.0.3
#Wrapper to run PAML 4.4e components on iPlant Collaborative's infrastructure
#Author: Naim Matasci <nmatasci@iplantcollaborative.org>
# 
# The contents of this file are subject to the terms listed in the LICENSE file you received with this code.
# Copyright (c) 2011, The Arizona Board of Regents on behalf of
# The University of Arizona
#
###############################################################################
use strict;

##Readonly only used in dev
#use Readonly;

#Readonly my $EXEC_PATH => '/home/nmatasci/paml44/bin';
#Readonly my $DATA_PATH => '/home/nmatasci/paml44/dat';
#Readonly my %WHITE_LIST => { codeml => 1 }; 
#Readonly my $CONTROL_KEY => 'ctrlfile';
#Readonly my $TOOL_KEY => 'program';
#Readonly my $COMMENT => '* Modified by iPlant';
#Readonly my $NFN => 'ipc-'.time.'-'.int(rand(10000)).'.ctl';


my $EXEC_PATH = '/usr/local2/nmatasci/paml44/bin';
my $DATA_PATH = '/usr/local2/nmatasci/paml44/dat';
my %WHITE_LIST = ( codeml => 1 ); 
my $CONTROL_KEY = 'ctrlfile';
my $TOOL_KEY = 'program';
my $COMMENT = '* Modified by iPlant';
my $NFN = 'ipc-'.time.'-'.int(rand(10000)).'.ctl';


main(@ARGV);

sub main{
my$arguments=get_arguments(@_);
my$ctl=delete($arguments->{$CONTROL_KEY});
my$command = $EXEC_PATH.'/'.delete($arguments->{$TOOL_KEY})." $NFN";

my$ctrl_in=read_ctrlfile($ctl);
my$new_ctrl=replace($arguments,$ctrl_in);
write_ctrl($NFN,$new_ctrl);
exec($command);
}


sub write_ctrl{
	my($name,$content)=@_;
	open my $OUTFH, ">$name" or die("Cannot write new control file: $!\n");
	print $OUTFH $content;
	close $OUTFH;
	return 1;	
}



sub replace {
	my($args,$ctrl)=@_;
	my@ctrl_lines=@{$ctrl};
	my$ctrl_text=join "", @ctrl_lines;	
	while( my($key,$val)  = each %{$args}){
		my$pattern=qr/(?:\s+|\*+|^)$key/;
		my@lines = grep /$pattern/, @ctrl_lines;
		
		if($key eq "seqfile" || $key eq "treefile"){
			$val=strip_path($val);
		}

		if(!@lines){
			$ctrl_text.="$key=$val\n";	
		}
		else{
			for (@lines){
				my$k=quotemeta($_);
				$ctrl_text =~ s/$k/$key = $val\t$COMMENT\n/;	
			}	
		}
		
	} 
	return fix_dat_path($ctrl_text);
}

sub strip_path{
	my$fullpath=shift;
	my@elements=split /\//, $fullpath;
	return pop @elements;
}

sub fix_dat_path{
	my$ctrl_text=shift;
	$ctrl_text =~ s/(?:\w+\/)*(\w+\.dat)/$DATA_PATH\/$1/g;
	return $ctrl_text;		
}

sub get_arguments {
	my@args=@_;
	my %arguments;
	for my $keyval (@args) {
		chomp $keyval;
		my ( $key, $val ) = split qr(=), $keyval;
		$arguments{$key} = $val;
	}
	if (!$arguments{$TOOL_KEY} || !$WHITE_LIST{$arguments{$TOOL_KEY}}){
		my$err_msg="Please specify a supported PAML program. Valid choices are: ".join (" ", keys %WHITE_LIST).".\n";
		die($err_msg);
	}
	elsif ( !$arguments{$CONTROL_KEY} ) {
		die("Please specify an input control file template\n");
	}
	else {
		return \%arguments;
	}
}

sub read_ctrlfile{
	my$path=shift;
	open my$CTRLF, "<$path" or die("Cannot load control file $path: $!\n");
	my@ctrl_entries=(<$CTRLF>);
	close $CTRLF;
	return \@ctrl_entries;
}


