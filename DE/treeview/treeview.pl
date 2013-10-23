#!/usr/bin/env perl
#treeview: A phylogenetic tree importer for iPlant tree viewer
#Author: Naim Matasci <nmatasci@iplantcollaborative.org>
#
# The contents of this file are subject to the terms listed in the LICENSE file you received with this code.
# Copyright (c) 2011, The Arizona Board of Regents on behalf of
# The University of Arizona
#
###############################################################################
use strict;
use Bio::Phylo;
use Bio::Phylo::IO qw(parse unparse);
use LWP::UserAgent;
use URI::Escape;



my $treefile=$ARGV[0];
my $tree_name=$treefile;
$tree_name=~ s/\.\w+$//m;
my$TEMPLATE='
<html>
<head>
<meta http-equiv="refresh" content="0;url=$URL">
</head>
</html>';
my $Server='http://estraven.iplantcollaborative.org';

my$tree=load_tree($treefile);
my$url=get_url($tree_name,$tree);
write_file($tree_name,$url,$TEMPLATE);
 
 
#get tree
sub load_tree{
	my$file=shift;
	my$project;

	open my $TREEFILE, "<$file" or die "Cannot load input tree $file: $!\n";
	my$input=join "",<$TREEFILE>;
	close $TREEFILE;

#	print "Parsing\n";

	if($input=~ /<nex:nexml/){
		$project=parse(
			-file => $file,
			-format => 'nexml',
			-as_project => 1
		);

	} 
	elsif($input=~ /#nexus/i){
		$project=parse(
			-file => $file,
			-format=> 'nexus',
			-as_project =>1
			);
	}elsif(check_newick($input)){
#		print "is newick\n";
		return $input;
	}
	else{
		die("Cannot recognize the input file $file.\nAccepted formats are newick, nexus and nexML.\n");	
	}
	
#	print "deparsing\n";
	my$forest=$project->get_forests();
	my$tree=@{$forest}[0]->first();
	return $tree->to_newick();
} 


#Send tree to service
sub get_url{
	my($name,$tree)=@_;
	$tree=uri_escape($tree);

  my $ua = LWP::UserAgent->new;
#  $ua->agent("Naim's treeview wrapper");

 

  
  # Create a request
  my $req = HTTP::Request->new(POST => join($Server,'/parseTree?') );
  $req->content_type('application/x-www-form-urlencoded');
  $req->content("newickData=$tree&name=$name");

  # Pass request to the user agent and get a response back
  my $res = $ua->request($req);

	
  if ($res->is_success) {
  	
      return $res->content;
  }
  else {

      exit("Cannot obtain tree address: $res->status_line\n");
  }
}


#Create file
sub write_file{
	my($name,$url,$template)=@_;
	$template=~ s/\$URL/$url/;
	$name=$name.'.html';
	open my $HTMLOUT, ">$name" or die("Cannot create output file $name: $!\n");
	print $HTMLOUT $template;
	close $HTMLOUT;
}

#Check Newick format
sub check_newick{
	my$string=shift;
	if($string!~/;$/m){ #Doesn't end with a semicolon.
		return 0;
	}
	
	$string=~ s/(\w)\s(\w)/$1_$2/g; #Whitespaces between words are transformed into underscores
	$string=~ s/(?:\n+|\s+)//g; #Removes newlines and whitespaces

#	Might not be necessary given the next test
	my@lpar=$string =~ m/\(/g;
	my@rpar=$string =~ m/\)/g;
	
	if(scalar(@lpar) == 0 || scalar(@rpar) == 0 || scalar(@lpar) != scalar(@rpar)){ #There are no parentheses or they are unbalanced. 
		return 0;	
	}


#What follows is a fast heuristic to decide wether a file is a valid newich string
#For large trees, parsing the string is time consuming and doesn't necessarily result in a valid string.
#This test identifies valid tokens that reresent a valid name, optional branch length and its bordering parentheses and commas
#These are then used to create a new newick string which is compared with the original one.

		my@species=$string=~ m/\(*(?:(?:\w+|\.|_|-|\'|\/|\+|\\|\")+(?::\d*\.?\d+)?)(?:,|\))+;?/g; #Splits the string into valid tokens
		my$species=join "",@species;
		if($species && length($species) eq length($string)){
			return 1;
		}
		return 0;
}
