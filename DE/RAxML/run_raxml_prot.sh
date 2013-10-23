#!/bin/sh

# raxml_prot.sh
# 
#
# Created by Naim Matasci on 7/5/11.
# Copyright 2011 The iPlant Collaborative. All rights reserved.

#Collapse
#!/bin/bash
 
while getopts ":s:m:n:" opt; do
  case $opt in
    s)
		seq=$OPTARG
#      echo "-s was triggered, Parameter: $OPTARG" >&2
      ;;
	m)
#      echo "-m was triggered, Parameter: $OPTARG" >&2
		mod=$OPTARG
      ;;
    n)
		nou=$OPTARG
#      echo "-n was triggered, Parameter: $OPTARG" >&2
      ;;

	
  esac
  
done
shift 6
for arg in "$@"
do
mod=$mod$arg
done

echo "$mod"
#/usr/local2/RAxML-7.2.8-ALPHA/raxmlHPC-SSE3 -s $seq -n $nou -m $mod