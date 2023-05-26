#!/bin/bash
# args:
# $1 - the fasta file
# $2 - the header of sequence to extract

# multi to single line 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $1 | grep -A 1 $2 