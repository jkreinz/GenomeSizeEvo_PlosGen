#!/bin/bash
set -e
set -u
set -o pipefail

#mkdir outinner and ltrdigest before running

# name the file stem based on suffixator index
SPECIES="tuberculatus"
GENOME=$1
DIR=$2

###########################################################################
## Run suffixerator to make a suffix array of the genome for genometools ##
###########################################################################

gt suffixerator -db $GENOME -indexname $SPECIES -tis -suf -lcp -des -ssp -sds -dna

#####################
## Run LTR harvest ##
#####################

#all defaults except maxlenltr, mindistltr, maxdistltr
mkdir -p ${DIR}/outinner

#this command was run1 which has TSD arguments
#gt ltrharvest -index $SPECIES  -gff3 ${DIR}/${SPECIES}.ltrharvest.gff3 -motif tgca -minlenltr 100 -maxlenltr 7000 -mindistltr 200 -maxdistltr 20000 -similar 85 -motifmis 1 -mintsd 5 -xdrop 5 -overlaps best -longoutput -outinner ${DIR}/outinner/${SPECIES}.ltrharvest.outinner.fa -out ${DIR}/${SPECIES}.ltrharvest.fa > ${DIR}/${SPECIES}.ltrharvest.out

#run2
gt ltrharvest -index $SPECIES  -gff3 ${DIR}/${SPECIES}.ltrharvest.gff3 -minlenltr 100 -maxlenltr 7000 -mindistltr 200 -maxdistltr 20000 -similar 85 -xdrop 5 -overlaps best -outinner ${DIR}/outinner/${SPECIES}.ltrharvest.outinner.fa -out ${DIR}/${SPECIES}.ltrharvest.fa > ${DIR}/${SPECIES}.ltrharvest.out

gt gff3 -sort ${DIR}/${SPECIES}.ltrharvest.gff3 > ${DIR}/${SPECIES}.ltrharvest.sorted.gff3

###################
## run ltrdigest ##
###################

gt -j 10 ltrdigest -outfileprefix ${DIR}/ltrdigest/${SPECIES}.ltrdigest -trnas ${DIR}/eukaryotic-tRNAs.fa -hmms ${DIR}/gydb_hmms/GyDB_collection/profiles/*.hmm -- ${DIR}/${SPECIES}.ltrharvest.sorted.gff3 $SPECIES > ${DIR}/${SPECIES}.ltrdigest.gff3
