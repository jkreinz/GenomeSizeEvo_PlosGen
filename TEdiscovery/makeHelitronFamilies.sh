#!/bin/bash
set -eou pipefail

GENOME=$1
CPU=$2


#########################
##   tab format output ##
#########################

cat ../${GENOME}.helitronsonly.hel.fa ../${GENOME}.helitronsonly.rc.hel.fa > ${GENOME}.HelitronScanner.draw.all.hel.fa

python helitron_scanner_out_to_tabout.py ${GENOME}.HelitronScanner.draw.all.hel.fa ${GENOME}.HelitronScanner.tabnames.fa > ${GENOME}.HelitronScanner.tabout

#########################
##  Make families      ##
#########################


### think about whether this should be the entire element or the earlier classification based on the terminal 30bp of the helitron.
### remember that the mtec helitrons have lots of N's in their internal regions, so this decision may have been due to data quality.

python get_last_30bp_fasta.py ${GENOME}.HelitronScanner.tabnames.fa > ${GENOME}.HelitronScanner.tabnames.terminal30bp.fa

vsearch -allpairs_global ${GENOME}.HelitronScanner.tabnames.terminal30bp.fa -blast6out ${GENOME}.terminal30bp.allvall.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads=$CPU

silix ${GENOME}.HelitronScanner.tabnames.terminal30bp.fa ${GENOME}.terminal30bp.allvall.out -f DHH -i 0.8 -r 0.8 > ${GENOME}.8080.fnodes
