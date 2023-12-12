#!/bin/bash

#########
#creating the nonoverlapping annotation file with the following precedence: rRNA>knownTE>simplerepeats>unknownrepeat
#rules for resolving the overlaps: just remove the section that was overlapping and at the end filter out anything that is shorter than 15 bp?
#####
#files
rDNAgff=/ohta1/solomiya.hnatovska/genome2/rnammer-1.2/amaranthusRun/AMATA_rdna_output.gff
RMoutput=/ohta1/solomiya.hnatovska/genome2/EDTA/noCDSrun/myRMrun/AMATA_finishedtohypo_renamedfinal.fasta.out
TEanno=/ohta1/solomiya.hnatovska/genome2/EDTA/noCDSrun/AMATA_finishedtohypo_renamedfinal.fasta.mod.EDTA.anno/AMATA_finishedtohypo_renamedfinal.fasta.mod.EDTA.TEanno.split.bed
#first make bed files for each element type
#rDNA
awk 'BEGIN {FS=" ";OFS="\t"} {if ($3~"rRNA") {print $1, $4, $5, $9, $7}}' $rDNAgff > tmp.rRNA.bed
bedtools sort -i tmp.rRNA.bed >tmp.rRNA.sorted.bed

#knownTE
awk 'BEGIN {FS=" ";OFS="\t"} {if ($13 !~ "repeat_region") {print $1, $2, $3, $13, $9}}' $TEanno > tmp.knownTEanno.bed
bedtools sort -i tmp.knownTEanno.bed >tmp.knownTEanno.sorted.bed

#simple repeats and low complexity regions
awk 'BEGIN {FS=" ";OFS="\t"} {if ($11~"Simple_repeat" || $11~"Low_complexity") {print $5, $6, $7, $11, $9}}' $RMoutput > tmp.simplerepeats.bed
bedtools sort -i tmp.simplerepeats.bed >tmp.simplerepeats.sorted.bed

#unknownTE
awk 'BEGIN {FS=" ";OFS="\t"} {if ($13 ~ "repeat_region") {print $1, $2, $3, $13, $9}}' $TEanno > tmp.unknownTEanno.bed
bedtools sort -i tmp.unknownTEanno.bed >tmp.unknownTEanno.sorted.bed

#then use bedtools subtract to remove overlaps
#first subtract TEregions that overlap with rRNA genes
bedtools subtract -a tmp.knownTEanno.sorted.bed -b tmp.rRNA.sorted.bed >tmp.knownTEanno_norRNA.bed
#next combine the file with rRNA and nonoverlapping TEs
cat tmp.rRNA.sorted.bed tmp.knownTEanno_norRNA.bed >tmp.rRNA_TEs.bed
bedtools sort -i tmp.rRNA_TEs.bed >tmp.rRNA_TEs.sorted.bed

#next get the simple repeats regions that dont overlap with  the TE/rRNA file:
bedtools subtract -a tmp.simplerepeats.sorted.bed -b tmp.rRNA_TEs.sorted.bed  >tmp.simplerepeats_norRNA_noTEs.bed
#next combine the non overlapping simple repeats with the rna and TEs
cat tmp.rRNA_TEs.bed tmp.simplerepeats_norRNA_noTEs.bed >tmp.simple_rRNA_TEs.bed
bedtools sort -i tmp.simple_rRNA_TEs.bed >tmp.simple_rRNA_TEs.sorted.bed

#next get unknown repeat regions that dont overlap with any other elements:
bedtools subtract -a tmp.unknownTEanno.sorted.bed -b tmp.simple_rRNA_TEs.bed >tmp.unknowns_nosimplerepeats_norRNA_noTEs.bed
#next combine non overlapping unknowns with the rest
cat tmp.simple_rRNA_TEs.bed tmp.unknowns_nosimplerepeats_norRNA_noTEs.bed >tmp.nonoverlapping_rRNA_TEs_simple_unknowns.bed
bedtools sort -i tmp.nonoverlapping_rRNA_TEs_simple_unknowns.bed >tmp.nonoverlapping_rRNA_TEs_simple_unknowns.sorted.bed

#next measure the lengths of the elements
awk '{if ($3-$2 >=20) {print $0, $3-$2}}' tmp.nonoverlapping_rRNA_TEs_simple_unknowns.sorted.bed >nonoverlapping_rRNA_TEs_simple_unknowns.sorted.minlength20.bed

#conversion to saf file format (elementID chrom start end strand)
awk 'BEGIN {FS=" ";OFS="\t"} {print $4, $1, $2, $3, $5}' nonoverlapping_rRNA_TEs_simple_unknowns.sorted.minlength20.bed >superfamilies_nonoverlapping_rRNA_TEs_simple_unknowns.sorted.minlength20.saf

#cleanup
rm tmp*
