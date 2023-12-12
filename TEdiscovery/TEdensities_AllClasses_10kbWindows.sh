#!/bin/bash

#Scriptfor getting the densities of different types of repeats in 10kb windows across the genome
#################
#get densities
#################
reference_fai=/ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/toshare/AMATA_finishedtohypo_renamedfinal.fasta.fai
#I am going to use the split gff file from EDTA for this:

RM_out=/ohta1/solomiya.hnatovska/genome2/EDTA/noCDSrun/AMATA_finishedtohypo_renamedfinal.fasta.mod.EDTA.anno/AMATA_finishedtohypo_renamedfinal.fasta.mod.EDTA.TEanno.split.bed
ref_anno=/ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/toshare/AMATA_finishedtohypo_renamedfinal.gff


#get start and end positions from fasta index
awk '{print $1 "\t" 0 "\t" $2}' $reference_fai > scaffolds_startend.fai.tmp
#chop scaffolds onto 10000 bp windows
bedops --chop 1000000 scaffolds_startend.fai.tmp > scaffolds_chopped.txt.tmp
#sort the scaffold windows
sort-bed scaffolds_chopped.txt.tmp > scaffolds_chopped_sorted.txt.tmp

#DNA (not including helitrons) only DNA TIR or MITE)
#pull out the scaffold names and ranges of lines with DNA including :DTA DTC DTH DTM DTT
awk '{if ($5 ~ /DT/) print $1, $2, $3}' $RM_out >DNA.bed.tmp
#sort bed files alphanumerically
sort-bed DNA.bed.tmp > DNA.sorted.bed.tmp

#Helitrons
#pull out the scaffold names and ranges of lines with helitron
awk '{if ($5 ~ /Helitron/) print $1, $2, $3}' $RM_out >Helitron.bed.tmp
#sort bed files alphanumerically
sort-bed Helitron.bed.tmp > Helitron.sorted.bed.tmp
#LTRs
#pull out the scaffold names and ranges of lines with LTR
awk '{if ($5 ~ /Copia/) print $1, $2, $3}' $RM_out >COPIA.bed.tmp
awk '{if ($5 ~ /Gypsy/) print $1, $2, $3}' $RM_out >GYPSY.bed.tmp
#sort bed files alphanumerically
sort-bed COPIA.bed.tmp > COPIA.sorted.bed.tmp
sort-bed GYPSY.bed.tmp > GYPSY.sorted.bed.tmp
#genes
#pull out the scaffold names and ranges of lines with GENE
awk '{if ($3 ~ /CDS/) print $1, $4, $5}' $ref_anno >gene.bed.tmp
#1086 in gene.bed.tmp genomic start and end are the same so it wont sort

#sort bed files alphanumerically
sort-bed gene.bed.tmp > gene.sorted.bed.tmp

#find overlap between reference window file and data position file
bedmap --delim '\t' --echo --count --bases-uniq-f scaffolds_chopped_sorted.txt.tmp COPIA.sorted.bed.tmp > COPIA_windowfraction.bed.tmp &
bedmap --delim '\t' --echo --count --bases-uniq-f scaffolds_chopped_sorted.txt.tmp GYPSY.sorted.bed.tmp > GYPSY_windowfraction.bed.tmp &
bedmap --delim '\t' --echo --count --bases-uniq-f scaffolds_chopped_sorted.txt.tmp Helitron.sorted.bed.tmp > Helitron_windowfraction.bed.tmp &
bedmap --delim '\t' --echo --count --bases-uniq-f scaffolds_chopped_sorted.txt.tmp DNA.sorted.bed.tmp > DNA_windowfraction.bed.tmp &
bedmap --delim '\t' --echo --count --bases-uniq-f scaffolds_chopped_sorted.txt.tmp gene.bed.tmp > gene_windowfraction.bed.tmp &

#for some reason this next section needs to be entered by hand does not work in the script
#combine all the densities into one file
#the order of columns should be: scaffold, start of window, end of window, genefraction, dnafraction, helitronfraction, copiafraction, gypsyfraction
paste gene_windowfraction.bed.tmp <(cut -f 5- DNA_windowfraction.bed.tmp) > genes_DNA_window_densities.bed.tmp
paste genes_DNA_window_densities.bed.tmp <(cut -f 5- Helitron_windowfraction.bed.tmp) > genes_DNA_Helitron_window_densities.bed.tmp
paste genes_DNA_Helitron_window_densities.bed.tmp <(cut -f 5- COPIA_windowfraction.bed.tmp) > genes_DNA_Helitron_COPIA_window_densities.bed.tmp
paste genes_DNA_Helitron_COPIA_window_densities.bed.tmp <(cut -f 5- GYPSY_windowfraction.bed.tmp) > genes_DNA_Helitron_COPIA_GYPSY_window_densities.bed.tmp
awk '{print $1,$2,$3,$5,$6,$7,$8,$9}' genes_DNA_Helitron_COPIA_GYPSY_window_densities.bed.tmp >AMATA.genes_DNA_Helitron_COPIA_GYPSY_window_densities.bed

#bedmap --delim '\t' --echo --count  --echo-map-id-uniq scaffolds_chopped_sorted.txt ltr_allpos_sorted.txt > ltr_window_counts.bed &

#cleanup
#rm *.tmp
