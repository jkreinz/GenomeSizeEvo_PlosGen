#############
#get densities for helitrons
#################
#######
#first change helitronscanner output fasta files into files with contig# and start and end position of helitron
#######

#take only lines with >HTG_ from helitron fasta file
awk '($1~">SCF_") {print $0}' ../amaranth.helitronsonly.hel.fa > helitronSCFonly.txt

#make new file with contig names and helitron ranges in columns
awk -F " " 'split ($1,a,"_") {print "SCF_"a[2] "\t" a[4]}' helitronSCFonly.txt >helitronranges.txt

#split the ranges into two columns
awk 'split($2,a,"-") {print $1 "\t" a[1] "\t" a[2]}' helitronranges.txt > helitronranges.split.txt

######reverse helitrons (same thing):

#take only lines with >HTG_ from helitron fasta file
awk '($1~">SCF_") {print $0}' ../amaranth.helitronsonly.rc.hel.fa > helitronSCFonly.rc.txt

#make new file with contig names and helitron ranges in columns
awk -F " " 'split ($1,a,"_") {print "SCF_"a[2] "\t" a[4]}' helitronSCFonly.rc.txt >helitronranges.rc.txt

#split the ranges into two columns
awk 'split($2,a,"-") {print $1 "\t" a[1] "\t" a[2]}' helitronranges.rc.txt > helitronranges.rc.split.txt

#switch the range start and end columns for the reverse strand helitrons
awk '{print $1 "\t" $3 "\t" $2}' helitronranges.rc.split.txt > helitronranges.rc.reversed.txt

#copy forward bed file into new file and add reversed reverse strand file to it
cp helitronranges.split.txt helitronranges.combined.txt
cat helitronranges.rc.reversed.txt >>helitronranges.combined.txt

#sort the rows alphanumerically
sort-bed helitronranges.combined.txt > helitronranges.combined.sorted.bed

#####################

#get start and end positions from fasta index
#awk '{print $1 "\t" 0 "\t" $2}' scaffolds.reduced.fa.fai > scaffolds_startend.fai

#chop scaffolds onto 10000 bp windows
#bedops --chop 100000 ../densities10kb/scaffolds_startend.fai > scaffolds_chopped100kb.txt

#make start of genes 0 based
awk '{print $1 "\t" $2-1 "\t" $3}' helitronranges.combined.sorted.bed > helitronranges.zerostart.bed
#add 1-based end position to allpos file and list every base of each gene/ltr in rows
#awk '{for (i=$2; i<=$NF ; i++) print $1 "\t" i}' ltrONLY_renamed_zerostart.gff3 | awk '{print $1 "\t" $2 "\t" $2+1}' >  ltr_allpos.txt

#sort bed files alphanumerically
#sort-bed scaffolds_chopped100kb.txt > scaffolds_chopped100kb_sorted.txt
#sort-bed ltrONLY_renamed_zerostart.gff3 > ltr_ranges_sorted.txt
#sort-bed ltr_allpos.txt > ltr_allpos_sorted.txt

#find overlap between reference window file and data position file

bedmap --delim '\t' --echo --count --bases-uniq-f ~/densities/densities100kb/scaffolds_chopped100kb_sorted.txt helitronranges.zerostart.bed > helitron.densities.100kb &

#bedmap --delim '\t' --echo --count  --echo-map-id-uniq scaffolds_chopped_sorted.txt ltr_allpos_sorted.txt > ltr_window_counts.bed &

#make column with ltr/gene densities
#awk '{print $0 "\t" $4/10000}' ltr_window_counts.bed > ltr_window_densities.bed



###################
#now run for genes
###################

#get just genes from annotation
#awk '($3=="gene") {print $0}' Amaranthus_tuberculatus_annos1-cds1_REDUCED.gff > amaranth_justgenes.gff

#get start and end positions from fasta index
#awk '{print $1 "\t" 0q "\t" $2}' scaffolds.reduced.fa.fai > scaffolds_startend.fai
#bedops --chop 10000 scaffolds_startend.fai > scaffolds_chopped.txt

#print every position between start and end of gene/LTR
#awk '{print $1 "\t" $4-1 "\t" $5}' amaranth_justgenes.gff  > amaranth_justgenes_zerostart.gff
#add 1-based end position to allpos file
#awk '{for (i=$2; i<=$NF ; i++) print $1 "\t" i}' amaranth_justgenes_zerostart.gff | awk '{print $1 "\t" $2 "\t" $2+1}' >  genes_allpos.txt

#sort bed files
#sort-bed scaffolds_chopped.txt > scaffolds_chopped_sorted.txt
#sort-bed amaranth_justgenes_zerostart.gff > genes_ranges_sorted.txt

#find overlap between reference window file and data position file
#bedmap --delim '\t' --echo --count --bases-uniq-f scaffolds_chopped100kb_sorted.txt ../densities10kb/genes_ranges_sorted.txt > genes_100kbwindow_fractions.bed &

#make column with ltr/gene densities
#awk '{print $0 "\t" $4/10000}' genes_window_counts_allchroms.bed > genes_window_densities.bed


###############
#now combine LTR and gene densities
###############

paste ~/densities/densities100kb/genes_100kbwindow_fractions.bed <(cut -f 4- ltr_100kbwindow_fractions.bed) > genes_LTR_100kbwindow_densities.bed
