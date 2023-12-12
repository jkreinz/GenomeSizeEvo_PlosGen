#!/bin/bash
hs=~/software/HelitronScanner_V1/HelitronScanner/HelitronScanner.jar
lcvs=~/software/HelitronScanner_V1/TrainingSet
genome=/ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/toshare/AMATA_finishedtohypo_renamedfinal.fasta
#heads
java -jar $hs scanHead -lf ${lcvs}/head.lcvs -g $genome -bs 0 -o tuberculatus.helitronscanner.heads -tl 10

#tails
java -jar $hs scanTail -lf ${lcvs}/tail.lcvs -g $genome -bs 0 -o tuberculatus.helitronscanner.tails -tl 10

#pairends
java -jar $hs pairends -hs tuberculatus.helitronscanner.heads -ts tuberculatus.helitronscanner.tails -o pairends

#draw
java -jar $hs draw -p pairends -g $genome -o tuberculatus.helitronsonly --pure

##reverse complementary mode

#heads
java -jar $hs scanHead -lf ${lcvs}/head.lcvs -g $genome -bs 0 -o tuberculatus.helitronscanner.heads.rc -tl 10 --rc

#tails
java -jar $hs scanTail -lf ${lcvs}/tail.lcvs -g $genome -bs 0 -o tuberculatus.helitronscanner.tails.rc -tl 10 --rc

#pairends
java -jar $hs pairends -hs tuberculatus.helitronscanner.heads.rc -ts tuberculatus.helitronscanner.tails.rc -o pairends.rc --rc

#draw
java -jar $hs draw -p pairends.rc -g $genome -o tuberculatus.helitronsonly.rc --pure
