###############
#This script contains analyses of repeat content and genome size by geographic, environmental, and organismal predictors
#later, it looks to analyze the relationship of flowering time and growth rate with genome size 

#Julia Kreiner, Dec 2023
###############

#FT = flowerting time, BM = biomass, LN = early leaf number, NN = node number at flowering, SW = stem width at flowering, PH = plant height at flowering, FC = flower color, SC=stem colour, Germ = time to germination, COTY = cotyledon width, HYPO= hypocotyl length
###library(data.table)
###library(dplyr)
library(lme4)
#install.packages("car")
library(car)
library(ggplot2)
#install.packages("ggpubr") #for gg arrange function
library("ggpubr")
###if(!require(devtools)) install.packages("devtools")
###devtools::install_github("kassambara/ggpubr")
setwd("~/Documents/UBC_Postdoc/Papers/Solomiya Paper/data/")
#install.packages("MetBrewer")
library(MetBrewer)
teip<-met.brewer("Tiepolo",18)
#install.packages("data.table")
#install.packages("tidyverse")
library(data.table)
library(tidyverse)

#read in dataframes
juliasdata<-read.table("commongarden_meta_familymeans.txt",header=T) #there are a lot of redundant column names in this so check names
tes<-fread("~/allinds_TEcov_mediangenescaled.txt",data.table = F)

#TE coverage by element
tes_copy <- tes %>% group_by(element) %>%
  summarise(across(starts_with("342"):ends_with("556"),mean, na.rm=T))
te_copy_check<-data.frame(element=tes_copy$element, var=apply(tes_copy[,-1], 1, var))
inds<-fread("~/sampleorder.txt")
dropcols<-c("scaf","start","end")

#scale TE coverage by element size
tes_bp <- tes %>%  mutate(sum_bp=end-start) %>%
  mutate(across(starts_with("342"):ends_with("556"), ~ . * sum_bp)) #ignore warning, validated
names(tes_bp)[4]<-"repeatclass"
less_tes<-tes_bp[,-c(1:3,192)]

#take the sum across classes, to get TE abundance by class
tesbp_byclass <- less_tes %>% group_by(repeatclass) %>% 
  summarise_all(sum)
tes_byclass_final<-tesbp_byclass

#do the same for genes
cov_bygene<-fread("~/allsamples_meandepth_bygene.txt",data.table = F)
head(cov_bygene)
ind_means<-apply(cov_bygene[,-c(1:3)],2,median)
names(cov_bygene)[1:3]<-c("scaf","start","end")
cov_bygene_bp <- cov_bygene %>%  mutate(sum_bp=end-start) %>%
  mutate(across(starts_with("V4"):ends_with("V190"), ~ . * sum_bp))  #ignore warning
cov_bygene_bpt<-cov_bygene_bp[,-c(1:3,191)]
cov_bygene_bpt_scaled <-  mapply('/', cov_bygene_bpt, ind_means)
inds<-read.table("~/sampleorder.txt")
inds$genic_content<-colSums(cov_bygene_bpt_scaled)
names(inds)[1]<-"Sample"

#do some reformatting for later
noclass<-t(tes_byclass_final)
tes_byclass_final_t<-as.data.frame((noclass[(-1),]))
tes_byclass_final$class<-gsub("/","_",tes_byclass_final$repeatclass)
names(tes_byclass_final_t)<-tes_byclass_final$class
tes_byclass_final_t <- tes_byclass_final_t %>% mutate_all( as.numeric)
sort(apply((tes_byclass_final_t/1000000),2,median,na.rm=T))

tes_byclass_final_t$Sample<-row.names(tes_byclass_final_t)
hist(tes_byclass_final_t$`18s_rRNA`,breaks=50)
hist(tes_byclass_final_t$nonTIR_helitron,breaks=50)

head(tes_byclass_final_t)

merged_wFlower_heteigen_x <- merge(tes_byclass_final_t, juliasdata, by = "Sample", all.x = TRUE, all.y = TRUE)
merged_wFlower_heteigen_x <- merge(merged_wFlower_heteigen_x, inds, by="Sample", all.x = TRUE, all.y = TRUE)
#remove the 5 samples with an estimated genome size that is bellow 300 and an abnormally high error rate 
merged_wFlower_heteigen <- merged_wFlower_heteigen_x[ which(merged_wFlower_heteigen_x$Error.Rate.x < 0.01040), ]

#get superfamiliy abundances and other heiarchecal levels of repeat abudances
merged_wFlower_heteigen$TEs<-merged_wFlower_heteigen$LTR_Gypsy + merged_wFlower_heteigen$LTR_unknown + merged_wFlower_heteigen$LTR_Copia +
  merged_wFlower_heteigen$TIR_CACTA + merged_wFlower_heteigen$TIR_hAT + merged_wFlower_heteigen$TIR_PIF_Harbinger + 
  merged_wFlower_heteigen$TIR_Mutator + merged_wFlower_heteigen$TIR_Tc1_Mariner +
  merged_wFlower_heteigen$nonLTR_LINE_element + merged_wFlower_heteigen$nonTIR_helitron

merged_wFlower_heteigen$LTRs<-merged_wFlower_heteigen$LTR_Gypsy + merged_wFlower_heteigen$LTR_unknown + merged_wFlower_heteigen$LTR_Copia
merged_wFlower_heteigen$TIRs<-merged_wFlower_heteigen$TIR_CACTA + merged_wFlower_heteigen$TIR_hAT + merged_wFlower_heteigen$TIR_PIF_Harbinger + 
  merged_wFlower_heteigen$TIR_Mutator + merged_wFlower_heteigen$TIR_Tc1_Mariner
merged_wFlower_heteigen$rRNAs<-merged_wFlower_heteigen$`8s_rRNA` +
  merged_wFlower_heteigen$`18s_rRNA` + merged_wFlower_heteigen$`28s_rRNA`

merged_wFlower_heteigen$allrepeats<-merged_wFlower_heteigen$`28s_rRNA`+
merged_wFlower_heteigen$`18s_rRNA`+
merged_wFlower_heteigen$`8s_rRNA`+
merged_wFlower_heteigen$repeat_region+
merged_wFlower_heteigen$Low_complexity+
merged_wFlower_heteigen$Simple_repeat+
merged_wFlower_heteigen$nonTIR_helitron+
merged_wFlower_heteigen$nonLTR_LINE_element+
merged_wFlower_heteigen$LTR_unknown+
merged_wFlower_heteigen$LTR_Gypsy+
merged_wFlower_heteigen$LTR_Copia+
merged_wFlower_heteigen$TIR_PIF_Harbinger+
merged_wFlower_heteigen$TIR_Mutator+
merged_wFlower_heteigen$TIR_hAT+
merged_wFlower_heteigen$TIR_Tc1_Mariner+
merged_wFlower_heteigen$TIR_CACTA

merged_wFlower_heteigen$simprepeats<-
  merged_wFlower_heteigen$repeat_region+
  merged_wFlower_heteigen$Low_complexity+
  merged_wFlower_heteigen$Simple_repeat

#while KMER coverage is already present in the dataset above (as estimated from genomescope2)
#lets also estimate depth based genome size

##
bpcov<-fread("~/column3_sum_and_averages.txt",sep=c(" ")) #per bp cumulative coverage across the genome
bpcov<-bpcov[,c(2,7,11,13)]
names(bpcov)<-c("Sample","cumcov","nsitesmap","avgcov")
bpcov$cumcov<-gsub(",","",bpcov$cumcov)
bpcov$nsitesmap<-gsub(",","",bpcov$nsitesmap)
bpcov$Sample<-gsub(".bpcov,","",bpcov$Sample)

#median genic coveerage
ind_means
inds$genic_mediancov<-ind_means
inds$Sample<-as.character(inds$Sample)
bpcov_wgeniccov<-inner_join(bpcov,inds,by="Sample")
bpcov_wgeniccov$cumcov<-as.numeric(bpcov_wgeniccov$cumcov)
bpcov_wgeniccov$genic_mediancov<-as.numeric(bpcov_wgeniccov$genic_mediancov)
bpcov_wgeniccov$GS_bpdepth<-bpcov_wgeniccov$cumcov/bpcov_wgeniccov$genic_mediancov

#dataframe with depth based genome size estimate
merged_wFlower_heteigen_bpcov<-inner_join(bpcov_wgeniccov,merged_wFlower_heteigen, by="Sample")


########
#some exploratory plots looking at the correlation between different genome size estimates
#KMER from genomescope that uses different upper cutoffs of depth in the SFS
#bp depth from mapped reads
########

ggplot(data=merged_wFlower_heteigen, aes((genic_content/1000000),Genome.Size..Mb..x)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw() +
  labs(y="KMER-based\nGenome Size Estimate (Mb)",
       x="Genic Content (bp)")

ggplot(data=merged_wFlower_heteigen, aes(allrepeats,Genome.Size..Mb..x)) +
  geom_point() +
  geom_smooth(method="loess") +
  theme_bw() +
  labs(y="KMER-based\nGenome Size Estimate (Mb)",
       x="All Repeat Content (Mb)")

P<-ggplot(data=merged_wFlower_heteigen_bpcov, aes(GS_bpdepth/1000000,Genome.Size..Mb..y)) +
  geom_smooth(method="lm") +
  geom_point(aes(color=Model.Fit.y)) +
  geom_abline(slope = 1,intercept = 0) +
  theme_classic() +
  lims(x=c(540,660),y=c(500,620)) +
  theme(legend.position = "none") + 
  scale_color_viridis_c() +
  labs(x="Scaled read-depth estimate of GS",y="Default GenomeScope2\nKMER-based estimate of GS") 

library(ggExtra)
ggExtra::ggMarginal(P, type = "histogram")

cor(merged_wFlower_heteigen_bpcov$GS_bpdepth, merged_wFlower_heteigen_bpcov$Genome.Size..Mb..y )


kmer_100k<-fread("~/round2_kmerGSests_100kcutoff.txt",sep="\t",fill=T)
names(kmer_100k)<-c("Sample","x","minGS","maxGS","y","minModelFit","maxModelFit")
kmer_100k$Sample<-as.character(kmer_100k$Sample)
merged_wFlower_heteigen_bpcov_kmer100<-inner_join(kmer_100k[,c(1,4,7)],merged_wFlower_heteigen_bpcov,by="Sample")
inliers<-merged_wFlower_heteigen_bpcov_kmer100[merged_wFlower_heteigen_bpcov_kmer100$maxGS> 500000000 & merged_wFlower_heteigen_bpcov_kmer100$maxGS< 700000000,]

q1<-ggplot(data=inliers, aes(maxGS/1000000,Genome.Size..Mb..y)) +
  geom_smooth(method="lm") +
  geom_point(aes(color=Model.Fit.y)) +
  geom_abline(slope = 1,intercept = 0) +
  theme_classic() +
 # lims(y=c(500,660),x=c(500,620)) +
  theme(legend.position = "none") + 
  scale_color_viridis_c() +
  labs(x="100k cov threshold GenomeScope2\nKMER-based estimate of GS",y="Default GenomeScope2\nKMER-based estimate of GS") 
#ggExtra::ggMarginal(q, type = "histogram")

summary(inliers$maxGS )

cor(inliers$Genome.Size..b..x, inliers$maxGS )

q2<-ggplot(data=inliers, aes(maxGS/1000000,GS_bpdepth/1000000)) +
  geom_smooth(method="lm") +
  geom_point(aes(color=Model.Fit.y)) +
  geom_abline(slope = 1,intercept = 0) +
  theme_classic() +
  # lims(y=c(540,660),x=c(500,620)) +
  theme(legend.position = "none") + 
  scale_color_viridis_c() +
  labs(y="Scaled read-depth estimate of GS",x="100k cov threshold GenomeScope2\nKMER-based estimate of GS") 

library(cowplot)
plot_grid(q1,q2,nrow=2,ncol=1)
cor(inliers$GS_bpdepth, inliers$maxGS )

##read depth based coverage gives values closer to expected based on 10.1016/j.ympev.2016.12.029
#it is also less biased from the 1:1 line compared to other estimators.
#read depth genome size used for all other anlayses!


###########
#how much does GS estimates vary by?
 
(summ_GSdepth<-summary(merged_wFlower_heteigen_bpcov$GS_bpdepth))

cor(merged_wFlower_heteigen_bpcov$GS_bpdepth, merged_wFlower_heteigen_bpcov$allrepeats)
cor(merged_wFlower_heteigen_bpcov$GS_bpdepth, merged_wFlower_heteigen_bpcov$genic_content.x)

hist(merged_wFlower_heteigen_bpcov$GS_bpdepth/1000000, breaks = 30, xlim=c(540,660),
     xlab="Read-depth based\nGenome Size Estimate",main="")
abline(v=median(merged_wFlower_heteigen_bpcov$GS_bpdepth/1000000),lty="dashed",cex=2)
mean(merged_wFlower_heteigen_bpcov$GS_bpdepth)

######
#coefficient of variation with GS updates

repeats <- data.frame("TIR_CACTA" = merged_wFlower_heteigen$TIR_CACTA, 
                      "TIR_hAT" = merged_wFlower_heteigen_bpcov$TIR_hAT, 
                      "TIR_PIF_Harbinger" = merged_wFlower_heteigen_bpcov$TIR_PIF_Harbinger, 
                      "TIR_Tc1_Mariner" = merged_wFlower_heteigen_bpcov$TIR_Tc1_Mariner,
                      "TIR_Mutator" = merged_wFlower_heteigen_bpcov$TIR_Mutator,
                      "LTR_Copia" = merged_wFlower_heteigen_bpcov$LTR_Copia,
                      "LTR_Gypsy" = merged_wFlower_heteigen_bpcov$LTR_Gypsy, 
                      "LTR_unknown" = merged_wFlower_heteigen_bpcov$LTR_unknown,
                      "nonTIR_helitron" = merged_wFlower_heteigen_bpcov$nonTIR_helitron, 
                      "nonLTR_LINE_element" = merged_wFlower_heteigen_bpcov$nonLTR_LINE_element, 
                      "Low_complexity" = merged_wFlower_heteigen_bpcov$Low_complexity,
                      "Simple_repeat" = merged_wFlower_heteigen_bpcov$Simple_repeat,
                      "repeat_region" = merged_wFlower_heteigen_bpcov$repeat_region,
                      "rDNA18s" = merged_wFlower_heteigen_bpcov$`18s_rRNA`, 
                      "rDNA8s" = merged_wFlower_heteigen_bpcov$`8s_rRNA`, 
                      "rDNA28s" = merged_wFlower_heteigen_bpcov$`28s_rRNA`,
                      "all repeats" = merged_wFlower_heteigen_bpcov$allrepeats,
                      "genic content"=merged_wFlower_heteigen_bpcov$genic_content.x,
                      "genome Size" = merged_wFlower_heteigen_bpcov$GS_bpdepth)

#calculate Coefficient of Variation (CV) for each column in data frame
x <-sapply(repeats, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100)
repeatVCs <- data.frame(x)

#change the column of row names into a categorical variable
repeat_VCs <- cbind(rownames(repeatVCs), data.frame(repeatVCs, row.names=NULL))
repeat_VCs$repeats <- factor(repeat_VCs$`rownames(repeatVCs)`,
                             levels = rev(c("nonTIR_helitron","TIR_hAT","TIR_Tc1_Mariner",
                                            "TIR_PIF_Harbinger","TIR_Mutator","TIR_CACTA",
                                            "LTR_Gypsy","LTR_Copia","LTR_unknown","nonLTR_LINE_element",
                                            "Simple_repeat","repeat_region","Low_complexity",
                                            "rDNA18s","rDNA8s","rDNA28s","all.repeats","genome.Size","genic.content")))

#plot a graph of the CVs
repeat_VCs$type<-c(rep("repeat",16),"other","other","other")
repeat_VCs$type<-factor(repeat_VCs$type, levels=c("repeat","other"))

p1<-ggplot(data=repeat_VCs, aes(x=reorder(repeats, desc(repeats)), y=x)) + 
  geom_bar(stat="identity", width=0.8) +
  # facet_grid(~type,scales="free_x",space="free_x") +
  xlab("Repeat Type") + ylab("Coefficient of variation in bps\nacross individuals") +
  theme_bw() +   theme(axis.text.x = element_text(angle = 90,hjust=1))

p1

#plot a graph of the variation
repeats2<-repeats
y <-sapply(repeats2, function(x) var(x, na.rm=T))
repeatvars <- data.frame(y)
repeatvars$repeats<-rownames(repeatvars)
repeatvars$repeats <- factor(repeatvars$repeats,
                             levels = c("nonTIR_helitron","TIR_hAT","TIR_Tc1_Mariner",
                                        "TIR_PIF_Harbinger","TIR_Mutator","TIR_CACTA",
                                        "LTR_Gypsy","LTR_Copia","LTR_unknown","nonLTR_LINE_element",
                                        "Simple_repeat","repeat_region","Low_complexity",
                                        "rDNA18s","rDNA8s","rDNA28s","all.repeats","genome.Size","genic.content"))
repeatvars$type<-c(rep("repeat",16),"other","other","other")
repeatvars$type<-factor(repeatvars$type, levels=c("repeat","other"))

p2<-ggplot(data=repeatvars[-c(17:19),], aes(x=repeats, y=y)) + 
  geom_bar(stat="identity", width=0.8) +
  xlab("Repeat Type") + ylab("Variance in bps\nacross individuals") +
  theme_bw() +
  # facet_grid(~type,scales="free_x",space="free_x") +
  theme(axis.text.x = element_text(angle = 90,hjust=1))

p3<-ggplot(data=repeatvars, aes(x=repeats, y=y)) + 
  geom_bar(stat="identity", width=0.8) +
  xlab("Repeat Type") + ylab("Variation") +
  theme_bw() +
  #facet_grid(~type,scales="free_x",space="free_x") +
  theme(axis.text.x = element_text(angle = 90,hjust=1))


###Figure 2 left side###
plot_grid(p1,p3,nrow=3,align = "v",rel_heights = c(.5,.5))
p2 #inset for p3
###

###Sup Figure 7###
merged_wFlower_heteigen$genomesize<-merged_wFlower_heteigen$Genome.Size..Mb..x
mean(merged_wFlower_heteigen$genomesize,na.rm=T)
quantile(merged_wFlower_heteigen$genomesize, na.rm = T)
ggplot(merged_wFlower_heteigen, aes(genomesize)) +
  geom_histogram(bins=30, fill="grey80",color="grey20") +
  theme_bw() +
  labs(y="Individual Frequency", x="KMER estimated genome size (Mb)") +
  ylim(0,18) +
  geom_vline(xintercept = mean(merged_wFlower_heteigen$genomesize,na.rm=T),lty="dashed")



##############
#investigation of predictors of TE content, along with FDR estimation for Sup Table 2
merged_wFlower_heteigen<-merged_wFlower_heteigen_bpcov

#relatedness matrix estimated from gemma
relmat<-read.table("commongarden_mergedSNPS_missing2_justCGinds.cXX.txt")
relmat_ID<-read.table("commongarden_mergedSNPS_missing2_justCGinds.fam")
row.names(relmat)<-colnames(relmat)<-relmat_ID$V2
class(relmat)

library(lme4qtl) #program to incorporate the relatedness matrix in linear regressions
#install.packages("lme4qtl")
library(lme4)

#checking alignment of the relatedness matrix with IDs in the dataset
names(merged_wFlower_heteigen)[1]<-"ID"
row.names(relmat)
class(merged_wFlower_heteigen$ID)
row.names(relmat)<-as.factor(row.names(relmat))
class(row.names(relmat))
merged_wFlower_heteigen$ID<-as.character(merged_wFlower_heteigen$ID)

merged_wFlower_heteigen<-merged_wFlower_heteigen[match(row.names(relmat),merged_wFlower_heteigen$ID),] #exclude individuals not found in tree
relmat<-as.matrix(relmat)

names(bpcov_wgeniccov)[1]<-"Sample.1.x"
merged_wFlower_heteigen$Sample.1.x<-as.character(merged_wFlower_heteigen$Sample.1.x)
merged_wFlower_heteigen$genomesize<-merged_wFlower_heteigen$GS_bpdepth/1000000


TIR_model<-relmatLmer(TIRs ~ Long + Lat + Sex.x * Env.x  + fam_K1+(1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
LTR_model<-relmatLmer(LTRs ~ Long + Lat + Sex.x * Env.x  + fam_K1+(1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
rRNAs_model<-relmatLmer(rRNAs ~ Long + Lat + Sex.x * Env.x  + fam_K1+(1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
CATCA_model<-relmatLmer(TIR_CACTA ~ Long + Lat + Sex.x * Env.x  + fam_K1+(1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
tc1_model<-relmatLmer(TIR_Tc1_Mariner ~ Long + Lat + Sex.x * Env.x  + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
pif_model<-relmatLmer(TIR_PIF_Harbinger ~ Long + Lat + Sex.x * Env.x + fam_K1+ (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
hat_model<-relmatLmer(TIR_hAT ~ Long + Lat + Sex.x * Env.x + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
mut_model<-relmatLmer(TIR_Mutator ~ Long + Lat + Sex.x * Env.x + fam_K1+ (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
heli_model<-relmatLmer(nonTIR_helitron ~ Long + Lat + Sex.x * Env.x + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
gyp_model<-relmatLmer(LTR_Gypsy ~ Long + Lat + Sex.x + Env.x  + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
copia_model<-relmatLmer(LTR_Copia ~ Long + Lat + Sex.x * Env.x  + fam_K1  + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
ltru_model<-relmatLmer(LTR_unknown ~  Long + Lat + Env.x + Sex.x *Env.x  + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
line_model<-relmatLmer(nonLTR_LINE_element ~ Long + Lat + Sex.x * Env.x  + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
rna8_model<-relmatLmer(`8s_rRNA` ~ Long + Lat + Env.x + Sex.x *Env.x  + fam_K1+(1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
rna18_model<-relmatLmer(`18s_rRNA` ~ Long + Lat + Sex.x * Env.x + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
rna28_model<-relmatLmer(`28s_rRNA` ~ Long + Lat + Sex.x * Env.x + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
sr_model<-relmatLmer(Simple_repeat ~ Long + Lat + Sex.x * Env.x + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
lowcop_model<-relmatLmer(Low_complexity ~ Long + Lat + Sex.x * Env.x  + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
rr_model<-relmatLmer(repeat_region ~ Long + Lat + Sex.x * Env.x + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))

#MULTIPLE TEST CORRECTION
#rm(list=matching_objects)

matching_objects <- ls(pattern = "model")
matching_objects
pvals<-list()
pred<-list()
yvar<-list()
pval_indajd<-list()
for (i in 1:length(matching_objects)){
  mod_tmp <- get(matching_objects[i])
  pvals[[i]]<-Anova(mod_tmp)$`Pr(>Chisq)`
  pval_indajd[[i]]<-p.adjust(Anova(mod_tmp)$`Pr(>Chisq)`,method="fdr")
  pred[[i]]<-row.names(Anova(mod_tmp))
  yvar[[i]]<-rep(as.character(formula(mod_tmp)[[2]]),each=length(row.names(Anova(mod_tmp))))
  
}

p_values <- unlist(pvals)
p_value_indajd <- unlist(pval_indajd)
y_vars<-unlist(yvar)
preds<-unlist(pred)

# Apply FDR correction
p_adjusted <- p.adjust(p_values, method = "fdr")

pval_df<-data.frame(yvars=y_vars, predictors=preds, uncor_p=p_values, indadj_pval=p_value_indajd)
pval_df$FDR_p<- p.adjust(pval_df$uncor_p, method="fdr")
#pval_df$bh_p<- p.adjust(pval_df$uncor_p, method="hochberg")

###this is Sup Table 2###
pval_df_p05<-pval_df[pval_df$uncor_p < 0.05,]
###

#Visualize plots for Figure 2 - correlates of repeat content by lat, ancestry, and habitat by sex

Anova(TIR_model,type = 3)
summary(TIR_model)

library(lsmeans)
lsmeans(TIR_model, specs=c("Lat","Sex.x"))
TIRout<-lsmip(TIR_model,  ~  Lat, ylab = "TIR Proportion",
              at=list(Lat=c(38,39,40,41,42)),
              xlab="Latitude",
              type="response",
              plotit=F)
head(TIRout)
lat1<-ggplot(data=TIRout, aes(Lat,yvar)) +
  geom_ribbon(aes(ymin=TIRout$yvar-TIRout$SE, ymax=TIRout$yvar+TIRout$SE), alpha=.5, linetype = 0) +
  # geom_point() +
  geom_smooth(method="lm", color="black") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  xlab("Latitude") +
  ylab("TIR (Mbs)") +
  theme(legend.position = "none")  +
  geom_point(data=merged_wFlower_heteigen, aes(Lat,TIRs),alpha=.2) +
  scale_y_continuous(labels = function(x) x/1000000)

Anova(heli_model,type = 3)
summary(heli_model)
summary(heli_model)
r.squaredGLMM(heli_model)
library(lsmeans)
lsmeans(heli_model, specs=c("Lat","Sex.x"))
heli_out<-lsmip(heli_model, ~  Lat, ylab = "Helitron Proportion",
                at=list(Lat=c(38,39,40,41,42)),
                xlab="Latitude",
                type="response",
                plotit=F)
head(test1)
lat2<-ggplot(data=heli_out, aes(Lat,yvar)) +
  geom_ribbon(aes(ymin=heli_out$yvar-heli_out$SE, ymax=heli_out$yvar+heli_out$SE), alpha=.5, linetype = 0) +
  # geom_point() +
  geom_smooth(method="lm",color="black") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  xlab("Latitude") +
  ylab("Helitron (Mbs)") +
  theme(legend.position = "none")  +
  geom_point(data=merged_wFlower_heteigen, aes(Lat,nonTIR_helitron),alpha=.2) +
  scale_y_continuous(labels = function(x) x/1000000)


Anova(copia_model,type = 3)
summary(copia_model)
copia_out<-lsmip(copia_model, ~  Lat, ylab = "Helitron Proportion",
                 at=list(Lat=c(38,39,40,41,42)),
                 xlab="Latitude",
                 type="response",
                 plotit=F)

copia_out2<-lsmip(copia_model, ~  fam_K1, ylab = "Helitron Proportion",
                  at=list(fam_K1=c(0,0.2,0.4,0.6,0.8,1)),
                  xlab="fam_K1",
                  type="response",
                  plotit=F)
head(test1)
lat3<-ggplot(data=copia_out, aes(Lat,yvar)) +
  geom_ribbon(aes(ymin=copia_out$yvar-copia_out$SE, ymax=copia_out$yvar+copia_out$SE), alpha=.5, linetype = 0) +
  # geom_point() +
  geom_smooth(method="lm",color="black") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  xlab("Latitude") +
  ylab("Copia (Mbs)") +
  theme(legend.position = "none")  +
  geom_point(data=merged_wFlower_heteigen, aes(Lat,LTR_Copia),alpha=.2) +
  scale_y_continuous(labels = function(x) x/1000000)

anc1<-ggplot(data=copia_out2, aes(fam_K1,yvar)) +
  geom_ribbon(aes(ymin=copia_out2$yvar-copia_out2$SE, ymax=copia_out2$yvar+copia_out2$SE), alpha=.5, linetype = 0) +
  # geom_point() +
  geom_smooth(method="lm",color="black") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  xlab("var. rudis Ancestry") +
  ylab("Copia (Mbs)") +
  theme(legend.position = "none")  +
  geom_point(data=merged_wFlower_heteigen, aes(fam_K1,LTR_Copia),alpha=.2) +
  scale_y_continuous(labels = function(x) x/1000000)

Anova(ltru_model,type = 3)
summary(ltru_model)
ltru_out2<-lsmip(ltru_model, ~  fam_K1, ylab = "Helitron Proportion",
                 at=list(fam_K1=c(0,0.2,0.4,0.6,0.8,1)),
                 xlab="fam_K1",
                 type="response",
                 plotit=F)

anc2<-ggplot(data=ltru_out2, aes(fam_K1,yvar)) +
  geom_ribbon(aes(ymin=ltru_out2$yvar-ltru_out2$SE, ymax=ltru_out2$yvar+ltru_out2$SE), alpha=.5, linetype = 0) +
  # geom_point() +
  geom_smooth(method="lm",color="black") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  xlab("var. rudis Ancestry") +
  ylab("unknown LTRs (Mbs)") +
  theme(legend.position = "none")  +
  geom_point(data=merged_wFlower_heteigen, aes(fam_K1,LTR_unknown),alpha=.2) +
  scale_y_continuous(labels = function(x) x/1000000)


rRNA_model<-relmatLmer(rRNAs ~ Long + Lat + Sex.x * Env.x   + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(rRNA_model,type = 3)
summary(rRNA_model)

library(lsmeans)
rRNAm<-lsmip(rRNA_model, Sex.x ~  Env.x, ylab = "rRNA Proportion",
             #at=list(Lat=c(38,39,40,41,42)),
             #xlab="Latitude",
             type="response",
             plotit=F)
head(rRNAm)

es2<-ggplot() +
  #geom_point(position=position_dodge(width=0.75)) +
  geom_errorbar(data=rRNAm, aes(x=interaction(Env.x, Sex.x), y=yvar, ymin=rRNAm$yvar-rRNAm$SE, ymax=rRNAm$yvar+rRNAm$SE, lty=Sex.x), position=position_dodge(width=0.75), alpha=.8,width=.6) +
  # geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  scale_linetype_manual(values=c("solid","twodash")) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  xlab("Environment by Sex") +
  ylab("rDNA (Mbs)") +
  theme(legend.position="none") +
  geom_jitter(data=merged_wFlower_heteigen, aes(interaction(Env.x, Sex.x),rRNAs), color="#a048a3",alpha=.4, cex=3,width = .1) +
  scale_y_continuous(labels = function(x) x/1000000) +  geom_errorbar(data=rRNAm, aes(x=interaction(Env.x, Sex.x), y=yvar, ymin=rRNAm$yvar-rRNAm$SE, ymax=rRNAm$yvar+rRNAm$SE, lty=Sex.x), position=position_dodge(width=0.75), alpha=.8,width=.6) +
  theme(axis.text.x = element_text(angle = -45, vjust=-.1))


rr_model<-relmatLmer(repeat_region ~ Long + Lat + Sex.x * Env.x + fam_K1 + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(rr_model,type = 3)
rr_m<-lsmip(rr_model, Sex.x ~  Env.x, ylab = "rRNA Proportion",
            #at=list(Lat=c(38,39,40,41,42)),
            #xlab="Latitude",
            type="response",
            plotit=F)
head(rr_m)
#es3<-


es1<-ggplot() +
  #geom_point(position=position_dodge(width=0.75)) +
  geom_errorbar(data=rr_m, aes(x=interaction(Env.x, Sex.x), y=yvar, ymin=rr_m$yvar-rr_m$SE, ymax=rr_m$yvar+rr_m$SE, lty=Sex.x), position=position_dodge(width=0.75), alpha=.8,width=.6) +
  # geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  scale_linetype_manual(values=c("solid","twodash")) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  xlab("Environment by Sex") +
  ylab("Unknown\nRepeats (Mbs)") +
  theme(legend.position="none") +
  geom_jitter(data=merged_wFlower_heteigen, aes(interaction(Env.x, Sex.x),repeat_region), color="pink",alpha=.4, cex=3,width = .1) +
  scale_y_continuous(labels = function(x) x/1000000) + 
  geom_errorbar(data=rr_m, aes(x=interaction(Env.x, Sex.x), y=yvar, ymin=rr_m$yvar-rr_m$SE, ymax=rr_m$yvar+rr_m$SE, lty=Sex.x), position=position_dodge(width=0.75), alpha=.8,width=.6) +
  theme(axis.text.x = element_text(angle = -45, vjust=-.1))


###Figure 2B####
library(patchwork)
library(cowplot)
space<-plot_spacer()
plot_grid(lat1, anc2, es1, lat3, anc1, es2, lat2, space, space, ncol =3, nrow=3,align = "hv") # sex*env
###


###############
#Correlates with genome size (analyses leading up to Fig 3)
library(MuMIn)
library(r2glmm)
FT_GS <- relmatLmer(FTmeans ~ genomesize +
                      Long + Lat  + Sex.x * Env.x + fam_K1 + (1|ID), 
                    data=merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(FT_GS,type=3)
summary(FT_GS)
r2beta(FT_GS)

PH_GS <- relmatLmer(PHmeans ~ genomesize +
                      Long + Lat  + Sex.x * Env.x + fam_K1 + (1|ID), 
                    data=merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(PH_GS,type=3)
summary(PH_GS)
r2beta(PH_GS)

GS <- relmatLmer(genomesize ~ 
                      Long + Lat  + Sex.x * Env.x + fam_K1 + (1|ID), 
                    data=merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(GS,type=3)
summary(GS)
r2beta(GS)

GS_mip<-lsmip(GS,~ Env.x * Sex.x, plotit=F )


GS1<-ggplot(data=GS_mip, aes(interaction(Env.x,Sex.x),yvar, lty=Sex.x)) +
  #geom_point(position=position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=GS_mip$yvar-GS_mip$SE, ymax=GS_mip$yvar+GS_mip$SE, lty=Sex.x), position=position_dodge(width=0.75),width=.6) +
  geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  theme_bw() +
  xlab("Environment") +
  ylab("Genome Size (Mb)") +
  theme(legend.position="none") +
  geom_point(data=merged_wFlower_heteigen, aes(interaction(Env.x,Sex.x),GS_bpdepth/1000000, fill=Sex.x), position=position_dodge(width=0.75),alpha=.5) +
  scale_linetype_manual(values=c("solid","dashed"))

GS_mip_lat<-lsmip(GS,~ Lat, plotit=F, at=list(Lat=c(38,39,40,41,42)) )

GS2<-ggplot(data=GS_mip_lat, aes(Lat,yvar)) +
  #geom_point(position=position_dodge(width=0.75)) +
  #geom_ribbon(aes(ymin=GS_mip_lat$yvar-GS_mip_lat$SE, ymax=GS_mip_lat$yvar+GS_mip_lat$SE), alpha=.8) +
  #geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  theme_bw() + xlim(38,42) +
  xlab("Latitude") +
  ylab("Genome Size (Mb)") +
  theme(legend.position="none") +
  geom_point(data=merged_wFlower_heteigen, aes(Lat,GS_bpdepth/1000000, fill=Sex.x),alpha=.8) +
  scale_linetype_manual(values=c("solid","dashed"))


###Figure 3###

library(lsmeans)
FT1<-lsmip(FT_GS, ~  Env.x + Sex.x,
           plotit=F)
head(FT1)

flowplot1<-ggplot(data=FT1, aes(interaction(Env.x,Sex.x),yvar, lty=Sex.x)) +
  #geom_point(position=position_dodge(width=0.75)) +
  geom_ribbon(aes(ymin=FT1$yvar-FT1$SE, ymax=FT1$yvar+FT1$SE, lty=Sex.x), position=position_dodge(width=0.75), alpha=.8,width=.6) +
  geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  theme_bw() +
  xlab("Environment") +
  ylab("Flowering Time\n(Family-mean, days)") +
  theme(legend.position="none") +
  geom_point(data=merged_wFlower_heteigen, aes(interaction(Env.x,Sex.x),FTmeans, fill=Sex.x ), position=position_dodge(width=0.75),alpha=.5) +
  scale_linetype_manual(values=c("solid","dashed"))


FT2<-lsmip(FT_GS, ~  genomesize,
           at=list(genomesize = c(540.000000,560.000000,580.000000,600.000000,620.000000,640.000000,650)),
           plotit=F)
head(FT2)
flowplot3<-ggplot(data=FT2, 
                  aes(genomesize, yvar)) +
  geom_ribbon(aes(ymin=FT2$yvar-FT2$SE, ymax=FT2$yvar+FT2$SE),alpha=.6,linetype=0) +
  geom_smooth(method="lm", color="black",alpha=.8) +
  theme_bw() +
  labs(y="Flowering Time\n(Family-mean, days)",x="Genome Size (Mb)") +
  geom_point(data=merged_wFlower_heteigen, aes(genomesize,FTmeans ),alpha=.8)


F3<-lsmip(FT_GS, ~  Lat,
          at=list(Lat=c(38,39,40,41,42)),
          plotit=F)
head(F3)

flowplot2<-ggplot(data=F3, aes(Lat,yvar)) +
  geom_ribbon(aes(ymin=F3$yvar-F3$SE, ymax=F3$yvar+F3$SE),alpha=.6,linetype=0) +
  geom_smooth(method="lm", color="black",alpha=.8) +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  xlim(38,42) +
  theme_bw() +
  ylab("Flowering Time\n(Family-mean, days)") +
  xlab("Latitude") +
  geom_point(data=merged_wFlower_heteigen, aes(Lat, FTmeans),alpha=.8) +
  theme(legend.position = "none")

plot_grid(flowplot2, flowplot1, flowplot3, ncol=3,nrow=1)

PH1<-lsmip(PH_GS, ~  Env.x + Sex.x,
           plotit=F)
head(PH1)

growplot1<-ggplot(data=PH1, aes(interaction(Sex.x,Env.x),yvar, lty=Sex.x)) +
  #geom_point(position=position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=PH1$yvar-PH1$SE, ymax=PH1$yvar+PH1$SE, lty=Sex.x), position=position_dodge(width=0.75), alpha=.8,width=.6) +
  geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  theme_bw() +
  xlab("Environment") +
  ylab("Growth Rate\n(Family-mean, cm/days)") +
  theme(legend.position="none") +
  geom_point(data=merged_wFlower_heteigen[merged_wFlower_heteigen$PHmeans>0.6,], aes(interaction(Sex.x,Env.x),PHmeans, fill=Sex.x ), position=position_dodge(width=0.75),alpha=.5) +
  scale_linetype_manual(values=c("solid","dashed"))

GR2<-lsmip(PH_GS, ~  genomesize,
           at=list(genomesize = c(540,560,580,600,620,640,650)),
           plotit=F)
head(GR2)
growplot3<-ggplot(data=GR2, 
                  aes(genomesize, yvar)) +
  geom_ribbon(aes(ymin=GR2$yvar-GR2$SE, ymax=GR2$yvar+GR2$SE),alpha=.6,linetype=0) +
  geom_smooth(method="lm", color="black",alpha=.8) +
  theme_bw() +
  labs(y="Growth Rate\n(Family-mean, cm/days)",x="Genome Size (Mb)") +
  geom_point(data=merged_wFlower_heteigen[merged_wFlower_heteigen$PHmeans > 0.6,], aes(genomesize,PHmeans ),alpha=.8)


PH3<-lsmip(PH_GS, ~  Lat,
           at=list(Lat=c(38,39,40,41,42)),
           plotit=F)
head(PH3)

growplot2<-ggplot(data=PH3, aes(Lat,yvar)) +
  #geom_ribbon(aes(ymin=PH3$yvar-PH3$SE, ymax=PH3$yvar+PH3$SE),alpha=.6,linetype=0) +
  #geom_smooth(method="lm", color="black",alpha=.8) +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  labs(y="Growth Rate\n(Family-mean, cm/days)",x="Latitude") +
  xlab("Latitude") +
  xlim(38,42) +
  geom_point(data=merged_wFlower_heteigen[merged_wFlower_heteigen$PHmeans > 0.6,], aes(Lat, PHmeans ),alpha=.8) +
  theme(legend.position = "none")

space<-plot_spacer()
###Figure 3!###
plot_grid(space, GS1,GS2, flowplot3, flowplot1, flowplot2, growplot3, growplot1,growplot2, ncol=3,nrow=3,align = "hv")
###




#Sup Figure 4 (map of TE abundance composition by pop)

repeats_fortom <- data.frame("TIR_CACTA" = merged_wFlower_heteigen$TIR_CACTA, 
                      "TIR_hAT" = merged_wFlower_heteigen$TIR_hAT, 
                      "TIR_PIF_Harbinger" = merged_wFlower_heteigen$TIR_PIF_Harbinger, 
                      "TIR_Tc1_Mariner" = merged_wFlower_heteigen$TIR_Tc1_Mariner,
                      "TIR_Mutator" = merged_wFlower_heteigen$TIR_Mutator,
                      "LTR_Copia" = merged_wFlower_heteigen$LTR_Copia,
                      "LTR_Gypsy" = merged_wFlower_heteigen$LTR_Gypsy, 
                      "LTR_unknown" = merged_wFlower_heteigen$LTR_unknown,
                      "nonTIR_helitron" = merged_wFlower_heteigen$nonTIR_helitron, 
                      "nonLTR_LINE_element" = merged_wFlower_heteigen$nonLTR_LINE_element, 
                      "Low_complexity" = merged_wFlower_heteigen$Low_complexity,
                      "Simple_repeat" = merged_wFlower_heteigen$Simple_repeat,
                      "repeat_region" = merged_wFlower_heteigen$repeat_region,
                      "rRNA18s" = merged_wFlower_heteigen$`18s_rRNA`, 
                      "rRNA8s" = merged_wFlower_heteigen$`8s_rRNA`, 
                      "rRNA28s" = merged_wFlower_heteigen$`28s_rRNA`,
                      "lat" = merged_wFlower_heteigen$Lat,
                      "long" = merged_wFlower_heteigen$Long,
                     # "sex" = merged_wFlower_heteigen$Sex.x,
                     # "sample"=merged_wFlower_heteigen$Sample.1.x,
                      #"env"=merged_wFlower_heteigen$Env,
                      #"anc"=merged_wFlower_heteigen$fam_K1,
                      "allrepeats"=merged_wFlower_heteigen$allrepeats
                     )

prop <- function(x, na.rm = FALSE) (x / repeats_fortom$allrepeats)

longrepeats<-melt(repeats_fortom, id.vars=c("allrepeats","lat","long","sex","sample","env","anc"))
vars<-names(repeats_fortom)[-c(17:19)]

repeats_fortom$lat<-as.factor(repeats_fortom$lat)
repeats_fortom$long<-as.factor(repeats_fortom$long)

column_sums_original <- repeats_fortom %>%
  group_by(lat,long) %>%
  summarise(mean_allrepeat=mean(allrepeats,na.rm=T))

formap<-repeats_fortom %>% group_by(lat,long) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>%
  mutate(across(where(is.numeric) , ~ . / allrepeats)) 

formap$allrepeats<-column_sums_original$mean_allrepeat
formap$long<-as.numeric(as.character(formap$long))
formap$lat<-as.numeric(as.character(formap$lat))
long_formap<-melt(formap, id.vars=c("lat","long","allrepeats"))

formap$lat_j<-jitter(formap$lat,amount = .4)
formap$long_j<-jitter(formap$long, amount = .2)

#get map base
library(raster)
states    <- c("New York","Pennsylvania","Maryland","West Virginia","Virginia","Kentucky","Ohio",
               "Michigan","Indiana","Illinois","Wisconsin","Minnesota","Iowa","Missouri",
               "Kansas","Nebraska","South Dakota","North Carolina","Tennessee","Mississippi", "Oklahoma",
               "Lake Michigan","Lake Ontario","Lake Superior")
provinces <- c("Ontario")

us <- getData("GADM",country="USA",level=1)
canada <- getData("GADM",country="CAN",level=1)

us.states <- us[us$NAME_1 %in% states,]
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

#lakes
lakes <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

#ocean
ocean <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'ocean', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

# rivers
rivers <- rnaturalearth::ne_download(scale = 10, 
                                     type = 'rivers_lake_centerlines', 
                                     category = 'physical')  %>% 
  sf::st_as_sf(lakes110, crs = 4269)

library("rnaturalearth")
#install.packages("rnaturalearthhires")
library("rnaturalearthhires")
library(rnaturalearth)
library(dplyr)
library(raster)
library(sf)
library(tidyverse)
library(ggrepel)

#get US states outlines
in_sf <- ne_states(geounit = "United States of America",
                   returnclass = "sf")

uss <- st_as_sf(us.states) %>% 
  mutate(
    lon = map_dbl(geometry, ~st_centroid(.x)[[1]]),
    lat = map_dbl(geometry, ~st_centroid(.x)[[2]]))

#move Michigan label and add Ontario label
uss[uss$NAME_1=="Michigan",]$lon<--85.0554
uss[uss$NAME_1=="Michigan",]$lat<-44.00902
uss<-uss %>% add_row(NAME_1 = "Ontario", lon = -78.5554, lat=45)


library(ggthemes)
plain1<- 
  ggplot()+
  geom_path(data=us.states,aes(x=long,y=lat,group=group))+
  geom_path(data=ca.provinces, aes(x=long,y=lat,group=group), fill="grey")+
  coord_map() +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "black")  +
  geom_sf(data = ocean,
          mapping = aes(geometry = geometry),
          color = "black")  +
  geom_sf(data = rivers,
          mapping = aes(geometry = geometry),
          color = "grey80",alpha=.75)  +
  #theme_nothing() +
  scale_x_continuous(limits=c(-97,-75)) +
  scale_y_continuous(limits=c(36,45)) +
  coord_sf(xlim=c(-96,-75)) +
  geom_text(
    data = uss,
    aes(x = lon,
        y = lat,
        label = NAME_1),cex=3, alpha=.8) 

purp<-met.brewer("VanGogh1",6)

##Sup Fig 4###
plain1 + 
  #geom_jitter(data=coord_byenv[coord_byenv$dataset == "Paired",],aes(y = lat, x = long, color=env, shape=env), size=5, alpha=.9, height=.2, width=.15) +
  geom_scatterpie(aes(x=long_j, lat_j, r=(allrepeats/1000000000)*1.2),
                  data=formap, alpha=.8, cols=vars) +
 # coord_equal() +
  scale_fill_manual(values=c(teip[10:6],teip[16:18],teip[1], "grey40", "#CC79A7","lightpink","pink", purp[1:3])) +
  #scale_shape_manual(values=c("circle","triangle")) +
  #labs(colour="Contemporary\nPaired Population\nEnvironment") +
  #new_scale_color() +
  #geom_jitter(data=cg[cg$dataset == "Herbarium",],aes(y = lat, x = long, color=year, shape=env), size=3, height=.1,width=.05, alpha=.9) +
  #scale_color_viridis_c() +
  xlab("Longitude") +
  ylab("Latitude")

head(longrepeats)

###Sup Fig 3####
repeats_fortom <- data.frame("TIR_CACTA" = merged_wFlower_heteigen$TIR_CACTA, 
                             "TIR_hAT" = merged_wFlower_heteigen$TIR_hAT, 
                             "TIR_PIF_Harbinger" = merged_wFlower_heteigen$TIR_PIF_Harbinger, 
                             "TIR_Tc1_Mariner" = merged_wFlower_heteigen$TIR_Tc1_Mariner,
                             "TIR_Mutator" = merged_wFlower_heteigen$TIR_Mutator,
                             "LTR_Copia" = merged_wFlower_heteigen$LTR_Copia,
                             "LTR_Gypsy" = merged_wFlower_heteigen$LTR_Gypsy, 
                             "LTR_unknown" = merged_wFlower_heteigen$LTR_unknown,
                             "nonTIR_helitron" = merged_wFlower_heteigen$nonTIR_helitron, 
                             "nonLTR_LINE_element" = merged_wFlower_heteigen$nonLTR_LINE_element, 
                             "Low_complexity" = merged_wFlower_heteigen$Low_complexity,
                             "Simple_repeat" = merged_wFlower_heteigen$Simple_repeat,
                             "repeat_region" = merged_wFlower_heteigen$repeat_region,
                             "rRNA18s" = merged_wFlower_heteigen$`18s_rRNA`, 
                             "rRNA8s" = merged_wFlower_heteigen$`8s_rRNA`, 
                             "rRNA28s" = merged_wFlower_heteigen$`28s_rRNA`,
                             "lat" = merged_wFlower_heteigen$Lat,
                             "long" = merged_wFlower_heteigen$Long,
                              "sex" = merged_wFlower_heteigen$Sex.x,
                              "sample"=merged_wFlower_heteigen$Sample.1.x,
                             #"env"=merged_wFlower_heteigen$Env,
                             #"anc"=merged_wFlower_heteigen$fam_K1,
                             "allrepeats"=merged_wFlower_heteigen$allrepeats
)

prop <- function(x, na.rm = FALSE) (x / repeats_fortom$allrepeats)
longrepeats<-melt(repeats_fortom, id.vars=c("allrepeats","lat","long","sex","sample"))
Frep<-ggplot(data=longrepeats[longrepeats$sex == "F",], aes(reorder(sample,lat), value, fill=variable)) + 
  geom_bar(stat="identity",width=1) +
  #xlab("Repeat Type") + ylab("Coefficient of Variation (CV)") +
  theme_bw() +
  scale_fill_manual(values=c(teip[10:6],teip[16:18],teip[1], "grey40", "#CC79A7","lightpink","pink", purp[1:3] )) +
  scale_color_manual(values=c("black","white")) +
  theme(axis.text.x = element_blank())

Mrep<-ggplot(data=longrepeats[longrepeats$sex == "M",], aes(reorder(sample,lat), value, fill=variable)) + 
    geom_bar(stat="identity",width=1) +
    #xlab("Repeat Type") + ylab("Coefficient of Variation (CV)") +
    theme_bw() +
  scale_color_manual(values=c("black","white")) +
    scale_fill_manual(values=c(teip[10:6],teip[16:18],teip[1], "grey40", "#CC79A7","lightpink","pink", purp[1:3] )) + 
  theme(axis.text.x = element_blank())

#install.packages("lemon")
library(lemon)
grid_arrange_shared_legend(Frep,Mrep,nrow = 1, ncol=2,position = "right")

















#####################################
m <- lm(formula = genomesize ~ fam_K1  + Long + Lat, data = merged_wFlower_heteigen)
summary(m)

#this comes from evo paper PC merged
merged_wFlower_heteigen$sample2
merged$sample2<-as.numeric(merged$sample)
merged$ancestry<-merged$V3
merged_wFlower_heteigen_more<-inner_join(merged_wFlower_heteigen,merged,by="sample2")

cor(merged_wFlower_heteigen_more$ancestry, merged_wFlower_heteigen_more$V3.x,use="complete.obs")
m <- lm(formula = genomesize ~ ancestry + fam_K1 + Long + Lat, data = merged_wFlower_heteigen_more)
summary(m)

library(lsmeans)
lsmeans(m , specs=c("fam_K1","V3"))
lsmip(m,  ~  fam_K1, ylab = "Genome Size",
              at=list(fam_K1=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)),
              #xlab="Latitude",
              #type="response",
              plotit=T)

lsmip(m,  ~  V3, ylab = "Genome Size",
      at=list(V3=c(seq(from=min(merged_wFlower_heteigen$V3,na.rm = T), to=max(merged_wFlower_heteigen$V3,na.rm = T),length.out = 8))),
      #xlab="Latitude",
      #type="response",
      plotit=T)

head(TIRout)
lat1<-ggplot(data=TIRout, aes(Lat,yvar, color=Sex.x)) +
  geom_ribbon(aes(ymin=TIRout$yvar-TIRout$SE, ymax=TIRout$yvar+TIRout$SE, fill=Sex.x), alpha=.15, linetype = 0) +
  # geom_point() +
  geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  xlab("Latitude") +
  ylab("TIR Mbs\nwithin Individuals") +
  theme(legend.position = "none")  +
  geom_point(data=merged_wFlower_heteigen, aes(Lat,TIRs),alpha=.2) +
  scale_y_continuous(labels = function(x) x/1000000)

merged_wFlower_heteigen_bpcov$rRNA<-merged_wFlower_heteigen_bpcov$`18s_rRNA` + merged_wFlower_heteigen_bpcov$`8s_rRNA` + merged_wFlower_heteigen_bpcov$`28s_rRNA`
###########################################  Multiple Linear Regressions:



######################################################################################  
relmat<-read.table("commongarden_mergedSNPS_missing2_justCGinds.cXX.txt")
relmat_ID<-read.table("commongarden_mergedSNPS_missing2_justCGinds.fam")

#relmat<-as.matrix(fread("~/commongarden_finalfilt_ldpruned.mibs"))
#relmat_ID<-read.table("~/commongarden_finalfilt_ldpruned.mibs.id")
#View(relmat_ID)
row.names(relmat)<-colnames(relmat)<-relmat_ID$V2
class(relmat)
#relmat<-as.matrix(relmat)

#devtools::install_github("variani/lme4qtl")
library(lme4qtl) #this is the program I found to incorporate the relatedness matrix
#install.packages("lme4qtl")
library(lme4)
?relmatLmer
?lme4
#install.packages("remotes")
#remotes::install_github("variani/lme4qtl")
#match order of relatedness matrix to df
names(merged_wFlower_heteigen)[1]<-"ID"
row.names(relmat)
class(merged_wFlower_heteigen$ID)
row.names(relmat)<-as.factor(row.names(relmat))
class(row.names(relmat))
merged_wFlower_heteigen$ID<-as.character(merged_wFlower_heteigen$ID)

merged_wFlower_heteigen<-merged_wFlower_heteigen[match(row.names(relmat),merged_wFlower_heteigen$ID),] #exclude individuals not found in tree
relmat<-as.matrix(relmat)

names(bpcov_wgeniccov)[1]<-"Sample.1.x"
merged_wFlower_heteigen$Sample.1.x<-as.character(merged_wFlower_heteigen$Sample.1.x)
merged_wFlower_heteigen_bpcovind<-inner_join(merged_wFlower_heteigen,bpcov_wgeniccov,by="Sample.1.x")
#Equations 6,7,8:
###########################################################################################

###########################################################################################

#########################################
#Everything in

nrow(merged_wFlower_heteigen)
length(cov_bygene_scaled_tr[,which(cov_bygene$V2 == 13152174)])
NADHcopy<-as.data.frame((cov_bygene_scaled_tr[,which(cov_bygene$V2 == 13152174)]))
names(NADHcopy)<-"copygeno"
NADHcopy$Sample.1.x<-as.integer(phenos_ordered$Sample)
library(dplyr)
mNADH<-inner_join(merged_wFlower_heteigen,NADHcopy, by="Sample.1.x")

row.names(relmat)
class(mNADH$ID)
row.names(relmat)<-as.factor(row.names(relmat))
class(row.names(relmat))
mNADH$ID<-as.character(mNADH$ID)

mNADH<-mNADH[match(row.names(relmat),mNADH$ID),] #exclude individuals not found in tree
head(mNADH)
#write.table(mNADH,"~/NADH_forrelatedness_BP.txt",quote=F, col.names = T, row.names = F)
#mNADH<-read.table("~/Downloads/Solomiya Paper/data/NADH_forrelatedness.txt",header=T)

ggplot(mNADH, aes(fam_K1, copygeno)) +
  geom_point() +
  geom_smooth(method="lm")

ggplot(mNADH, aes(Lat, copygeno)) +
  geom_point() +
  geom_smooth(method="lm")



some_mNADH<-mNADH[,c("FTmeans","PHmeans","rRNA","repeat_region","genomesize","copygeno","Sex.x","ID")]
some_mNADH$rRNA_scale<-scale(some_mNADH$rRNA)
some_mNADH$genomesize_scale<-scale(some_mNADH$genomesize)
some_mNADH$copygeno_scale<-scale(some_mNADH$copygeno)
some_mNADH$rr_scale<-scale(some_mNADH$repeat_region)
FT_long<-melt(some_mNADH, id.vars=c("FTmeans","PHmeans","Sex.x","ID"))


FT_long$variable <- factor(FT_long$variable, rev(levels(FT_long$variable)[c(1,5,2,8,3,6,4,7)]))

FT_long %>% filter(variable != "rRNA" & variable != "genomesize" & variable != "copygeno" & variable != "repeat_region" ) %>%
  ggplot( aes( value, FTmeans, group=variable, color=variable)) +
  geom_point(alpha=.8,size=2) +
  geom_smooth(method="lm",alpha=.2) +
  labs(x="Scaled Predictor Value", y="Flowering Time Family Mean") +
  theme_bw() +
  theme(legend.position = "right") +
  scale_color_manual(labels=c("NADH Cluster Copy Number","Genome Size","Unknown Repeat Regions","rRNA"),
                     values=c("forestgreen", "grey60", "pink",purp[2]))

FT_long %>% filter(variable != "rRNA" & variable != "genomesize" & variable != "copygeno" & variable != "repeat_region" & variable != "copygeno_scale") %>%
  ggplot( aes( value, PHmeans, group=variable, color=variable)) +
  geom_point(alpha=.8,size=2) +
  geom_smooth(method="lm",alpha=.2) +
  labs(x="Scaled Predictor Value", y="Plant Height Family Mean") +
  theme_bw() +
  ylim(0.75,2) +
  theme(legend.position = "right") +
  scale_color_manual(labels=c("Genome Size","Unknown Repeat Regions","rRNA"),
                     values=c( "grey60", "pink",purp[2]))

model_full_geno2 <- relmatLmer(FTmeans ~ GS_scale  + 
                                 + Lat + Sex.x+ Env.x +  (1|ID), 
                               mNADH, relmat = list(ID=relmat))
model_full_geno3 <- relmatLmer(FTmeans ~ copygeno  + 
                                 + Lat + Sex.x+ Env.x +  (1|ID), 
                               mNADH, relmat = list(ID=relmat))

model_full <- relmatLmer(FTmeans ~ Genome.Size..Mb..y + Long + Lat + Sex.x + Env.x + (1|ID), 
                         mNADH, relmat = list(ID=relmat))
model_some <- relmatLmer(FTmeans ~ copygeno + Genome.Size..Mb..y + (1|ID), 
                         mNADH, relmat = list(ID=relmat))
model_geno <- relmatLmer(FTmeans ~ copygeno + (1|ID), 
                         mNADH, relmat = list(ID=relmat))
model_size <- relmatLmer(FTmeans ~ Genome.Size..Mb..y + (1|ID), 
                         mNADH, relmat = list(ID=relmat))
Anova(model_full_geno,type = 3)
r.squaredGLMM(model_full_geno)

summary(model)
r.squaredGLMM(model_full_geno)
r.squaredGLMM(model_full_geno2)
Anova(model_full_geno,type = 3)
AIC(model_full_geno, model_full_geno2, model_full_geno3)

r.squaredGLMM(model_full)
r.squaredGLMM(model_some)
r.squaredGLMM(model_geno)
r.squaredGLMM(model_size)



library(MuMIn)
r.squaredGLMM(model_full)
r.squaredGLMM(model_red)

r.squaredGLMM(model_red2)


model_full <- relmatLmer(PHmeans ~ Genome.Size..Mb..y + Long + Lat + Sex.x + Env.x + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
model_red <- relmatLmer(PHmeans ~  Long + Lat + Sex.x + Env.x + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
model_red2 <- relmatLmer(PHmeans ~  Genome.Size..Mb..y + Long + Lat  + Env.x + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(model_full)
Anova(model_red)
Anova(model_red2)


library(MuMIn)
r.squaredGLMM(model_full)
r.squaredGLMM(model_red)
r.squaredGLMM(model_red2)


FT_rRNA_model <- relmatLmer(FTmeans ~ genomesize + TEs + rRNA + Simple_repeat + Low_complexity + repeat_region + Long + Lat + Sex.x + Env.x + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(FT_rRNA_model,type = 3)
summary(FT_rRNA_model)

library(lsmeans)
test1<-lsmip(FT_rRNA_model, ~  rRNA, ylab = "Flowering Time",
             at=list(rRNA=c(0.01,0.02, 0.03, 0.04, 0.05, 0.06, 0.07)),
             xlab="rRNA",
             type="response",
             plotit=F)
head(test1)

ggplot(data=test1, aes(rRNA,yvar)) +
  geom_ribbon(aes(ymin=test1$yvar-test1$SE, ymax=test1$yvar+test1$SE), alpha=.15, linetype = 0) +
  # geom_point() +
  geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  ylab("Flowering Time") +
  xlab("rRNA proportion") +
  theme(legend.position = "none") 

ggplot(data=merged_wFlower_heteigen, aes(rRNA,yvar)) +
  geom_ribbon(aes(ymin=test1$yvar-test1$SE, ymax=test1$yvar+test1$SE), alpha=.15, linetype = 0) +
  # geom_point() +
  geom_smooth(method="lm") +
  #geom_smooth(data=test2, lty="dashed",aes(alpha=.8)) +
  # ylim(0.03, 0.085) +
  theme_bw() +
  ylab("Flowering Time") +
  xlab("rRNA proportion") +
  theme(legend.position = "none") 

model <- relmatLmer(FTmeans ~ TIRs + nonTIR_helitron + LTRs + nonLTR_LINE_element + rRNA + Simple_repeat + Low_complexity + repeat_region + Long + Lat + Sex.x + Env.x + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(model,type = 3)
summary(model)


#same models as above but explaining biomass:
merged_wFlower_heteigen_PH<-merged_wFlower_heteigen[merged_wFlower_heteigen$PHmeans > 0.55,]

relmat_PH<-relmat[na.omit(merged_wFlower_heteigen_PH$ID),na.omit(merged_wFlower_heteigen_PH$ID)]

model <- relmatLmer(PHmeans ~ scale(Genome.Size..Mb..y) + scale(rRNA) + Long + Lat + Sex.x + Env.x + (1|ID), merged_wFlower_heteigen_PH, relmat = list(ID=relmat_PH))
Anova(model,type = 3)

model <- relmatLmer(PHmeans ~ scale(Genome.Size..Mb..y) + scale(rRNA) + scale(repeat_region) 
                    + Long + Lat + Sex.x + Env.x + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(model,type = 3)


cor(merged_wFlower_heteigen$PHmeans, merged_wFlower_heteigen$rRNA,use = "complete.obs")
cor(merged_wFlower_heteigen$PHmeans, merged_wFlower_heteigen$genomesize,use = "complete.obs")

summary(model)


ggplot(data=merged_wFlower_heteigen, aes(Genome.Size..Mb..y, PHmeans)) +
  geom_point() +
  geom_smooth(method="lm") +
  ylim(0.6,2)

model <- relmatLmer(PHmeans ~ TEs + rRNA + Simple_repeat + Low_complexity + repeat_region + Long + Lat + Sex.x + Env.x + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(model,type = 3)
summary(model)

model <- relmatLmer(PHmeans ~ TIRs + nonTIR_helitron + LTRs + nonLTR_LINE_element + rRNA + Simple_repeat + Low_complexity + repeat_region + Long + Lat + Sex.x + Env.x + (1|ID), merged_wFlower_heteigen, relmat = list(ID=relmat))
Anova(model,type = 3)
summary(model)



