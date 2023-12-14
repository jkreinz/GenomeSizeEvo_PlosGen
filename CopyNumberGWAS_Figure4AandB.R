#Script for performing a copy number GWAS by gene and TE annotation
#Julia Kreiner, Dec 2023

#load in data and clean
library(data.table)
library(tidyverse)
cov_bygene<-fread("~/allsamples_meandepth_bygene.txt",data.table = F)
cov_byte<-fread("~/allinds_TEcov_mediangenescaled.txt",data.table = F)
head(cov_bygene)

#get median gene and TE copy number for scaling, by ind
ind_means<-apply(cov_bygene[,-c(1:3)],2,median)
ind_means_te<-apply(cov_byte[,-c(1:4)],2,median)

#rescale gene coverage based on median cov per sample to get copy number
cov_bygene_tr <- as.data.frame(t(cov_bygene[,-c(1:3)]))
cov_bygene_scaled <- cov_bygene_tr %>% 
  mutate(across(everything(), ~ . / ind_means))

#rounding to assign copy number (dont use this)
presabs <- ifelse(cov_bygene_scaled < 0.75, 0, ifelse(cov_bygene_scaled <= 1.5, 1, 2))

#assign sample names
inds<-read.table("~/sampleorder.txt")
indsmedians<-data.frame(cov_median=ind_means, sample=inds$V1)
#write.table(indsmedians, "median_geniccov_bysample.txt",quote=F,row.names = F, col.names = F)

row.names(cov_bygene_scaled)<- inds$V1
row.names(presabs)<- inds$V1

#read in phenos
phenos<-read.table("~/metadata.txt",header=T)

#drop coverage data for samples without phenotypes
todrop<-setdiff(row.names(cov_bygene_scaled), as.character(phenos$Sample))
cov_bygene_scaled_tr <- cov_bygene_scaled[!(rownames(cov_bygene_scaled) %in% todrop), ]
presabs <- presabs[!(rownames(presabs) %in% todrop), ]

cov_byte_t<-t(cov_byte[,-c(1:4)])
todrop<-setdiff(row.names(cov_byte_t), as.character(phenos$Sample))
cov_byte_t_matched <- cov_byte_t[!(rownames(cov_byte_t) %in% todrop), ]

#reorder phenos to match row order of cov data
phenos_ordered <- phenos[order(match(phenos$Sample, row.names(cov_bygene_scaled))), ]
phenos_ordered2 <- phenos[order(match(phenos$Sample, row.names(cov_byte_t_matched))), ]

#run flowering time GWAS based on gene copy number
pval<-list()
tval<-list()
slope<-list()
slope_error<-list()
temp<-data.frame("geno"=NA, "pheno"="NA")
b0<-list()

#here you could sub the cov_by_gen_scaled_tr dataframe, which represents a matrix of individual copy number at each gene...
#for cov_byte_t_matched, if you wanted to perform a TE copy number GWAS

for (i in 1:(ncol(cov_bygene_scaled_tr))) {
  temp<-data.frame(geno=as.numeric(cov_bygene_scaled_tr[,i]), pheno=phenos_ordered$FTmeans)
  
  lm1<-lm(data=temp, unlist(pheno) ~ temp[,1])
  temp2<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-temp2$`Pr(>|t|)`[2]
  tval[[i]]<-temp2$`t value`[2]
  slope[[i]]<-temp2$Estimate[2]
  slope_error[[i]]<-temp2$`Std. Error`[2]
  b0[[i]]<-temp2$Estimate[1]

}

names(cov_bygene)[1:3]<-c("scaf","start","end")
#names(cov_byte)[1:4]<-c("scaf","start","end")

#results<-data.frame(scaf=cov_byte$scaf, start=cov_byte$start, end=cov_byte$end, pval=unlist(pval),tval=unlist(tval),slope=unlist(slope),slope_error=unlist(slope_error))
results<-data.frame(scaf=cov_bygene$scaf, start=cov_bygene$start, end=cov_bygene$end, pval=unlist(pval),tval=unlist(tval),slope=unlist(slope),slope_error=unlist(slope_error))
#write.table(results,"FT_copynumberGWAS_results.txt",quote=F)
#results<-read.table("~/FT_copynumberGWAS_results.txt")

hist(-log10(results$pval))
results$scaf <- gsub("Scaffold_","",results$scaf)

results$fdr<-p.adjust(results$pval, method = "fdr") 
results_fdrsig<-results[results$fdr < 0.05,]

results$bon<-p.adjust(results$pval, method = "bonferroni") 
results_bonsig<-results[results$bon < 0.05,]

###Sup Figure 9 Manhattan plots
results$scaf<-as.numeric(results$scaf)
ggplot(data=results, aes(start/1000000, -log10(pval))) +
  geom_point(alpha=.15) +
  facet_grid(~scaf, scales = "free_x",space = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab("Genomic Position (Mb)") +
  geom_hline(yintercept = -log10(max(results_fdrsig$pval)),lty="dashed") +
  geom_hline(yintercept = -log10(max(results_bonsig$pval)),lty="dashed")

scaf10_results<-results[results$scaf==10,]
library(tidyverse)


###Figure 4 A###
p1<-scaf10_results %>%
  filter(start > 12500000 & start < 14000000) %>%
  ggplot( aes(start/1000000, -log10(pval))) +
  geom_point(alpha=.5) +
  facet_grid(~scaf, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab("Genomic Position (Mb)") +
  geom_hline(yintercept = -log10(max(results_fdrsig$pval,na.rm=T)),lty="dashed",color="grey60") +
  geom_hline(yintercept = -log10(max(results_bonsig$pval,na.rm=T)),lty="dashed",color="grey80")+
  geom_vline(xintercept = 13152174/1000000,lty="dashed")


cov_bygene_scaled_tr2<-as.data.frame(t(cov_bygene_scaled))
cov_bygene_scaled_tr2$scaf<-cov_bygene[,1]
#cov_bygene_scaled_tr2$scaf <- as.numeric(cov_bygene_scaled_wpos$scaf)
cov_bygene_scaled_tr2$start<-cov_bygene[,2]
cov_bygene_scaled_tr2$end<-cov_bygene[,3]
cov_bygene_scaled_tr2$midpos<-cov_bygene_scaled_tr2$start+((cov_bygene_scaled_tr2$end-cov_bygene_scaled_tr2$start)/2)

library(reshape2)
names(cov_bygene_scaled_tr2)
head(cov_bygene_scaled_tr2)
some_long<-melt(cov_bygene_scaled_tr2, id.vars = c("scaf","start","end","midpos"))

some_long$scaf<-gsub("Scaffold_","",some_long$scaf)
some_long$scaf<-as.numeric(some_long$scaf)

#some_long %>% 
#  filter(scaf==11) %>%
#  #filter(start > 12500000 & start < 14000000) %>%
#  ggplot( aes(end/1000000, value,color=variable)) +
#  geom_point(alpha=.1) +
#  facet_grid(~scaf, scales = "free_x") + theme_bw() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#        legend.position = "none") + 
#  xlab("Genomic Position (Mb)") +
#  ylab("Copy Number") +
#  ylim(0,10) + 
#  geom_vline(xintercept = 21401220/1000000,lty="dashed") +
#  xlim(19.401220,23.401220)


standard_error <- function(x) sd(x) / sqrt(length(x)) # Create own function

p2<-some_long %>% 
  filter(scaf==10) %>%
  filter(start > 12500000 & start < 14000000) %>%
  group_by(scaf,midpos) %>% summarise(meancopy=mean(value),secopy=standard_error(value)) %>%
  ggplot( aes(midpos/1000000, meancopy)) +
  geom_errorbar(aes(ymin=meancopy-secopy, ymax=meancopy+secopy)) +
  geom_point(alpha=.5) +
  facet_grid(~scaf, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") + 
  xlab("Genomic Position (Mb)") +
  ylab("Mean Copy\nNumber") +
  geom_vline(xintercept = 13152174/1000000,lty="dashed")

some_long %>% 
  #filter(scaf==10) %>%
  #filter(start > 12500000 & start < 14000000) %>%
  group_by(scaf,midpos) %>% summarise(meancopy=mean(value),secopy=standard_error(value)) %>%
  ggplot( aes(midpos/1000000, meancopy)) +
  geom_errorbar(aes(ymin=meancopy-secopy, ymax=meancopy+secopy)) +
  geom_point(alpha=.5) +
  facet_grid(~scaf, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") + 
  xlab("Genomic Position (Mb)") +
  ylab("Mean Copy\nNumber") +
  ylim(0,10)

some_long %>% 
  #filter(scaf==10) %>%
  #filter(start > 12500000 & start < 14000000) %>%
  group_by(scaf,midpos) %>% summarise(secopy=standard_error(value),secopy=standard_error(value)) %>%
  ggplot( aes(midpos/1000000, secopy)) +
  #geom_errorbar(aes(ymin=meancopy-secopy, ymax=meancopy+secopy)) +
  geom_point(alpha=.5) +
  facet_grid(~scaf, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") + 
  xlab("Genomic Position (Mb)") +
  ylab("Mean Copy\nNumber") +
  ylim(0,10)

###Figure 4A###
library(cowplot)
plot_grid(p1,p2, nrow=2,align = "hv")
###

some_long %>% filter(scaf==10) %>%
  filter(start > 13000000 & start < 13500000) %>%
  group_by(scaf,midpos) %>% summarise(meancopy=mean(value),secopy=standard_error(value)) %>%
  ggplot( aes(midpos/1000000, meancopy)) +
  geom_errorbar(aes(ymin=meancopy-secopy, ymax=meancopy+secopy)) +
  geom_point(alpha=.5) +
  facet_grid(~scaf, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") + 
  xlab("Genomic Position (Mb)") +
  ylab("Mean Copy\nNumber") +
  geom_vline(xintercept = 13152174/1000000)

#write.table(results_sig, "floweringtime_copynumerGWAS.bed",quote=F, row.names = F, col.names = F)
results_sig<-fread("~/floweringtime_copynumerGWAS.bed")


###
FTGWA<-fread("~/commongarden_finalfilt_FTmean.assoc.txt")
FTGWA$fdr<-p.adjust(FTGWA$p_lrt, method = "fdr") 
results_fdrsig<-FTGWA[FTGWA$fdr < 0.05,]
nrow(results_fdrsig)

FTGWA$bon<-p.adjust(FTGWA$p_lrt, method = "bonferroni") 
results_bonsig<-FTGWA[FTGWA$bon < 0.05,]
nrow(results_fdrsig)

FTGWA$fdrten<-p.adjust(FTGWA$p_lrt, method = "fdr") 
results_fdrsig_ten<-FTGWA[FTGWA$fdrten < 0.1,]
nrow(results_fdrsig_ten)

FTGWA$chr <- gsub("Scaffold_","",FTGWA$chr)
FTGWA$chr <- as.numeric(FTGWA$chr)

###Figure 4B, top###
ggplot(data=FTGWA, aes(ps/1000000, -log10(p_lrt))) +
  geom_point(alpha=.1) +
  facet_grid(~chr, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8),
        panel.background = element_blank(),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16)) + 
  labs(x="Genomic Position (Mb)",y="-log10(pval)") +
 # geom_hline(yintercept = -log10(max(results_fdrsig$p_lrt)),lty="dashed",color="#52514f") +
  geom_hline(yintercept = -log10(max(results_bonsig$p_lrt)),lty="dashed",color="#828282") +
  geom_hline(yintercept = -log10(max(results_fdrsig_ten$p_lrt)),lty="dashed",color="black")
####

which(cov_bygene$start == 13152174) #for ATP locus
lm1<-lm( phenos_ordered$FTmeans ~ cov_bygene_scaled_tr[,which(cov_bygene$start == 13152174)])
summary(lm1)

phenos_ordered$ATP_cov<-cov_bygene_scaled_tr[,which(cov_bygene$start == 13152174)]

#which(cov_bygene$start == 21395774) #for FTD
#FTD<-cov_bygene_scaled_tr[,which(cov_bygene$start == 21395774)]
#FTD_df<-data.frame(FTD=FTD, ID=phenos_ordered$Sample)
#lm1<-lm( phenos_ordered$FTmeans ~ FTD)
#summary(lm1,type=3)


###Figure 4A, right side###
plot(phenos_ordered$FTmeans ~ cov_bygene_scaled_tr[,which(cov_bygene$start == 13152174)])
phenos_ordered$FThit_scaled<-cov_bygene_scaled_tr[,which(cov_bygene$start == 13152174)]
ggplot(data=phenos_ordered,
       aes(FThit_scaled,FTmeans)) +
  geom_point(alpha=.8) +
  geom_smooth(method="lm",color="black") +
  theme_bw() +
  labs(y="Family Mean\nFlowering Time", 
       x="Copy Number at\nATP Synthesis Locus") +
  geom_text(aes(label="r2 = 0.20", x=12.5,y=33))

