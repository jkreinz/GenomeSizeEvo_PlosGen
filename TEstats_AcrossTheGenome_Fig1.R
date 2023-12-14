#Scripts for investigation TE content across the genome, along with recombination rate and gene density
#Julia Kreiner, December 2023

library(data.table)
library(reshape2)
library(ggplot2)
library(MetBrewer)

setwd("~/Documents/UBC_Postdoc/Papers/Solomiya Paper/data/")
beds<-read.table("AMATA.genes_DNA_Helitron_COPIA_GYPSY_window_densities_1MB_CDS.bed") #windowed density of repeat content, 1MB
head(beds)
names(beds)<-c("scaf","start","end","dna","helitron","copia","gypsy","CDS")

#convert to long formatt for plotting
bed_long<-melt(beds, id.vars=c("scaf","start","end","CDS") )
teip<-met.brewer("Tiepolo",18)

#look within directory for files with stats on length of TE annotations from bed, and distances from nearest gene
myFiles <- list.files(pattern=".tmp")
myFiles
ldf <- list() # creates a list

for (k in 1:length(myFiles)){
  ldf[[k]] <- read.table(myFiles[k])
}

class<-list()
meanlength<-list()
medianlength<-list()
selength<-list()
standard_error <- function(x) sd(x) / sqrt(length(x)) # Create own function

for (k in 1:length(myFiles)){
  class[[k]]<-gsub(".sorted.bed.tmp","",x=myFiles[k])
  meanlength[[k]]<-mean(ldf[[k]]$V3-ldf[[k]]$V2)
  medianlength[[k]]<-median(ldf[[k]]$V3-ldf[[k]]$V2)
  selength[[k]]<-standard_error(ldf[[k]]$V3-ldf[[k]]$V2)
}

repeat_lengths<-data.frame(unlist(class),unlist(meanlength),unlist(medianlength),unlist(selength))
repeat_lengths<-repeat_lengths[-1,]
repeat_lengths$unlist.class.[repeat_lengths$unlist.class.=="Gypsy"]<-"Ty3"
repeat_lengths$unlist.class.[repeat_lengths$unlist.class.=="DTA"]<-"hAT"
repeat_lengths$unlist.class.[repeat_lengths$unlist.class.=="DTC"]<-"CACTA"
repeat_lengths$unlist.class.[repeat_lengths$unlist.class.=="DTH"]<-"Pif/Harbinger"
repeat_lengths$unlist.class.[repeat_lengths$unlist.class.=="DTM"]<-"Mutator"
repeat_lengths$unlist.class.[repeat_lengths$unlist.class.=="DTT"]<-"Tc1/Mariner"
repeat_lengths$unlist.class.[repeat_lengths$unlist.class.=="LTRunknown"]<-"Unknown LTR"

repeat_lengths$unlist.class. <- factor(repeat_lengths$unlist.class.,
                             levels = c("Helitron","hAT","Tc1/Mariner",
                                        "Pif/Harbinger","Mutator","CACTA",
                                        "Ty3","Copia","Unknown LTR","LINE","gene","LTR"))
repeat_lengths$unlist.class.
#counts of TE annotations by variable/class (taken from tail +n2 | wc -l of bed files )
repeat_lengths$number<-c(98935,NA,76336, 47548,9175,63730, 168077, NA, 106574,182443,1841,314717, 109205 )

library(tidyverse)
#first row of Fig 1B
p1A<-repeat_lengths %>% filter(unlist.class.!="gene" & unlist.class.!="LTR") %>%
  ggplot(aes(unlist.class., number, color=unlist.class.)) +
  geom_point() + 
  #geom_errorbar(aes(ymin=unlist.meanlength.-unlist.selength., 
  #                  ymax=unlist.meanlength.+unlist.selength.)) +
  labs(y="Number",x="TE Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c(teip[1], teip[6:10],teip[16:19],"grey40"))

library(tidyverse)
#second row of Fig 1B
repeat_lengths %>% filter(unlist.class.!="gene" & unlist.class.!="LTR")  %>% dplyr::summarise(total=sum(number))
p1<-repeat_lengths %>% filter(unlist.class.!="gene" & unlist.class.!="LTR") %>%
  ggplot(aes(unlist.class., unlist.meanlength., color=unlist.class.)) +
  geom_point() + 
  geom_errorbar(aes(ymin=unlist.meanlength.-unlist.selength., 
                ymax=unlist.meanlength.+unlist.selength.)) +
  labs(y="Size (bp)",x="TE Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c(teip[1], teip[6:10],teip[16:19],"grey40"))


#distance to gene (third row of Fig 1B)
myFiles <- list.files(pattern=".distance")
myFiles
longdf<-data.frame()

ldf2 <- list() 
for (k in 1:length(myFiles)){
  ldf2[[k]] <- read.table(myFiles[k],)
  names(ldf2[[k]])<-c("Scaffold","testart","teend","scaf2","genestart",
                      "gendend","distance")
  ldf2[[k]]$class<-rep(gsub(".distance","",myFiles[k]),nrow(ldf2[[k]]))
  longdf<-rbind(longdf,ldf2[[k]])
}

#just to visualize
ggplot(data=longdf, aes(distance, color=class)) +
  geom_density() +
  xlim(0,5000) +
  theme_bw()

class<-list()
meandistance<-list()
sedist<-list()
standard_error <- function(x) sd(x) / sqrt(length(x))

#do some stats on TE distance from gene by class
for (k in 1:length(myFiles)){
  class[[k]]<-gsub(".distance","",x=myFiles[k])
  meandistance[[k]]<-mean(ldf2[[k]]$distance)
  sedist[[k]]<-standard_error(ldf2[[k]]$distance)
}

repeat_dist<-data.frame(unlist(class),unlist(meandistance),unlist(sedist))
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="Gypsy"]<-"Ty3"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTA"]<-"hAT"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTC"]<-"CACTA"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTH"]<-"Pif/Harbinger"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTM"]<-"Mutator"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTT"]<-"Tc1/Mariner"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="LTRunknown"]<-"Unknown LTR"

repeat_dist$unlist.class. <- factor(repeat_dist$unlist.class.,
                                       levels = c("Helitron","hAT","Tc1/Mariner",
                                                  "Pif/Harbinger","Mutator","CACTA",
                                                  "Ty3","Copia","Unknown LTR","LTR","LINE"))

p2<-repeat_dist %>% filter(unlist.class.!="gene" & unlist.class.!="LTR")  %>%
  ggplot(aes(unlist.class., unlist.meandistance., color=unlist.class.)) +
  geom_point() + 
  geom_errorbar(aes(ymin=unlist.meandistance.-unlist.sedist., 
                    ymax=unlist.meandistance.+unlist.sedist.)) +
  labs(y="Distance to\nnearest gene (bp)",x="TE Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c(teip[1], teip[6:10],teip[16:19],"grey40"))

library(cowplot)
#this is part of figure 1B, need more stats below
plot_grid(p1A,p1,p2,align = "hv",nrow=3)


###########
#visual proportions of the genome by repeat class in 1Mb windows
prop<-read.table("allrepeatclasses_1Mb_prop.txt")
repeat_dist<-data.frame(unlist(class),unlist(meandistance),unlist(sedist))
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="Gypsy"]<-"Ty3"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTA"]<-"hAT"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTC"]<-"CACTA"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTH"]<-"Pif/Harbinger"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTM"]<-"Mutator"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTT"]<-"Tc1/Mariner"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="LTRunknown"]<-"Unknown LTR"
names(prop)<-c("scaf","start","end",repeat_dist$unlist.class[c(1:9,11)],"CDS")
colMeans(prop[4:14])
prop$`Unknown LTR`<-prop$LTR - (prop$Copia + prop$Ty3)
LINEp<-read.table("LINE.prop1MB")
head(prop)
prop$LINE<-LINEp$V5
rnaprop<-read.table("AMATA_rdna_output.prop100kb")
#prop$rRNA<-rnaprop$V5

teip<-met.brewer("Tiepolo",18)

prop_long<-melt(prop, id.vars = c("scaf","start","end"))

levels(prop_long$variable)
prop_long$variable <- factor(prop_long$variable,
                             levels = c("Helitron","hAT","Tc1/Mariner",
                                        "Pif/Harbinger","Mutator","CACTA",
                                        "Ty3","Copia","Unknown LTR","CDS","LTR","DNA",
                                        "LINE","Unknown"))

library(patchwork)
prop_long$scaf<-gsub("Scaffold_","",prop_long$scaf)
prop_long$scaf<-as.numeric(prop_long$scaf)

###part of Figure 1C ###
allscaf_TEs <-  prop_long %>% #filter(scaf=="Scaffold_1") %>%
  filter(variable!="LTR" & variable!="DNA") %>%
  ggplot(aes(x = start, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~scaf,scales="free_x",space = "free_x"
  ) +
  labs(y= "Fraction of 1MB window", x = "") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #   axis.line.x = element_blank(),
        panel.spacing = unit(0.15, "lines")) +
  scale_fill_manual(values=c(teip[1], teip[6:10],teip[16:18], "black","grey"))


#get total percent of the genome composed by each repeat class
prop<-read.table("allrepeatclasses_100kb_prop.txt")
repeat_dist<-data.frame(unlist(class),unlist(meandistance),unlist(sedist))
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="Gypsy"]<-"Ty3"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTA"]<-"hAT"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTC"]<-"CACTA"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTH"]<-"Pif/Harbinger"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTM"]<-"Mutator"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="DTT"]<-"Tc1/Mariner"
repeat_dist$unlist.class.[repeat_dist$unlist.class.=="LTRunknown"]<-"Unknown LTR"
names(prop)<-c("scaf","start","end",repeat_dist$unlist.class[c(1:9,11)],"CDS")
colMeans(prop[4:14])
prop$`Unknown LTR`<-prop$LTR - (prop$Copia + prop$Ty3)
Linep100k<-read.table("LINE.prop")
prop$LINE<-Linep100k$V5
prop$scaf<-gsub("Scaffold_","",prop$scaf)
prop$scaf<-as.numeric(prop$scaf)

byfam<-prop_long %>% group_by(variable) %>% summarise(mean_value=round(mean(value)*100,digits = 1))

byfam<-rbind(byfam, c("Unknown",2.4))
byfam<-rbind(byfam, c("LINE",0.1))
byfam$mean_value<-as.numeric(byfam$mean_value)

###Figure 1A ###
view(byfam)

byfam %>% filter(variable != "LTR" & variable != "DNA") %>% 
  ggplot(aes(x = mean_value, y = variable, fill=variable)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label=paste(mean_value,"%",sep = "")), nudge_x = 1.2, size=3,color="grey40") +
  theme_classic() +
  scale_fill_manual(values=c(teip[1], teip[6:10],teip[16:18], "black","grey40","grey60")) +
  labs(x="Percentage of Genome",y="Superfamily") 
###
#Investigations of TE content by recombination rate and gene density
#note that intially I looked at diversity and differentiation as well, not included in manuscript
#############

allchr_recomb<-fread("allscafs_recombrate.txt")
head(allchr_recomb)
allchr_recomb$V1<-gsub("tuberculatus_final_lcl\\|","",allchr_recomb$V1)
names(allchr_recomb)<-c("scaf","Pos","GenPos")

head(prop)
recomb_bywin<- 
  allchr_recomb %>% 
  mutate(win=floor(Pos/100000)) %>%
  group_by(scaf, win) %>%
  summarise(mean_genpos=mean(GenPos)) %>%
  mutate(start=win*100000, end=((win+1)*100000))


pi<-fread("~/4fold_commongarden_invar_minfilt_wmultiallelic_pi.txt")
fst<-fread("~/4fold_commongarden_invar_minfilt_wmultiallelic_fst.txt")

head(pi)
pi_bywin<-pi %>% 
  mutate(win=floor(window_pos_1/100000)) %>%
  group_by(chromosome, win) %>%
  summarise(mean_pi=mean(avg_pi,na.rm=T)) %>%
  mutate(start=win*100000, end=((win+1)*100000))
names(pi_bywin)[1]<-c("scaf")
prop$scaf<-paste("Scaffold_",prop$scaf,sep = "")
prop_pi<-inner_join(prop,pi_bywin,by=c("scaf","start","end"))  

prop_recomb<-inner_join(prop,recomb_bywin,by=c("scaf","start","end"))  
names(prop_recomb)
library(reshape2)
prop_recomb_long<-melt(prop_recomb, id.vars=c("win","mean_genpos","CDS","scaf","start","end"))
library(ggrepel)
#install.packages("directlabels")
library(directlabels)

prop_recomb_long$variable <- factor(prop_recomb_long$variable,
                             levels = c("Helitron","hAT","Tc1/Mariner",
                                        "Pif/Harbinger","Mutator","CACTA",
                                        "Ty3","Copia","Unknown LTR","CDS","LTR","DNA",
                                        "LINE","Unknown"))

###Sup Figure 1###
prop_recomb_long$value<-as.numeric(as.character(prop_recomb_long$value))
prop_recomb_long %>% filter(variable != "DNA") %>% filter(variable != "LTR") %>%
  ggplot(aes(value, CDS,color=variable)) +
  geom_point(alpha=.05) + 
  geom_smooth(method="lm") +
  #xlim(0,0.25) +
  #geom_dl(label=prop_recomb_long$variable, method="maxvar.points", inherit.aes=T) +
  theme_bw() +
  labs(x="Proportion of 100kb Window", color="TE superfamily") + 
  scale_color_manual(values=c(teip[1], teip[6:10],teip[16:18], "black","grey40","grey60"))
###

##Sup Figure 2###
prop_recomb_long %>% filter(variable != "DNA") %>% filter(variable != "LTR") %>%
  ggplot(aes(value, mean_genpos,color=variable)) +
  geom_point(alpha=.05) + 
  geom_smooth(method="lm") +
  #xlim(0,0.25) +
  #geom_dl(label=prop_recomb_long$variable, method="maxvar.points", inherit.aes=T) +
  theme_bw() +
  labs(y="Mean 4NeR",x="Proportion of 100kb Window", color="TE superfamily") + 
  scale_color_manual(values=c(teip[1], teip[6:10],teip[16:18], "black","grey40","grey60"))
###


prop_pi$scaf<-gsub("Scaffold_","",prop_pi$scaf)
prop_pi$scaf<-as.numeric(prop_pi$scaf)
recomb_bywin$scaf<-gsub("Scaffold_","",recomb_bywin$scaf)
recomb_bywin$scaf<-as.numeric(recomb_bywin$scaf)
recomb_pi_bywin<-inner_join(recomb_bywin, prop_pi, by=c("scaf","win"))

############


recomb_pi_bywin$allrepeats<-recomb_pi_bywin$Copia + recomb_pi_bywin$DNA + recomb_pi_bywin$hAT +
  recomb_pi_bywin$`Pif/Harbinger` + recomb_pi_bywin$Mutator + recomb_pi_bywin$`Tc1/Mariner` +
  recomb_pi_bywin$`Ty3` + recomb_pi_bywin$Helitron + recomb_pi_bywin$LTR+
  recomb_pi_bywin$`Unknown LTR` 

#LD based recombination rate across the genome (calculated in Kreiner et al., PNAS 2019)
recomb_allscaf<-recomb_bywin %>% #filter(scaf=="Scaffold_1") %>% 
  ggplot(aes(win, mean_genpos)) +
  #geom_point(alpha=.2) + 
  geom_smooth(method="loess",span=.2,color="black",se = F,size=.75) +
  facet_grid(~scaf, scales = "free_x", space = "free_x")+
  #xlim(0,30) +
  labs(y="Mean 4Ner", x="Chromsome 1 Position") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.15, "lines"))

#pi_bywin$scaf<-gsub("Scaffold_","",pi_bywin$scaf)
#pi_bywin$scaf<-as.numeric(pi_bywin$scaf)
#
#pi_allscaf<-pi_bywin %>% #filter(scaf=="Scaffold_1") %>% 
#  ggplot(aes(win, mean_pi)) +
#  #geom_point(alpha=.2) + 
#  geom_smooth(method="loess",span=.2,color="black",se = F,cex=.75) +
#  facet_grid(~scaf, scales = "free_x", space = "free_x")+
#  #xlim(0,30) +
#  labs(y="Mean 4fold pi", x="") +
#  theme_classic() +
#  theme(strip.background = element_blank(),
#        strip.text.x = element_blank(),
#        panel.spacing = unit(0.15, "lines"),
#        axis.title.x =element_blank(),
#        axis.ticks.x = element_blank(),
#        axis.text.x = element_blank())


names(rnaprop)<-c("scaf","start","end","x","prop")
rnaprop$scaf<-gsub("Scaffold_","",rnaprop$scaf)
rnaprop$scaf<-as.numeric(rnaprop$scaf)

rrna_allscafs<-rnaprop %>% 
  ggplot(aes(start, prop)) +
  geom_point(alpha=.7) + 
  geom_path(alpha=.5) +
  # geom_smooth(method="loess",span=.color="black",se = F,cex=.75) +
  facet_grid(~scaf, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        axis.title.x =element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("rRNA prop")


####Figure 1C!!!###
allscaf_TEs + plot_spacer() + rrna_allscafs +
  plot_spacer() + recomb_allscaf +
  plot_layout(
    ncol = 2, 
    nrow = 3, 
    widths = c(6, 0,0),
    heights = c(6, 1,1)
  ) 

####

#OK finally get the rest of figure 1B (correlations of repeat content with recombination rate and gene density)


Helitron<-lm(data=recomb_pi_bywin, Helitron ~ mean_genpos + CDS  + scaf )
hAT<-lm(data=recomb_pi_bywin, hAT ~ mean_genpos + CDS  + scaf )
TM<-lm(data=recomb_pi_bywin, `Tc1/Mariner` ~ mean_genpos + CDS  + scaf )
PH<-lm(data=recomb_pi_bywin, `Pif/Harbinger` ~ mean_genpos + CDS  + scaf )
Mutator<-lm(data=recomb_pi_bywin, Mutator ~ mean_genpos + CDS  + scaf )
CACTA<-lm(data=recomb_pi_bywin, CACTA ~ mean_genpos + CDS  + scaf )
Ty3<-lm(data=recomb_pi_bywin, Ty3 ~ mean_genpos + CDS  + scaf )
Copia<-lm(data=recomb_pi_bywin, Copia ~ mean_genpos + CDS  + scaf )
UnknownLTR <-lm(data=recomb_pi_bywin, `Unknown LTR` ~ mean_genpos + CDS  + scaf )
LINE <-lm(data=recomb_pi_bywin, LINE ~ mean_genpos + CDS  + scaf )


CDS_slopes<-data.frame(rbind(summary(Helitron)$coefficients[3,1:2],
                             summary(hAT)$coefficients[3,1:2],
                             summary(TM)$coefficients[3,1:2],
                             summary(PH)$coefficients[3,1:2],
                             summary(Mutator)$coefficients[3,1:2],
                             summary(CACTA)$coefficients[3,1:2],
                             summary(Ty3)$coefficients[3,1:2],
                             summary(Copia)$coefficients[3,1:2],
                             summary(UnknownLTR)$coefficients[3,1:2],
                             summary(LINE)$coefficients[3,1:2]
))

genpos_slopes<-data.frame(rbind(summary(Helitron)$coefficients[2,1:2],
                                summary(hAT)$coefficients[2,1:2],
                                summary(TM)$coefficients[2,1:2],
                                summary(PH)$coefficients[2,1:2],
                                summary(Mutator)$coefficients[2,1:2],
                                summary(CACTA)$coefficients[2,1:2],
                                summary(Ty3)$coefficients[2,1:2],
                                summary(Copia)$coefficients[2,1:2],
                                summary(UnknownLTR)$coefficients[2,1:2],
                                summary(LINE)$coefficients[2,1:2]
))

CDS_slopes$class<-c("Helitron","hAT","Tc1/Mariner","Pif/Harbinger","Mutator","CACTA",
                    "Ty3","Copia","Unknown LTR","LINE")
genpos_slopes$class<-c("Helitron","hAT","Tc1/Mariner","Pif/Harbinger","Mutator","CACTA",
                       "Ty3","Copia","Unknown LTR","LINE")

genpos_slopes$class <- factor(genpos_slopes$class,
                              levels = c("Helitron","hAT","Tc1/Mariner",
                                         "Pif/Harbinger","Mutator","CACTA",
                                         "Ty3","Copia","Unknown LTR","LINE"))
CDS_slopes$class <- factor(CDS_slopes$class,
                           levels = c("Helitron","hAT","Tc1/Mariner",
                                      "Pif/Harbinger","Mutator","CACTA",
                                      "Ty3","Copia","Unknown LTR","LINE"))


gen_p<-genpos_slopes %>% 
  ggplot(aes(class, Estimate, color=class)) +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_point() + 
  geom_errorbar(aes(ymin=Estimate-Std..Error, 
                    ymax=Estimate+Std..Error)) +
  labs(y="Effect of Recombination Rate",x="TE Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c(teip[1], teip[6:10],teip[16:19],"grey30"))


cds_p<-CDS_slopes %>% 
  ggplot(aes(class, Estimate, color=class)) +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_point() + 
  geom_errorbar(aes(ymin=Estimate-Std..Error, 
                    ymax=Estimate+Std..Error)) +
  labs(y="Effect of CDS",x="TE Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c(teip[1], teip[6:10],teip[16:19],"grey30"))

###Figure 1C###
library(cowplot)
plot_grid(p1A,p1,p2,cds_p, gen_p,align = "hv",nrow=5)
###
  