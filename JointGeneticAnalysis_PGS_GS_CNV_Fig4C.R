#Script for estimating PGS and comparative analysis of genetic features underlying Flowering Time
#Julia Kreiner, Dec 2023

library(data.table)
library(dplyr)
library(ggplot2)
library(r2glmm)

setwd("~/Documents/UBC_Postdoc/Papers/Solomiya Paper/data/")
  
library(data.table)
  PH_gwas<-fread("FTnocovs_FDR10p_sigsnps.txt")
  #PH_gwas<-fread("~/commongarden_finalfilt_FTmean_recheck_nokinship.sigFDR10")
  names_gwas<-fread("FTnocovs_FDR10p_sigsnps.txt",nrows = 1,header = T)
  #names(PH_gwas)<-c("chr","rs","ps","n_mis","n_obs","allele1","allele0","af","beta","se","p_wald")
  
  #need genotype matrix here to get PRS
  PH_genos<-read.table("FTnocovs_FDR10p_sigsnps.raw",header=T, sep=" ") #genotype matrix of loci signficant for its association with FT past the FDR10% threshold
 # PH_genos<-read.table("~/commongarden_finalfilt_FTmean_recheck_nokinship_sigFDR10.raw",header=T, sep=" ")
  names(PH_genos)
  PH_genos_t<-t(PH_genos[7:ncol(PH_genos)])
  colnames(PH_genos_t)<-PH_genos$IID
  PH_genos_t<-as.data.frame(PH_genos_t)
  
  #loop for calculating polygenic score, first at each signficant locus, then summing in a later step
  ind_list<-list()
  snp_list<-list()
  for (i in 1:176) {
    for (k in 1:(ncol(PH_genos)-6)) {
      snp_list[[k]]<-2*(PH_gwas$beta[k])*PH_genos_t[k,i]
    } 
    ind_list[[i]]<-sum(unlist(snp_list),na.rm=T)
  }
  
  PH_PRS<-data.frame(ID=PH_genos$IID,PRS=unlist(ind_list))
  hist(PH_PRS$PRS)
  
  #dataframe containing copy number at NADH gene, along w depth based genome size estimate
  mNADH<-read.table("~/NADH_forrelatedness_BP.txt",header=T)
  names(mNADH)

  #merge PRS with predictor dataset
  mNADH$ID<-as.integer(mNADH$ID)
  test<-inner_join(mNADH, PH_PRS,by="ID")
  names(test)
  #rescale predictors as Z-values
  test$PRS_scale<-scale(test$PRS)
  test$cg_scale<-scale(test$copygeno) #NADH locus copy number
  test$ID<-as.integer(test$ID)
  #merge with all metadata
  merged_wFlower_heteigen<-read.table("merged_wFlower_heteigen.txt",header=T)
  test2<-inner_join(test, merged_wFlower_heteigen, by="ID")
    hist(test2$PRS_scale)
  hist(test2$FTmeans.x)
  test2$genomesize_scale<-(scale(test2$GS_bpdepth))
  
  #note the linear fit isn't great here, do logistic below
  ggplot(test, aes(PRS_scale, FTmeans)) +
    geom_point() +
    geom_smooth(method="lm") +
   # labs(x="mean centered PGRS for FT",y="var. rudis Ancestry") +
    theme_bw() 

###############################
#joint analysis of genetic predictors of flowering time
library(lme4qtl) #this is the program I found to incorporate the relatedness matrix
library(lme4)
library(car)
    
relmat<-read.table("commongarden_mergedSNPS_missing2_justCGinds.cXX.txt")
relmat_ID<-read.table("commongarden_mergedSNPS_missing2_justCGinds.fam")
row.names(relmat)<-colnames(relmat)<-relmat_ID$V2
class(relmat)
relmat<-as.matrix(relmat)

qqnorm(log(test2$PRS)) #looks great!
min(test2$PRS) #some negative values, add to each before log transformation

test2$PRS_adj<- scale(log(test2$PRS+67.28941)) #z transformation of log transformed values (scaled to be positive by adding the min value)

#key result!
PH_PRS_model <- relmatLmer(FTmeans.x ~  genomesize_scale  + PRS_adj + cg_scale + 
                               + (1|ID) , test2, relmat = list(ID=relmat))
Anova(PH_PRS_model,type=3)
library(r2glmm)
r.squaredGLMM(PH_PRS_model)
  
#Sup Figure 5
library(viridis)
library(heatmap3)

summary(lm(data=test2, FTmeans.x~PRS_adj))

ATP_model <- relmatLmer(copygeno ~   Long.x + Lat.x + Sex.x.x * Env.x.x + fam_K1.x 
                        + (1|ID)  , test2, relmat = list(ID=relmat))
Anova(ATP_model,type=3)

test2$PRS_l<- scale(log(test2$PRS+67.29))

PGRS_model <- relmatLmer(log(PRS+67.29) ~   Long.x + Lat.x + Sex.x.x * Env.x.x + fam_K1.x 
                         + (1|ID) , test2, relmat = list(ID=relmat))
Anova(PGRS_model,type=3)

###Figure 4B###
ggplot(data=test2, aes(PRS_adj,FTmeans.x)) +
  geom_point() +
  lims(x=c(-1.5,3)) +
  #xlim(-1,0.5) +
  geom_smooth(method="lm",color="black") +
  theme_bw() +
  labs(y="Flowering Time\nFamily Mean", x="PVGFT")
####
library(r2glmm)
library(lsmeans)
gs_lsmeans<-lsmip(PH_PRS_model, ~ genomesize_scale ,
                   at=list(genomesize_scale=c(-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)),
                   type="response",
                   plotit=F)

 
PRS_lsmeans<-lsmip(PH_PRS_model, ~  PRS_adj,
                  at=list(PRS_adj=c(-2,-1.5,-1,-.5,0,0.5,1,1.5,2,2.5)),
                  type="response",
                  plotit=F)

copygeno_lsmeans<-lsmip(PH_PRS_model, ~ cg_scale ,
                       at=list(cg_scale=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3)),
                       type="response",
                       plotit=F)



FT_long<-melt(test2, id.vars=c("FTmeans.x","PHmeans.x","Lat.x","ID"),
             measure.vars =c("genomesize_scale","cg_scale","PRS_adj") )
# FT_long$variable <- factor(FT_long$variable, rev(levels(FT_long$variable)[c(1,5,2,8,3,6,4,7,8,9)]))

levels(FT_long$variable)

purp<-met.brewer("VanGogh1",6)

###Figure 4C!!###  
    ggplot() +
      geom_point(data=FT_long, aes(as.numeric(value), FTmeans.x, group=variable, color=variable),alpha=.4,size=2) +
      #geom_smooth(alpha=.2, method="lm") +
      #geom_ribbon(aes(x=x_PRS$PRS, ymin=x_PRS$lower, ymax=x_PRS$upper), alpha= 0.3, fill="blue") +
      labs(x="Scaled Predictor Value", y="Flowering Time Family Mean Value") +
      theme_bw() +
      xlim(-3,3.5)  + ylim(25,80) +
      theme(legend.position = "right") +
      scale_color_manual(labels=c("Genome Size","NADH Cluster Copy Number","PRS"),
                         values=c("grey60", "forestgreen","#7dba14")) +
      geom_ribbon(data=copygeno_lsmeans, aes(y=yvar,x=cg_scale,ymin=yvar-SE,ymax=yvar+SE), fill="forestgreen",alpha=.9) +
      geom_ribbon(data=gs_lsmeans, aes(y=yvar,x=genomesize_scale,ymin=yvar-SE,ymax=yvar+SE), fill="grey60",alpha=.9) +
      geom_ribbon(data=PRS_lsmeans, aes(y=yvar,x=PRS_adj,ymin=yvar-SE,ymax=yvar+SE), fill="#7dba14",alpha=.9) 
    
    FT_fixcor<-summary(PH_PRS_model)
    fit_cov<-as.matrix( FT_fixcor$vcov[c("PRS_adj","genomesize_scale","cg_scale"),
                                       c("PRS_adj","genomesize_scale","cg_scale")])
    #names(fit_cov)<-c("Low complexity repeats","Genome size","PRS","Genome Size")
    #row.names(fit_cov)<-c("Low complexity repeats","Genome size","PRS","Genome Size")
    
    #heatmap(fit_cov,  symm = F,)
    library(viridis)
    greys<-RColorBrewer::brewer.pal(n = 9, "Greys")
    heatmap3::heatmap3(x = cor(test2[,c(125,253,251)]),scale = "none",margins = c(10,10),cexRow =1, cexCol = 1,symm = T,
                       #labCol = c("ATP Synthase CN","Low complexity repeats","Latitude","Genome size","PRS"),
                       #labCol = c("PRS","Genome size","Latitude","Low complexity repeats","ATP Synthase CN"),
                       # labRow = c("PRS","Genome size","Latitude","Low complexity repeats","ATP Synthase CN"),
                       col=viridis(n=100))
    
    cor(test2[,c(125,253,251)])
###
