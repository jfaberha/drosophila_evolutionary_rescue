##########################################################################################
### This script is associated with a study on the evolution of spinosad-resistance in  ###
### field populations of Drosophila melanogaster run in 2023 by the Rudman lab (WSUV). ###
### This particular script covers a portion of the bioinformatic analysis that examines###
### two general questions; (1) is there enrichment of an a priori list of insecticide  ###
### resistance-associated genes in the various lists of outlier genes in our study,    ###
### determined through GLM analysis on all SNPs, and (2) is there a pleiotropic cost to###
### empirically predicted spinosyn adaptation? 										   ###																	   
###																					   ###
### The first enrichment question uses two approaches; 1) hypergeometric overlap tests ###
### and 2) Wilcoxon Rank Sum analysis. The first approach simply tests whether the a   ###
### priori candidate genes are overrepresented in lists of genes containing outlier    ###
### SNPs from various GLM contrasts. The second approach is a non-parametric analysis  ###
### that tests whether the significance-ranks of these candidate genes are higher than ###
### random matched sets of genes across the various contrasts. 						   ###
###																					   ###
### As an add-on, we also wanted to test whether genes containing outlier SNPs in      ###
### several groups of spinosad-exposed populations (vs. controls) were enriched in     ###
### other contrasts in our experiment, using the same methods. These analyses seek to  ###
### answer the question of whether we could detect genomic convergence between flies   ###
### that were exposed to spinosad at the start of the experiment (SE.T1 & SP.T4) to    ###
### those that were previously adapted to the insecticide (PA).                        ###
###																					   ###
### The second cost question uses a variety of approaches as well. We identified lists ###
### of putative spinosyn-adaptive loci using either the PAvE or PAvE.T1 contrasts      ###
### (mostly the latter), and queried whether allele frequencies in temporal contrasts  ###
### (responding to the environment) we changing in the direction of PA-biased or       ###
### E-biased allele frequencies. If changing in the PA-biased direction, it would seem ###
### spinosyn-adaptive alleles have no detectable environmental costs, and vice versa.  ###
### We use binomial sign tests and t-tests to look for E-biased changes.			   ###
###																					   ###
### Next we attempted to answer this question using hypergeometric overlap tests. Here ###
### we asked whether PAvE or PAvE.T1 outlier SNPs are overrepresented within temporal  ###
### contrast lists in either the PA-biased or E-biased direction. Because temporal     ###
### contrast GLMs failed to detect any substantial evolutionary signal, we used an     ###
###  alternate method to define lists based on strictly parallel AF differences across ###
### all or nearly all cages in a particular treatment/time-point.					   ###
##########################################################################################
 

### In R ###
#configure r environment
setwd("/scratch/user/jfaberha/20260902_101357/admera/gp_analysis/rudflies_2023_redo/r")
library(emmeans)
library(matrixStats)
#install.packages("~/Downloads/ACER-master", repos=NULL, type="source")
library(ACER)
library(poolSeq)
library(ggplot2)
library(ggpubr)
library(stats)
library(ggfortify)
library(dplyr)
library(zoo)
library(rstatix)
#install.packages("slider")
library(slider)
library(stringr)
library(RColorBrewer)


#################################
### Load required input files ###
#################################

## Variant Effect Predictor (VEP) annotation file with appended FLYCADD scores
vep <- read.table("filtered-all.annot.vcf.FLYCADD.tsv", header=TRUE) 
## More detailed gene information from gff
snp.gff.overlap <- read.delim("rudflies_2023_hafpipe_loci_genic_overlap_info.bed", sep="\t", header=TRUE)
## Sample metadata table
haf.meta <- read.table("rudflies_2023_meta.tsv", header=TRUE)
## Hafpipe imputed allele frequency table
haf.freq <- read.delim("rudflies_2023_hafpipe.csv", header=TRUE, sep = ",")

## Do a bit of reformatting for the frequency table header
names(haf.freq) <- gsub("[.]af","",names(haf.freq))
names(haf.freq) <- gsub("X","",names(haf.freq))
names(haf.freq)[1]="CHROM"
names(haf.freq)[2]="POS"

## Filter for founder samples for PA, S, and E
haf.meta.filt <- haf.meta[haf.meta$batch == "a" & haf.meta$experiment == "spino" & (haf.meta$treat.fix == "PA" | haf.meta$treat.fix == "S" | haf.meta$treat.fix == "E"),]
## Use this next line if you want to run aCMH with balanced pops, since all pops represented in TPT3
#haf.meta.filt <- haf.meta.filt[which((haf.meta.filt$cage.fix %in% haf.meta.filt[haf.meta.filt$tpt==3,]$cage.fix)==TRUE),]
## Find matching sample names in the frequency table, post metadata filtering
haf.freq.filt <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt$samp)==TRUE)]
## And run the reciprocal filter
haf.meta.filt <- haf.meta.filt[which((haf.meta.filt$samp %in% names(haf.freq.filt))==TRUE),]

## additional filter to remove low variance loci
haf.sites.filt <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.filt))>0.001,c(1:2)])
haf.freq.filt <- na.omit(haf.freq.filt[rowVars(as.matrix(haf.freq.filt))>0.001,])

## tables need a bit of reformatting for some downstream plotting
haf.meta.filt$treat.fix <- as.factor(haf.meta.filt$treat.fix)
haf.meta.filt$tpt <- as.factor(haf.meta.filt$tpt)

######################################
### For downstream analysis create ###
### a TPT1-only metadata table     ###
######################################

## Create a new column, "condition", which contains the population outcomes for S samples
## For example; SE - extinct S populations, SP - persistent S populations
filt.table <- cbind(cage=c("3","7","11","15","21","27","33","37","41","45","2","6","10","14","20","24","31","36","40","44","1","5","9","13","19","23","29","35","39","43","1","5","9","13","19","23","29","35","39","43"),tpt=c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"),condition=c("SP","SP","SE","SP","SE","SE","SP","SP","SE","SE","PA","PA","PA","PA","PA","PA","PA","PA","PA","PA","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E","E"))
haf.meta.T1filt <- merge(haf.meta, filt.table, by=c("cage","tpt"))

haf.meta.T1filt <- unique(haf.meta.T1filt[haf.meta.T1filt$batch == "a" & (haf.meta.T1filt$treat == "S" | haf.meta.T1filt$treat == "PA" | haf.meta.T1filt$treat == "E"),])
haf.meta.T1filt <- haf.meta.T1filt[which((haf.meta.T1filt$samp %in% c(1:9,30:171))==TRUE),]
haf.freq.T1filt <- haf.freq[, which((names(haf.freq) %in% haf.meta.T1filt$samp)==TRUE)]

#additional filter to remove extreme low variance loci
haf.sites.T1filt <- haf.freq[rowVars(as.matrix(haf.freq.T1filt))>0.001,c(1:2)]
haf.freq.T1filt <- haf.freq.T1filt[rowVars(as.matrix(haf.freq.T1filt))>0.001,]

#tables need a bit of reformatting for glm to run
haf.meta.T1filt$treat <- as.factor(haf.meta.T1filt$treat)
haf.meta.T1filt$tpt <- as.factor(haf.meta.T1filt$tpt)
haf.meta.T1filt$condition <- as.factor(haf.meta.T1filt$condition)
haf.meta.T1filt <- haf.meta.T1filt[order(haf.meta.T1filt$samp),]


##########################################
### Create master table of GLM results ###
##########################################

### Start with PA and S contrasts ###
## Load PA vs S, all samples combined
contrast.PAvS.treat.all <- read.table("rudflies_2023_PAvS_treat.all.table.GLMcontrast.txt", header=TRUE)
contrast.PAvS.treat.all <- contrast.PAvS.treat.all[,-3] #remove AF mean column
names(contrast.PAvS.treat.all)[3] <- "PAvS" #label glm results column
## Load PA vs S for each timepoint
contrast.PAvS.treat <- read.table("rudflies_2023_PAvS_treat.GLMcontrast.txt", header=TRUE)
#names(contrast.PAvS.treat) <- c("CHROM","POS","PAvS.T1","PAvS.T2","PAvS.T3","PAvS.T4")
## Load timepoint contrasts using PA and S combined
contrast.PAvS.tpt <- read.table("rudflies_2023_PAvS_tpt.GLMcontrast.txt", header=TRUE)
#names(contrast.PAvS.tpt) <- c("CHROM","POS","PA.S.T1vT2","PA.S.T1.T3","PA.S.T1vT4","PA.S.T2vT3","PA.S.T2vT4","PA.S.T3vT4") #label glm results columns
## Load PA and S "timepoint:treatment" interaction
contrast.PAvS.int <- read.table("rudflies_2023_PAvS_int.GLMcontrast.txt", header=TRUE)
#names(contrast.PAvS.int) <- c("CHROM","POS","int.PAvS.T1vT2","int.PAvS.T1vT3","int.PAvS.T1vT4","int.PAvS.T2vT3","int.PAvS.T2vT4","int.PAvS.T3vT4") #label glm results columns
contrast.PAvS <- cbind(contrast.PAvS.treat, contrast.PAvS.tpt[,-c(1,2)], contrast.PAvS.int[,-c(1,2)])
contrast.PAvS <- merge(contrast.PAvS.treat.all, contrast.PAvS, by=c("CHROM","POS"))
contrast.PAvS <- contrast.PAvS[order(contrast.PAvS[,1], contrast.PAvS[,2]), ] #sort by locus
## Save all PA and S contrasts for posterity
#write.table(contrast.PAvS, file="rudflies_2023_redo.PAvS.multiGLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Look at PA-only, E-only, and S-only time-point contrast GLM results ##
# Load S-only results
contrast.tpt.S.table <- read.table("rudflies_2023_S_tpt.table.GLMcontrast.txt", header=TRUE)
#names(contrast.tpt.S.table) <- c("CHROM","POS","S.af.mean","S.T1vT2","S.T1vT3","S.T1vT4","S.T2vT3","S.T2vT4","S.T3vT4") #label glm results columns
# Load PA-only results
contrast.tpt.PA.table <- read.table("rudflies_2023_PA_tpt.table.GLMcontrast.txt", header=TRUE)
#names(contrast.tpt.PA.table) <- c("CHROM","POS","PA.af.mean","PA.T1vT2","PA.T1vT3","PA.T1vT4","PA.T2vT3","PA.T2vT4","PA.T3vT4") #label glm results columns
# Merge the two
# Load E-only results
contrast.tpt.E.table <- read.table("rudflies_2023_E_tpt.table.GLMcontrast.txt", header=TRUE)
names(contrast.tpt.E.table) <- c("CHROM","POS","E.af.mean","E.T1vT2","E.T1vT3","E.T1vT4","E.T2vT3","E.T2vT4","E.T3vT4") #label glm results columns
# Merge them
contrast.tpt.S_only.PA_only.E_only <- merge(merge(contrast.tpt.S.table, contrast.tpt.PA.table, by=c("CHROM","POS")),contrast.tpt.E.table, by=c("CHROM","POS"))
contrast.tpt.S_only.PA_only.E_only <- contrast.tpt.S_only.PA_only.E_only[,-c(3,10,17)] #remove AF mean cols
# Save relabeled tables for posterity
#write.table(contrast.tpt.S.table, file="rudflies_2023_S_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
#write.table(contrast.tpt.PA.table, file="rudflies_2023_PA_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
#write.table(contrast.tpt.E.table, file="rudflies_2023_E_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## TPT1-only GLM contrasts: PA vs S vs SE vs E - all pairwise ##
contrast.PAvSvSEvE.table <- read.table("rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", header=TRUE)
#names(contrast.PAvSvSEvE.table) <- c("CHROM","POS","EvPA.T1","EvSE.T1","EvSP.T1","PAvSE.T1","PAvSP.T1","SEvSP.T1") #label glm results columns
# Save relabeled table for posterity
# write.table(contrast.PAvSvSEvE.table, file="rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Load PA and E contrasts ##
# PA vs E, all time-points combined
contrastout.PAvE.table <- read.table("rudflies_2023_PAvE_treat.all.table.GLMcontrast.txt", header=TRUE)
# PA vs E, founders
contrastout.PAvE.founder.table <- read.table("rudflies_2023_PAvE_treat.table.GLMcontrast.founders.txt", header=TRUE)
# PA vs E, at individual time-points
contrastout.treat.PA.E.table <- read.table("rudflies_2023_PAvE_treat.table.GLMcontrast2.txt", header=TRUE)
# Combine these three tables
contrastout.PAvE.table <- merge(contrastout.PAvE.table[,-3], contrastout.PAvE.founder.table[,-3], by=c("CHROM","POS"))
contrastout.PAvE.table <- merge(contrastout.PAvE.table, contrastout.treat.PA.E.table[,-3], by=c("CHROM","POS"))
names(contrastout.PAvE.table) <- c("CHROM","POS","PAvE","PAvE.F","PAvE.T1","PAvE.T2","PAvE.T3","PAvE.T4") #label glm results columns
# Save relabeled table for posterity
#write.table(contrastout.PAvE.table, file="rudflies_2023_PAvE_treat.all.table.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Load SE and E contrasts ##
# S vs E, all time-points combined
contrastout.SvE.table <- read.table("rudflies_2023_SvE_treat.all.table.GLMcontrast.txt", header=TRUE)[,c(1,2,4)]
# S vs E, at individual time-points
contrastout.treat.S.E.table <- read.table("rudflies_2023_SvE_treat.GLMcontrast.table.txt", header=TRUE)
# Combine these two tables
contrastout.SvE.table <- merge(contrastout.SvE.table, contrastout.treat.S.E.table, by=c("CHROM","POS"))
# S and E: "treatment:time-point" interaction
contrastout.int.S.E.table <- read.table("rudflies_2023_SvE_int.GLMcontrast.table.txt", header=TRUE)
# Add these results
contrastout.SvE.table <- merge(contrastout.SvE.table, contrastout.int.S.E.table, by=c("CHROM","POS"))
names(contrastout.SvE.table) <- c("CHROM","POS","SvE","SvE.T1","SvE.T2","SvE.T3","SvE.T4","int.SvE.T1vT2","int.SvE.T1vT3","int.SvE.T1vT4","int.SvE.T2vT3","int.SvE.T2vT4","int.SvE.T3vT4") #label glm results columns
# Save full table for posterity
#write.table(contrastout.SvE.table, file="rudflies_2023_SvE_treat.all.table.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)


### Merge all glm results ###
glm.all <- merge(contrast.PAvS, contrast.tpt.S_only.PA_only.E_only, by=c("CHROM","POS"))
glm.all <- merge(glm.all, contrast.PAvSvSEvE.table, by=c("CHROM","POS"))
glm.all <- merge(glm.all, contrastout.PAvE.table, by=c("CHROM","POS"))
glm.all <- merge(glm.all, contrastout.SvE.table, by=c("CHROM","POS"))
glm.all <- glm.all[glm.all$CHROM != "4",] #remove chromsome 4
#glm.all <- glm.all[is.finite(rowSums(glm.all[,-c(1,2)])), ]
glm.all <- glm.all %>%
  		arrange(CHROM,POS) #sort by locus

## How many sites remain after merging all filtered GLM results tables?
dim(glm.all)
#[1] 1183949      60

## P-value correction
x=ncol(glm.all) #use number of columns
for(i in c(3:x)) { #start after loci columns and append FDR cols at the end of the table
  glm.all <- cbind(glm.all, p.adjust(glm.all[,i], method = "fdr"))
  colnames(glm.all)[i+(x-2)] <- paste(names(glm.all)[i], ".fdr", sep="")
}

## Double check number of cols after adding FDR cols
dim(glm.all)
#[1] 1183949     118

## Calculate -log10 for p-values
y=ncol(glm.all) #use number of columns
for(i in 3:x) { #start after loci columns and append logp cols at the end of the table
  #find minimum non-zero p-value and divide by 2 to reassign to zero values
  min.temp <- min(glm.all[glm.all[,i]!=0,i])/2 
  col.temp <- glm.all[,i]
  col.temp[col.temp==0] <- min.temp #reassign zero p-values for plotting only
  #perform -log10(p) calculations
  glm.all <- cbind(glm.all, as.numeric(-log10(glm.all[,i])))
  colnames(glm.all)[i+(y-2)] <- paste(names(glm.all)[i], ".logp", sep="")
}

## Only write table if needed, these files are huge!
#write.table(glm.all, file=paste("rudflies_2023_redo.glm.masterfile.txt", sep=""), quote = FALSE, row.names = F)

## Merge raw GLM results table (pre-rolling window average) with annotation info 
glm.all.annot <- merge(glm.all, vep, by=c("CHROM","POS"))

##################################################
### For downstream matched set analysis create ###
### a non-redundant SNP-consequence VEP table  ###
##################################################

## Calculate E and PA founder means for all loci. E was the founder for S treatments too.
## This can be used to filter for intermediate starting allele frequencies, since 
## significant allele frequencies at these loci are more likely to have population-level
## phenotypic consequences.
haf.freq.Fmean <- cbind(haf.freq[,c(1:2)], E=rowMeans(haf.freq[,c(3:5)]), PA=rowMeans(haf.freq[,c(6:8)]))

## Calculating mean allele frequencies for each treat by TPT and across TPTs
freq_means <- cbind(haf.sites.filt, 
  E.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"]),
  S.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"]),
  PA.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"]),
  E.T1.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"&haf.meta.filt$tpt=="1"]),
  E.T2.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"&haf.meta.filt$tpt=="2"]),
  E.T3.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"&haf.meta.filt$tpt=="3"]),
  E.T4.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"&haf.meta.filt$tpt=="4"]),
  S.T1.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"&haf.meta.filt$tpt=="1"]),
  S.T2.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"&haf.meta.filt$tpt=="2"]),
  S.T3.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"&haf.meta.filt$tpt=="3"]),
  S.T4.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="S"&haf.meta.filt$tpt=="4"]),
  PA.T1.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"&haf.meta.filt$tpt=="1"]),
  PA.T2.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"&haf.meta.filt$tpt=="2"]),
  PA.T3.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"&haf.meta.filt$tpt=="3"]),
  PA.T4.af=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"&haf.meta.filt$tpt=="4"]))

## Merge mean frequency tables with GLM results for filtering
glm.all.wMeans <- merge(glm.all,merge(haf.freq.Fmean,freq_means,by=c("CHROM","POS")),by=c("CHROM","POS"))

## Create a non-redundant locus list that has priority allele consequence impacts
vep_cons <- unique(cbind(glm.all.annot[,c(1,2)],Consequence=glm.all.annot$Consequence))
vep_score <- as.data.frame(cbind(Consequence=vep$Consequence,FLYCADD=vep$FLYCADD,Extra=vep$Extra))

## Fix "Extra" column syntax
vep_score$Extra <- sub(";.*", "", vep_score$Extra)
vep_score$Extra <- sub(".*=", "", vep_score$Extra)
## Make FlyCADD scores numeric
vep_score$FLYCADD <- as.numeric(vep_score$FLYCADD)

## Summarise unique SNP consequences by FlyCADD score, used to rank consequences
vep_score_mean <- vep_score %>%
  group_by(Consequence,Extra) %>% 
  summarise_all(mean)

## Reformat mean table  
vep_score_mean <- as.data.frame(vep_score_mean)
vep_score_mean$Impact <- vep_score_mean$Extra

## Here, we'll convert "Impact" to numbers for primary ranking
vep_score_mean$Impact <- sub("MODIFIER", "1", vep_score_mean$Impact)
vep_score_mean$Impact <- sub("LOW", "2", vep_score_mean$Impact)
vep_score_mean$Impact <- sub("MODERATE", "3", vep_score_mean$Impact)
vep_score_mean$Impact <- sub("HIGH", "4", vep_score_mean$Impact)
vep_score_mean$Impact <- as.numeric(vep_score_mean$Impact)

## Rank in order of "Impact" then "FLYCADD" to find most impactful consqeuence per SNP
vep_score_mean <- vep_score_mean[order(-vep_score_mean[,4], -vep_score_mean[,3]), ]

## The following loop will select unique vep entries per locus prioritized in order of
## the variant consequence's impact and mean flycadd score
vep_priority <- c()
vep_unique <- merge(vep_cons,haf.freq.Fmean,by=c("CHROM","POS"))
for(r in 1:nrow(vep_score_mean)) {
	temp_impact <- unique(vep_unique[vep_unique$Consequence == vep_score_mean[r,1],])
	vep_priority <- rbind(vep_priority,temp_impact)
	vep_unique <- anti_join(vep_unique, temp_impact, by = c("CHROM", "POS"))
}


#################################
### Calculate rolling windows ###
#################################

### 21 consecutive SNP window ###

## Create indexes for each chromosome to slide over and calculate mean over 21-SNP windows
## in order to reduce noise and detect real signals from from regions under strong 
## selection (evident from likely selective sweeps).

## If we don't do this per chromosome, the average will be calculated across consecutive 
## gaps between chromosomes.

## 2L
glm_2L <- unique(glm.all[glm.all$CHROM=="2L",]) 
glm_2L <- glm_2L %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2L <- glm_2L[is.finite(rowSums(glm_2L[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2L$POS <- as.integer(glm_2L$POS)
glm_2L.rolwin21 <- glm_2L[,c(1,2)]

for(i in 3:ncol(glm_2L)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2L[is.finite(glm_2L[,i])=="TRUE",i])+100
	glm_2L[is.infinite(glm_2L[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2L.rolwin21 <- cbind(glm_2L.rolwin21,slide_mean(glm_2L[,i], before=10, after=10))
    colnames(glm_2L.rolwin21)[i] <- paste(names(glm_2L)[i], ".rolwin21", sep="")
}

## 2R
glm_2R <- unique(glm.all[glm.all$CHROM=="2R",])
glm_2R <- glm_2R %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2R <- glm_2R[is.finite(rowSums(glm_2R[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2R$POS <- as.integer(glm_2R$POS)
glm_2R.rolwin21 <- glm_2R[,c(1,2)]

for(i in 3:ncol(glm_2R)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2R[is.finite(glm_2R[,i])=="TRUE",i])+100
	glm_2R[is.infinite(glm_2R[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2R.rolwin21 <- cbind(glm_2R.rolwin21,slide_mean(glm_2R[,i], before=10, after=10))
    colnames(glm_2R.rolwin21)[i] <- paste(names(glm_2R)[i], ".rolwin21", sep="")
}

## 3L
glm_3L <- unique(glm.all[glm.all$CHROM=="3L",])
glm_3L <- glm_3L %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3L <- glm_3L[is.finite(rowSums(glm_3L[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3L$POS <- as.integer(glm_3L$POS)
glm_3L.rolwin21 <- glm_3L[,c(1,2)]

for(i in 3:ncol(glm_3L)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3L[is.finite(glm_3L[,i])=="TRUE",i])+100
	glm_3L[is.infinite(glm_3L[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3L.rolwin21 <- cbind(glm_3L.rolwin21,slide_mean(glm_3L[,i], before=10, after=10))
    colnames(glm_3L.rolwin21)[i] <- paste(names(glm_3L)[i], ".rolwin21", sep="")
}

## 3R
glm_3R <- unique(glm.all[glm.all$CHROM=="3R",])
glm_3R <- glm_3R %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3R <- glm_3R[is.finite(rowSums(glm_3R[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3R$POS <- as.integer(glm_3R$POS)
glm_3R.rolwin21 <- glm_3R[,c(1,2)]

for(i in 3:ncol(glm_3R)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3R[is.finite(glm_3R[,i])=="TRUE",i])+100
	glm_3R[is.infinite(glm_3R[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3R.rolwin21 <- cbind(glm_3R.rolwin21,slide_mean(glm_3R[,i], before=10, after=10))
    colnames(glm_3R.rolwin21)[i] <- paste(names(glm_3R)[i], ".rolwin21", sep="")
}

## X
glm_X <- unique(glm.all[glm.all$CHROM=="X",])
glm_X <- glm_X %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_X <- glm_X[is.finite(rowSums(glm_X[,-c(1:2)])),]
## create new empty file for filtered rows
glm_X$POS <- as.integer(glm_X$POS) #reformat
glm_X.rolwin21 <- glm_X[,c(1,2)]

for(i in 3:ncol(glm_X)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_X[is.finite(glm_X[,i])=="TRUE",i])+100
	glm_X[is.infinite(glm_X[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_X.rolwin21 <- cbind(glm_X.rolwin21,slide_mean(glm_X[,i], before=10, after=10))
    colnames(glm_X.rolwin21)[i] <- paste(names(glm_X)[i], ".rolwin21", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm.all.rolwin21 <- rbind(glm_2L.rolwin21,glm_2R.rolwin21,glm_3L.rolwin21,glm_3R.rolwin21,glm_X.rolwin21)

## Remove NAs before proceeding
glm.all.rolwin21 <- na.omit(glm.all.rolwin21)

## Merge annotations for filtering and highlighting
glm.all.rolwin21.annot <- merge(glm.all.rolwin21, vep, by=c("CHROM","POS"))

## Make a rolling window table with just p-value, remove logp and fdr columns 
glm.all.rolwin21.pval <- select(glm.all.rolwin21,-contains("logp"),-contains("fdr"))

## Make a rolling window table with just logp, selects them and append loci
glm.all.rolwin21.logp <- cbind(glm.all.rolwin21.pval[,c(1:2)],select(glm.all.rolwin21,contains("logp")))

## Make a rolling window table with just fdr, selects them and append loci
glm.all.rolwin21.fdr <- cbind(glm.all.rolwin21.pval[,c(1:2)],select(glm.all.rolwin21,contains("fdr")))

## Load spinosad-resistance candidate gene list
spino.cand <- read.table("spino.cand.list.txt", header=FALSE)
names(spino.cand) <- "Gene"

## Now merge to find All SNPs in and around candidate genes
glm.all.rolwin21.annot.spino <- merge(spino.cand, glm.all.rolwin21.annot, by="Gene")


####################################
### Hypergeometric overlap tests ###
####################################

## Calculate minimum p-values and FDR values per annotated spino candidate gene. We will 
## assign the most significant FDR value overlapping the gene, plus its immediate upstream 
## and downstream regions, to represent the gene in comparison with background gene sets.

## Pull FDR columns only for spino candidates
glm.all.fdr.rolwin21.annot.spino <- cbind(glm.all.rolwin21.annot.spino[,c(1:3)],glm.all.rolwin21.annot.spino$Consequence,select(glm.all.rolwin21.annot.spino,contains("fdr")))

## Remove duplicates
glm.all.fdr.rolwin21.annot.spino <- unique(glm.all.fdr.rolwin21.annot.spino)

## Loop through each GLM-results FDR column to calculate spino candidate minimum values
glm.all.fdr.rolwin21.annot.spino.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin21.annot.spino)) {
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin21.annot.spino$Gene,glm.all.fdr.rolwin21.annot.spino[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"
    #assign minimum FDR values
    glm.all.fdr.rolwin21.annot.spino.min <- cbind(glm.all.fdr.rolwin21.annot.spino.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin21.annot.spino.min) <- names(glm.all.fdr.rolwin21.annot.spino)[c(5:ncol(glm.all.fdr.rolwin21.annot.spino))]
## Add gene names to the table
glm.all.fdr.rolwin21.annot.spino.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin21.annot.spino.min)
        	
## Save minimum FDR results for downstream use
write.table(glm.all.fdr.rolwin21.annot.spino.min, file="rudflies_2023_redo.glm.all.rolwin21.minfdr.spino.txt", quote = FALSE, sep = "\t", row.names = F)
        

## Pull columns for all genes/features and SNP consequences
glm.all.fdr.rolwin21.annot <- cbind(Gene=glm.all.rolwin21.annot$Gene,glm.all.rolwin21.annot[,c(1,2)],Consequence=glm.all.rolwin21.annot$Consequence)

## Remove duplicates
glm.all.fdr.rolwin21.annot <- unique(glm.all.fdr.rolwin21.annot)

## Finally, append FDR values
glm.all.fdr.rolwin21.annot <- merge(glm.all.fdr.rolwin21.annot, glm.all.rolwin21.fdr, by=c("CHROM","POS"))

## Loop through each GLM-results FDR column to calculate minimum values for all genes
glm.all.fdr.rolwin21.annot.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin21.annot)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin21.annot$Gene,glm.all.fdr.rolwin21.annot[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.rolwin21.annot.min <- cbind(glm.all.fdr.rolwin21.annot.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin21.annot.min) <- names(glm.all.fdr.rolwin21.annot)[c(5:ncol(glm.all.fdr.rolwin21.annot))]
## Add gene names to the table
glm.all.fdr.rolwin21.annot.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin21.annot.min)


### Redo above with missense SNPs in candidates only ###

## Merge to find missense SNPs in and around candidate genes
glm.all.rolwin21.annot.spino.missense <- merge(spino.cand, glm.all.rolwin21.annot[glm.all.rolwin21.annot$Consequence=="missense_variant",], by="Gene")

## Pull FDR columns only for spino candidate missense SNPs
glm.all.fdr.rolwin21.annot.spino.missense <- unique(cbind(glm.all.rolwin21.annot.spino.missense[,c(1:3)],glm.all.rolwin21.annot.spino.missense$Consequence,select(glm.all.rolwin21.annot.spino.missense,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate spino candidate minimum values
glm.all.fdr.rolwin21.annot.spino.missense.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin21.annot.spino.missense)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin21.annot.spino.missense$Gene,glm.all.fdr.rolwin21.annot.spino.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.rolwin21.annot.spino.missense.min <- cbind(glm.all.fdr.rolwin21.annot.spino.missense.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin21.annot.spino.missense.min) <- names(glm.all.fdr.rolwin21.annot.spino.missense)[c(5:ncol(glm.all.fdr.rolwin21.annot.spino.missense))]
## Add gene names to the table
glm.all.fdr.rolwin21.annot.spino.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin21.annot.spino.missense.min)

## Save minimum FDR results from missense SNPs for downstream use
write.table(glm.all.fdr.rolwin21.annot.spino.missense.min, file="rudflies_2023_redo.glm.all.rolwin21.minfdr.spino.missense.txt", quote = FALSE, sep = "\t", row.names = F)


## Pull FDR columns only for all genes/features, filtered for missense SNPs only
glm.all.fdr.rolwin21.annot.missense <- glm.all.fdr.rolwin21.annot[glm.all.fdr.rolwin21.annot$Consequence == "missense_variant",]

## Loop through each GLM-results FDR column to calculate minimum values for all genes
glm.all.fdr.rolwin21.annot.missense.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin21.annot.missense)) {
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin21.annot.missense$Gene,glm.all.fdr.rolwin21.annot.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.rolwin21.annot.missense.min <- cbind(glm.all.fdr.rolwin21.annot.missense.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin21.annot.missense.min) <- names(glm.all.fdr.rolwin21.annot.missense)[c(5:ncol(glm.all.fdr.rolwin21.annot.missense))]
## Add gene names to the table
glm.all.fdr.rolwin21.annot.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin21.annot.missense.min)


## Now we've found minimum FDR values for candidate gene and background gene lists, 
## let's run hypergeometric tests on the overlap of candidates and GLM outliers.

## run hypergeometric overlap test per GLM contrast, 21-SNP rolling windows
## First with all gene-related/adjacent SNPS
## Significance threshold FDR<0.05
phyper_all_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.spino.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.spino.min[glm.all.fdr.rolwin21.annot.spino.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.spino.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.min[glm.all.fdr.rolwin21.annot.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_contrasts05 <- append(phyper_all_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)    
}

## Create results summary table
phyper_spinoCand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.spino.min[,2:ncol(glm.all.fdr.rolwin21.annot.spino.min)]),phyper_pval=phyper_all_contrasts05,phyper_p_sig=phyper_all_contrasts05 < 0.05,phyper_fdr=p.adjust(phyper_all_contrasts05, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_all_contrasts05, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_all_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.spino.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.spino.min[glm.all.fdr.rolwin21.annot.spino.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.spino.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.min[glm.all.fdr.rolwin21.annot.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_contrasts01 <- append(phyper_all_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_spinoCand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.spino.min[,2:ncol(glm.all.fdr.rolwin21.annot.spino.min)]),phyper_pval=phyper_all_contrasts01,phyper_p_sig=phyper_all_contrasts01 < 0.05,phyper_fdr=p.adjust(phyper_all_contrasts01, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_all_contrasts01, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)


### Next with only missense SNPS ###

## Significance threshold FDR<0.05
phyper_missense_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.spino.missense.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.spino.missense.min[glm.all.fdr.rolwin21.annot.spino.missense.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.spino.missense.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[glm.all.fdr.rolwin21.annot.missense.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_contrasts05 <- append(phyper_missense_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_spinoCand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.spino.missense.min[,2:ncol(glm.all.fdr.rolwin21.annot.spino.missense.min)]),phyper_pval=phyper_missense_contrasts05,phyper_p_sig=phyper_missense_contrasts05 < 0.05,phyper_fdr=p.adjust(phyper_missense_contrasts05, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_missense_contrasts05, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_missense_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.spino.missense.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.spino.missense.min[glm.all.fdr.rolwin21.annot.spino.missense.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.spino.missense.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[glm.all.fdr.rolwin21.annot.missense.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_contrasts01 <- append(phyper_missense_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_spinoCand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.spino.missense.min[,2:ncol(glm.all.fdr.rolwin21.annot.spino.missense.min)]),phyper_pval=phyper_missense_contrasts01,phyper_p_sig=phyper_missense_contrasts01 < 0.05,phyper_fdr=p.adjust(phyper_missense_contrasts01, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_missense_contrasts01, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.05)
write.table(phyper_spinoCand_FDR05, file="rudflies_2023_redo.phyper_allSNPs_spinoCand_FDR05.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.01)
write.table(phyper_spinoCand_FDR01, file="rudflies_2023_redo.phyper_allSNPs_spinoCand_FDR01.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.05), missense only
write.table(phyper_missense_spinoCand_FDR05, file="rudflies_2023_redo.phyper_missenseSNPs_spinoCand_FDR05.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.01), missense only
write.table(phyper_missense_spinoCand_FDR01, file="rudflies_2023_redo.phyper_missenseSNPs_spinoCand_FDR01.txt", quote = FALSE, sep = "\t", row.names = F)

### Which contrasts show significant overlap, prior to p-value correction

## spino candidates vs GLM outliers (FDR < 0.05), all SNPs
phyper_spinoCand_FDR05[phyper_spinoCand_FDR05[,2]<0.05,]
## spino candidates vs GLM outliers (FDR < 0.01), all SNPs
phyper_spinoCand_FDR01[phyper_spinoCand_FDR01[,2]<0.05,]
## spino candidates vs GLM outliers (FDR < 0.05), missense SNPs
phyper_missense_spinoCand_FDR05[phyper_missense_spinoCand_FDR05[,2]<0.05,]
## spino candidates vs GLM outliers (FDR < 0.01), missense SNPs
phyper_missense_spinoCand_FDR01[phyper_missense_spinoCand_FDR01[,2]<0.05,]


#############################################################
### Examine candidate gene p-values via Wilcoxon rank sum ###
#############################################################

## Hypergeometric overlap analysis relies on GLM outlier cutoff decisions, which can 
## impact results. Wilcoxon Rank Sum tests are non-parametric and rely on rank order of
## all genes in the dataset. Therefore, this test is more sensitive to more modest 
## AF changes across a broader gene network or pathway.

### Remake annotation stat summary lists for pvalues instead of fdr values ###

## Pull raw p-value columns only for spino candidates
glm.all.pval.rolwin21.annot.spino <- cbind(glm.all.rolwin21.annot.spino[,c(1:3)],Consequence=glm.all.rolwin21.annot.spino$Consequence,select(glm.all.rolwin21.annot.spino,contains("rolwin21"),-contains("fdr"),-contains("logp")))

## Remove duplicates
glm.all.pval.rolwin21.annot.spino <- unique(glm.all.pval.rolwin21.annot.spino)

## Remove na's 
glm.all.pval.rolwin21.annot.spino <- na.omit(glm.all.pval.rolwin21.annot.spino)

## Loop through each GLM-results pval column to calculate spino candidate minimum values
glm.all.pval.rolwin21.annot.spino.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin21.annot.spino)) {
  tryCatch({
  	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin21.annot.spino$Gene,glm.all.pval.rolwin21.annot.spino[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum p-values    
    glm.all.pval.rolwin21.annot.spino.min <- cbind(glm.all.pval.rolwin21.annot.spino.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin21.annot.spino.min) <- names(glm.all.pval.rolwin21.annot.spino)[c(5:ncol(glm.all.pval.rolwin21.annot.spino))]
## Add gene names to the table
glm.all.pval.rolwin21.annot.spino.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin21.annot.spino.min)
        	
## Save minimum pval results for downstream use
write.table(glm.all.pval.rolwin21.annot.spino.min, file="rudflies_2023_redo.glm.all.rolwin21.minpval.spino.txt", quote = FALSE, sep = "\t", row.names = F)
        	
        	
## Pull columns for all genes/features and SNP consequences
glm.all.pval.rolwin21.annot <- cbind(Gene=glm.all.rolwin21.annot$Gene,glm.all.rolwin21.annot[,c(1,2)],Consequence=glm.all.rolwin21.annot$Consequence)

## Remove duplicates 
glm.all.pval.rolwin21.annot <- unique(glm.all.pval.rolwin21.annot)

## Remove na's 
glm.all.pval.rolwin21.annot <- na.omit(glm.all.pval.rolwin21.annot)

## Finally, append p-values
glm.all.pval.rolwin21.annot <- merge(glm.all.pval.rolwin21.annot, glm.all.rolwin21.pval, by=c("CHROM","POS"))

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.rolwin21.annot.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin21.annot)) {
  tryCatch({
    #create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin21.annot$Gene,glm.all.pval.rolwin21.annot[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"  
    #assign minimum p-values  
    glm.all.pval.rolwin21.annot.min <- cbind(glm.all.pval.rolwin21.annot.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin21.annot.min) <- names(glm.all.pval.rolwin21.annot)[c(5:ncol(glm.all.pval.rolwin21.annot))]
## Add gene names to the table
glm.all.pval.rolwin21.annot.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin21.annot.min) 	
        	
## Pull pval columns only for all genes/features, filtered for missense SNPs only
glm.all.pval.rolwin21.annot.spino.missense <- cbind(glm.all.rolwin21.annot.spino.missense[,c(1:3)],glm.all.rolwin21.annot.spino.missense$Consequence,select(glm.all.rolwin21.annot.spino.missense,contains("rolwin21"),-contains("fdr"),-contains("logp")))

## Remove duplicates
glm.all.pval.rolwin21.annot.spino.missense <- unique(glm.all.pval.rolwin21.annot.spino.missense)

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.rolwin21.annot.spino.missense.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin21.annot.spino.missense)) { #nrow(haf.freq.cand)
  tryCatch({
	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin21.annot.spino.missense$Gene,glm.all.pval.rolwin21.annot.spino.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum FDR values    
    glm.all.pval.rolwin21.annot.spino.missense.min <- cbind(glm.all.pval.rolwin21.annot.spino.missense.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin21.annot.spino.missense.min) <- names(glm.all.pval.rolwin21.annot.spino.missense)[c(5:ncol(glm.all.pval.rolwin21.annot.spino.missense))]
## Add gene names to the table
glm.all.pval.rolwin21.annot.spino.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin21.annot.spino.missense.min)


## Save minimum pval results from missense SNPs for downstream use
write.table(glm.all.pval.rolwin21.annot.spino.missense.min, file="rudflies_2023_redo.glm.all.rolwin21.minpval.spino.missense.txt", quote = FALSE, sep = "\t", row.names = F)


## Pull pval columns only for all genes/features, filtered for missense SNPs only
glm.all.pval.rolwin21.annot.missense <- glm.all.pval.rolwin21.annot[glm.all.pval.rolwin21.annot$Consequence == "missense_variant",]

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.rolwin21.annot.missense.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin21.annot.missense)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin21.annot.missense$Gene,glm.all.pval.rolwin21.annot.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"   
    #assign minimum p-values 
    glm.all.pval.rolwin21.annot.missense.min <- cbind(glm.all.pval.rolwin21.annot.missense.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin21.annot.missense.min) <- names(glm.all.pval.rolwin21.annot.missense)[c(5:ncol(glm.all.pval.rolwin21.annot.missense))]
## Add gene names to the table
glm.all.pval.rolwin21.annot.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin21.annot.missense.min)
        

## Time to run the Wicoxon Rank Sum analyses. To do so, we will select one matched 
## background gene for every spinosad candidate gene. This matched set selection is 
## based on the following criteria:
## 	 -Gene must be similar length from start to stop (+/-25%)
## 	 -Gene must have similar number of SNPs (+/-25%)
## 	 -Gene must be the same type (protein-coding, ncRNA, etc)
## 	 -Gene must be on the same chromosome

## After all matching genes are discovered, we randomly select one match for each 
## candidate gene, then perform the Wilcoxon rank sum test between the candidate and 
## background genes to see if candidates show significantly higher p-value ranks (lower 
## p-values), or greater parallel allele frequency differences in the respective GLM 
## contrast. This process is repeated for 1K iterations of matched set selection to 
## mitigate any potential random sampling bias. 


##Find gene characteristics of candidates from GFF for matching
glm.all.pval.rolwin21.annot.spino.min.gff <- merge(snp.gff.overlap, glm.all.pval.rolwin21.annot.spino.min, by.x="ID", by.y="Gene")

##Find gene characteristics of background genes from GFF for matching
glm.all.pval.rolwin21.annot.min.gff <- merge(snp.gff.overlap, glm.all.pval.rolwin21.annot.min, by.x="ID", by.y="Gene")

## Pull spino candidate gene IDs
spino_IDs <- glm.all.pval.rolwin21.annot.spino.min.gff$ID

## Remove spino candidates from potential background gene list
glm.all.pval.rolwin21.annot.bg.min.gff <- glm.all.pval.rolwin21.annot.min.gff[!(glm.all.pval.rolwin21.annot.min.gff$ID %in% spino_IDs),]

## Create empty results table with 1K rows (for each iteration) and one column for 
## each GLM contrast.
wilcox.rolwin21.p <- data.frame(matrix(NA, nrow = 1000, ncol = ncol(glm.all.pval.rolwin21.annot.spino.min.gff)-ncol(snp.gff.overlap)))
## Modify column names to reflect GLM contrasts
names(wilcox.rolwin21.p) <- names(glm.all.pval.rolwin21.annot.spino.min.gff)[c(9:ncol(glm.all.pval.rolwin21.annot.spino.min.gff))]

#generate candidate gene filtering tables
for(k in 1:1000) { #set number of iterations
  #Establish background non-candidate set
  bg.samp.list <- c()
  #This loop finds a random set of matched background genes
  for(i in c(1:nrow(glm.all.pval.rolwin21.annot.spino.min.gff))) {
  	#Select matched lists based on following criteria
	length <- glm.all.pval.rolwin21.annot.spino.min.gff[i,]$LENGTH
	count <- glm.all.pval.rolwin21.annot.spino.min.gff[i,]$SNP_COUNT
	type <- glm.all.pval.rolwin21.annot.spino.min.gff[i,]$TYPE
	chrom <- glm.all.pval.rolwin21.annot.spino.min.gff[i,]$CHROM
	#Pull one match per candidate gene
	bg.samp.temp <- sample_n(glm.all.pval.rolwin21.annot.bg.min.gff %>%
	 filter( 
	    TYPE==type,
		CHROM==chrom,
		LENGTH < length*1.25,
		LENGTH > length*0.75,
		SNP_COUNT < count*1.25,
		SNP_COUNT > count*0.75
	  ),1)
	#append to background list
	bg.samp.list <- rbind(bg.samp.list,bg.samp.temp)
  }
  #This loop will run Wilcoxon Rank Sum tests for each of "j" GLM contrasts
  for(j in 1:(ncol(glm.all.pval.rolwin21.annot.spino.min.gff)-ncol(snp.gff.overlap))) { #set number of iterations to the number of contrasts
    #Establish Rank Sum test sets for p values
    #candidate set
    candidates <- glm.all.pval.rolwin21.annot.spino.min.gff[,8+j]  
    #background non-candidate set
    bg.samp <- as.numeric(bg.samp.list[,8+j])
    candidates.samp <- as.numeric(candidates)
    #this will subsample the SNPs contained in one list to the number of SNPs in the shorter list
    if(length(candidates.samp) > length(bg.samp)){
      candidates.samp <- sample(candidates.samp, length(bg.samp))
    } else {
      bg.samp <- sample(bg.samp, length(candidates.samp))
    }
    #finally, run the tests and append results
    wilcox.rolwin21.p[k,j] <- wilcox.test(candidates.samp, bg.samp, alternative = "less")$p.value
}
}

## Save all results
# Wilcoxon rank sum p-values for all iterations
write.table(wilcox.rolwin21.p, file="rudflies_2023_redo.glm.all.wilcox.rolwin21.p.txt", quote = FALSE, sep = "\t", row.names = F)

# Median Wilcoxon rank sum p-values across iterations
write.table(colMedians(as.matrix(wilcox.rolwin21.p)), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.p.medians.txt", quote = FALSE, sep = "\t", row.names = F)

# Mean Wilcoxon rank sum p-values across iterations
write.table(colMeans(as.matrix(wilcox.rolwin21.p)), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.p.means.txt", quote = FALSE, sep = "\t", row.names = F)

## Calculate the percent of p-values significant at a raw p<0.05
colMedians(as.matrix(wilcox.rolwin21.p))
colMeans(as.matrix(wilcox.rolwin21.p))
wilcox.rolwin21.contrast <- c()
wilcox.rolwin21.contrast.p05 <- wilcox.rolwin21.p<0.05
for(l in 1:ncol(wilcox.rolwin21.p)) {
   wilcox.rolwin21.contrast <- append(wilcox.rolwin21.contrast,sum(wilcox.rolwin21.contrast.p05[,l], na.rm=TRUE)*0.1)
}
cbind(contrast=names(wilcox.rolwin21.p),perc_sig=wilcox.rolwin21.contrast)

# Save percent of significant Wilcoxon rank sum p-values across iterations
write.table(cbind(contrast=names(wilcox.rolwin21.p),perc_sig=wilcox.rolwin21.contrast), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.p.percsig.txt", quote = FALSE, sep = "\t", row.names = F)





##############################################
###   Use EvSE.T1 candidates as test set   ###
##############################################

## Find EvSE.T1 candidate gene list (exclude flanking features)
EvSE.T1.cand <- data.frame(unique(na.omit(glm.all.rolwin21.annot[glm.all.rolwin21.annot$EvSE.T1.fdr < 0.05 & glm.all.rolwin21.annot$Consequence != "downstream_variant" & glm.all.rolwin21.annot$Consequence != "upstream_variant" & glm.all.rolwin21.annot$Feature_type == "Transcript",])$Gene))
names(EvSE.T1.cand) <- "Gene"

## Now merge to find All SNPs in and around candidate genes
glm.all.rolwin21.annot.EvSE.T1 <- merge(EvSE.T1.cand, glm.all.rolwin21.annot, by="Gene")


####################################
### Hypergeometric overlap tests ###
####################################

## Calculate minimum p-values and FDR values per annotated EvSE.T1 candidate genes. We 
## will assign the most significant FDR value overlapping the gene, plus its immediate 
## upstream and downstream regions, to represent the gene in comparison with background 
## gene sets.

## Pull gene and SNP consequence columns for EvSE.T1 candidates
glm.all.fdr.rolwin21.annot.EvSE.T1 <- cbind(glm.all.rolwin21.annot.EvSE.T1[,c(1:3)],glm.all.rolwin21.annot.EvSE.T1$Consequence)

## Remove duplicates
glm.all.fdr.rolwin21.annot.EvSE.T1 <- unique(glm.all.fdr.rolwin21.annot.EvSE.T1)

## Finally, add FDR values
glm.all.fdr.rolwin21.annot.EvSE.T1 <- merge(glm.all.fdr.rolwin21.annot.EvSE.T1, glm.all.rolwin21.fdr, by=c("CHROM","POS"))

## Loop through each GLM-results FDR column to calculate EvSE.T1 candidate minimum values
glm.all.fdr.rolwin21.annot.EvSE.T1.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1)) {
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin21.annot.EvSE.T1$Gene,glm.all.fdr.rolwin21.annot.EvSE.T1[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"
    #assign minimum FDR values
    glm.all.fdr.rolwin21.annot.EvSE.T1.min <- cbind(glm.all.fdr.rolwin21.annot.EvSE.T1.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin21.annot.EvSE.T1.min) <- names(glm.all.fdr.rolwin21.annot.EvSE.T1)[c(5:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1))]
## Add gene names to the table
glm.all.fdr.rolwin21.annot.EvSE.T1.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin21.annot.EvSE.T1.min)
        	
## Save minimum FDR results for downstream use
write.table(glm.all.fdr.rolwin21.annot.EvSE.T1.min, file="rudflies_2023_redo.glm.all.rolwin21.minfdr.EvSE.T1.txt", quote = FALSE, sep = "\t", row.names = F)
        	

### Redo above with missense SNPs in candidates only ###

## Merge to find missense SNPs in and around candidate genes
glm.all.rolwin21.annot.EvSE.T1.missense <- merge(EvSE.T1.cand, glm.all.rolwin21.annot[glm.all.rolwin21.annot$Consequence=="missense_variant",], by="Gene")

## Pull FDR columns only for EvSE.T1 candidate missense SNPs
glm.all.fdr.rolwin21.annot.EvSE.T1.missense <- unique(cbind(glm.all.rolwin21.annot.EvSE.T1.missense[,c(1:3)],glm.all.rolwin21.annot.EvSE.T1.missense$Consequence,select(glm.all.rolwin21.annot.EvSE.T1.missense,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate EvSE.T1 candidate minimum values
glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.missense)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin21.annot.EvSE.T1.missense$Gene,glm.all.fdr.rolwin21.annot.EvSE.T1.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min <- cbind(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min) <- names(glm.all.fdr.rolwin21.annot.EvSE.T1.missense)[c(5:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.missense))]
## Add gene names to the table
glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min)

## Save minimum FDR results from missense SNPs for downstream use
write.table(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min, file="rudflies_2023_redo.glm.all.rolwin21.minfdr.EvSE.T1.missense.txt", quote = FALSE, sep = "\t", row.names = F)



## Now we've found minimum FDR values for candidate gene and background gene lists, 
## let's run hypergeometric tests on the overlap of candidates and GLM outliers.

## run hypergeometric overlap test per GLM contrast, 21-SNP rolling windows
## First with all gene-related/adjacent SNPS
## Significance threshold FDR<0.05
phyper_all_EvSE_T1_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.EvSE.T1.min[glm.all.fdr.rolwin21.annot.EvSE.T1.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.EvSE.T1.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.min[glm.all.fdr.rolwin21.annot.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_EvSE_T1_contrasts05 <- append(phyper_all_EvSE_T1_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)    
}

## Create results summary table
phyper_EvSE.T1Cand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.EvSE.T1.min[,2:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.min)]),phyper_pval=phyper_all_EvSE_T1_contrasts05,phyper_p_sig=phyper_all_EvSE_T1_contrasts05 < 0.05,phyper_fdr=p.adjust(phyper_all_EvSE_T1_contrasts05, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_all_EvSE_T1_contrasts05, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_all_EvSE_T1_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.EvSE.T1.min[glm.all.fdr.rolwin21.annot.EvSE.T1.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.EvSE.T1.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.min[glm.all.fdr.rolwin21.annot.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_EvSE_T1_contrasts01 <- append(phyper_all_EvSE_T1_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_EvSE.T1Cand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.EvSE.T1.min[,2:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.min)]),phyper_pval=phyper_all_EvSE_T1_contrasts01,phyper_p_sig=phyper_all_EvSE_T1_contrasts01 < 0.05,phyper_fdr=p.adjust(phyper_all_EvSE_T1_contrasts01, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_all_EvSE_T1_contrasts01, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)


### Next with only missense SNPS ###

## Significance threshold FDR<0.05
phyper_missense_EvSE_T1_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min[glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[glm.all.fdr.rolwin21.annot.missense.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_EvSE_T1_contrasts05 <- append(phyper_missense_EvSE_T1_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_EvSE.T1Cand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min[,2:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min)]),phyper_pval=phyper_missense_EvSE_T1_contrasts05,phyper_p_sig=phyper_missense_EvSE_T1_contrasts05 < 0.05,phyper_fdr=p.adjust(phyper_missense_EvSE_T1_contrasts05, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_missense_EvSE_T1_contrasts05, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_missense_EvSE_T1_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min[glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[glm.all.fdr.rolwin21.annot.missense.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_EvSE_T1_contrasts01 <- append(phyper_missense_EvSE_T1_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_EvSE.T1Cand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min[,2:ncol(glm.all.fdr.rolwin21.annot.EvSE.T1.missense.min)]),phyper_pval=phyper_missense_EvSE_T1_contrasts01,phyper_p_sig=phyper_missense_EvSE_T1_contrasts01 < 0.05,phyper_fdr=p.adjust(phyper_missense_EvSE_T1_contrasts01, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_missense_EvSE_T1_contrasts01, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Save results table for EvSE.T1 candidates vs GLM outliers (FDR < 0.05)
write.table(phyper_EvSE.T1Cand_FDR05, file="rudflies_2023_redo.phyper_allSNPs_EvSE.T1Cand_FDR05.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for EvSE.T1 candidates vs GLM outliers (FDR < 0.01)
write.table(phyper_EvSE.T1Cand_FDR01, file="rudflies_2023_redo.phyper_allSNPs_EvSE.T1Cand_FDR01.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for EvSE.T1 candidates vs GLM outliers (FDR < 0.05), missense only
write.table(phyper_missense_EvSE.T1Cand_FDR05, file="rudflies_2023_redo.phyper_missenseSNPs_EvSE.T1Cand_FDR05.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for EvSE.T1 candidates vs GLM outliers (FDR < 0.01), missense only
write.table(phyper_missense_EvSE.T1Cand_FDR01, file="rudflies_2023_redo.phyper_missenseSNPs_EvSE.T1Cand_FDR01.txt", quote = FALSE, sep = "\t", row.names = F)

### Which contrasts show significant overlap, prior to p-value correction

## EvSE.T1 candidates vs GLM outliers (FDR < 0.05), all SNPs
phyper_EvSE.T1Cand_FDR05[as.numeric(phyper_EvSE.T1Cand_FDR05[,2])<0.05,]
## EvSE.T1 candidates vs GLM outliers (FDR < 0.01), all SNPs
phyper_EvSE.T1Cand_FDR01[as.numeric(phyper_EvSE.T1Cand_FDR01[,2])<0.05,]
## EvSE.T1 candidates vs GLM outliers (FDR < 0.05), missense SNPs
phyper_missense_EvSE.T1Cand_FDR05[as.numeric(phyper_missense_EvSE.T1Cand_FDR05[,2])<0.05,]
## EvSE.T1 candidates vs GLM outliers (FDR < 0.01), missense SNPs
phyper_missense_EvSE.T1Cand_FDR01[as.numeric(phyper_missense_EvSE.T1Cand_FDR01[,2])<0.05,]


#############################################################
### Examine candidate gene p-values via Wilcoxon rank sum ###
#############################################################

## Hypergeometric overlap analysis relies on GLM outlier cutoff decisions, which can 
## impact results. Wilcoxon Rank Sum tests are non-parametric and rely on rank order of
## all genes in the dataset. Therefore, this test is more sensitive to more modest 
## AF changes across a broader gene network or pathway.

### Remake annotation stat summary lists for pvalues instead of fdr values ###

## Pull raw p-value columns only for EvSE.T1 candidates
glm.all.pval.rolwin21.annot.EvSE.T1 <- cbind(glm.all.rolwin21.annot.EvSE.T1[,c(1:3)],Consequence=glm.all.rolwin21.annot.EvSE.T1$Consequence,select(glm.all.rolwin21.annot.EvSE.T1,contains("rolwin21"),-contains("fdr"),-contains("logp")))

## Remove duplicates
glm.all.pval.rolwin21.annot.EvSE.T1 <- unique(glm.all.pval.rolwin21.annot.EvSE.T1)

## Loop through each GLM-results pval column to calculate EvSE.T1 candidate minimum values
glm.all.pval.rolwin21.annot.EvSE.T1.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin21.annot.EvSE.T1)) {
  tryCatch({
  	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin21.annot.EvSE.T1$Gene,glm.all.pval.rolwin21.annot.EvSE.T1[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum p-values    
    glm.all.pval.rolwin21.annot.EvSE.T1.min <- cbind(glm.all.pval.rolwin21.annot.EvSE.T1.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin21.annot.EvSE.T1.min) <- names(glm.all.pval.rolwin21.annot.EvSE.T1)[c(5:ncol(glm.all.pval.rolwin21.annot.EvSE.T1))]
## Add gene names to the table
glm.all.pval.rolwin21.annot.EvSE.T1.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin21.annot.EvSE.T1.min)
        	
## Save minimum pval results for downstream use
write.table(glm.all.pval.rolwin21.annot.EvSE.T1.min, file="rudflies_2023_redo.glm.all.rolwin21.minpval.EvSE.T1.txt", quote = FALSE, sep = "\t", row.names = F)
        		
        	
## Pull pval columns only for all genes/features, filtered for missense SNPs only
glm.all.pval.rolwin21.annot.EvSE.T1.missense <- unique(cbind(glm.all.rolwin21.annot.EvSE.T1.missense[,c(1:3)],glm.all.rolwin21.annot.EvSE.T1.missense$Consequence,select(glm.all.rolwin21.annot.EvSE.T1.missense,contains("rolwin21"),-contains("fdr"),-contains("logp"))))

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.rolwin21.annot.EvSE.T1.missense.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin21.annot.EvSE.T1.missense)) { #nrow(haf.freq.cand)
  tryCatch({
	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin21.annot.EvSE.T1.missense$Gene,glm.all.pval.rolwin21.annot.EvSE.T1.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum FDR values    
    glm.all.pval.rolwin21.annot.EvSE.T1.missense.min <- cbind(glm.all.pval.rolwin21.annot.EvSE.T1.missense.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin21.annot.EvSE.T1.missense.min) <- names(glm.all.pval.rolwin21.annot.EvSE.T1.missense)[c(5:ncol(glm.all.pval.rolwin21.annot.EvSE.T1.missense))]
## Add gene names to the table
glm.all.pval.rolwin21.annot.EvSE.T1.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin21.annot.EvSE.T1.missense.min)


## Save minimum pval results from missense SNPs for downstream use
write.table(glm.all.pval.rolwin21.annot.EvSE.T1.missense.min, file="rudflies_2023_redo.glm.all.rolwin21.minpval.EvSE.T1.missense.txt", quote = FALSE, sep = "\t", row.names = F)        	
        	
        	

## Time to run the Wicoxon Rank Sum analyses. To do so, we will select one matched 
## background gene for every EvSE.T1 candidate gene. This matched set selection is 
## based on the following criteria:
## 	 -Gene must be similar length from start to stop (+/-25%)
## 	 -Gene must have similar number of SNPs (+/-25%)
## 	 -Gene must be the same type (protein-coding, ncRNA, etc)
## 	 -Gene must be on the same chromosome

## After all matching genes are discovered, we randomly select one match for each 
## candidate gene, then perform the Wilcoxon rank sum test between the candidate and 
## background genes to see if candidates show significantly higher p-value ranks (lower 
## p-values), or greater parallel allele frequency differences in the respective GLM 
## contrast. This process is repeated for 1K iterations of matched set selection to 
## mitigate any potential random sampling bias. 


##Find gene characteristics of candidates from GFF for matching
glm.all.pval.rolwin21.annot.EvSE.T1.min.gff <- merge(snp.gff.overlap, glm.all.pval.rolwin21.annot.EvSE.T1.min, by.x="ID", by.y="Gene")

##Find gene characteristics of background genes from GFF for matching
glm.all.pval.rolwin21.annot.min.gff <- merge(snp.gff.overlap, glm.all.pval.rolwin21.annot.min, by.x="ID", by.y="Gene")

## Pull EvSE.T1 candidate gene IDs
EvSE.T1_IDs <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff$ID

## Remove EvSE.T1 candidates from potential background gene list
glm.all.pval.rolwin21.annot.bg.min.gff <- glm.all.pval.rolwin21.annot.min.gff[!(glm.all.pval.rolwin21.annot.min.gff$ID %in% EvSE.T1_IDs),]

## Create empty results table with 1K rows (for each iteration) and one column for 
## each GLM contrast.
wilcox.rolwin21.EvSE.T1.p <- data.frame(matrix(NA, nrow = 1000, ncol = ncol(glm.all.pval.rolwin21.annot.EvSE.T1.min.gff)-ncol(snp.gff.overlap)))
## Modify column names to reflect GLM contrasts
names(wilcox.rolwin21.EvSE.T1.p) <- names(glm.all.pval.rolwin21.annot.EvSE.T1.min.gff)[c(9:ncol(glm.all.pval.rolwin21.annot.EvSE.T1.min.gff))]

###First, determine how many matches so we can exclude those without any
bg.samp.counts <- c()
  #This loop finds a random set of matched background genes
  for(i in c(1:nrow(glm.all.pval.rolwin21.annot.EvSE.T1.min.gff))) {
  	#Select matched lists based on following criteria
	length <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[i,]$LENGTH
	count <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[i,]$SNP_COUNT
	type <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[i,]$TYPE
	chrom <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[i,]$CHROM
	#Pull one match per candidate gene
	bg.samp.temp <- nrow(glm.all.pval.rolwin21.annot.bg.min.gff %>%
	 filter( 
	    TYPE==type,
		CHROM==chrom,
		LENGTH < length*1.25,
		LENGTH > length*0.75,
		SNP_COUNT < count*1.25,
		SNP_COUNT > count*0.75))
	#append to background list
	bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
  }

#generate candidate gene filtering tables
for(k in 1:1000) { #set number of iterations
  #Establish background non-candidate set
  bg.samp.list <- c()
  #This loop finds a random set of matched background genes
  for(i in as.numeric(rownames(glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[bg.samp.counts[,1]>4,]))) {
  	#Select matched lists based on following criteria
	length <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[i,]$LENGTH
	count <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[i,]$SNP_COUNT
	type <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[i,]$TYPE
	chrom <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[i,]$CHROM
	#Pull one match per candidate gene
	bg.samp.temp <- sample_n(glm.all.pval.rolwin21.annot.bg.min.gff %>%
	 filter( 
	    TYPE==type,
		CHROM==chrom,
		LENGTH < length*1.25,
		LENGTH > length*0.75,
		SNP_COUNT < count*1.25,
		SNP_COUNT > count*0.75
	  ),1)
	#append to background list
	bg.samp.list <- rbind(bg.samp.list,bg.samp.temp)
  }
  #This loop will run Wilcoxon Rank Sum tests for each of "j" GLM contrasts
  for(j in 1:(ncol(glm.all.pval.rolwin21.annot.EvSE.T1.min.gff)-ncol(snp.gff.overlap))) { #set number of iterations to the number of contrasts
    #Establish Rank Sum test sets for p values
    #candidate set
    candidates <- glm.all.pval.rolwin21.annot.EvSE.T1.min.gff[,8+j]  
    #background non-candidate set
    bg.samp <- as.numeric(bg.samp.list[,8+j])
    candidates.samp <- as.numeric(candidates)
    #this will subsample the SNPs contained in one list to the number of SNPs in the shorter list
    if(length(candidates.samp) > length(bg.samp)){
      candidates.samp <- sample(candidates.samp, length(bg.samp))
    } else {
      bg.samp <- sample(bg.samp, length(candidates.samp))
    }
    #finally, run the tests and append results
    wilcox.rolwin21.EvSE.T1.p[k,j] <- wilcox.test(candidates.samp, bg.samp, alternative = "less")$p.value
}
}

## Save all results
# Wilcoxon rank sum p-values for all iterations
write.table(wilcox.rolwin21.EvSE.T1.p, file="rudflies_2023_redo.glm.all.wilcox.rolwin21.EvSE.T1.p.txt", quote = FALSE, sep = "\t", row.names = F)

# Median Wilcoxon rank sum p-values across iterations
write.table(colMedians(as.matrix(wilcox.rolwin21.EvSE.T1.p)), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.EvSE.T1.p.medians.txt", quote = FALSE, sep = "\t", row.names = F)

# Mean Wilcoxon rank sum p-values across iterations
write.table(colMeans(as.matrix(wilcox.rolwin21.EvSE.T1.p)), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.EvSE.T1.p.means.txt", quote = FALSE, sep = "\t", row.names = F)

## Calculate the percent of p-values significant at a raw p<0.05
colMedians(as.matrix(wilcox.rolwin21.EvSE.T1.p))
colMeans(as.matrix(wilcox.rolwin21.EvSE.T1.p))
wilcox.rolwin21.EvSE.T1.contrast <- c()
wilcox.rolwin21.EvSE.T1.contrast.p05 <- wilcox.rolwin21.EvSE.T1.p<0.05
for(l in 1:ncol(wilcox.rolwin21.EvSE.T1.p)) {
   wilcox.rolwin21.EvSE.T1.contrast <- append(wilcox.rolwin21.EvSE.T1.contrast,sum(wilcox.rolwin21.EvSE.T1.contrast.p05[,l], na.rm=TRUE)*0.1)
}
cbind(contrast=names(wilcox.rolwin21.EvSE.T1.p),perc_sig=wilcox.rolwin21.EvSE.T1.contrast)

# Save percent of significant Wilcoxon rank sum p-values across iterations
write.table(cbind(contrast=names(wilcox.rolwin21.EvSE.T1.p),perc_sig=wilcox.rolwin21.EvSE.T1.contrast), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.EvSE.T1.p.percsig.txt", quote = FALSE, sep = "\t", row.names = F)


##############################################
###   Use EvSP.T4 candidates as test set   ###
##############################################

## Find EvSP.T4 candidate gene list (excluding flanking features)
EvSP.T4.cand <- data.frame(unique(na.omit(glm.all.rolwin21.annot[glm.all.rolwin21.annot$SvE.T4.fdr < 0.05 & glm.all.rolwin21.annot$Consequence != "downstream_variant" & glm.all.rolwin21.annot$Consequence != "upstream_variant" & glm.all.rolwin21.annot$Feature_type == "Transcript",])$Gene))
names(EvSP.T4.cand) <- "Gene"

## Now merge to find All SNPs in and around candidate genes
glm.all.rolwin21.annot.EvSP.T4 <- merge(EvSP.T4.cand, glm.all.rolwin21.annot, by="Gene")


####################################
### Hypergeometric overlap tests ###
####################################

## Calculate minimum p-values and FDR values per annotated EvSP.T4 candidate genes. We 
## will assign the most significant FDR value overlapping the gene, plus its immediate 
## upstream and downstream regions, to represent the gene in comparison with background 
## gene sets.

## Pull gene and SNP consequence columns for EvSP.T4 candidates
glm.all.fdr.rolwin21.annot.EvSP.T4 <- cbind(glm.all.rolwin21.annot.EvSP.T4[,c(1:3)],glm.all.rolwin21.annot.EvSP.T4$Consequence)

## Remove duplicates
glm.all.fdr.rolwin21.annot.EvSP.T4 <- unique(glm.all.fdr.rolwin21.annot.EvSP.T4 )

## Finally, add FDR values
glm.all.fdr.rolwin21.annot.EvSP.T4 <- merge(glm.all.fdr.rolwin21.annot.EvSP.T4, glm.all.rolwin21.fdr, by=c("CHROM","POS"))

## Loop through each GLM-results FDR column to calculate EvSP.T4 candidate minimum values
glm.all.fdr.rolwin21.annot.EvSP.T4.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4)) {
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin21.annot.EvSP.T4$Gene,glm.all.fdr.rolwin21.annot.EvSP.T4[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"
    #assign minimum FDR values
    glm.all.fdr.rolwin21.annot.EvSP.T4.min <- cbind(glm.all.fdr.rolwin21.annot.EvSP.T4.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin21.annot.EvSP.T4.min) <- names(glm.all.fdr.rolwin21.annot.EvSP.T4)[c(5:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4))]
## Add gene names to the table
glm.all.fdr.rolwin21.annot.EvSP.T4.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin21.annot.EvSP.T4.min)
        	
## Save minimum FDR results for downstream use
write.table(glm.all.fdr.rolwin21.annot.EvSP.T4.min, file="rudflies_2023_redo.glm.all.rolwin21.minfdr.EvSP.T4.txt", quote = FALSE, sep = "\t", row.names = F)


### Redo above with missense SNPs in candidates only ###

## Merge to find missense SNPs in and around candidate genes
glm.all.rolwin21.annot.EvSP.T4.missense <- merge(EvSP.T4.cand, glm.all.rolwin21.annot[glm.all.rolwin21.annot$Consequence=="missense_variant",], by="Gene")

## Pull FDR columns only for EvSP.T4 candidate missense SNPs
glm.all.fdr.rolwin21.annot.EvSP.T4.missense <- unique(cbind(glm.all.rolwin21.annot.EvSP.T4.missense[,c(1:3)],glm.all.rolwin21.annot.EvSP.T4.missense$Consequence,select(glm.all.rolwin21.annot.EvSP.T4.missense,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate EvSP.T4 candidate minimum values
glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.missense)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.rolwin21.annot.EvSP.T4.missense$Gene,glm.all.fdr.rolwin21.annot.EvSP.T4.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min <- cbind(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min) <- names(glm.all.fdr.rolwin21.annot.EvSP.T4.missense)[c(5:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.missense))]
## Add gene names to the table
glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min)

## Save minimum FDR results from missense SNPs for downstream use
write.table(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min, file="rudflies_2023_redo.glm.all.rolwin21.minfdr.EvSP.T4.missense.txt", quote = FALSE, sep = "\t", row.names = F)


## Now we've found minimum FDR values for candidate gene and background gene lists, 
## let's run hypergeometric tests on the overlap of candidates and GLM outliers.

## run hypergeometric overlap test per GLM contrast, 21-SNP rolling windows
## First with all gene-related/adjacent SNPS
## Significance threshold FDR<0.05
phyper_all_EvSP_T4_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.EvSP.T4.min[glm.all.fdr.rolwin21.annot.EvSP.T4.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.EvSP.T4.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.min[glm.all.fdr.rolwin21.annot.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_EvSP_T4_contrasts05 <- append(phyper_all_EvSP_T4_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)    
}

## Create results summary table
phyper_EvSP.T4Cand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.EvSP.T4.min[,2:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.min)]),phyper_pval=phyper_all_EvSP_T4_contrasts05,phyper_p_sig=phyper_all_EvSP_T4_contrasts05 < 0.05,phyper_fdr=p.adjust(phyper_all_EvSP_T4_contrasts05, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_all_EvSP_T4_contrasts05, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_all_EvSP_T4_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.EvSP.T4.min[glm.all.fdr.rolwin21.annot.EvSP.T4.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.EvSP.T4.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.min[glm.all.fdr.rolwin21.annot.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_EvSP_T4_contrasts01 <- append(phyper_all_EvSP_T4_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_EvSP.T4Cand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.EvSP.T4.min[,2:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.min)]),phyper_pval=phyper_all_EvSP_T4_contrasts01,phyper_p_sig=phyper_all_EvSP_T4_contrasts01 < 0.05,phyper_fdr=p.adjust(phyper_all_EvSP_T4_contrasts01, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_all_EvSP_T4_contrasts01, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)


### Next with only missense SNPS ###

## Significance threshold FDR<0.05
phyper_missense_EvSP_T4_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min[glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[glm.all.fdr.rolwin21.annot.missense.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_EvSP_T4_contrasts05 <- append(phyper_missense_EvSP_T4_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_EvSP.T4Cand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min[,2:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min)]),phyper_pval=phyper_missense_EvSP_T4_contrasts05,phyper_p_sig=phyper_missense_EvSP_T4_contrasts05 < 0.05,phyper_fdr=p.adjust(phyper_missense_EvSP_T4_contrasts05, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_missense_EvSP_T4_contrasts05, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_missense_EvSP_T4_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min)) {
  q <- length(unique(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min[glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min[,1]))
  n <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.rolwin21.annot.missense.min[glm.all.fdr.rolwin21.annot.missense.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_EvSP_T4_contrasts01 <- append(phyper_missense_EvSP_T4_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_EvSP.T4Cand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min[,2:ncol(glm.all.fdr.rolwin21.annot.EvSP.T4.missense.min)]),phyper_pval=phyper_missense_EvSP_T4_contrasts01,phyper_p_sig=phyper_missense_EvSP_T4_contrasts01 < 0.05,phyper_fdr=p.adjust(phyper_missense_EvSP_T4_contrasts01, method = "fdr"),phyper_fdr_sig=p.adjust(phyper_missense_EvSP_T4_contrasts01, method = "fdr") < 0.05,q_all,m_all,k_all,n_all)

## Save results table for EvSP.T4T4idates vs GLM outliers (FDR < 0.05)
write.table(phyper_EvSP.T4Cand_FDR05, file="rudflies_2023_redo.phyper_allSNPs_EvSP.T4Cand_FDR05.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for EvSP.T4 candidates vs GLM outliers (FDR < 0.01)
write.table(phyper_EvSP.T4Cand_FDR01, file="rudflies_2023_redo.phyper_allSNPs_EvSP.T4Cand_FDR01.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for EvSP.T4 candidates vs GLM outliers (FDR < 0.05), missense only
write.table(phyper_missense_EvSP.T4Cand_FDR05, file="rudflies_2023_redo.phyper_missenseSNPs_EvSP.T4Cand_FDR05.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for EvSP.T4 candidates vs GLM outliers (FDR < 0.01), missense only
write.table(phyper_missense_EvSP.T4Cand_FDR01, file="rudflies_2023_redo.phyper_missenseSNPs_EvSP.T4Cand_FDR01.txt", quote = FALSE, sep = "\t", row.names = F)

### Which contrasts show significant overlap, prior to p-value correction

## EvSP.T4 candidates vs GLM outliers (FDR < 0.05), all SNPs
phyper_EvSP.T4Cand_FDR05[as.numeric(phyper_EvSP.T4Cand_FDR05[,2])<0.05,]
## EvSP.T4 candidates vs GLM outliers (FDR < 0.01), all SNPs
phyper_EvSP.T4Cand_FDR01[as.numeric(phyper_EvSP.T4Cand_FDR01[,2])<0.05,]
## EvSP.T4 candidates vs GLM outliers (FDR < 0.05), missense SNPs
phyper_missense_EvSP.T4Cand_FDR05[as.numeric(phyper_missense_EvSP.T4Cand_FDR05[,2])<0.05,]
## EvSP.T4 candidates vs GLM outliers (FDR < 0.01), missense SNPs
phyper_missense_EvSP.T4Cand_FDR01[as.numeric(phyper_missense_EvSP.T4Cand_FDR01[,2])<0.05,]


#############################################################
### Examine candidate gene p-values via Wilcoxon rank sum ###
#############################################################

## Hypergeometric overlap analysis relies on GLM outlier cutoff decisions, which can 
## impact results. Wilcoxon Rank Sum tests are non-parametric and rely on rank order of
## all genes in the dataset. Therefore, this test is more sensitive to more modest 
## AF changes across a broader gene network or pathway.

### Remake annotation stat summary lists for pvalues instead of fdr values ###

## Pull gene and SNP consequence columns for EvSP.T4 candidates
glm.all.pval.rolwin21.annot.EvSP.T4 <- cbind(glm.all.rolwin21.annot.EvSP.T4[,c(1:3)],Consequence=glm.all.rolwin21.annot.EvSP.T4$Consequence)

## Remove duplicates
glm.all.pval.rolwin21.annot.EvSP.T4 <- unique(glm.all.pval.rolwin21.annot.EvSP.T4)

## Finally, add p-values 
glm.all.pval.rolwin21.annot.EvSP.T4 <- merge(glm.all.pval.rolwin21.annot.EvSP.T4, glm.all.rolwin21.pval, by=c("CHROM","POS"))

## Loop through each GLM-results pval column to calculate EvSP.T4 candidate minimum values
glm.all.pval.rolwin21.annot.EvSP.T4.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin21.annot.EvSP.T4)) {
  tryCatch({
  	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin21.annot.EvSP.T4$Gene,glm.all.pval.rolwin21.annot.EvSP.T4[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum p-values    
    glm.all.pval.rolwin21.annot.EvSP.T4.min <- cbind(glm.all.pval.rolwin21.annot.EvSP.T4.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin21.annot.EvSP.T4.min) <- names(glm.all.pval.rolwin21.annot.EvSP.T4)[c(5:ncol(glm.all.pval.rolwin21.annot.EvSP.T4))]
## Add gene names to the table
glm.all.pval.rolwin21.annot.EvSP.T4.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin21.annot.EvSP.T4.min)
        	
## Save minimum pval results for downstream use
write.table(glm.all.pval.rolwin21.annot.EvSP.T4.min, file="rudflies_2023_redo.glm.all.rolwin21.minpval.EvSP.T4.txt", quote = FALSE, sep = "\t", row.names = F)
        	
        	
## Pull pval columns only for all genes/features, filtered for missense SNPs only
glm.all.pval.rolwin21.annot.EvSP.T4.missense <- unique(cbind(glm.all.rolwin21.annot.EvSP.T4.missense[,c(1:3)],glm.all.rolwin21.annot.EvSP.T4.missense$Consequence,select(glm.all.rolwin21.annot.EvSP.T4.missense,contains("rolwin21"),-contains("fdr"),-contains("logp"))))

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.rolwin21.annot.EvSP.T4.missense.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.rolwin21.annot.EvSP.T4.missense)) { #nrow(haf.freq.cand)
  tryCatch({
	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.rolwin21.annot.EvSP.T4.missense$Gene,glm.all.pval.rolwin21.annot.EvSP.T4.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum FDR values    
    glm.all.pval.rolwin21.annot.EvSP.T4.missense.min <- cbind(glm.all.pval.rolwin21.annot.EvSP.T4.missense.min,
    	temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.rolwin21.annot.EvSP.T4.missense.min) <- names(glm.all.pval.rolwin21.annot.EvSP.T4.missense)[c(5:ncol(glm.all.pval.rolwin21.annot.EvSP.T4.missense))]
## Add gene names to the table
glm.all.pval.rolwin21.annot.EvSP.T4.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.rolwin21.annot.EvSP.T4.missense.min)


## Save minimum pval results from missense SNPs for downstream use
write.table(glm.all.pval.rolwin21.annot.EvSP.T4.missense.min, file="rudflies_2023_redo.glm.all.rolwin21.minpval.EvSP.T4.missense.txt", quote = FALSE, sep = "\t", row.names = F)
        	


## Time to run the Wicoxon Rank Sum analyses. To do so, we will select one matched 
## background gene for every EvSP.T4 candidate gene. This matched set selection is 
## based on the following criteria:
## 	 -Gene must be similar length from start to stop (+/-25%)
## 	 -Gene must have similar number of SNPs (+/-25%)
## 	 -Gene must be the same type (protein-coding, ncRNA, etc)
## 	 -Gene must be on the same chromosome

## After all matching genes are discovered, we randomly select one match for each 
## candidate gene, then perform the Wilcoxon rank sum test between the candidate and 
## background genes to see if candidates show significantly higher p-value ranks (lower 
## p-values), or greater parallel allele frequency differences in the respective GLM 
## contrast. This process is repeated for 1K iterations of matched set selection to 
## mitigate any potential random sampling bias. 


##Find gene characteristics of candidates from GFF for matching
glm.all.pval.rolwin21.annot.EvSP.T4.min.gff <- merge(snp.gff.overlap, glm.all.pval.rolwin21.annot.EvSP.T4.min, by.x="ID", by.y="Gene")

##Find gene characteristics of background genes from GFF for matching
glm.all.pval.rolwin21.annot.min.gff <- merge(snp.gff.overlap, glm.all.pval.rolwin21.annot.min, by.x="ID", by.y="Gene")

## Pull EvSP.T4 candidate gene IDs
EvSP.T4_IDs <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff$ID

## Remove EvSP.T4 candidates from potential background gene list
glm.all.pval.rolwin21.annot.bg.min.gff <- glm.all.pval.rolwin21.annot.min.gff[!(glm.all.pval.rolwin21.annot.min.gff$ID %in% EvSP.T4_IDs),]

## Create empty results table with 1K rows (for each iteration) and one column for 
## each GLM contrast.
wilcox.rolwin21.EvSP.T4.p <- data.frame(matrix(NA, nrow = 1000, ncol = ncol(glm.all.pval.rolwin21.annot.EvSP.T4.min.gff)-ncol(snp.gff.overlap)))
## Modify column names to reflect GLM contrasts
names(wilcox.rolwin21.EvSP.T4.p) <- names(glm.all.pval.rolwin21.annot.EvSP.T4.min.gff)[c(9:ncol(glm.all.pval.rolwin21.annot.EvSP.T4.min.gff))]

###First, determine how many matches so we can exclude those without any
bg.samp.counts <- c()
  #This loop finds a random set of matched background genes
  for(i in c(1:nrow(glm.all.pval.rolwin21.annot.EvSP.T4.min.gff))) {
  	#Select matched lists based on following criteria
	length <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[i,]$LENGTH
	count <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[i,]$SNP_COUNT
	type <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[i,]$TYPE
	chrom <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[i,]$CHROM
	#Pull one match per candidate gene
	bg.samp.temp <- nrow(glm.all.pval.rolwin21.annot.bg.min.gff %>%
	 filter( 
	    TYPE==type,
		CHROM==chrom,
		LENGTH < length*1.25,
		LENGTH > length*0.75,
		SNP_COUNT < count*1.25,
		SNP_COUNT > count*0.75))
	#append to background list
	bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
  }

#generate candidate gene filtering tables
for(k in 1:1000) { #set number of iterations
  #Establish background non-candidate set
  bg.samp.list <- c()
  #This loop finds a random set of matched background genes
  for(i in as.numeric(rownames(glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[bg.samp.counts[,1]>4,]))) {
  	#Select matched lists based on following criteria
	length <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[i,]$LENGTH
	count <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[i,]$SNP_COUNT
	type <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[i,]$TYPE
	chrom <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[i,]$CHROM
	#Pull one match per candidate gene
	bg.samp.temp <- sample_n(glm.all.pval.rolwin21.annot.bg.min.gff %>%
	 filter( 
	    TYPE==type,
		CHROM==chrom,
		LENGTH < length*1.25,
		LENGTH > length*0.75,
		SNP_COUNT < count*1.25,
		SNP_COUNT > count*0.75
	  ),1)
	#append to background list
	bg.samp.list <- rbind(bg.samp.list,bg.samp.temp)
  }
  #This loop will run Wilcoxon Rank Sum tests for each of "j" GLM contrasts
  for(j in 1:(ncol(glm.all.pval.rolwin21.annot.EvSP.T4.min.gff)-ncol(snp.gff.overlap))) { #set number of iterations to the number of contrasts
    #Establish Rank Sum test sets for p values
    #candidate set
    candidates <- glm.all.pval.rolwin21.annot.EvSP.T4.min.gff[,8+j]  
    #background non-candidate set
    bg.samp <- as.numeric(bg.samp.list[,8+j])
    candidates.samp <- as.numeric(candidates)
    #this will subsample the SNPs contained in one list to the number of SNPs in the shorter list
    if(length(candidates.samp) > length(bg.samp)){
      candidates.samp <- sample(candidates.samp, length(bg.samp))
    } else {
      bg.samp <- sample(bg.samp, length(candidates.samp))
    }
    #finally, run the tests and append results
    wilcox.rolwin21.EvSP.T4.p[k,j] <- wilcox.test(candidates.samp, bg.samp, alternative = "less")$p.value
}
}

## Save all results
# Wilcoxon rank sum p-values for all iterations
write.table(wilcox.rolwin21.EvSP.T4.p, file="rudflies_2023_redo.glm.all.wilcox.rolwin21.EvSP.T4.p.txt", quote = FALSE, sep = "\t", row.names = F)

# Median Wilcoxon rank sum p-values across iterations
write.table(colMedians(as.matrix(wilcox.rolwin21.EvSP.T4.p)), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.EvSP.T4.p.medians.txt", quote = FALSE, sep = "\t", row.names = F)

# Mean Wilcoxon rank sum p-values across iterations
write.table(colMeans(as.matrix(wilcox.rolwin21.EvSP.T4.p)), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.EvSP.T4.p.means.txt", quote = FALSE, sep = "\t", row.names = F)

## Calculate the percent of p-values significant at a raw p<0.05
colMedians(as.matrix(wilcox.rolwin21.EvSP.T4.p))
colMeans(as.matrix(wilcox.rolwin21.EvSP.T4.p))
wilcox.rolwin21.EvSP.T4.contrast <- c()
wilcox.rolwin21.EvSP.T4.contrast.p05 <- wilcox.rolwin21.EvSP.T4.p<0.05
for(l in 1:ncol(wilcox.rolwin21.EvSP.T4.p)) {
   wilcox.rolwin21.EvSP.T4.contrast <- append(wilcox.rolwin21.EvSP.T4.contrast,sum(wilcox.rolwin21.EvSP.T4.contrast.p05[,l], na.rm=TRUE)*0.1)
}
cbind(contrast=names(wilcox.rolwin21.EvSP.T4.p),perc_sig=wilcox.rolwin21.EvSP.T4.contrast)

# Save percent of significant Wilcoxon rank sum p-values across iterations
write.table(cbind(contrast=names(wilcox.rolwin21.EvSP.T4.p),perc_sig=wilcox.rolwin21.EvSP.T4.contrast), file="rudflies_2023_redo.glm.all.wilcox.rolwin21.EvSP.T4.p.percsig.txt", quote = FALSE, sep = "\t", row.names = F)



##########################################################################################
###						  Pleitropic cost of spinosyn adaptation 					   ###
##########################################################################################

### In this section we will run various enrichment tests that evaluate whether alleles 
### that we predict to be associated with spinosyn-adaptation are under negative 
### directional selection in natural conditions. These lists of predicted adaptive alleles 
### come from and various PAvE.T1 and PAvE population contrasts. 


#### Create p-value correlation tables to further examine the cost of adaptation ####

### As a first pass, these correlations do not look at directionality, just significance 
### to see if putatively spinosyn-adaptive loci are show similar levels of significance in 
### temporal contrasts. These are pretty broad tests, but informative.

## Define temporal contrasts
temporal_contrast <- c("PA.T1vT2","PA.T1vT3","PA.T1vT4","PA.T2vT3","PA.T2vT4","PA.T3vT4","E.T1vT2","E.T1vT3","E.T1vT4","E.T2vT3","E.T2vT4","E.T3vT4")

pval_corr_table <- c()

## Use PAvE.T1 as the reference adaptive contrast
for(x in temporal_contrast) { #cycle through all temporal contrast
		temp_row <- t(cor.test(as.numeric(glm.all$PAvE.T1), as.numeric(glm.all[[x]]), method="pearson", alternative="greater")[c(1,3,4)])
		pval_corr_table <- rbind(pval_corr_table, cbind(set1="PAvE.T1",set2=x,temp_row))
}

## Use PAvE as the reference adaptive contrast
for(x in temporal_contrast) { #cycle through all temporal contrast
		temp_row <- t(cor.test(as.numeric(glm.all$PAvE), as.numeric(glm.all[[x]]), method="pearson", alternative="greater")[c(1,3,4)])
		pval_corr_table <- rbind(pval_corr_table, cbind(set1="PAvE",set2=x,temp_row))
}

## Save table
write.table(pval_corr_table, file="rudflies_2023_redo.spinosyn_vs_temporal_pval_corr_table.txt", sep = "\t", quote = FALSE, row.names = F)


### From now on, we'll consider directionality in addition to significance level ###

## Let's calculate mean frequencies, then calculate pairwise differences
## We'll use these values to assign positive or negative signs to -log10p values
## The downstream values will reflect both the magnitude and direction of AF change
freq_diff_conv <- cbind(haf.sites.filt[,c(1:2)], #locus info
	PAvE.T1.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1"]), #PAvE.T1
	PAvE.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E"]), #PAvE
	PA.T1vT2.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="2"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1"]), #PA.T1vT2
	PA.T1vT3.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="3"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1"]), #PA.T1vT2
	PA.T1vT4.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="4"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1"]), #PA.T1vT2
	PA.T2vT3.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="3"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="2"]), #PA.T1vT2
	PA.T2vT4.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="4"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="2"]), #PA.T1vT2
	PA.T3vT4.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="4"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="3"]), #PA.T1vT2
	E.T1vT2.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="2"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1"]), #E.T1vT2
	E.T1vT3.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="3"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1"]), #E.T1vT2
	E.T1vT4.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="4"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1"]), #E.T1vT2
	E.T2vT3.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="3"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="2"]), #E.T1vT2
	E.T2vT4.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="4"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="2"]), #E.T1vT2
	E.T3vT4.diff=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="4"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="3"])) #E.T1vT2

## Convert all frequency differences to +1/-1 to reflect positive or negative AF diff
freq_diff_sign <- cbind(freq_diff_conv[,c(1:2)],freq_diff_conv[,c(3:ncol(freq_diff_conv))]/abs(freq_diff_conv[,c(3:ncol(freq_diff_conv))]))

## Join PAvE and temporal GLM results (-log10p) with new freq_diff_sign table
freq_diff_sign_temp <- merge(freq_diff_sign, cbind(glm.all[,c(1,2)], 
	PAvE.T1.logp=glm.all$PAvE.T1.logp, 
	PAvE.logp=glm.all$PAvE.logp, 
	PA.T1vT2.logp=glm.all$PA.T1vT2.logp, 
	PA.T1vT3.logp=glm.all$PA.T1vT3.logp, 
	PA.T1vT4.logp=glm.all$PA.T1vT4.logp, 
	PA.T2vT3.logp=glm.all$PA.T2vT3.logp, 
	PA.T2vT4.logp=glm.all$PA.T2vT4.logp, 
	PA.T3vT4.logp=glm.all$PA.T3vT4.logp, 
	E.T1vT2.logp=glm.all$E.T1vT2.logp, 
	E.T1vT3.logp=glm.all$E.T1vT3.logp, 
	E.T1vT4.logp=glm.all$E.T1vT4.logp, 
	E.T2vT3.logp=glm.all$E.T2vT3.logp, 
	E.T2vT4.logp=glm.all$E.T2vT4.logp, 
	E.T3vT4.logp=glm.all$E.T3vT4.logp), 
	by=c("CHROM","POS"))
	

## Make new table to store sign-corrected -log10(p) values
freq_diff_sign_logp <- freq_diff_sign_temp[,c(1,2)]

## Conduct sign-correction through multiplication
for(i in c(3:16)) { #start after loci columns and append FDR cols at the end of the table
  freq_diff_sign_logp <- cbind(freq_diff_sign_logp, freq_diff_sign_temp[,i] * freq_diff_sign_temp[,i+14])
  colnames(freq_diff_sign_logp)[i] <- paste(names(freq_diff_sign_temp)[i+14], ".sign", sep="")
}

## Sort by locus
freq_diff_sign_logp <- freq_diff_sign_logp[order(freq_diff_sign_logp[,1], freq_diff_sign_logp[,2]), ]

## Now calculate AF means across treatments for additional filtering of table
haf.freq.T1filt.af_mean <- cbind(haf.sites.T1filt,af_mean=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix=="E" | haf.meta.T1filt$treat.fix=="S" | haf.meta.T1filt$treat.fix=="PA"]))

## Merge new means with GLM/freq_diff table
freq_diff_sign_logp <- merge(freq_diff_sign_logp, haf.freq.T1filt.af_mean, by=c("CHROM","POS"))

## Filter if AF means are less than 0.15 or greater than 0.85.
## This preserves ranking patterns for biologically relevant loci while removing those 
## that appear to have major treatment differences due to REF or ALT alleles being rare.
freq_diff_sign_logp <- freq_diff_sign_logp[freq_diff_sign_logp$af_mean>0.15 & freq_diff_sign_logp$af_mean<0.85,]


### Directional stacked bar plots for top 10K PAvE.T1 SNPs, ranked by GLM FDR ###

## First, find top 10K SNPs
PAvE.T1.10K <- head(unique(glm.all[order(glm.all$PAvE.T1.fdr),]),10000)

## Merge to retrieve sign info
freq_diff_sign_PAvE_T1_top10K <- merge(PAvE.T1.10K[,c(1:2)], freq_diff_sign, by=c("CHROM","POS"))

## Here we recalibrate sign-correction factor for all temporal comparison to match 
## direction of PAvE.T1 contrast for easier comparison
freq_diff_sign_PAvE_T1_top10K_rel <- freq_diff_sign_PAvE_T1_top10K[,c(5:ncol(freq_diff_sign_PAvE_T1_top10K))] * freq_diff_sign_PAvE_T1_top10K[,3]

## Get counts of positive (PA-biased) & negative (E-biased) AF directions per contrast
freq_diff_sign_counts_PAvE_T1_top10K_rel <- lapply(freq_diff_sign_PAvE_T1_top10K_rel, table)

## Make data frame
freq_diff_sign_counts_PAvE_T1_top10K_rel <- do.call(rbind, lapply(freq_diff_sign_counts_PAvE_T1_top10K_rel, as.data.frame))

## Make contrast field to keep track and filter
freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast <- row.names(freq_diff_sign_counts_PAvE_T1_top10K_rel)

## Remove number suffixes so like-contrasts can be plotted together
freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast <- str_sub(freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast, end = -8)

##Give Var1 field informative values and header
freq_diff_sign_counts_PAvE_T1_top10K_rel$Var1 <- gsub("-1","E-biased",freq_diff_sign_counts_PAvE_T1_top10K_rel$Var1)
freq_diff_sign_counts_PAvE_T1_top10K_rel$Var1 <- gsub("1","PA-biased",freq_diff_sign_counts_PAvE_T1_top10K_rel$Var1)
names(freq_diff_sign_counts_PAvE_T1_top10K_rel)[1] <- "PAvE.T1_contrast"


### Plot all p-value correlations to check for cost of adaptation
pdf(file = "rudflies_2023_redo.PAvE.T1_vs_temporal.signed.stackbarplot_top10K.pdf", width=12, height=4)
	ggplot(freq_diff_sign_counts_PAvE_T1_top10K_rel) +
  		geom_bar(aes(x = contrast, y = Freq, fill = PAvE.T1_contrast), 
           position = "stack", stat = "identity") +
      	geom_hline(yintercept=5000,color="black",linetype="dashed",linewidth=.25) +
      	ggtitle("Direction of temporal AF change for top 10k outlier SNPs in PAvE.T1 contrast") +
      	theme_classic()    
dev.off()

## Save plot as object for multi-panel plot
freq_diff_sign_counts_PAvE_T1_top10K_rel_bar <-	ggplot(freq_diff_sign_counts_PAvE_T1_top10K_rel) +
  		geom_bar(aes(x = contrast, y = Freq, fill = PAvE.T1_contrast), 
           position = "stack", stat = "identity") +
      	geom_hline(yintercept=5000,color="black",linetype="dashed",linewidth=.25) +
      	ggtitle("Direction of temporal AF change for top 10k outlier SNPs in PAvE.T1 contrast") +
      	theme_classic()  
      	
## Save a version with only T2vT4 contrasts for multi-panel plot
freq_diff_sign_counts_PAvE_T1_top10K_rel_T2vT4 <- freq_diff_sign_counts_PAvE_T1_top10K_rel[freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast=="PA.T2vT4" | freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast=="E.T2vT4",]

freq_diff_sign_counts_PAvE_T1_top10K_rel_T2vT4_bar <-	ggplot(freq_diff_sign_counts_PAvE_T1_top10K_rel_T2vT4) +
  		geom_bar(aes(x = contrast, y = Freq, fill = PAvE.T1_contrast), 
           position = "stack", stat = "identity") +
      	geom_hline(yintercept=5000,color="black",linetype="dashed",linewidth=.25) +
      	ggtitle("Top 10k outliers (PAvE.T1)") +
      	theme_classic()  

### Statistical test associated with these stacked barplots would be binomial sign test
## This will test if there are more E-biased AF differences than expected than 50/50

## list contrasts to test
temporal_contrast <- c("PA.T1vT2","PA.T1vT3","PA.T1vT4","PA.T2vT3","PA.T2vT4","PA.T3vT4","E.T1vT2","E.T1vT3","E.T1vT4","E.T2vT3","E.T2vT4","E.T3vT4")

## Initialize results table
binom_res_10kSNP_PAvE.T1 <- c()

## loop through all contrasts
for (contrast in temporal_contrast) { #cycle through all contrasts
	Ebias <- freq_diff_sign_counts_PAvE_T1_top10K_rel[freq_diff_sign_counts_PAvE_T1_top10K_rel$PAvE.T1_contrast=="E-biased" & freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast==contrast,]$Freq #save E-biased counts
	PAbias <- freq_diff_sign_counts_PAvE_T1_top10K_rel[freq_diff_sign_counts_PAvE_T1_top10K_rel$PAvE.T1_contrast=="PA-biased" & freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast==contrast,]$Freq #save PA-biased counts
	##Run binomial test
	binom_res_10kSNP_PAvE.T1 <- rbind(binom_res_10kSNP_PAvE.T1,cbind(contrast,(data.frame(t(binom.test(x = Ebias, n = Ebias + PAbias, p = 0.5, alternative = "greater")[c(1,3,5)])))))
}

## Save table
fwrite(data.frame(binom_res_10kSNP_PAvE.T1), file="rudflies_2023_redo.binom_res_10kSNP_PAvE.T1.txt")


### Directional stacked bar plots for all SNPs ###

## Here we recalibrate sign-correction factor for all temporal comparison to match 
## direction of PAvE.T1 contrast for easier comparison
freq_diff_sign_rel <- freq_diff_sign[,c(5:ncol(freq_diff_sign))] * freq_diff_sign[,3]

## Get counts of positive (PA-biased) & negative (E-biased) AF directions per contrast
freq_diff_sign_counts_rel <- lapply(freq_diff_sign_rel, table)

## Make data frame
freq_diff_sign_counts_rel <- do.call(rbind, lapply(freq_diff_sign_counts_rel, as.data.frame))

## Make contrast field to keep track and filter
freq_diff_sign_counts_rel$contrast <- row.names(freq_diff_sign_counts_rel)

## Remove number suffixes so like-contrasts can be plotted together
freq_diff_sign_counts_rel$contrast <- str_sub(freq_diff_sign_counts_rel$contrast, end = -8)

##Give Var1 field informative values and header
freq_diff_sign_counts_rel$Var1 <- gsub("-1","E-biased",freq_diff_sign_counts_rel$Var1)
freq_diff_sign_counts_rel$Var1 <- gsub("1","PA-biased",freq_diff_sign_counts_rel$Var1)
names(freq_diff_sign_counts_rel)[1] <- "PAvE.T1_contrast"

### Plot all p-value correlations to check for cost of adaptation
pdf(file = "rudflies_2023_redo.PAvE.T1_vs_temporal.signed.stackbarplot_allSNPs.pdf", width=12, height=4)
	ggplot(freq_diff_sign_counts_rel) +
  		geom_bar(aes(x = contrast, y = Freq, fill = PAvE.T1_contrast), 
           position = "stack", stat = "identity") +
      	geom_hline(yintercept=nrow(freq_diff_sign)/2, color="black", linetype="dashed", linewidth=.25) +
      	ggtitle("Direction of temporal AF change for all SNPs") +
      	theme_classic()    
dev.off()

## Save plot as object for multi-panel plot
freq_diff_sign_counts_rel_bar <- ggplot(freq_diff_sign_counts_rel) +
  		geom_bar(aes(x = contrast, y = Freq, fill = PAvE.T1_contrast), 
           position = "stack", stat = "identity") +
      	geom_hline(yintercept=nrow(freq_diff_sign)/2,color="black",linetype="dashed",linewidth=.25) +
      	ggtitle("Direction of temporal AF change for all SNPs") +
      	theme_classic()    
      	
## Save a version with only T2vT4 contrasts for multi-panel plot
freq_diff_sign_counts_rel_T2vT4 <- freq_diff_sign_counts_rel[freq_diff_sign_counts_rel$contrast=="PA.T2vT4" | freq_diff_sign_counts_rel$contrast=="E.T2vT4",]

freq_diff_sign_counts_rel_T2vT4_bar <-	ggplot(freq_diff_sign_counts_rel_T2vT4) +
  		geom_bar(aes(x = contrast, y = Freq, fill = PAvE.T1_contrast), 
           position = "stack", stat = "identity") +
      	geom_hline(yintercept=nrow(freq_diff_sign)/2, color="black", linetype="dashed", linewidth=.25) +
      	ggtitle("All SNPs") +
      	theme_classic()  

### Plot all p-value correlations to check for cost of adaptation
pdf(file = "rudflies_2023_redo.PAvE.T1_vs_temporal.signed.multibarplot.pdf", width=12, height=8)
	ggarrange(freq_diff_sign_counts_PAvE_T1_top10K_rel_bar, freq_diff_sign_counts_rel_bar,
		ncol = 1, nrow = 2)
dev.off()


### Plot all p-value correlations to check for cost of adaptation
pdf(file = "rudflies_2023_redo.PAvE.T1_vs_T2vT4.signed.multibarplot.pdf", width=8, height=4)
	ggarrange(freq_diff_sign_counts_PAvE_T1_top10K_rel_T2vT4_bar, freq_diff_sign_counts_rel_T2vT4_bar,
		ncol = 2, nrow = 1)
dev.off()


### Statistical test associated with these stacked barplots would be binomial sign test
## This will test if there are more E-biased AF differences than expected than 50/50

## Initialize results table
binom_res_allSNP_PAvE.T1 <- c()

## loop through all contrasts
for (contrast in temporal_contrast) { #cycle through all contrasts
	Ebias <- freq_diff_sign_counts_PAvE_T1_top10K_rel[freq_diff_sign_counts_PAvE_T1_top10K_rel$PAvE.T1_contrast=="E-biased" & freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast==contrast,]$Freq #save E-biased counts
	PAbias <- freq_diff_sign_counts_PAvE_T1_top10K_rel[freq_diff_sign_counts_PAvE_T1_top10K_rel$PAvE.T1_contrast=="PA-biased" & freq_diff_sign_counts_PAvE_T1_top10K_rel$contrast==contrast,]$Freq #save PA-biased counts
	##Run binomial test
	binom_res_allSNP_PAvE.T1 <- rbind(binom_res_allSNP_PAvE.T1,cbind(contrast,(data.frame(t(binom.test(x = Ebias, n = Ebias + PAbias, p = 0.5, alternative = "greater")[c(1,3,5)])))))
}

## Save table
fwrite(data.frame(binom_res_allSNP_PAvE.T1), file ="rudflies_2023_redo.binom_res_allSNP_PAvE.T1.txt")


### Directional histograms for top 10k PAvE.T1 SNPs ###

#PAvE.T1.10K <- head(unique(glm.all[order(glm.all$PAvE.T1.fdr),]),10000)
freq_diff_sign_logp_PAvE_T1_top10K <- merge(PAvE.T1.10K[,c(1:2)], freq_diff_sign_logp, by=c("CHROM","POS"))

## Loop over temporal contrasts to make histograms of top 10k sign-corrected -log10p vals
for (x in c(5:16)) { #cycle through all contrasts
	templist <- freq_diff_sign_logp_PAvE_T1_top10K[,c(1,2,x)]
	names(templist)[3] <- "value"
	assign(paste("PAvE.T1_vs_",temporal_contrast[x-4],"_signed_hist_top10K",sep=""), 
		ggplot(templist, 
	  		aes(value)) + 
	  		geom_histogram(bins = 100) +
	  		xlab(paste(temporal_contrast[x-4]," p-values",sep="")) +
	  		theme_classic())
}


### Plot all p-value correlations to check for cost of adaptation
pdf(file = "rudflies_2023_redo.PAvE.T1_vs_temporal.signed_logp.hist_top10K.pdf", width=8, height=10)
	ggarrange(PAvE.T1_vs_PA.T1vT2_signed_hist_top10K, 
		PAvE.T1_vs_E.T1vT2_signed_hist_top10K, 
		PAvE.T1_vs_PA.T1vT3_signed_hist_top10K, 
		PAvE.T1_vs_E.T1vT3_signed_hist_top10K, 
		PAvE.T1_vs_PA.T1vT4_signed_hist_top10K, 
		PAvE.T1_vs_E.T1vT4_signed_hist_top10K, 
		PAvE.T1_vs_PA.T2vT3_signed_hist_top10K, 
		PAvE.T1_vs_E.T2vT3_signed_hist_top10K, 
		PAvE.T1_vs_PA.T2vT4_signed_hist_top10K, 
		PAvE.T1_vs_E.T2vT4_signed_hist_top10K, 
		PAvE.T1_vs_PA.T3vT4_signed_hist_top10K, 
		PAvE.T1_vs_E.T3vT4_signed_hist_top10K, 
        ncol = 2, nrow = 6)
dev.off()


## Run t-tests to see if mean differences are negative (E-biased) for all contrasts
## Specifically, we will test if the mean of the focal sets is less than 0
ttest_res_10kSNP <- c()
for (x in c(5:16)) { #cycle through all contrasts
	contrast <- temporal_contrast[x-4]
	templist <- freq_diff_sign_logp_PAvE_T1_top10K[,c(1,2,x)]
	names(templist)[3] <- "value"
	ttest_res_10kSNP <- rbind(ttest_res_10kSNP, cbind(contrast, t(t.test(templist$value, mu = 0, alternative = "less")[c(1,3,5)])))
}

## Save
write.table(ttest_res_10KPAvE_T1, file="rudflies_2023_redo.ttest_res_10KPAvE_T1.txt", sep = "\t", quote = FALSE, row.names = F)


### Directional histograms for top 1k PAvE.T1 SNPs ###

#PAvE.T1.1K <- head(unique(glm.all[order(glm.all$PAvE.T1.fdr),]),10000)
freq_diff_sign_logp_PAvE_T1_top1K <- merge(PAvE.T1.1K[,c(1:2)], freq_diff_sign_logp, by=c("CHROM","POS"))

## Loop over temporal contrasts to make histograms of top 1k sign-corrected -log10p vals
for (x in c(5:16)) { #cycle through all contrasts
	templist <- freq_diff_sign_logp_PAvE_T1_top1K[,c(1,2,x)]
	names(templist)[3] <- "value"
	assign(paste("PAvE.T1_vs_",temporal_contrast[x-4],"_signed_hist_top1K",sep=""), 
		ggplot(templist, 
	  		aes(value)) + 
	  		geom_histogram(bins = 100) +
	  		xlab(paste(temporal_contrast[x-4]," p-values",sep="")) +
	  		theme_classic())
}


### Plot all p-value correlations to check for cost of adaptation
pdf(file = "rudflies_2023_redo.PAvE.T1_vs_temporal.signed_logp.hist_top1K.pdf", width=8, height=10)
	ggarrange(PAvE.T1_vs_PA.T1vT2_signed_hist_top1K, 
		PAvE.T1_vs_E.T1vT2_signed_hist_top1K, 
		PAvE.T1_vs_PA.T1vT3_signed_hist_top1K, 
		PAvE.T1_vs_E.T1vT3_signed_hist_top1K, 
		PAvE.T1_vs_PA.T1vT4_signed_hist_top1K, 
		PAvE.T1_vs_E.T1vT4_signed_hist_top1K, 
		PAvE.T1_vs_PA.T2vT3_signed_hist_top1K, 
		PAvE.T1_vs_E.T2vT3_signed_hist_top1K, 
		PAvE.T1_vs_PA.T2vT4_signed_hist_top1K, 
		PAvE.T1_vs_E.T2vT4_signed_hist_top1K, 
		PAvE.T1_vs_PA.T3vT4_signed_hist_top1K, 
		PAvE.T1_vs_E.T3vT4_signed_hist_top1K, 
        ncol = 2, nrow = 6)
dev.off()

## Run t-tests to see if mean differences are more negative (E-biased) than matched sets 
## for all contrasts. 

## Matched set selection is more robust than checking for significantly negative values, 
## since many loci across the genome, unrelated to spinosyn adaption, may be experiencing  
## negative directional selection in response to outdoor adaption. Unfortunately, matched  
## set selection on 10k SNPs for many iterarations for each contrast is memory and time  
## intensive, so here we'll focus on only the top 1k SNPs from the PAvE.T1 contrast. While  
## this is a smaller set than we tested above, these SNPs show a consistently stronger  
## signal of adaptation in PA populations, therefore the tests of pleiotropic costs of  
## spinosyn adaption become even more robust. 

## Find top 1K SNPs
PAvE.T1.1K <- head(unique(glm.all[order(glm.all$PAvE.T1.fdr),]),1000)

## Merge to retrieve sign-corrected -log10(p) values
freq_diff_sign_logp_PAvE_T1_top1K <- merge(PAvE.T1.1K[,c(1:2)], freq_diff_sign_logp, by=c("CHROM","POS"))

## Find sign-corrected -log10(p) values for top 1k PAvE.T1 SNPs
freq_diff_vep_sig <- merge(freq_diff_sign_logp_PAvE_T1_top1K, vep_priority, by=c("CHROM","POS"))
## Find sign-corrected -log10(p) values for non-top 1k PAvE.T1 SNPs
freq_diff_vep <- merge(freq_diff_sign_logp, vep_priority, by=c("CHROM","POS"))
freq_diff_vep_nonsig <- anti_join(freq_diff_vep, freq_diff_vep_sig, by = c("CHROM","POS"))

## Prior to running t-tests, determine how many matches so we can exclude those with <5
bg.samp.counts <- c()
#This loop finds a random set of matched background genes
for(j in c(1:nrow(freq_diff_vep_sig))) {
	#Select matched lists based on following criteria
	type <- freq_diff_vep_sig[j,]$Consequence
	chrom <- freq_diff_vep_sig[j,]$CHROM
	found <- freq_diff_vep_sig[j,]$E
	pos <- freq_diff_vep_sig[j,]$POS
	#Pull one match per candidate gene
	bg.samp.temp <- nrow(freq_diff_vep_nonsig %>%
		filter( 
		Consequence==type,
		CHROM==chrom,
		((E < found*1.25) & (E > found*.75)),
		abs(POS - pos) > 50000))
	#append to background list
	bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
}

## Initialize results table
ttest_res_1kSNPs_vs_match <- c()

## Run t-tests for all temporal contrast sign-corrected -log10(p) values
for(i in c(5:16)) { #cycle through all contrasts
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(freq_diff_vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- freq_diff_vep_sig[j,]$Consequence
			chrom <- freq_diff_vep_sig[j,]$CHROM
			found <- freq_diff_vep_sig[j,]$E
			pos <- freq_diff_vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(freq_diff_vep_nonsig %>%
		 		filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000),1)
			#append to background list
			bg.samp.list <- rbind(bg.samp.list,bg.samp.temp)
	  	}
  		
  		## Correct sign of delta AF from left-out sample using sign of left-in samples
  		## This is done by multiplying positive or negative 1 
		focal.afdiff <- freq_diff_vep_sig[,i] * (abs(freq_diff_vep_sig[,3])/freq_diff_vep_sig[,3])
		bg.afdiff <- bg.samp.list[,i] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		ttest_res_1kSNPs_vs_match <- rbind(ttest_res_1kSNPs_vs_match,cbind(contrast=temporal_contrast[i-4],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "less")[c(1,2,3,7)])))
	}
}

#########################################################################################
####### Test for costs using parallel AF differences to filter temporal contrasts ####### 
#########################################################################################

### New approach to look at cost of spinosyn adaptation. Since signal of temporal 
### selection appears weak, particularly in E populations (with very few outliers) 
### according to GLM results, we will instead use an parallel allele frequency difference 
### filtering approach. First, we will consider SNPs subject to parallel seasonal/temporal 
### adaptation if they show >0.03 AF difference in the same direction in all cages per 
### treatment. For a slightly less conservative approach, we will alternatively consider 
### SNPs subject to temporal adaptation if they show the same minimum AF difference in all 
### but one cage per treatment. 

### This will actually allow us to compare these temporal outlier SNPs lists with similar 
### lists generated from PAvE.T1 contrasts, which contain SNPs putatively subject to 
### parallel adaption. The lists will be compared through hypergeometric overlap tests, 
### and we can run separate tests for temporal SNP lists that are either parallel in the 
### E-biased or PA-biased direction.


############ Look for concordance of temporal selection in all PA cages ############

## Let's calculate mean frequencies, then calculate pairwise differences
## We'll use these values to assign positive or negative signs to -log10p values
## The resulting values will reflect both the magnitude and direction of AF change

PA_cages <- unique(haf.meta.filt[haf.meta.filt$treat.fix=="PA",]$cage.fix)
#[1] "2"  "6"  "10" "14" "20" "24" "31" "36" "40" "44"

for(x in PA_cages) { #cycle through all PA cages
	## left-out cage
	assign(paste("freq_diff_only",x,sep=""),
		cbind(haf.sites.filt,
			PA.T1vT2.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="2" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1" & haf.meta.filt$cage.fix==x],
			PA.T1vT3.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="3" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1" & haf.meta.filt$cage.fix==x],
			PA.T1vT4.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="4" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1" & haf.meta.filt$cage.fix==x],
			PA.T2vT3.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="3" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="2" & haf.meta.filt$cage.fix==x],
			PA.T2vT4.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="4" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="2" & haf.meta.filt$cage.fix==x],
			PA.T3vT4.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="4" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="3" & haf.meta.filt$cage.fix==x],
			PAvE.T1.diff=(rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1"]))/abs(rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1"]))))
}


### Reorganize to By temporal contrast ###
## Tedious way to do this, I admit, but easier to keep track as opposed to nested looping 
## with multiple variables in object/field names.

## T1vT2
PA.T1vT2.diff_ind <- cbind(freq_diff_only2[,c(1,2)],PAvE.T1.sign=freq_diff_only2$PAvE.T1.diff,cage2=freq_diff_only2$PA.T1vT2.diff, cage6=freq_diff_only6$PA.T1vT2.diff, cage10=freq_diff_only10$PA.T1vT2.diff, cage14=freq_diff_only14$PA.T1vT2.diff, cage20=freq_diff_only20$PA.T1vT2.diff, cage24=freq_diff_only24$PA.T1vT2.diff, cage31=freq_diff_only31$PA.T1vT2.diff, cage36=freq_diff_only36$PA.T1vT2.diff, cage40=freq_diff_only40$PA.T1vT2.diff, cage44=freq_diff_only44$PA.T1vT2.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
PA.T1vT2.diff_signed <- PA.T1vT2.diff_ind[,c(1:2)]
for(x in 4:ncol(PA.T1vT2.diff_ind)) {
	PA.T1vT2.diff_signed <- cbind(PA.T1vT2.diff_signed, PA.T1vT2.diff_ind[,x] * PA.T1vT2.diff_ind[,3])
	names(PA.T1vT2.diff_signed)[x-1] <- names(PA.T1vT2.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
PA.T1vT2.diff_signed <- cbind(PA.T1vT2.diff_signed, maxes=rowMaxs(as.matrix(PA.T1vT2.diff_signed[,c(3:12)])))

## Add column containing min sign-corredted -log10(p) values per row
PA.T1vT2.diff_signed <- cbind(PA.T1vT2.diff_signed, mins=rowMins(as.matrix(PA.T1vT2.diff_signed[,c(3:12)])))

## Add column containing PAvE.T1 GLM FDR values
PA.T1vT2.diff_signed <- merge(PA.T1vT2.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
PA.T1vT2.diff_signed <- merge(PA.T1vT2.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$PA.T1vT2.fdr), by=c("CHROM","POS"))

## T1vT3
PA.T1vT3.diff_ind <- cbind(freq_diff_only2[,c(1,2)],PAvE.T1.sign=freq_diff_only2$PAvE.T1.diff,cage2=freq_diff_only2$PA.T1vT3.diff, cage6=freq_diff_only6$PA.T1vT3.diff, cage10=freq_diff_only10$PA.T1vT3.diff, cage14=freq_diff_only14$PA.T1vT3.diff, cage20=freq_diff_only20$PA.T1vT3.diff, cage24=freq_diff_only24$PA.T1vT3.diff, cage31=freq_diff_only31$PA.T1vT3.diff, cage36=freq_diff_only36$PA.T1vT3.diff, cage40=freq_diff_only40$PA.T1vT3.diff, cage44=freq_diff_only44$PA.T1vT3.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
PA.T1vT3.diff_signed <- PA.T1vT3.diff_ind[,c(1:2)]
for(x in 4:ncol(PA.T1vT3.diff_ind)) {
	PA.T1vT3.diff_signed <- cbind(PA.T1vT3.diff_signed, PA.T1vT3.diff_ind[,x] * PA.T1vT3.diff_ind[,3])
	names(PA.T1vT3.diff_signed)[x-1] <- names(PA.T1vT3.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
PA.T1vT3.diff_signed <- cbind(PA.T1vT3.diff_signed, maxes=rowMaxs(as.matrix(PA.T1vT3.diff_signed[,c(3:12)])))

## Add column containing min sign-corredted -log10(p) values per row
PA.T1vT3.diff_signed <- cbind(PA.T1vT3.diff_signed, mins=rowMins(as.matrix(PA.T1vT3.diff_signed[,c(3:12)])))

## Add column containing PAvE.T1 GLM FDR values
PA.T1vT3.diff_signed <- merge(PA.T1vT3.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
PA.T1vT3.diff_signed <- merge(PA.T1vT3.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$PA.T1vT3.fdr), by=c("CHROM","POS"))

## T1vT4
PA.T1vT4.diff_ind <- cbind(freq_diff_only2[,c(1,2)],PAvE.T1.sign=freq_diff_only2$PAvE.T1.diff,cage2=freq_diff_only2$PA.T1vT4.diff, cage6=freq_diff_only6$PA.T1vT4.diff, cage10=freq_diff_only10$PA.T1vT4.diff, cage14=freq_diff_only14$PA.T1vT4.diff, cage20=freq_diff_only20$PA.T1vT4.diff, cage24=freq_diff_only24$PA.T1vT4.diff, cage31=freq_diff_only31$PA.T1vT4.diff, cage36=freq_diff_only36$PA.T1vT4.diff, cage40=freq_diff_only40$PA.T1vT4.diff, cage44=freq_diff_only44$PA.T1vT4.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
PA.T1vT4.diff_signed <- PA.T1vT4.diff_ind[,c(1:2)]
for(x in 4:ncol(PA.T1vT4.diff_ind)) {
	PA.T1vT4.diff_signed <- cbind(PA.T1vT4.diff_signed, PA.T1vT4.diff_ind[,x] * PA.T1vT4.diff_ind[,3])
	names(PA.T1vT4.diff_signed)[x-1] <- names(PA.T1vT4.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
PA.T1vT4.diff_signed <- cbind(PA.T1vT4.diff_signed, maxes=rowMaxs(as.matrix(PA.T1vT4.diff_signed[,c(3:12)])))

## Add column containing min sign-corredted -log10(p) values per row
PA.T1vT4.diff_signed <- cbind(PA.T1vT4.diff_signed, mins=rowMins(as.matrix(PA.T1vT4.diff_signed[,c(3:12)])))

## Add column containing PAvE.T1 GLM FDR values
PA.T1vT4.diff_signed <- merge(PA.T1vT4.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
PA.T1vT4.diff_signed <- merge(PA.T1vT4.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$PA.T1vT4.fdr), by=c("CHROM","POS"))

## T2vT3
PA.T2vT3.diff_ind <- cbind(freq_diff_only2[,c(1,2)],PAvE.T1.sign=freq_diff_only2$PAvE.T1.diff,cage2=freq_diff_only2$PA.T2vT3.diff, cage6=freq_diff_only6$PA.T2vT3.diff, cage10=freq_diff_only10$PA.T2vT3.diff, cage14=freq_diff_only14$PA.T2vT3.diff, cage20=freq_diff_only20$PA.T2vT3.diff, cage24=freq_diff_only24$PA.T2vT3.diff, cage31=freq_diff_only31$PA.T2vT3.diff, cage36=freq_diff_only36$PA.T2vT3.diff, cage40=freq_diff_only40$PA.T2vT3.diff, cage44=freq_diff_only44$PA.T2vT3.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
PA.T2vT3.diff_signed <- PA.T2vT3.diff_ind[,c(1:2)]
for(x in 4:ncol(PA.T2vT3.diff_ind)) {
	PA.T2vT3.diff_signed <- cbind(PA.T2vT3.diff_signed, PA.T2vT3.diff_ind[,x] * PA.T2vT3.diff_ind[,3])
	names(PA.T2vT3.diff_signed)[x-1] <- names(PA.T2vT3.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
PA.T2vT3.diff_signed <- cbind(PA.T2vT3.diff_signed, maxes=rowMaxs(as.matrix(PA.T2vT3.diff_signed[,c(3:12)])))

## Add column containing min sign-corredted -log10(p) values per row
PA.T2vT3.diff_signed <- cbind(PA.T2vT3.diff_signed, mins=rowMins(as.matrix(PA.T2vT3.diff_signed[,c(3:12)])))

## Add column containing PAvE.T1 GLM FDR values
PA.T2vT3.diff_signed <- merge(PA.T2vT3.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
PA.T2vT3.diff_signed <- merge(PA.T2vT3.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$PA.T2vT3.fdr), by=c("CHROM","POS"))

## T2vT4
PA.T2vT4.diff_ind <- cbind(freq_diff_only2[,c(1,2)],PAvE.T1.sign=freq_diff_only2$PAvE.T1.diff,cage2=freq_diff_only2$PA.T2vT4.diff, cage6=freq_diff_only6$PA.T2vT4.diff, cage10=freq_diff_only10$PA.T2vT4.diff, cage14=freq_diff_only14$PA.T2vT4.diff, cage20=freq_diff_only20$PA.T2vT4.diff, cage24=freq_diff_only24$PA.T2vT4.diff, cage31=freq_diff_only31$PA.T2vT4.diff, cage36=freq_diff_only36$PA.T2vT4.diff, cage40=freq_diff_only40$PA.T2vT4.diff, cage44=freq_diff_only44$PA.T2vT4.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
PA.T2vT4.diff_signed <- PA.T2vT4.diff_ind[,c(1:2)]
for(x in 4:ncol(PA.T2vT4.diff_ind)) {
	PA.T2vT4.diff_signed <- cbind(PA.T2vT4.diff_signed, PA.T2vT4.diff_ind[,x] * PA.T2vT4.diff_ind[,3])
	names(PA.T2vT4.diff_signed)[x-1] <- names(PA.T2vT4.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
PA.T2vT4.diff_signed <- cbind(PA.T2vT4.diff_signed, maxes=rowMaxs(as.matrix(PA.T2vT4.diff_signed[,c(3:12)])))

## Add column containing min sign-corredted -log10(p) values per row
PA.T2vT4.diff_signed <- cbind(PA.T2vT4.diff_signed, mins=rowMins(as.matrix(PA.T2vT4.diff_signed[,c(3:12)])))

## Add column containing PAvE.T1 GLM FDR values
PA.T2vT4.diff_signed <- merge(PA.T2vT4.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
PA.T2vT4.diff_signed <- merge(PA.T2vT4.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$PA.T2vT4.fdr), by=c("CHROM","POS"))

## T3vT4
PA.T3vT4.diff_ind <- cbind(freq_diff_only2[,c(1,2)],PAvE.T1.sign=freq_diff_only2$PAvE.T1.diff,cage2=freq_diff_only2$PA.T3vT4.diff, cage6=freq_diff_only6$PA.T3vT4.diff, cage10=freq_diff_only10$PA.T3vT4.diff, cage14=freq_diff_only14$PA.T3vT4.diff, cage20=freq_diff_only20$PA.T3vT4.diff, cage24=freq_diff_only24$PA.T3vT4.diff, cage31=freq_diff_only31$PA.T3vT4.diff, cage36=freq_diff_only36$PA.T3vT4.diff, cage40=freq_diff_only40$PA.T3vT4.diff, cage44=freq_diff_only44$PA.T3vT4.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
PA.T3vT4.diff_signed <- PA.T3vT4.diff_ind[,c(1:2)]
for(x in 4:ncol(PA.T3vT4.diff_ind)) {
	PA.T3vT4.diff_signed <- cbind(PA.T3vT4.diff_signed, PA.T3vT4.diff_ind[,x] * PA.T3vT4.diff_ind[,3])
	names(PA.T3vT4.diff_signed)[x-1] <- names(PA.T3vT4.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
PA.T3vT4.diff_signed <- cbind(PA.T3vT4.diff_signed, maxes=rowMaxs(as.matrix(PA.T3vT4.diff_signed[,c(3:12)])))

## Add column containing min sign-corredted -log10(p) values per row
PA.T3vT4.diff_signed <- cbind(PA.T3vT4.diff_signed, mins=rowMins(as.matrix(PA.T3vT4.diff_signed[,c(3:12)])))

## Add column containing PAvE.T1 GLM FDR values
PA.T3vT4.diff_signed <- merge(PA.T3vT4.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
PA.T3vT4.diff_signed <- merge(PA.T3vT4.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$PA.T3vT4.fdr), by=c("CHROM","POS"))

####################################
### Hypergeometric overlap tests ###
####################################

### Parallel temporal SNP overlap with positive (PA-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed") #object names
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$mins > 0.03 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$mins > 0.03,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_temporal_positive_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_temporal_positive_results, file="rudflies_2023_redo.PAvE_PA_temporal_positive_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$maxes < -0.03 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$maxes < -0.03,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_temporal_negative_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_temporal_negative_results, file="rudflies_2023_redo.PAvE_PA_temporal_negative_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Significant parallel temporal SNP overlap with positive (PA-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$mins > 0.03 & PA_temporal_input$temporal.fdr < 0.05 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$mins > 0.03 & PA_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_sig_temporal_positive_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_sig_temporal_positive_results, file="rudflies_2023_redo.PAvE_PA_sig_temporal_positive_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Significant parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$maxes < -0.03 & PA_temporal_input$temporal.fdr < 0.05 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$maxes < -0.03 & PA_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_sig_temporal_negative_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_sig_temporal_negative_results, file="rudflies_2023_redo.PAvE_PA_sig_temporal_negative_results.txt", sep = "\t", quote = FALSE, row.names = F)



##############################################
### Redo with condordance in 9/10 PA cages ###
##############################################

## For this strategy, we will sort the sign-corrected -log10(p) values per row and choose 
## the 2nd lowest and highest values to make filtering less conservative prior to phyper.
## We will add these new columns to the existing data-frames used for phyper seen above.

## T1vT2
temp_sort <- t(apply(PA.T1vT2.diff_signed[,c(3:12)], 1, sort))
PA.T1vT2.diff_signed <- cbind(PA.T1vT2.diff_signed, maxes9=temp_sort[,9])
PA.T1vT2.diff_signed <- cbind(PA.T1vT2.diff_signed, mins9=temp_sort[,2])

## T1vT3
temp_sort <- t(apply(PA.T1vT3.diff_signed[,c(3:12)], 1, sort))
PA.T1vT3.diff_signed <- cbind(PA.T1vT3.diff_signed, maxes9=temp_sort[,9])
PA.T1vT3.diff_signed <- cbind(PA.T1vT3.diff_signed, mins9=temp_sort[,2])

## T1vT4
temp_sort <- t(apply(PA.T1vT4.diff_signed[,c(3:12)], 1, sort))
PA.T1vT4.diff_signed <- cbind(PA.T1vT4.diff_signed, maxes9=temp_sort[,9])
PA.T1vT4.diff_signed <- cbind(PA.T1vT4.diff_signed, mins9=temp_sort[,2])

## T2vT3
temp_sort <- t(apply(PA.T2vT3.diff_signed[,c(3:12)], 1, sort))
PA.T2vT3.diff_signed <- cbind(PA.T2vT3.diff_signed, maxes9=temp_sort[,9])
PA.T2vT3.diff_signed <- cbind(PA.T2vT3.diff_signed, mins9=temp_sort[,2])

## T2vT4
temp_sort <- t(apply(PA.T2vT4.diff_signed[,c(3:12)], 1, sort))
PA.T2vT4.diff_signed <- cbind(PA.T2vT4.diff_signed, maxes9=temp_sort[,9])
PA.T2vT4.diff_signed <- cbind(PA.T2vT4.diff_signed, mins9=temp_sort[,2])

## T3vT4
temp_sort <- t(apply(PA.T3vT4.diff_signed[,c(3:12)], 1, sort))
PA.T3vT4.diff_signed <- cbind(PA.T3vT4.diff_signed, maxes9=temp_sort[,9])
PA.T3vT4.diff_signed <- cbind(PA.T3vT4.diff_signed, mins9=temp_sort[,2])


### Parallel temporal SNP overlap with positive (PA-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$mins9 > 0.03 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$mins9 > 0.03,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_temporal_9samp_positive_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_temporal_9samp_positive_results, file="rudflies_2023_redo.PAvE_PA_9samp_temporal_positive_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$maxes9 < -0.03 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$maxes9 < -0.03,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_temporal_9samp_negative_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_temporal_9samp_negative_results, file="rudflies_2023_redo.PAvE_PA_9samp_temporal_negative_results.txt", sep = "\t", quote = FALSE, row.names = F)



### Significant parallel temporal SNP overlap with positive (PA-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$mins9 > 0.03 & PA_temporal_input$temporal.fdr < 0.05 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$mins9 > 0.03 & PA_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_sig_temporal_9samp_positive_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_sig_temporal_9samp_positive_results, file="rudflies_2023_redo.PAvE_PA_9samp_sig_temporal_positive_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Significant parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$maxes9 < -0.03 & PA_temporal_input$temporal.fdr < 0.05 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$maxes9 < -0.03 & PA_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_sig_temporal_9samp_negative_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_sig_temporal_9samp_negative_results, file="rudflies_2023_redo.PAvE_PA_9samp_sig_temporal_negative_results.txt", sep = "\t", quote = FALSE, row.names = F)


### GLM FDR outlier hypergeomtric test only ###

### Significant parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("PA.T1vT2.diff_signed", "PA.T1vT3.diff_signed", "PA.T1vT4.diff_signed", "PA.T2vT3.diff_signed", "PA.T2vT4.diff_signed", "PA.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	PA_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$temporal.fdr < 0.05 & PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(PA_temporal_input) - m)
  	k <- as.numeric(nrow(PA_temporal_input[PA_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
PA_sig_temporal_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(PA_sig_temporal_results, file="rudflies_2023_redo.PA_sig_temporal_results.txt", sep = "\t", quote = FALSE, row.names = F)


############ Look for concordance of temporal selection in all E cages ############

## Let's calculate mean frequencies, then calculate pairwise differences
## We'll use these values to assign positive or negative signs to -log10p values
## The resulting values will reflect both the magnitude and direction of AF change

E_cages <- unique(haf.meta.filt[haf.meta.filt$treat.fix=="E",]$cage.fix)
#[1] "1"  "5"  "9"  "13" "19" "23" "29" "35" "39" "43"

for(x in E_cages) { #cycle through all E cages
	## calculate 
	assign(paste("freq_diff_only",x,sep=""),
		cbind(haf.sites.filt,
			E.T1vT2.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="2" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1" & haf.meta.filt$cage.fix==x],
			E.T1vT3.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="3" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1" & haf.meta.filt$cage.fix==x],
			E.T1vT4.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="4" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1" & haf.meta.filt$cage.fix==x],
			E.T2vT3.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="3" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="2" & haf.meta.filt$cage.fix==x],
			E.T2vT4.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="4" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="2" & haf.meta.filt$cage.fix==x],
			E.T3vT4.diff=haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="4" & haf.meta.filt$cage.fix==x] - haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="3" & haf.meta.filt$cage.fix==x],
			PAvE.T1.diff=(rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1"]))/abs(rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="PA" & haf.meta.filt$tpt=="1"]) - rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" & haf.meta.filt$tpt=="1"]))))
}



### By temporal contrast ###

## T1vT2
E.T1vT2.diff_ind <- cbind(freq_diff_only1[,c(1,2)],PAvE.T1.sign=freq_diff_only1$PAvE.T1.diff,cage1=freq_diff_only1$E.T1vT2.diff, cage5=freq_diff_only5$E.T1vT2.diff, cage9=freq_diff_only9$E.T1vT2.diff, cage13=freq_diff_only13$E.T1vT2.diff, cage23=freq_diff_only23$E.T1vT2.diff, cage29=freq_diff_only29$E.T1vT2.diff, cage35=freq_diff_only35$E.T1vT2.diff, cage39=freq_diff_only39$E.T1vT2.diff, cage43=freq_diff_only43$E.T1vT2.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
E.T1vT2.diff_signed <- E.T1vT2.diff_ind[,c(1:2)]
for(x in 4:ncol(E.T1vT2.diff_ind)) {
	E.T1vT2.diff_signed <- cbind(E.T1vT2.diff_signed, E.T1vT2.diff_ind[,x] * E.T1vT2.diff_ind[,3])
	names(E.T1vT2.diff_signed)[x-1] <- names(E.T1vT2.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
E.T1vT2.diff_signed <- cbind(E.T1vT2.diff_signed, maxes=rowMaxs(as.matrix(E.T1vT2.diff_signed[,c(3:11)])))

## Add column containing min sign-corredted -log10(p) values per row
E.T1vT2.diff_signed <- cbind(E.T1vT2.diff_signed, mins=rowMins(as.matrix(E.T1vT2.diff_signed[,c(3:11)])))

## Add column containing PAvE.T1 GLM FDR values
E.T1vT2.diff_signed <- merge(E.T1vT2.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
E.T1vT2.diff_signed <- merge(E.T1vT2.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$E.T1vT2.fdr), by=c("CHROM","POS"))

## T1vT3
E.T1vT3.diff_ind <- cbind(freq_diff_only1[,c(1,2)],PAvE.T1.sign=freq_diff_only1$PAvE.T1.diff,cage1=freq_diff_only1$E.T1vT3.diff, cage5=freq_diff_only5$E.T1vT3.diff, cage9=freq_diff_only9$E.T1vT3.diff, cage13=freq_diff_only13$E.T1vT3.diff, cage23=freq_diff_only23$E.T1vT3.diff, cage29=freq_diff_only29$E.T1vT3.diff, cage35=freq_diff_only35$E.T1vT3.diff, cage39=freq_diff_only39$E.T1vT3.diff, cage43=freq_diff_only43$E.T1vT3.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
E.T1vT3.diff_signed <- E.T1vT3.diff_ind[,c(1:2)]
for(x in 4:ncol(E.T1vT3.diff_ind)) {
	E.T1vT3.diff_signed <- cbind(E.T1vT3.diff_signed, E.T1vT3.diff_ind[,x] * E.T1vT3.diff_ind[,3])
	names(E.T1vT3.diff_signed)[x-1] <- names(E.T1vT3.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
E.T1vT3.diff_signed <- cbind(E.T1vT3.diff_signed, maxes=rowMaxs(as.matrix(E.T1vT3.diff_signed[,c(3:11)])))

## Add column containing min sign-corredted -log10(p) values per row
E.T1vT3.diff_signed <- cbind(E.T1vT3.diff_signed, mins=rowMins(as.matrix(E.T1vT3.diff_signed[,c(3:11)])))

## Add column containing PAvE.T1 GLM FDR values
E.T1vT3.diff_signed <- merge(E.T1vT3.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
E.T1vT3.diff_signed <- merge(E.T1vT3.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$E.T1vT3.fdr), by=c("CHROM","POS"))

## T1vT4
E.T1vT4.diff_ind <- cbind(freq_diff_only1[,c(1,2)],PAvE.T1.sign=freq_diff_only1$PAvE.T1.diff,cage1=freq_diff_only1$E.T1vT4.diff, cage5=freq_diff_only5$E.T1vT4.diff, cage9=freq_diff_only9$E.T1vT4.diff, cage13=freq_diff_only13$E.T1vT4.diff, cage23=freq_diff_only23$E.T1vT4.diff, cage29=freq_diff_only29$E.T1vT4.diff, cage35=freq_diff_only35$E.T1vT4.diff, cage39=freq_diff_only39$E.T1vT4.diff, cage43=freq_diff_only43$E.T1vT4.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
E.T1vT4.diff_signed <- E.T1vT4.diff_ind[,c(1:2)]
for(x in 4:ncol(E.T1vT4.diff_ind)) {
	E.T1vT4.diff_signed <- cbind(E.T1vT4.diff_signed, E.T1vT4.diff_ind[,x] * E.T1vT4.diff_ind[,3])
	names(E.T1vT4.diff_signed)[x-1] <- names(E.T1vT4.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
E.T1vT4.diff_signed <- cbind(E.T1vT4.diff_signed, maxes=rowMaxs(as.matrix(E.T1vT4.diff_signed[,c(3:11)])))

## Add column containing min sign-corredted -log10(p) values per row
E.T1vT4.diff_signed <- cbind(E.T1vT4.diff_signed, mins=rowMins(as.matrix(E.T1vT4.diff_signed[,c(3:11)])))

## Add column containing PAvE.T1 GLM FDR values
E.T1vT4.diff_signed <- merge(E.T1vT4.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
E.T1vT4.diff_signed <- merge(E.T1vT4.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$E.T1vT4.fdr), by=c("CHROM","POS"))

## T2vT3
E.T2vT3.diff_ind <- cbind(freq_diff_only1[,c(1,2)],PAvE.T1.sign=freq_diff_only1$PAvE.T1.diff,cage1=freq_diff_only1$E.T2vT3.diff, cage5=freq_diff_only5$E.T2vT3.diff, cage9=freq_diff_only9$E.T2vT3.diff, cage13=freq_diff_only13$E.T2vT3.diff, cage23=freq_diff_only23$E.T2vT3.diff, cage29=freq_diff_only29$E.T2vT3.diff, cage35=freq_diff_only35$E.T2vT3.diff, cage39=freq_diff_only39$E.T2vT3.diff, cage43=freq_diff_only43$E.T2vT3.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
E.T2vT3.diff_signed <- E.T2vT3.diff_ind[,c(1:2)]
for(x in 4:ncol(E.T2vT3.diff_ind)) {
	E.T2vT3.diff_signed <- cbind(E.T2vT3.diff_signed, E.T2vT3.diff_ind[,x] * E.T2vT3.diff_ind[,3])
	names(E.T2vT3.diff_signed)[x-1] <- names(E.T2vT3.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
E.T2vT3.diff_signed <- cbind(E.T2vT3.diff_signed, maxes=rowMaxs(as.matrix(E.T2vT3.diff_signed[,c(3:11)])))

## Add column containing min sign-corredted -log10(p) values per row
E.T2vT3.diff_signed <- cbind(E.T2vT3.diff_signed, mins=rowMins(as.matrix(E.T2vT3.diff_signed[,c(3:11)])))

## Add column containing PAvE.T1 GLM FDR values
E.T2vT3.diff_signed <- merge(E.T2vT3.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
E.T2vT3.diff_signed <- merge(E.T2vT3.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$E.T2vT3.fdr), by=c("CHROM","POS"))

## T2vT4
E.T2vT4.diff_ind <- cbind(freq_diff_only1[,c(1,2)],PAvE.T1.sign=freq_diff_only1$PAvE.T1.diff,cage1=freq_diff_only1$E.T2vT4.diff, cage5=freq_diff_only5$E.T2vT4.diff, cage9=freq_diff_only9$E.T2vT4.diff, cage13=freq_diff_only13$E.T2vT4.diff, cage23=freq_diff_only23$E.T2vT4.diff, cage29=freq_diff_only29$E.T2vT4.diff, cage35=freq_diff_only35$E.T2vT4.diff, cage39=freq_diff_only39$E.T2vT4.diff, cage43=freq_diff_only43$E.T2vT4.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
E.T2vT4.diff_signed <- E.T2vT4.diff_ind[,c(1:2)]
for(x in 4:ncol(E.T2vT4.diff_ind)) {
	E.T2vT4.diff_signed <- cbind(E.T2vT4.diff_signed, E.T2vT4.diff_ind[,x] * E.T2vT4.diff_ind[,3])
	names(E.T2vT4.diff_signed)[x-1] <- names(E.T2vT4.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
E.T2vT4.diff_signed <- cbind(E.T2vT4.diff_signed, maxes=rowMaxs(as.matrix(E.T2vT4.diff_signed[,c(3:11)])))

## Add column containing min sign-corredted -log10(p) values per row
E.T2vT4.diff_signed <- cbind(E.T2vT4.diff_signed, mins=rowMins(as.matrix(E.T2vT4.diff_signed[,c(3:11)])))

## Add column containing PAvE.T1 GLM FDR values
E.T2vT4.diff_signed <- merge(E.T2vT4.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
E.T2vT4.diff_signed <- merge(E.T2vT4.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$E.T2vT4.fdr), by=c("CHROM","POS"))

## T3vT4
E.T3vT4.diff_ind <- cbind(freq_diff_only1[,c(1,2)],PAvE.T1.sign=freq_diff_only1$PAvE.T1.diff,cage1=freq_diff_only1$E.T3vT4.diff, cage5=freq_diff_only5$E.T3vT4.diff, cage9=freq_diff_only9$E.T3vT4.diff, cage13=freq_diff_only13$E.T3vT4.diff, cage23=freq_diff_only23$E.T3vT4.diff, cage29=freq_diff_only29$E.T3vT4.diff, cage35=freq_diff_only35$E.T3vT4.diff, cage39=freq_diff_only39$E.T3vT4.diff, cage43=freq_diff_only43$E.T3vT4.diff)

#This loop corrects the signs of AF differences in the direction of PA vs E
E.T3vT4.diff_signed <- E.T3vT4.diff_ind[,c(1:2)]
for(x in 4:ncol(E.T3vT4.diff_ind)) {
	E.T3vT4.diff_signed <- cbind(E.T3vT4.diff_signed, E.T3vT4.diff_ind[,x] * E.T3vT4.diff_ind[,3])
	names(E.T3vT4.diff_signed)[x-1] <- names(E.T3vT4.diff_ind)[x]
}

## Add column containing max sign-corredted -log10(p) values per row
E.T3vT4.diff_signed <- cbind(E.T3vT4.diff_signed, maxes=rowMaxs(as.matrix(E.T3vT4.diff_signed[,c(3:11)])))

## Add column containing min sign-corredted -log10(p) values per row
E.T3vT4.diff_signed <- cbind(E.T3vT4.diff_signed, mins=rowMins(as.matrix(E.T3vT4.diff_signed[,c(3:11)])))

## Add column containing PAvE.T1 GLM FDR values
E.T3vT4.diff_signed <- merge(E.T3vT4.diff_signed, cbind(glm.all[,c(1,2)], PAvE.T1.fdr=glm.all$PAvE.T1.fdr), by=c("CHROM","POS"))

## Add column containing temporal contrast GLM FDR values
E.T3vT4.diff_signed <- merge(E.T3vT4.diff_signed, cbind(glm.all[,c(1,2)], temporal.fdr=glm.all$E.T3vT4.fdr), by=c("CHROM","POS"))

### Parallel temporal SNP overlap with positive (PA-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$mins > 0.03 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$mins > 0.03,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_temporal_positive_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_temporal_positive_results, file="rudflies_2023_redo.PAvE_E_temporal_positive_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$maxes < -0.03 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$maxes < -0.03,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_temporal_negative_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_temporal_negative_results, file="rudflies_2023_redo.PAvE_E_temporal_negative_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Significant parallel temporal SNP overlap with positive (PA-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$mins > 0.03 & E_temporal_input$temporal.fdr < 0.05 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$mins > 0.03 & E_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_sig_temporal_positive_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_sig_temporal_positive_results, file="rudflies_2023_redo.PAvE_E_sig_temporal_positive_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Significant parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$maxes < -0.03 & E_temporal_input$temporal.fdr < 0.05 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$maxes < -0.03 & E_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_sig_temporal_negative_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_sig_temporal_negative_results, file="rudflies_2023_redo.PAvE_E_sig_temporal_negative_results.txt", sep = "\t", quote = FALSE, row.names = F)



############################################
### Redo with condordance in 8/9 E cages ###
############################################

## For this strategy, we will sort the sign-corrected -log10(p) values per row and choose 
## the 2nd lowest and highest values to make filtering less conservative prior to phyper.
## We will add these new columns to the existing data-frames used for phyper seen above.

## T1vT2
temp_sort <- t(apply(E.T1vT2.diff_signed[,c(3:11)], 1, sort))
E.T1vT2.diff_signed <- cbind(E.T1vT2.diff_signed, maxes8=temp_sort[,8])
E.T1vT2.diff_signed <- cbind(E.T1vT2.diff_signed, mins8=temp_sort[,2])

## T1vT3
temp_sort <- t(apply(E.T1vT3.diff_signed[,c(3:11)], 1, sort))
E.T1vT3.diff_signed <- cbind(E.T1vT3.diff_signed, maxes8=temp_sort[,8])
E.T1vT3.diff_signed <- cbind(E.T1vT3.diff_signed, mins8=temp_sort[,2])

## T1vT4
temp_sort <- t(apply(E.T1vT4.diff_signed[,c(3:11)], 1, sort))
E.T1vT4.diff_signed <- cbind(E.T1vT4.diff_signed, maxes8=temp_sort[,8])
E.T1vT4.diff_signed <- cbind(E.T1vT4.diff_signed, mins8=temp_sort[,2])

## T2vT3
temp_sort <- t(apply(E.T2vT3.diff_signed[,c(3:11)], 1, sort))
E.T2vT3.diff_signed <- cbind(E.T2vT3.diff_signed, maxes8=temp_sort[,8])
E.T2vT3.diff_signed <- cbind(E.T2vT3.diff_signed, mins8=temp_sort[,2])

## T2vT4
temp_sort <- t(apply(E.T2vT4.diff_signed[,c(3:11)], 1, sort))
E.T2vT4.diff_signed <- cbind(E.T2vT4.diff_signed, maxes8=temp_sort[,8])
E.T2vT4.diff_signed <- cbind(E.T2vT4.diff_signed, mins8=temp_sort[,2])

## T3vT4
temp_sort <- t(apply(E.T3vT4.diff_signed[,c(3:11)], 1, sort))
E.T3vT4.diff_signed <- cbind(E.T3vT4.diff_signed, maxes8=temp_sort[,8])
E.T3vT4.diff_signed <- cbind(E.T3vT4.diff_signed, mins8=temp_sort[,2])


### Parallel temporal SNP overlap with positive (PA-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$mins8 > 0.03 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$mins8 > 0.03,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_temporal_8samp_positive_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_temporal_8samp_positive_results, file="rudflies_2023_redo.PAvE_E_8samp_temporal_positive_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$maxes8 < -0.03 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$maxes8 < -0.03,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_temporal_8samp_negative_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_temporal_8samp_negative_results, file="rudflies_2023_redo.PAvE_E_8samp_temporal_negative_results.txt", sep = "\t", quote = FALSE, row.names = F)



### Significant parallel temporal SNP overlap with positive (PA-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$mins8 > 0.03 & E_temporal_input$temporal.fdr < 0.05 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$mins8 > 0.03 & E_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_sig_temporal_8samp_positive_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_sig_temporal_8samp_positive_results, file="rudflies_2023_redo.PAvE_E_8samp_sig_temporal_positive_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Significant parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$maxes8 < -0.03 & E_temporal_input$temporal.fdr < 0.05 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$maxes8 < -0.03 & E_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_sig_temporal_8samp_negative_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_sig_temporal_8samp_negative_results, file="rudflies_2023_redo.PAvE_E_8samp_sig_temporal_negative_results.txt", sep = "\t", quote = FALSE, row.names = F)


### GLM FDR outlier hypergeomtric test only ###

### Significant parallel temporal SNP overlap with negative (E-biased) PAvE.T1 SNPs ###
temporal <- c("E.T1vT2.diff_signed", "E.T1vT3.diff_signed", "E.T1vT4.diff_signed", "E.T2vT3.diff_signed", "E.T2vT4.diff_signed", "E.T3vT4.diff_signed")
temp_res <- c()
#create empty object
q <- c() #candidate - outlier overlap SNP count
m <- c() #candidate SNP count
n <- c() #non-candidate SNP count
k <- c() #outlier SNP count
## Loop through all temporal contrasts
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
  	q <- as.numeric(nrow(E_temporal_input[E_temporal_input$temporal.fdr < 0.05 & E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	m <- as.numeric(nrow(E_temporal_input[E_temporal_input$temporal.fdr < 0.05,]))
  	n <- as.numeric(nrow(E_temporal_input) - m)
  	k <- as.numeric(nrow(E_temporal_input[E_temporal_input$PAvE.T1.fdr < 0.05,]))
  	## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  	temp_res <- rbind(temp_res, c(q,m,n,k,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)))
}

## Reformat to dataframe
E_sig_temporal_results <- cbind(temporal,as.data.frame(temp_res))

## Save table
write.table(E_sig_temporal_results, file="rudflies_2023_redo.E_sig_temporal_results.txt", sep = "\t", quote = FALSE, row.names = F)


### Check correlation between the max AF changes and the second highest maxes for QC ###
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
	assign(paste0(temporal[y],".max.corr"), ggplot(E_temporal_input, 
	  	aes(x=maxes, y=maxes8)) + 
	  	geom_bin2d(bins = 50) +
	  	scale_fill_continuous(type = "viridis") +
	  	xlab("maxes") + 
	  	ylab("maxes8") +
	  	theme_classic())
}

### Check correlation between the min AF changes and the second lowest minsfor QC ###
for(y in 1:length(temporal)) {
	E_temporal_input <- get(temporal[y])
	assign(paste0(temporal[y],".min.corr"), ggplot(E_temporal_input, 
	  	aes(x=mins, y=mins8)) + 
	  	geom_bin2d(bins = 50) +
	  	scale_fill_continuous(type = "viridis") +
	  	xlab("mins") + 
	  	ylab("mins8") +
	  	theme_classic())
}

## Save all correlation plots together
pdf(file = "rudflies_2023_redo.E_temporal_min_max_corr.pdf", 
	width=8, height=10)
	ggarrange(E.T1vT2.diff_signed.min.corr, E.T1vT3.diff_signed.min.corr, E.T1vT4.diff_signed.min.corr, E.T2vT3.diff_signed.min.corr, E.T2vT4.diff_signed.min.corr, E.T3vT4.diff_signed.min.corr, E.T1vT2.diff_signed.max.corr, E.T1vT3.diff_signed.max.corr, E.T1vT4.diff_signed.max.corr, E.T2vT3.diff_signed.max.corr, E.T2vT4.diff_signed.max.corr, E.T3vT4.diff_signed.max.corr, 
              ncol = 2, nrow = 6)
dev.off()

