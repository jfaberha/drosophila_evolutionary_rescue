##########################################################################################
### This script is associated with a study on the evolution of spinosad-resistance in  ###
### field populations of Drosophila melanogaster run in 2023 by the Rudman lab (WSUV). ###
### This particular script covers a portion of the bioinformatic analysis that examines###
### levels of parallelism in the adaptive outlier loci. More specifically, we expected ###
### GLM outliers to be parallel, but significant GLM results can result from either    ###
### consistent AF changes across all samples in a group or extreme changes in multiple,###
### but not all, samples. This analysis seeks to quantify parallelism, as that is the  ###
### key metric to discern selection from drift.                                        ###
###                                                                                    ###
### There are several supplemental analyses included here as well. Notably we sought to###
### quantify whether putatively selected-on loci in key time-point/treatment group     ###
### combinations showed parallel selection in other key groups. We saw very little     ###
### overlap when comparing EvSE.T1 and EvSP.T1 contrast outliers, for example, but     ###
### perhaps we can identify reciprocal trends of parallelism in loci that show lower,  ###
### but consistent, levels of AF differentiation. In addition to EvSE.T1 and           ###
### EvSP.T1 SNP lists, we compared EvSP.T4 outlier SNPs as well, since that is the     ###
### experimental endpoint for adapted SP flies.                                        ###
##########################################################################################

### In R ###
#configure r environment
setwd("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/parallel_analysis")
library(emmeans)
library(matrixStats)
library(ACER)
library(poolSeq)
library(ggplot2)
library(gplots)
library(ggpubr)
library(stats)
library(ggfortify)
library(dplyr)
library(zoo)
library(rstatix)
library(slider)
library(RColorBrewer)
library(geiger)
library(RRHO)
library(RRHO2)
library(stringr)
library(scales)
library(slider)
library(GenomicRanges)
library(gridExtra)
library(moments)

## set treatment group color scheme upfront for plotting
#treatment group order: E, PA, SE, SP
sample_cols <- c("#D26183","#495184","#848556","#D9B851")

#################################
### Load required input files ###
#################################

## Variant Effect Predictor (VEP) annotation file with appended FLYCADD scores
vep <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/filtered-all.annot.vcf.FLYCADD.tsv", header=TRUE) 
## Sample metadata table
haf.meta <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_meta.tsv", header=TRUE)
## Hafpipe imputed allele frequency table
haf.freq <- read.delim("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_hafpipe.csv", header=TRUE, sep = ",")

## Do a bit of reformatting for the frequency table header
names(haf.freq) <- gsub("[.]af","",names(haf.freq))
names(haf.freq) <- gsub("X","",names(haf.freq))
names(haf.freq)[1]="CHROM"
names(haf.freq)[2]="POS"

## Filter for PA, S, and E samples
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

## Calculate founder means for all loci (E was the founder for E and S treatments)
haf.freq.Fmean <- cbind(haf.freq[,c(1:2)], E=rowMeans(haf.freq[,c(3:5)]), PA=rowMeans(haf.freq[,c(6:8)]))

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

## For additional filtering, calculate mean allele frequencies for each treatment at TPT1
freq_means_t1 <- cbind(haf.sites.T1filt,
	E.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "E" & haf.meta.T1filt$tpt=="1"]),
	S.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "S" & haf.meta.T1filt$tpt=="1"]),
	SE.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition == "SE" & haf.meta.T1filt$tpt=="1"]),
	SP.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition == "SP" & haf.meta.T1filt$tpt=="1"]),
	PA.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "PA"&haf.meta.T1filt$tpt=="1"]))

## Let's calculate mean frequencies, then calculate pairwise differences
## We'll use these values to assign positive or negative signs to -log10p values
## The resulting values will reflect both the magnitude and direction of AF change
freq_diff_t1 <- cbind(EvSE.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]),
	EvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]),
	EvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SEvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SPvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SEvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]))

## For use in bedtools, we'll create a bed file from AF differences
freq_diff_t1_bed <- cbind(haf.sites.T1filt,STOP=haf.sites.T1filt$POS+1, freq_diff_t1)
freq_diff_t1_bed <- freq_diff_t1_bed[order(freq_diff_t1_bed[,1], freq_diff_t1_bed[,2]), ]
## save
options(scipen=20)
write.table(freq_diff_t1_bed, file="rudflies_2023_redo.freq_diff.T1.bed",sep = "\t", quote = FALSE, row.names = F)
options(scipen=0)


######################################
### For downstream analysis create ###
### a TPT4-only metadata table     ###
######################################

## Make TPT4 AF table subset
## Filter metadata for TPT4 samples for PA, S, and E
haf.meta.T4filt <- haf.meta[haf.meta$tpt == "4" & haf.meta$batch == "a" & (haf.meta$treat.fix == "PA" | haf.meta$treat.fix == "S" | haf.meta$treat.fix == "E"),]

## Filter AF table for TPT4 samples for PA, S, and E
haf.freq.T4filt <- haf.freq[, which((names(haf.freq) %in% haf.meta.T4filt$samp)==TRUE)]
haf.meta.T4filt <- haf.meta.T4filt[which((haf.meta.T4filt$samp %in% names(haf.freq.T4filt))==TRUE),]

#additional filter to remove extreme low variance loci
haf.sites.T4filt <- haf.freq[rowVars(as.matrix(haf.freq.T4filt))>0.001,c(1:2)]
haf.freq.T4filt <- haf.freq.T4filt[rowVars(as.matrix(haf.freq.T4filt))>0.001,]

#tables need a bit of reformatting for plotting
haf.meta.T4filt$treat.fix <- as.factor(haf.meta.T4filt$treat.fix)
haf.meta.T4filt$tpt <- as.factor(haf.meta.T4filt$tpt)


## Let's calculate mean frequencies, then calculate pairwise differences
## We'll use these values to assign positive or negative signs to -log10p values
## The resulting values will reflect both the magnitude and direction of AF change
freq_diff_t4 <- cbind(EvSP.T4.diff=rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="E"]) - rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="S"]),
	EvPA.T4.diff=rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="E"]) - rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="PA"]),
	PAvSP.T4.diff=rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="PA"]) - rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="S"]))
	
## For use in bedtools, we'll create a bed file from AF differences
freq_diff_t4_bed <- cbind(haf.sites.T4filt,STOP=haf.sites.T4filt$POS+1, freq_diff_t4)
freq_diff_t4_bed <- freq_diff_t4_bed[order(freq_diff_t4_bed[,1], freq_diff_t4_bed[,2]), ]

## save
options(scipen=20)
write.table(freq_diff_t4_bed, file="rudflies_2023_redo.freq_diff.T4.bed",sep = "\t", quote = FALSE, row.names = F)
options(scipen=0)


##########################################
### Create master table of GLM results ###
##########################################

### Start with PA and S contrasts ###
## Load PA vs S, all samples combined
contrast.PAvS.treat.all <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvS_treat.all.table.GLMcontrast.txt", header=TRUE)
contrast.PAvS.treat.all <- contrast.PAvS.treat.all[,-3] #remove AF mean column
names(contrast.PAvS.treat.all)[3] <- "PAvS" #label glm results column
## Load PA vs S for each timepoint
contrast.PAvS.treat <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvS_treat.GLMcontrast.txt", header=TRUE)
#names(contrast.PAvS.treat) <- c("CHROM","POS","PAvS.T1","PAvS.T2","PAvS.T3","PAvS.T4")
## Load timepoint contrasts using PA and S combined
contrast.PAvS.tpt <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvS_tpt.GLMcontrast.txt", header=TRUE)
#names(contrast.PAvS.tpt) <- c("CHROM","POS","PA.S.T1vT2","PA.S.T1.T3","PA.S.T1vT4","PA.S.T2vT3","PA.S.T2vT4","PA.S.T3vT4") #label glm results columns
## Load PA and S "timepoint:treatment" interaction
contrast.PAvS.int <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvS_int.GLMcontrast.txt", header=TRUE)
#names(contrast.PAvS.int) <- c("CHROM","POS","int.PAvS.T1vT2","int.PAvS.T1vT3","int.PAvS.T1vT4","int.PAvS.T2vT3","int.PAvS.T2vT4","int.PAvS.T3vT4") #label glm results columns
contrast.PAvS <- cbind(contrast.PAvS.treat, contrast.PAvS.tpt[,-c(1,2)], contrast.PAvS.int[,-c(1,2)])
contrast.PAvS <- merge(contrast.PAvS.treat.all, contrast.PAvS, by=c("CHROM","POS"))
contrast.PAvS <- contrast.PAvS[order(contrast.PAvS[,1], contrast.PAvS[,2]), ] #sort by locus
## Save all PA and S contrasts for posterity
#write.table(contrast.PAvS, file="rudflies_2023_redo.PAvS.multiGLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Look at PA-only, E-only, and S-only time-point contrast GLM results ##
# Load S-only results
contrast.tpt.S.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_S_tpt.table.GLMcontrast.txt", header=TRUE)
#names(contrast.tpt.S.table) <- c("CHROM","POS","S.af.mean","S.T1vT2","S.T1vT3","S.T1vT4","S.T2vT3","S.T2vT4","S.T3vT4") #label glm results columns
# Load PA-only results
contrast.tpt.PA.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PA_tpt.table.GLMcontrast.txt", header=TRUE)
#names(contrast.tpt.PA.table) <- c("CHROM","POS","PA.af.mean","PA.T1vT2","PA.T1vT3","PA.T1vT4","PA.T2vT3","PA.T2vT4","PA.T3vT4") #label glm results columns
# Load E-only results
contrast.tpt.E.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_E_tpt.table.GLMcontrast.txt", header=TRUE)
names(contrast.tpt.E.table) <- c("CHROM","POS","E.af.mean","E.T1vT2","E.T1vT3","E.T1vT4","E.T2vT3","E.T2vT4","E.T3vT4") #label glm results columns
# Merge them
contrast.tpt.S_only.PA_only.E_only <- merge(merge(contrast.tpt.S.table, contrast.tpt.PA.table, by=c("CHROM","POS")),contrast.tpt.E.table, by=c("CHROM","POS"))
contrast.tpt.S_only.PA_only.E_only <- contrast.tpt.S_only.PA_only.E_only[,-c(3,10,17)] #remove AF mean cols
# Save relabeled tables for posterity
#write.table(contrast.tpt.S.table, file="rudflies_2023_S_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
#write.table(contrast.tpt.PA.table, file="rudflies_2023_PA_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)
#write.table(contrast.tpt.E.table, file="rudflies_2023_E_tpt.table.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## TPT1-only GLM contrasts: PA vs S vs SE vs E - all pairwise ##
contrast.PAvSvSEvE.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", header=TRUE)
#names(contrast.PAvSvSEvE.table) <- c("CHROM","POS","EvPA.T1","EvSE.T1","EvSP.T1","PAvSE.T1","PAvSP.T1","SEvSP.T1") #label glm results columns
# Save relabeled table for posterity
# write.table(contrast.PAvSvSEvE.table, file="rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Load PA and E contrasts ##
# PA vs E, all time-points combined
contrastout.PAvE.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvE_treat.all.table.GLMcontrast.txt", header=TRUE)
# PA vs E, founders
contrastout.PAvE.founder.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvE_treat.table.GLMcontrast.founders.txt", header=TRUE)
# PA vs E, at individual time-points
contrastout.treat.PA.E.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvE_treat.table.GLMcontrast2.txt", header=TRUE)
# Combine these three tables
contrastout.PAvE.table <- merge(contrastout.PAvE.table[,-3], contrastout.PAvE.founder.table[,-3], by=c("CHROM","POS"))
contrastout.PAvE.table <- merge(contrastout.PAvE.table, contrastout.treat.PA.E.table[,-3], by=c("CHROM","POS"))
names(contrastout.PAvE.table) <- c("CHROM","POS","PAvE","PAvE.F","PAvE.T1","PAvE.T2","PAvE.T3","PAvE.T4") #label glm results columns
# Save relabeled table for posterity
#write.table(contrastout.PAvE.table, file="rudflies_2023_PAvE_treat.all.table.wLoci.GLMcontrast.txt", sep="\t", quote = FALSE, row.names = F)

## Load SE and E contrasts ##
# S vs E, all time-points combined
contrastout.SvE.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_SvE_treat.all.table.GLMcontrast.txt", header=TRUE)[,c(1,2,4)]
# S vs E, at individual time-points
contrastout.treat.S.E.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_SvE_treat.GLMcontrast.table.txt", header=TRUE)
# Combine these two tables
contrastout.SvE.table <- merge(contrastout.SvE.table, contrastout.treat.S.E.table, by=c("CHROM","POS"))
# S and E: "treatment:time-point" interaction
contrastout.int.S.E.table <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_SvE_int.GLMcontrast.table.txt", header=TRUE)
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
  glm.all[,i] <- as.numeric(glm.all[,i])
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

## Merge annotations for filtering and highlighting
glm.all.annot <- merge(glm.all, vep, by=c("CHROM","POS"))

## Load spinosad-resistance candidate gene list
spino.cand <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/spino.cand.list.txt", header=FALSE)
names(spino.cand) <- "Gene"

## Now merge to find All SNPs in and around candidate genes
glm.all.annot.spino <- merge(spino.cand, glm.all.annot, by="Gene")


##################################################################################
### Create a non-redundant list of SNP consequences ranked by predicted impact ###
##################################################################################

### This will be used for selecting matched loci for downstream statistical tests

## Create a list of loci and each SNP consequence
vep_cons <- unique(cbind(glm.all.annot[,c(1,2)],Consequence=glm.all.annot$Consequence))

## Create a list of SNP consequences and each associated FlyCADD score and VEP impact
vep_score <- as.data.frame(cbind(Consequence=vep$Consequence,FLYCADD=vep$FLYCADD,Extra=vep$Extra))

## Modify column formats for ranking
vep_score$Extra <- sub(";.*", "", vep_score$Extra)
vep_score$Extra <- sub(".*=", "", vep_score$Extra)
vep_score$FLYCADD <- as.numeric(vep_score$FLYCADD)

## Calculate mean FlyCADD score for each unique consequence/impact
vep_score_mean <- vep_score %>%
  group_by(Consequence,Extra) %>% 
  summarise_all(mean)

## Reformat and assign values to SNP impacts for sorting
vep_score_mean <- as.data.frame(vep_score_mean)
vep_score_mean$Impact <- vep_score_mean$Extra
vep_score_mean$Impact <- sub("MODIFIER", "1", vep_score_mean$Impact)
vep_score_mean$Impact <- sub("LOW", "2", vep_score_mean$Impact)
vep_score_mean$Impact <- sub("MODERATE", "3", vep_score_mean$Impact)
vep_score_mean$Impact <- sub("HIGH", "4", vep_score_mean$Impact)
vep_score_mean$Impact <- as.numeric(vep_score_mean$Impact)

## Sort SNP consequences by VEP impact then FlyCADD score, high to low
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

## Make a rolling window table with just p-value, remove logp and fdr columns 
glm.all.rolwin21.pval <- select(glm.all.rolwin21,-contains("logp"),-contains("fdr"))
## Make a rolling window table with just logp, selects them and append loci
glm.all.logp.rolwin21 <- cbind(glm.all.rolwin21.pval[,c(1:2)],select(glm.all.rolwin21,contains("logp")))
## Make a rolling window table with just fdr, selects them and append loci
glm.all.fdr.rolwin21 <- cbind(glm.all.rolwin21.pval[,c(1:2)],select(glm.all.rolwin21,contains("fdr")))

## Merge annotations for filtering and highlighting
glm.all.rolwin21.annot <- merge(glm.all.rolwin21, vep, by=c("CHROM","POS"))

## Now merge to find All SNPs in and around candidate genes
glm.all.rolwin21.annot.spino <- merge(spino.cand, glm.all.rolwin21.annot, by="Gene")


#################################################
### Parallelism vs drift from founder to TPT1 ###
#################################################

## Here, we want to calculate population allele frequency differentiation from founders to 
## TPT1, then see whether within-condition (SE, SP, or S combined) are better correlated 
## (or more parallel) than between-condition correlations. Parallelism differentiates 
## directional selection from genetic drift, and by looking at loci with high divergence 
## from founder, correlation metrics between all population pairs will give us a sense of 
## how strong directional selection is (via high correlation metrics 
## within-condition/treatment).

## Make a new metadata table filtered to include only TPT1 S and E samples
## We only use these samples because they were both derived from the same "E" founder pop
haf.meta.T1.noPAfilt <- unique(haf.meta.T1filt[haf.meta.T1filt$batch == "a" & (haf.meta.T1filt$treat == "S" | haf.meta.T1filt$treat == "E"),])
haf.meta.T1filt <- haf.meta.T1filt[which((haf.meta.T1filt$samp %in% c(1:9,30:171))==TRUE),] #only spinosad experiment samples
haf.freq.T1.noPAfilt <- haf.freq.T1filt[, which((names(haf.freq.T1filt) %in% haf.meta.T1.noPAfilt$samp)==TRUE)] #match samples across tables


## Calculate E founder mean for all loci (E was the founder for E and S treatments)
head(haf.freq[,c(1:5)]) #samples 1, 2, and 3 are the E founder samples
haf.freq.Fmean <- cbind(haf.freq[,c(1:2)], F=rowMeans(haf.freq[,c(3:5)]))

## Use this color scheme for "condition" side colors
sidecols_EvS <- haf.meta.T1.noPAfilt$condition 
sidecols_EvS <- gsub("SP","#D9B851",sidecols_EvS)
sidecols_EvS <- gsub("SE","#848556",sidecols_EvS)
sidecols_EvS <- gsub("E","#D26183",sidecols_EvS)

### We want to filter based on the highest AF changes in any direction in any sample ###
### First, choose top 2% of AF changes in any SE, SP, or E samples at TPT1
haf.freq.Fmean.T1filt <- merge(haf.freq.Fmean, cbind(haf.sites.T1filt, haf.freq.T1.noPAfilt), by=c("CHROM","POS")) 

## Calculate absolute difference from founder mean
haf.freq.Fdiff.T1filt <- apply(abs(haf.freq.Fmean.T1filt[,c(4:ncol(haf.freq.Fmean.T1filt))] - haf.freq.Fmean.T1filt$F), 2, rank)

## Find the minimum difference from founder across all samples and apply ranking
min_T1_AFrank <- cbind(haf.freq.Fmean.T1filt[,c(1:2)], minRank =rowMins(haf.freq.Fdiff.T1filt))

## This filter will pull loci where any sample falls within the top 2% of AF changes
top_T1_AFdiff <- min_T1_AFrank[min_T1_AFrank$minRank < (0.02 * nrow(min_T1_AFrank)),]

## Retrieve the actual allele frequencies for top ranked loci
top_T1_AFdiff_table <- merge(top_T1_AFdiff[,-3], merge(haf.sites.T1filt,haf.freq.Fmean.T1filt,by=c("CHROM","POS")), by=c("CHROM","POS"))

## We want to calculate AF change from TPT0 and TPT1 for all cages and all filtered loci
top_T1_AFdiff_distance <- c()#create empty object
for(n in c(4:ncol(top_T1_AFdiff_table))) { #loop through each TPT1 cage
	temp_dist <- top_T1_AFdiff_table[,n] - top_T1_AFdiff_table[,3] #sample "n" - E T0 mean
  	#add row to table
  	top_T1_AFdiff_distance <- cbind(top_T1_AFdiff_distance,temp_dist)
  	#name new column
  	colnames(top_T1_AFdiff_distance)[n-3] <- colnames(top_T1_AFdiff_table)[n]
}

## Calculate Spearman correlation metrics between cages for all loci
top_T1_AFdiff_cormat <- cor(top_T1_AFdiff_distance, method = "spearman", use = "pairwise.complete.obs")
diag(top_T1_AFdiff_cormat) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_AFdiff_distance.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_AFdiff_cormat, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


### Second, choose top 2% of AF changes in any SE sample at TPT1
## Filter for SE samples only
haf.meta.T1.SEfilt <- unique(haf.meta.T1filt[haf.meta.T1filt$batch == "a" & (haf.meta.T1filt$condition == "SE"),]) #filter for SE samples only
haf.freq.T1.SEfilt <- haf.freq.T1filt[, which((names(haf.freq.T1filt) %in% haf.meta.T1.SEfilt$samp)==TRUE)] #match samples across tables

## Merge founder mean AF frequencies
haf.freq.Fmean.T1.SEfilt <- merge(haf.freq.Fmean, cbind(haf.sites.T1filt, haf.freq.T1.SEfilt), by=c("CHROM","POS")) 

## Calculate absolute difference from founder mean
haf.freq.Fdiff.T1.SEfilt <- apply(abs(haf.freq.Fmean.T1.SEfilt[,c(4:ncol(haf.freq.Fmean.T1.SEfilt))]-haf.freq.Fmean.T1.SEfilt$F), 2, rank)

## Find the minimum difference from founder across all SE samples and apply ranking
min_T1_SE_AFrank <- cbind(haf.freq.Fmean.T1.SEfilt[,c(1:2)], minRank =rowMins(haf.freq.Fdiff.T1filt))

## This filter will pull loci where any SE sample falls within the top 2% of AF changes
top_T1_SE_AFdiff <- min_T1_SE_AFrank[min_T1_SE_AFrank$minRank < (0.02 * nrow(min_T1_SE_AFrank)),]

## Retrieve the actual allele frequencies for top ranked loci
top_T1_SE_AFdiff_table <- merge(top_T1_SE_AFdiff[,-3], merge(haf.sites.T1filt,haf.freq.Fmean.T1filt,by=c("CHROM","POS")), by=c("CHROM","POS"))

## We want to calculate AF change from TPT0 and TPT1 for all cages and all filtered loci
top_T1_SE_AFdiff_distance <- c()#create empty object
for(n in c(4:ncol(top_T1_SE_AFdiff_table))) { #loop through each TPT1 SE cage
	temp_dist <- top_T1_SE_AFdiff_table[,n]-top_T1_SE_AFdiff_table[,3] #sample "n"-E T0 mean
  	#add row to table
  	top_T1_SE_AFdiff_distance <- cbind(top_T1_SE_AFdiff_distance,temp_dist)
  	#name new column
  	colnames(top_T1_SE_AFdiff_distance)[n-3] <- colnames(top_T1_SE_AFdiff_table)[n]
}

## Calculate Spearman correlation metrics between cages for all loci
top_T1_SE_AFdiff_cormat <- cor(top_T1_SE_AFdiff_distance, method = "spearman", use = "pairwise.complete.obs")
diag(top_T1_SE_AFdiff_cormat) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_SE_AFdiff_distance.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_SE_AFdiff_cormat, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


### Third, choose top 2% of AF changes in any SP sample at TPT1
## Filter for SP samples only
haf.meta.T1.SPfilt <- unique(haf.meta.T1filt[haf.meta.T1filt$batch == "a" & (haf.meta.T1filt$condition == "SP"),]) #filter for SP samples only
haf.freq.T1.SPfilt <- haf.freq.T1filt[, which((names(haf.freq.T1filt) %in% haf.meta.T1.SPfilt$samp)==TRUE)] #match samples across tables

## Merge founder mean AF frequencies
haf.freq.Fmean.T1.SPfilt <- merge(haf.freq.Fmean, cbind(haf.sites.T1filt, haf.freq.T1.SPfilt), by=c("CHROM","POS")) 

## Calculate absolute difference from founder mean
haf.freq.Fdiff.T1.SPfilt <- apply(abs(haf.freq.Fmean.T1.SPfilt[,c(4:ncol(haf.freq.Fmean.T1.SPfilt))]-haf.freq.Fmean.T1.SPfilt$F), 2, rank)

## Find the minimum difference from founder across all SP samples and apply ranking
min_T1_SP_AFrank <- cbind(haf.freq.Fmean.T1.SPfilt[,c(1:2)], minRank =rowMins(haf.freq.Fdiff.T1filt))

## This filter will pull loci where any SP sample falls within the top 2% of AF changes
top_T1_SP_AFdiff <- min_T1_SP_AFrank[min_T1_SP_AFrank$minRank < (0.02 * nrow(min_T1_SP_AFrank)),]

## Retrieve the actual allele frequencies for top ranked loci
top_T1_SP_AFdiff_table <- merge(top_T1_SP_AFdiff[,-3], merge(haf.sites.T1filt,haf.freq.Fmean.T1filt,by=c("CHROM","POS")), by=c("CHROM","POS"))

## We want to calculate AF change from TPT0 and TPT1 for all cages and all filtered loci
top_T1_SP_AFdiff_distance <- c()#create empty object
for(n in c(4:ncol(top_T1_SP_AFdiff_table))) { #loop through each TPT1 SP cage
	temp_dist <- top_T1_SP_AFdiff_table[,n]-top_T1_SP_AFdiff_table[,3] #sample "n"-E T0 mean
  	#add row to table
  	top_T1_SP_AFdiff_distance <- cbind(top_T1_SP_AFdiff_distance,temp_dist)
  	#name new column
  	colnames(top_T1_SP_AFdiff_distance)[n-3] <- colnames(top_T1_SP_AFdiff_table)[n]
}

## Calculate Spearman correlation metrics between cages for all loci
top_T1_SP_AFdiff_cormat <- cor(top_T1_SP_AFdiff_distance, method = "spearman", use = "pairwise.complete.obs")
diag(top_T1_SP_AFdiff_cormat) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_SP_AFdiff_distance.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_SP_AFdiff_cormat, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


### Fourth, choose top 2% of AF changes in any S sample (SE or SP) at TPT1
## Filter for SE or SP samples, no E
haf.meta.T1.Sfilt <- unique(haf.meta.T1filt[haf.meta.T1filt$batch == "a" & (haf.meta.T1filt$condition == "SE" | haf.meta.T1filt$condition == "SP"),])
haf.freq.T1.Sfilt <- haf.freq.T1filt[, which((names(haf.freq.T1filt) %in% haf.meta.T1.Sfilt$samp)==TRUE)] #match samples across tables

## Merge founder mean AF frequencies
haf.freq.Fmean.T1.Sfilt <- merge(haf.freq.Fmean, cbind(haf.sites.T1filt, haf.freq.T1.Sfilt), by=c("CHROM","POS")) 

## Calculate absolute difference from founder mean
haf.freq.Fdiff.T1.Sfilt <- apply(abs(haf.freq.Fmean.T1.Sfilt[,c(4:ncol(haf.freq.Fmean.T1.Sfilt))]-haf.freq.Fmean.T1.Sfilt$F), 2, rank)

## Find the minimum difference from founder across all S samples and apply ranking
min_T1_S_AFrank <- cbind(haf.freq.Fmean.T1.Sfilt[,c(1:2)], minRank =rowMins(haf.freq.Fdiff.T1filt))

## This filter will pull loci where any S sample falls within the top 2% of AF changes
top_T1_S_AFdiff <- min_T1_S_AFrank[min_T1_S_AFrank$minRank < (0.02 * nrow(min_T1_S_AFrank)),]

## Retrieve the actual allele frequencies for top ranked loci
top_T1_S_AFdiff_table <- merge(top_T1_S_AFdiff[,-3], merge(haf.sites.T1filt,haf.freq.Fmean.T1filt,by=c("CHROM","POS")), by=c("CHROM","POS"))

## We want to calculate AF change from TPT0 and TPT1 for all cages and all filtered loci
top_T1_S_AFdiff_distance <- c()#create empty object
for(n in c(4:ncol(top_T1_S_AFdiff_table))) { #loop through each TPT1 S cage
	temp_dist <- top_T1_S_AFdiff_table[,n] - top_T1_S_AFdiff_table[,3] #sample "n"-E T0 mean
  	#add row to table
  	top_T1_S_AFdiff_distance <- cbind(top_T1_S_AFdiff_distance,temp_dist)
  	#name new column
  	colnames(top_T1_S_AFdiff_distance)[n-3] <- colnames(top_T1_S_AFdiff_table)[n]
}

## Calculate Spearman correlation metrics between cages for all loci
top_T1_S_AFdiff_cormat <- cor(top_T1_S_AFdiff_distance, method = "spearman", use = "pairwise.complete.obs")
diag(top_T1_S_AFdiff_cormat) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_S_AFdiff_distance.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_S_AFdiff_cormat, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()

### Finally, use AF changes from TPT0 to TPT1 for all loci
## We want to calculate AF change from TPT0 and TPT1 for all cages and all loci
all_T1_AFdiff_distance <- c()#create empty object
for(n in c(4:ncol(haf.freq.Fmean.T1filt))) { #loop through each TPT1 cage
	temp_dist <- haf.freq.Fmean.T1filt[,n] - haf.freq.Fmean.T1filt[,3] #sample "n" - E T0 mean
  	#add row to table
  	all_T1_AFdiff_distance <- cbind(all_T1_AFdiff_distance,temp_dist)
  	#name new column
  	colnames(all_T1_AFdiff_distance)[n-3] <- colnames(haf.freq.Fmean.T1filt)[n]
}

## Calculate Spearman correlation metrics between cages for all loci
all_T1_AFdiff_cormat <- cor(all_T1_AFdiff_distance, method = "spearman", use = "pairwise.complete.obs")
diag(all_T1_AFdiff_cormat) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.all_T1_AFdiff_distance.heatmap.pdf", width=10, height=10)
	heatmap.2(all_T1_AFdiff_cormat, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


#########################################################################
### Look for non-parallel AF change via Chi-Square and Binomial tests ###
#########################################################################

## This time, we don't care as much about the level of AF change from TPT0 founders, aside 
## from the initial filtering of loci, just the direction of change.

## Grab allele frequencies for top 2% of AF changes in any SE, SP, or E samples at TPT1
temp_table <- merge(top_T1_AFdiff[,c(1:2)],haf.freq.Fmean.T1filt, by=c("CHROM","POS"))

## Set results matrix sizes as N samples X N samples
matsize <- ncol(haf.freq.Fmean.T1filt)-3
top_T1_chi_table <- matrix(nrow = matsize, ncol = matsize)#create empty object
top_T1_binomial_table <- matrix(nrow = matsize, ncol = matsize)#create empty object
top_T1_cor_table <- matrix(nrow = matsize, ncol = matsize)#create empty object

## Outer loop will cycle through all samples for matrix columns
for(n in c(4:ncol(haf.freq.Fmean.T1filt))) { 
	## Inner loop will cycle through all samples for matrix rows
	for(m in c(4:ncol(haf.freq.Fmean.T1filt))) { 
			x1 <- (temp_table[,n] - temp_table[,3]) < 0 #is col sample AF down from TPT0
			x2 <- (temp_table[,m] - temp_table[,3]) < 0 #is row sample AF down from TPT0
  			temp_df <- cbind(x1,x2) #join outcomes for comparison
  			temp_mat <- matrix(nrow = 2, ncol = 2) #make outcome table
  			temp_mat[1,1] <- nrow(temp_df[temp_df[,1]== "FALSE" & temp_df[,2]== "FALSE",])
  			temp_mat[1,2] <- nrow(temp_df[temp_df[,1]== "TRUE" & temp_df[,2]== "FALSE",])
  			temp_mat[2,1] <- nrow(temp_df[temp_df[,1]== "FALSE" & temp_df[,2]== "TRUE",])
  			temp_mat[2,2] <- nrow(temp_df[temp_df[,1]== "TRUE" & temp_df[,2]== "TRUE",])
			# Run Chi-square test
			test_results1 <- chisq.test(na.omit(temp_mat))
			# Save Chi-square statistic
			top_T1_chi_table[n-3,m-3] <- test_results1$statistic
			# Run Binomial test
			test_results2 <- binom.test(c(temp_mat[1,1] + temp_mat[2,2], 
								temp_mat[1,2] + temp_mat[2,1]), p = 1/2)
			# Save -log10(p)-value
			top_T1_binomial_table[n-3,m-3] <- test_results2$statistic
			#top_T1_binomial_table[n-3,m-3] <- -log10(test_results2$p.value)
			# Run binary correlation test
			top_T1_cor_table[n-3,m-3] <- cor(x1, x2)
	}
}

### Process Chi-square test results ###
## name rows and columns in chi-square test results table
rownames(top_T1_chi_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
colnames(top_T1_chi_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
diag(top_T1_chi_table) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_chi_table.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_chi_table, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()

## Experimental heatmap to see if we can see patterns in ranks of pairwise comparisons.
## This one is non-parametric and clusters agnostically to the overall range of values in 
## the results table.
top_T1_chi_rank <- apply(t(top_T1_chi_table), 2, rank) #find ranks within each column

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_chi_rank.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_chi_rank, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


### Process Binomial results ###
## name rows and columns in binomial test results table
rownames(top_T1_binomial_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
colnames(top_T1_binomial_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
diag(top_T1_binomial_table) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_binomial_table.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_binomial_table, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()

## Experimental heatmap to see if we can see patterns in ranks of pairwise comparisons.
## This one is non-parametric and clusters agnostically to the overall range of values in 
## the results table.
top_T1_binomial_rank <- apply(t(top_T1_binomial_table), 2, rank) #find ranks within each column

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_binomial_rank.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_binomial_rank, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


### Process binary correlation test results ###
## name rows and columns in correlation results table
rownames(top_T1_cor_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
colnames(top_T1_cor_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
diag(top_T1_cor_table) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_cor_table.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_cor_table, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()

## Experimental heatmap to see if we can see patterns in ranks of pairwise comparisons.
## This one is non-parametric and clusters agnostically to the overtop range of values in 
## the results table.
top_T1_cor_rank <- apply(t(top_T1_cor_table), 2, rank) #find ranks within each column

## Make the heatmap
pdf(file = "rudflies_2023_redo.top_T1_cor_rank.heatmap.pdf", width=10, height=10)
	heatmap.2(top_T1_cor_rank, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


### Instead of subsetting values by top AF change, what happens if we use all loci?
## Set results matrix sizes as N samples X N samples
matsize <- ncol(haf.freq.Fmean.T1filt)-3
all_T1_chi_table <- matrix(nrow = matsize, ncol = matsize)#create empty object
all_T1_binomial_table <- matrix(nrow = matsize, ncol = matsize)#create empty object
all_T1_cor_table <- matrix(nrow = matsize, ncol = matsize)#create empty object

## Outer loop will cycle through all samples for matrix columns
for(n in c(4:ncol(haf.freq.Fmean.T1filt))) { 
	## Inner loop will cycle through all samples for matrix rows
	for(m in c(4:ncol(haf.freq.Fmean.T1filt))) { 
			#is col sample AF down from TPT0, pull from table with all loci
			x1 <- (haf.freq.Fmean.T1filt[,n] - haf.freq.Fmean.T1filt[,3]) < 0
			#is row sample AF down from TPT0, pull from table with all loci
			x2 <- (haf.freq.Fmean.T1filt[,m] - haf.freq.Fmean.T1filt[,3]) < 0
  			temp_df <- cbind(x1,x2) #join outcomes for comparison
  			temp_mat <- matrix(nrow = 2, ncol = 2) #make outcome table
  			temp_mat[1,1] <- nrow(temp_df[temp_df[,1]== "FALSE" & temp_df[,2]== "FALSE",])
  			temp_mat[1,2] <- nrow(temp_df[temp_df[,1]== "TRUE" & temp_df[,2]== "FALSE",])
  			temp_mat[2,1] <- nrow(temp_df[temp_df[,1]== "FALSE" & temp_df[,2]== "TRUE",])
  			temp_mat[2,2] <- nrow(temp_df[temp_df[,1]== "TRUE" & temp_df[,2]== "TRUE",])
			# Run Chi-square test
			test_results1 <- chisq.test(na.omit(temp_mat))
			# Save Chi-square statistic
			all_T1_chi_table[n-3,m-3] <- test_results1$statistic
			# Run Binomial test
			test_results2 <- binom.test(c(temp_mat[1,1] + temp_mat[2,2], 
								temp_mat[1,2] + temp_mat[2,1]), p = 1/2)
			# Save -log10(p)-value
			all_T1_binomial_table[n-3,m-3] <- -log10(test_results2$p.value)
			# Run binary correlation test
			all_T1_cor_table[n-3,m-3] <- cor(x1, x2)
	}
}

### Process Chi-square test results ###
## name rows and columns in chi-square test results table
rownames(all_T1_chi_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
colnames(all_T1_chi_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
diag(all_T1_chi_table) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.all_T1_chi_table.heatmap.pdf", width=10, height=10)
	heatmap.2(all_T1_chi_table, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()

## Experimental heatmap to see if we can see patterns in ranks of pairwise comparisons.
## This one is non-parametric and clusters agnostically to the overall range of values in 
## the results table.
all_T1_chi_rank <- apply(t(all_T1_chi_table), 2, rank) #find ranks within each column

## Make the heatmap
pdf(file = "rudflies_2023_redo.all_T1_chi_rank.heatmap.pdf", width=10, height=10)
	heatmap.2(all_T1_chi_rank, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


### Process Binomial results ###
## name rows and columns in binomial test results table
rownames(all_T1_binomial_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
colnames(all_T1_binomial_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
diag(all_T1_binomial_table) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.all_T1_binomial_table.heatmap.pdf", width=10, height=10)
	heatmap.2(all_T1_binomial_table, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()

## Experimental heatmap to see if we can see patterns in ranks of pairwise comparisons.
## This one is non-parametric and clusters agnostically to the overall range of values in 
## the results table.
all_T1_binomial_rank <- apply(t(all_T1_binomial_table), 2, rank) #find ranks within each column

## Make the heatmap
pdf(file = "rudflies_2023_redo.all_T1_binomial_rank.heatmap.pdf", width=10, height=10)
	heatmap.2(all_T1_binomial_rank, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


### Process binary correlation test results ###
## name rows and columns in correlation results table
rownames(all_T1_cor_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
colnames(all_T1_cor_table) <- names(haf.freq.Fmean.T1filt)[-c(1:3)]
diag(all_T1_cor_table) <- NA #make the self-self correlations NA

## Make the heatmap
pdf(file = "rudflies_2023_redo.all_T1_cor_table.heatmap.pdf", width=10, height=10)
	heatmap.2(all_T1_cor_table, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()

## Experimental heatmap to see if we can see patterns in ranks of pairwise comparisons.
## This one is non-parametric and clusters agnostically to the overall range of values in 
## the results table.
all_T1_cor_rank <- apply(t(all_T1_cor_table), 2, rank) #find ranks within each column

## Make the heatmap
pdf(file = "rudflies_2023_redo.all_T1_cor_rank.heatmap.pdf", width=10, height=10)
	heatmap.2(all_T1_cor_rank, dendrogram="both", ColSideColors=sidecols_EvS, RowSideColors=sidecols_EvS, trace = "none")
dev.off()


######################################
### Load leave-one-out GLM results ###
######################################

## These GLM analyses were performed while iteratively dropping each S population and 
## running TPT1 filtered pairwise treatment contrasts. Due to the low N for both SE 
## (extinct) and SP (persistent) populations, we can assess whether any indivual S sample 
## had an outsized impact on results, compared to the GLM results utilizing all samples.

## Load all GLM results generated via the separate R script:
## "glm.rudflies2023.PAvSvSEvE.leave1out.r"

# Dropped cage 3 (SP)
glm_no3 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no3.txt", header=TRUE)
names(glm_no3) <- c("CHROM","POS","EvPA.T1_no3","EvSE.T1_no3","EvSP.T1_no3","PAvSE.T1_no3","PAvSP.T1_no3","SEvSP.T1_no3")

# Dropped cage 7 (SP)
glm_no7 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no7.txt", header=TRUE)
names(glm_no7) <- c("CHROM","POS","EvPA.T1_no7","EvSE.T1_no7","EvSP.T1_no7","PAvSE.T1_no7","PAvSP.T1_no7","SEvSP.T1_no7")

# Dropped cage 11 (SE)
glm_no11 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no11.txt", header=TRUE)
names(glm_no11) <- c("CHROM","POS","EvPA.T1_no11","EvSE.T1_no11","EvSP.T1_no11","PAvSE.T1_no11","PAvSP.T1_no11","SEvSP.T1_no11")

# Dropped cage 15 (SP)
glm_no15 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no15.txt", header=TRUE)
names(glm_no15) <- c("CHROM","POS","EvPA.T1_no15","EvSE.T1_no15","EvSP.T1_no15","PAvSE.T1_no15","PAvSP.T1_no15","SEvSP.T1_no15")

# Dropped cage 21 (SE)
glm_no21 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no21.txt", header=TRUE)
names(glm_no21) <- c("CHROM","POS","EvPA.T1_no21","EvSE.T1_no21","EvSP.T1_no21","PAvSE.T1_no21","PAvSP.T1_no21","SEvSP.T1_no21")

# Dropped cage 27 (SE)
glm_no27 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no27.txt", header=TRUE)
names(glm_no27) <- c("CHROM","POS","EvPA.T1_no27","EvSE.T1_no27","EvSP.T1_no27","PAvSE.T1_no27","PAvSP.T1_no27","SEvSP.T1_no27")

# Dropped cage 33 (SP)
glm_no33 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no33.txt", header=TRUE)
names(glm_no33) <- c("CHROM","POS","EvPA.T1_no33","EvSE.T1_no33","EvSP.T1_no33","PAvSE.T1_no33","PAvSP.T1_no33","SEvSP.T1_no33")

# Dropped cage 37 (SP)
glm_no37 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no37.txt", header=TRUE)
names(glm_no37) <- c("CHROM","POS","EvPA.T1_no37","EvSE.T1_no37","EvSP.T1_no37","PAvSE.T1_no37","PAvSP.T1_no37","SEvSP.T1_no37")

# Dropped cage 41 (SE)
glm_no41 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no41.txt", header=TRUE)
names(glm_no41) <- c("CHROM","POS","EvPA.T1_no41","EvSE.T1_no41","EvSP.T1_no41","PAvSE.T1_no41","PAvSP.T1_no41","SEvSP.T1_no41")

# Dropped cage 45 (SE)
glm_no45 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no45.txt", header=TRUE)
names(glm_no45) <- c("CHROM","POS","EvPA.T1_no45","EvSE.T1_no45","EvSP.T1_no45","PAvSE.T1_no45","PAvSP.T1_no45","SEvSP.T1_no45")

# No dropped cages for reference
glm_PAvSvSEvE_all <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt", header=TRUE)
glm_PAvSvSEvE_all <- na.omit(glm_PAvSvSEvE_all)
names(glm_PAvSvSEvE_all) <- c("CHROM","POS","EvPA.T1","EvSE.T1","EvSP.T1","PAvSE.T1","PAvSP.T1","SEvSP.T1")

## Make a master file of all results by merging each table by locus
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no3, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no7, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no11, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no15, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no21, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no27, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no33, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no37, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no41, by=c("CHROM","POS"))
glm_PAvSvSEvE_all <- merge(glm_PAvSvSEvE_all, glm_no45, by=c("CHROM","POS"))

# filter out chromosome 4
glm_PAvSvSEvE_all <- glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$CHROM != "4",]

## P-value correction
x=ncol(glm_PAvSvSEvE_all) #use number of columns
for(i in c(3:x)) { #start after loci columns and append FDR cols at the end of the table
  glm_PAvSvSEvE_all[,i] <- as.numeric(glm_PAvSvSEvE_all[,i])
  glm_PAvSvSEvE_all <- cbind(glm_PAvSvSEvE_all, p.adjust(glm_PAvSvSEvE_all[,i], method = "fdr"))
  colnames(glm_PAvSvSEvE_all)[i+(x-2)] <- paste(names(glm_PAvSvSEvE_all)[i], ".fdr", sep="")
}

## Calculate -log10 for p-values
y=ncol(glm_PAvSvSEvE_all) #use number of columns
for(i in 3:x) { #start after loci columns and append logp cols at the end of the table
  #find minimum non-zero p-value and divide by 2 to reassign to zero values
  min.temp <- min(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all[,i]!=0,i])/2 
  col.temp <- glm_PAvSvSEvE_all[,i]
  col.temp[col.temp==0] <- min.temp #reassign zero p-values for plotting only
  #perform -log10(p) calculations
  glm_PAvSvSEvE_all <- cbind(glm_PAvSvSEvE_all, as.numeric(-log10(glm_PAvSvSEvE_all[,i])))
  colnames(glm_PAvSvSEvE_all)[i+(y-2)] <- paste(names(glm_PAvSvSEvE_all)[i], ".logp", sep="")
}


#################################################################
### Check correlations of full glm results vs. leave-one-out  ###
### results for the "SEvSP" contrast.						  ###
#################################################################

## All samples vs. no cage 3 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1),
		 as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no3), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no3)
#t = 2373.5, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8829463 0.8836302
#sample estimates:
#      cor 
#0.8832887 

## All samples vs. no cage 7 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no7), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no7)
#t = 3141.9, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9279494 0.9283804
#sample estimates:
#      cor 
#0.9281652 

## All samples vs. no cage 11 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no11), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no11)
#t = 3863.7, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9505866 0.9508857
#sample estimates:
#      cor 
#0.9507364 

## All samples vs. no cage 15 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no15), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no15)
#t = 2097.7, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8568608 0.8576856
#sample estimates:
#      cor 
#0.8572738 

## All samples vs. no cage 21 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no21), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no21)
#t = 2430.5, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8874921 0.8881510
#sample estimates:
#     cor 
#0.887822 

## All samples vs. no cage 27 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no27), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no27)
#t = 3111.8, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9266982 0.9271364
#sample estimates:
#      cor 
#0.9269176 

## All samples vs. no cage 33 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no33), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no33)
#t = 2918.1, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9178477 0.9183365
#sample estimates:
#      cor 
#0.9180925 

## All samples vs. no cage 37 (SP)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no37), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no37)
#t = 3101.3, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9262571 0.9266978
#sample estimates:
#      cor 
#0.9264777

## All samples vs. no cage 41 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no41), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no41)
#t = 2214.7, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8688292 0.8695899
#sample estimates:
#      cor 
#0.8692101 

## All samples vs. no cage 45 (SE)
cor.test(as.numeric(glm_PAvSvSEvE_all$SEvSP.T1), as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no45), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvSEvE_all$SEvSP.T1) and as.numeric(glm_PAvSvSEvE_all$SEvSP.T1_no45)
#t = 3046.1, df = 1587113, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9238606 0.9243150
#sample estimates:
#      cor 
#0.9240881 



### Build Manhattan Plots ###

# No dropped cages for reference
manh_SEvSP_all <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") + 
	 # highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]$SEvSP.T1.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: all samples")

# Dropped cage 3 (SP)
manh_SEvSP_no3 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no3.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no3.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no3.fdr < 0.05,]$SEvSP.T1_no3.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 3 (SP)")

# Dropped cage 7 (SP)
manh_SEvSP_no7 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no7.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no7.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no7.fdr < 0.05,]$SEvSP.T1_no7.logp)) ,color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 7 (SP)")

# Dropped cage 11 (SE)
manh_SEvSP_no11 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no11.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no11.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no11.fdr < 0.05,]$SEvSP.T1_no11.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 11 (SE)")

# Dropped cage 15 (SP)
manh_SEvSP_no15 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no15.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no15.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no15.fdr < 0.05,]$SEvSP.T1_no15.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 15 (SP)")
  	
# Dropped cage 21 (SE)
manh_SEvSP_no21 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no21.logp)) +
	geom_line(alpha = 0.5, colour = "grey") + 
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no21.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no21.fdr < 0.05,]$SEvSP.T1_no21.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 21 (SE)")

# Dropped cage 27 (SE)
manh_SEvSP_no27 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no27.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no27.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no27.fdr < 0.05,]$SEvSP.T1_no27.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 27 (SE)")

# Dropped cage 33 (SP)
manh_SEvSP_no33 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no33.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no33.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no33.fdr < 0.05,]$SEvSP.T1_no33.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 33 (SP)")
  
# Dropped cage 37 (SP)
manh_SEvSP_no37 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no37.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no37.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no37.fdr < 0.05,]$SEvSP.T1_no37.logp)), color="black",linetype="dashed",linewidth=.25) +
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 37 (SP)")
  	
# Dropped cage 41 (SE)
manh_SEvSP_no41 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no41.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") + 
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no41.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no41.fdr < 0.05,]$SEvSP.T1_no41.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 41 (SE)")

# Dropped cage 45 (SE)
manh_SEvSP_no45 <- ggplot(glm_PAvSvSEvE_all, aes(POS, SEvSP.T1_no45.logp)) +
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	geom_point(data=na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1.fdr < 0.05,]), aes(POS, SEvSP.T1_no45.logp, color = "red"), size = 1) +
	# draw a dotted line at the minimum value of all significant loci (FDR < 0.05) 
	geom_hline(yintercept = min(na.omit(glm_PAvSvSEvE_all[glm_PAvSvSEvE_all$SEvSP.T1_no45.fdr < 0.05,]$SEvSP.T1_no45.logp)), color="black",linetype="dashed",linewidth=.25) +
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$SEvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("SE vs SP: no 45 (SE)")

 
### Figure S3 - Plot all leave-one-out manhattans together 
pdf(file = "rudflies_2023_redo.SEvSP.T1.LOO.manh.pdf", width=10, height=12)
	ggarrange(manh_SEvSP_all, manh_SEvSP_no3, manh_SEvSP_no7, manh_SEvSP_no33, manh_SEvSP_no37, manh_SEvSP_no11, manh_SEvSP_no15, manh_SEvSP_no21, manh_SEvSP_no27, manh_SEvSP_no41, manh_SEvSP_no45,
              ncol = 2, nrow = 6)
dev.off()


######################################################
### Build Manhattan plots for the "EvSP" contrast. ###
######################################################

# No dropped cages for reference
manh_EvSP_all <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: all samples")

# Dropped cage 3 (SP)
manh_EvSP_no3 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no3.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  	
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 3")

# Dropped cage 7 (SP)
manh_EvSP_no7 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no7.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# highlight significant loci (FDR < 0.05) from GLM contrast using all samples
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 7")
  	
# Dropped cage 15 (SP)
manh_EvSP_no15 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no15.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") + 
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 15")

# Dropped cage 33 (SP)
manh_EvSP_no33 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no33.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 33")

# Dropped cage 37 (SP)
manh_EvSP_no37 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSP.T1_no37.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSP.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SP: no 37")
  
### Unused in paper - Plot all leave-one-out manhattans together 
pdf(file = "rudflies_2023_redo.EvSP.T1.LOO.manh.pdf", width=10, height=16)
	ggarrange(manh_EvSP_all, manh_EvSP_no3, manh_EvSP_no7, manh_EvSP_no15, manh_EvSP_no33, manh_EvSP_no37,
              ncol = 1, nrow = 6)
dev.off()


############################################################
### Build bonus Manhattan plots for the "EvSE" contrast. ###
############################################################

# No dropped cages for reference
manh_EvSE_all <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: all samples")

# Dropped cage 11 (SE)
manh_EvSE_no11 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no11.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 11")

# Dropped cage 21 (SE)
manh_EvSE_no21 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no21.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 21")
  
# Dropped cage 27 (SE)
manh_EvSE_no27 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no27.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 27")
  	
# Dropped cage 41 (SE)
manh_EvSE_no41 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no41.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 41")
  	
# Dropped cage 45 (SE)
manh_EvSE_no45 <- ggplot(glm_PAvSvSEvE_all, aes(POS, EvSE.T1_no45.logp)) + 
	geom_line(alpha = 0.5, colour = "grey") +  
	# visual options
  	scale_fill_discrete(guide="none") +
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvSEvE_all$EvSE.T1.logp)))) +
  	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col="candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE: no 45")

### Unused in paper - Plot all leave-one-out manhattans together 
pdf(file = "rudflies_2023_redo.EvSE.T1.LOO.manh.pdf", width=10, height=16)
	ggarrange(manh_EvSE_all, manh_EvSE_no11, manh_EvSE_no21, manh_EvSE_no27, manh_EvSE_no41, manh_EvSE_no45,
              ncol = 1, nrow = 6)
dev.off()



#################
### Figure S9 ###
#################

## Goal here is to calculate pairwise distances of TPT4 SP vs E outlier SNPs among T1 E 
## populations, then and compare with distances between E and other TPT1 treatments.
## Mainly, we want to see if T4 adaptive convergent candidates are within the range of 
## negative control population variation to confirm they did not start divergent through 
## drift or founder effect.

## First, let's grab the allele frequencies for SNPs in the filtered list of GLM results
haf.freq.glm.filt <- merge(glm.all.rolwin21[,c(1:2)], haf.freq, by=c("CHROM","POS"))

## Select TPT4 SP vs E outliers (FDR < 0.05) from frequency table
haf.freq.glm.SvE.T4.filt <- na.omit(haf.freq.glm.filt[glm.all.rolwin21$SvE.T4.fdr.rolwin21 < 0.05,])

## How many are there?
dim(haf.freq.glm.SvE.T4.filt)
#[1] 1944  191

## Alternatively, just pick the TPT4 convergent candidates used for figure S12

## Retrieve "TPT4 convergent candidate loci" using the following criteria:
#	TPT4 SP vs E GLM contrast significant at FDR < 0.05 (1st target threshold)
#	TPT4 PA vs E GLM contrast significant at FDR < 0.01 (2nd target threshold)
#	TPT1 PA vs E GLM contrast significant at FDR < 0.01 (ensures PA & E start divergent)
#	TPT1 E allele frequency < 0.9 (not rare alleles with possible inflated AF diff)
#	TPT1 E allele frequency > 0.1 (not rare alleles with possible inflated AF diff)
#	TPT1 absolute AF difference between PA and E > 0.1 (sizable AF diff between PA & E)
glm.all.rolwin21.TPT4.convergent_candidates <- unique(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep, by=c("CHROM","POS")))

## Grab sites for TPT4 convergent candidates
haf.freq.glm.SvE.T4.filt <- merge(unique(glm.all.rolwin21.TPT4.convergent_candidates[,c(1:2)]), haf.freq, by=c("CHROM","POS"))

## How many are there?
dim(haf.freq.glm.SvE.T4.filt)
#[1] 179 191

### USE ONE OR THE OTHER SET (haf.freq.glm.SvE.T4.filt) FROM ABOVE ###

## Now we want to grab some key metadata columns for TPT1 treatment comparisons 
## Specifically, we'll grab sample ID and condition
haf.meta.T1filt.slim <- haf.meta.T1filt[,c(3,13)]

## Pull only TPT1 samples from the TPT4 convergent candidate AF table
haf.freq.glm.SvE.T4.cand.T1filt <- unique(subset(haf.freq.glm.SvE.T4.filt, select = as.character(haf.meta.T1filt$samp)))

## Make a Euclidean distance matrix from the filtered frequency table
haf.freq.glm.SvE.T4.cand.T1.dist <- as.matrix(dist(t(haf.freq.glm.SvE.T4.cand.T1filt), method="euclidean", diag=TRUE, upper=FALSE))

## Convert 0's to NA (since they are self-self comparisons and not useful)
haf.freq.glm.SvE.T4.cand.T1.dist[haf.freq.glm.SvE.T4.cand.T1.dist == 0] <- NA

## Melt the distance matrix into a dataframe that we can use for plotting
haf.freq.glm.SvE.T4.cand.T1.dist.df <- na.omit(reshape2::melt(as.matrix(haf.freq.glm.SvE.T4.cand.T1.dist), varnames = c("row", "col")))

## Now we want to grab some key metadata columns for TPT1 treatment comparisons 
## Specifically, we'll grab sample ID and condition
haf.meta.T1filt.slim <- haf.meta.T1filt[,c(3,13)]

## Merge them with the above distance dataframe by sample A
haf.freq.glm.SvE.T4.cand.T1.dist.df <- merge(haf.freq.glm.SvE.T4.cand.T1.dist.df, haf.meta.T1filt.slim, by.x="row", by.y="samp")

## Merge them with the above distance dataframe by sample B
haf.freq.glm.SvE.T4.cand.T1.dist.df <- merge(haf.freq.glm.SvE.T4.cand.T1.dist.df, haf.meta.T1filt.slim, by.x="col", by.y="samp")

## Now, using the metadata for samples A and B, we'll label treatment contrasts 
## First, we'll make a black treatment contrast
haf.freq.glm.SvE.T4.cand.T1.dist.df$contrast <- "NA"

## Label SEvSE
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SE" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SE",]$contrast <- "SEvSE"

## Label SPvSP
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SP" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SP",]$contrast <- "SPvSP"

## Label EvE
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "E" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "E",]$contrast <- "EvE"

## Label PAvPA
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "PA" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "PA",]$contrast <- "PAvPA"

## Label SEvSP
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SP" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SE",]$contrast <- "SEvSP"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SE" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SP",]$contrast <- "SEvSP"

## Label EvSE
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "E" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SE",]$contrast <- "EvSE"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SE" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "E",]$contrast <- "EvSE"

## Label EvSP
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "E" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SP",]$contrast <- "EvSP"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SP" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "E",]$contrast <- "EvSP"

## Label EvPA
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "E" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "PA",]$contrast <- "EvPA"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "PA" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "E",]$contrast <- "EvPA"

## Label PAvSE
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SE" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "PA",]$contrast <- "PAvSE"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "PA" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SE",]$contrast <- "PAvSE"

## Label PAvSP
haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "SP" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "PA",]$contrast <- "PAvSP"

haf.freq.glm.SvE.T4.cand.T1.dist.df[haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.x == "PA" & haf.freq.glm.SvE.T4.cand.T1.dist.df$condition.y == "SP",]$contrast <- "PAvSP"

## Format from plotting
haf.freq.glm.SvE.T4.dist.df <- as.data.frame(unique(haf.freq.glm.SvE.T4.cand.T1.dist.df[,c(6,3)]))
haf.freq.glm.SvE.T4.dist.df$value <- as.numeric(haf.freq.glm.SvE.T4.dist.df$value)


### Plot all contrasts
## Boxplot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.boxplot.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df, aes(x=contrast, y=value)) + 
	geom_boxplot(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()

## Jitter plot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.jitter.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df, aes(x=contrast, y=value)) + 
	geom_jitter(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()

### Plot without PA
## Take out PA-related contrasts
haf.freq.glm.SvE.T4.dist.df.noPA <- haf.freq.glm.SvE.T4.dist.df[haf.freq.glm.SvE.T4.dist.df$contrast != "EvPA" & haf.freq.glm.SvE.T4.dist.df$contrast != "PAvPA" & haf.freq.glm.SvE.T4.dist.df$contrast != "PAvSE" & haf.freq.glm.SvE.T4.dist.df$contrast != "PAvSP",]

## Boxplot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.noPA.boxplot.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df.noPA, aes(x=contrast, y=value)) + 
	geom_boxplot(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()

## Jitter plot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.noPA.jitter.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df.noPA, aes(x=contrast, y=value)) + 
	geom_jitter(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()


### Figure S9 ###
## Only within E and between SE and SP
haf.freq.glm.SvE.T4.dist.df.figS9 <- haf.freq.glm.SvE.T4.dist.df[haf.freq.glm.SvE.T4.dist.df$contrast == "EvE" | haf.freq.glm.SvE.T4.dist.df$contrast == "SEvSP",]

## Jitter plot
pdf(file = "haf.freq.glm.SvE.T4.T1dist.df.figS9.jitter.pdf")
	ggplot(data = haf.freq.glm.SvE.T4.dist.df.figS9, aes(x=contrast, y=value)) + 
	geom_jitter(aes(colour=contrast)) + 
	ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + 
	theme_classic()
dev.off()



########################################################################
### Mean AF differences	for leave-one-out and leave-one-in contrasts ###
########################################################################

## Let's calculate mean frequencies with and without each SE cage, then calculate pairwise 
## differences from other conditions. The resulting values will reflect both the magnitude 
## and direction of AF change. Later, we'll use these values look at concordance of 
## left-out and left-in populations, as well as assign positive or negative signs to 
## -log10p values. 

## Make list of all SE cages in T1
SE_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SE",]$cage
## Cycle through all SE T1 cages
for(i in SE_loo_cages) { 
	## Calculate pairwise differences of mean frequencies leaving out one SE cage
	assign(paste("freq_diff_no",i,sep=""),
		cbind(EvSE.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE" & haf.meta.T1filt$cage!=i]),
			EvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]),
			EvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SEvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE" & haf.meta.T1filt$cage!=i]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SPvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SEvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE" & haf.meta.T1filt$cage!=i]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"])))
	
	## For use in bedtools, we'll create a bed file from AF differences
	tempdf <- get(paste("freq_diff_no",i,sep="")) #grab above data frame
	tempdf2 <- cbind(haf.sites.T1filt,STOP=haf.sites.T1filt$POS+1, tempdf) #grab loci
	tempdf2 <- tempdf2[order(tempdf2[,1], tempdf2[,2]), ] #sort
	assign(paste("freq_diff_no",i,"_bed",sep=""),tempdf2) 

	## Calculate pairwise differences of mean frequencies leaving in only one SE cage
	assign(paste("freq_diff_only",i,sep=""),
		cbind(EvSE.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - haf.freq.T1filt[,haf.meta.T1filt$condition=="SE" & haf.meta.T1filt$cage==i],
			EvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]),
			EvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SEvPA.T1.diff=haf.freq.T1filt[,haf.meta.T1filt$condition=="SE" & haf.meta.T1filt$cage==i] - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SPvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SEvSP.T1.diff=haf.freq.T1filt[,haf.meta.T1filt$condition=="SE" & haf.meta.T1filt$cage==i] - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"])))

	## For use in bedtools, we'll create a bed file from AF differences
	tempdf <- get(paste("freq_diff_only",i,sep="")) #grab above data frame
	tempdf2 <- cbind(haf.sites.T1filt,STOP=haf.sites.T1filt$POS+1, tempdf) #grab loci
	tempdf2 <- tempdf2[order(tempdf2[,1], tempdf2[,2]), ] #sort
	assign(paste("freq_diff_only",i,"_bed",sep=""),tempdf2) 
}

## Let's calculate mean frequencies with and without each SP cage, then calculate pairwise 
## differences from other conditions. The resulting values will reflect both the magnitude 
## and direction of AF change. Later, we'll use these values look at concordance of 
## left-out and left-in populations, as well as assign positive or negative signs to 
## -log10p values. 

## Make list of all SP cages in T1
SP_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SP",]$cage
## Cycle through all SP T1 cages
for(i in SP_loo_cages) { 
	## Calculate pairwise differences of mean frequencies leaving out one SE cage
	assign(paste("freq_diff_no",i,sep=""),
		cbind(EvSE.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]),
			EvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP" & haf.meta.T1filt$cage!=i]),
			EvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SEvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SPvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP" & haf.meta.T1filt$cage!=i]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SEvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP" & haf.meta.T1filt$cage!=i])))
	
	## For use in bedtools, we'll create a bed file from AF differences
	tempdf <- get(paste("freq_diff_no",i,sep="")) #grab above data frame
	tempdf2 <- cbind(haf.sites.T1filt,STOP=haf.sites.T1filt$POS+1, tempdf) #grab loci
	tempdf2 <- tempdf2[order(tempdf2[,1], tempdf2[,2]), ] #sort
	assign(paste("freq_diff_no",i,"_bed",sep=""),tempdf2) 

	## Calculate pairwise differences of mean frequencies leaving in only one SP cage
	assign(paste("freq_diff_only",i,sep=""),
		cbind(EvSE.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]),
			EvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - haf.freq.T1filt[,haf.meta.T1filt$condition=="SP" & haf.meta.T1filt$cage==i],
			EvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SEvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SPvPA.T1.diff=haf.freq.T1filt[,haf.meta.T1filt$condition=="SP" & haf.meta.T1filt$cage==i] - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
			SEvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - haf.freq.T1filt[,haf.meta.T1filt$condition=="SP" & haf.meta.T1filt$cage==i]))

	## For use in bedtools, we'll create a bed file from AF differences
	tempdf <- get(paste("freq_diff_only",i,sep="")) #grab above data frame
	tempdf2 <- cbind(haf.sites.T1filt,STOP=haf.sites.T1filt$POS+1, tempdf) #grab loci
	tempdf2 <- tempdf2[order(tempdf2[,1], tempdf2[,2]), ] #sort
	assign(paste("freq_diff_only",i,"_bed",sep=""),tempdf2) 
}


###############################################
### Try to find clusters for data reduction ###
###############################################

## Start with raw T1 GLM contrast table and filter to remove chr4
contrast.PAvSvSEvE.table <- contrast.PAvSvSEvE.table[contrast.PAvSvSEvE.table$CHROM != "4",]

## Instead of the "EvPA.T1" column, let's merge and use the "PAvE.T1" from the time series 
## time-point GLM table since it has more p-value precision at lower p-values. 
contrast.PAvSvSEvE.table2 <- merge(contrast.PAvSvSEvE.table, contrastout.PAvE.table[,c(1,2,5)],by=c("CHROM","POS"))

## Now replace "EvPA.T1" column
contrast.PAvSvSEvE.table2$EvPA.T1 <- contrast.PAvSvSEvE.table2$PAvE.T1

## Remove original "PAvE.T1" column
contrast.PAvSvSEvE.table2 <- contrast.PAvSvSEvE.table2[,-9]

## Duplicate table as template to store corrected p-values
contrast.PAvSvSEvE.table.fdr <- contrast.PAvSvSEvE.table2
## P-value correction
for(i in c(3:8)) { #start after loci columns and append FDR cols at the end of the table
  contrast.PAvSvSEvE.table.fdr[,i] <- as.numeric(contrast.PAvSvSEvE.table.fdr[,i])
  contrast.PAvSvSEvE.table.fdr[,i] <- p.adjust(contrast.PAvSvSEvE.table.fdr[,i], method = "fdr")
}

## Initialize another table to store -log10(p) values with loci columns
contrast.PAvSvSEvE.table.logp <- contrast.PAvSvSEvE.table2[,c(1:2)]
## Calculate -log10 for p-values
for(i in 3:8) { #start after loci columns and append logp cols at the end of the table
  contrast.PAvSvSEvE.table2[,i] <- as.numeric(contrast.PAvSvSEvE.table2[,i])
  #find minimum non-zero p-value and divide by 2 to reassign to zero values
  min.temp <- min(contrast.PAvSvSEvE.table2[contrast.PAvSvSEvE.table2[,i]!=0,i])/2 
  col.temp <- contrast.PAvSvSEvE.table2[,i]
  col.temp[col.temp==0] <- min.temp #reassign zero p-values for plotting 
  #perform -log10(p) calculations
  contrast.PAvSvSEvE.table.logp <- cbind(contrast.PAvSvSEvE.table.logp, as.numeric(-log10(contrast.PAvSvSEvE.table2[,i])))
  colnames(contrast.PAvSvSEvE.table.logp)[i] <- paste(names(contrast.PAvSvSEvE.table2)[i], ".logp", sep="")
}

## Join FDR and delta frequency cols that will be used to calculate scores for clustering
glm_PAvSvSEvE_score <- merge(contrast.PAvSvSEvE.table.fdr, cbind(haf.sites.T1filt, abs(freq_diff_t1)), by=c("CHROM","POS"))
## Sort by locus
glm_PAvSvSEvE_score <- glm_PAvSvSEvE_score[order(glm_PAvSvSEvE_score[,1], glm_PAvSvSEvE_score[,2]), ]

## Assign scores based on a combination of significance level and mean frequency diffs
## These criteria were borrowed and modified from Rudman et. al, 2022:
## "Direct observation of adaptive tracking on ecological time scales in Drosophila"
for(f in c(3:8)) {
	glm_PAvSvSEvE_score[glm_PAvSvSEvE_score[,f] > 0.2,f+12] <- 0
	glm_PAvSvSEvE_score[glm_PAvSvSEvE_score[,f] < 0.2 | (glm_PAvSvSEvE_score[,f] > 0.2 & glm_PAvSvSEvE_score[,f+6] > 0.02),f+12] <- 1
	glm_PAvSvSEvE_score[glm_PAvSvSEvE_score[,f] < 0.05 & glm_PAvSvSEvE_score[,f+6] > 0.02,f+12] <- 2
	glm_PAvSvSEvE_score[glm_PAvSvSEvE_score[,f] < 0.01 & glm_PAvSvSEvE_score[,f+6] > 0.02,f+12] <- 3
	glm_PAvSvSEvE_score[glm_PAvSvSEvE_score[,f] < 0.001 & glm_PAvSvSEvE_score[,f+6] > 0.02,f+12] <- 4
	glm_PAvSvSEvE_score[glm_PAvSvSEvE_score[,f] < 0.0001 & glm_PAvSvSEvE_score[,f+6] > 0.02,f+12] <- 5
	glm_PAvSvSEvE_score[glm_PAvSvSEvE_score[,f] < 0.00001 & glm_PAvSvSEvE_score[,f+6] > 0.02,f+12] <- 6
	names(glm_PAvSvSEvE_score)[f+12] <- paste(names(glm_PAvSvSEvE_score[f]),".score",sep="")
}

## Main idea is testing whether window-based scores are significantly higher than randomly 
## assigned scores. Here, let's shuffle scores per chromosome to control for observed 
## variable signals of selective sweeps across chromosomes.
set.seed(42) 
glm_PAvSvSEvE_score_rand <- c() #initialize data frame

## Loop through all chromosomes
for(i in unique(glm_PAvSvSEvE_score$CHROM)) {
	shuffle_idx <- c()
	tempCHROM <- glm_PAvSvSEvE_score[glm_PAvSvSEvE_score$CHROM==i,]
	shuffle_idx <- append(shuffle_idx, sample(1:nrow(tempCHROM)))
	tempCHROM[, c("CHROM", "POS")] <- tempCHROM[shuffle_idx, c("CHROM", "POS")]
	glm_PAvSvSEvE_score_rand <- rbind(glm_PAvSvSEvE_score_rand,tempCHROM)
}

## Sort based on randomized loci
glm_PAvSvSEvE_score_rand <- glm_PAvSvSEvE_score_rand[order(glm_PAvSvSEvE_score_rand[,1], glm_PAvSvSEvE_score_rand[,2]), ]


### True loci window scores ###
## Create indexes for each chromosome to slide over and calculate mean scores over 501-SNP 
## windows. If we don't do this per chromosome, the average will be calculated across 
## consecutive chromosomes.

## 2L
glm_2L_true <- unique(glm_PAvSvSEvE_score[glm_PAvSvSEvE_score$CHROM=="2L",]) 
glm_2L_true <- glm_2L_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2L_true <- glm_2L_true[is.finite(rowSums(glm_2L_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2L_true$POS <- as.integer(glm_2L_true$POS)
glm_2L_true.rolwin501 <- glm_2L_true[,c(1,2)]

for(i in 3:ncol(glm_2L_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2L_true[is.finite(glm_2L_true[,i])=="TRUE",i])+100
	glm_2L_true[is.infinite(glm_2L_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2L_true.rolwin501 <- cbind(glm_2L_true.rolwin501,slide_mean(glm_2L_true[,i], before=250, after=250, step = 100))
    colnames(glm_2L_true.rolwin501)[i] <- paste(names(glm_2L_true)[i], ".rolwin501", sep="")
}

## 2R
glm_2R_true <- unique(glm_PAvSvSEvE_score[glm_PAvSvSEvE_score$CHROM=="2R",])
glm_2R_true <- glm_2R_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2R_true <- glm_2R_true[is.finite(rowSums(glm_2R_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2R_true$POS <- as.integer(glm_2R_true$POS)
glm_2R_true.rolwin501 <- glm_2R_true[,c(1,2)]

for(i in 3:ncol(glm_2R_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2R_true[is.finite(glm_2R_true[,i])=="TRUE",i])+100
	glm_2R_true[is.infinite(glm_2R_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2R_true.rolwin501 <- cbind(glm_2R_true.rolwin501,slide_mean(glm_2R_true[,i], before=250, after=250, step = 100))
    colnames(glm_2R_true.rolwin501)[i] <- paste(names(glm_2R_true)[i], ".rolwin501", sep="")
}

## 3L
glm_3L_true <- unique(glm_PAvSvSEvE_score[glm_PAvSvSEvE_score$CHROM=="3L",])
glm_3L_true <- glm_3L_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3L_true <- glm_3L_true[is.finite(rowSums(glm_3L_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3L_true$POS <- as.integer(glm_3L_true$POS)
glm_3L_true.rolwin501 <- glm_3L_true[,c(1,2)]

for(i in 3:ncol(glm_3L_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3L_true[is.finite(glm_3L_true[,i])=="TRUE",i])+100
	glm_3L_true[is.infinite(glm_3L_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3L_true.rolwin501 <- cbind(glm_3L_true.rolwin501,slide_mean(glm_3L_true[,i], before=250, after=250, step = 100))
    colnames(glm_3L_true.rolwin501)[i] <- paste(names(glm_3L_true)[i], ".rolwin501", sep="")
}

## 3R
glm_3R_true <- unique(glm_PAvSvSEvE_score[glm_PAvSvSEvE_score$CHROM=="3R",])
glm_3R_true <- glm_3R_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3R_true <- glm_3R_true[is.finite(rowSums(glm_3R_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3R_true$POS <- as.integer(glm_3R_true$POS)
glm_3R_true.rolwin501 <- glm_3R_true[,c(1,2)]

for(i in 3:ncol(glm_3R_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3R_true[is.finite(glm_3R_true[,i])=="TRUE",i])+100
	glm_3R_true[is.infinite(glm_3R_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3R_true.rolwin501 <- cbind(glm_3R_true.rolwin501,slide_mean(glm_3R_true[,i], before=250, after=250, step = 100))
    colnames(glm_3R_true.rolwin501)[i] <- paste(names(glm_3R_true)[i], ".rolwin501", sep="")
}

## X
glm_X_true <- unique(glm_PAvSvSEvE_score[glm_PAvSvSEvE_score$CHROM=="X",])
glm_X_true <- glm_X_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_X_true <- glm_X_true[is.finite(rowSums(glm_X_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_X_true$POS <- as.integer(glm_X_true$POS) #reformat
glm_X_true.rolwin501 <- glm_X_true[,c(1,2)]

for(i in 3:ncol(glm_X_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_X_true[is.finite(glm_X_true[,i])=="TRUE",i])+100
	glm_X_true[is.infinite(glm_X_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_X_true.rolwin501 <- cbind(glm_X_true.rolwin501,slide_mean(glm_X_true[,i], before=250, after=250, step = 100))
    colnames(glm_X_true.rolwin501)[i] <- paste(names(glm_X_true)[i], ".rolwin501", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm_PAvSvSEvE_score.rolwin501 <- na.omit(rbind(glm_2L_true.rolwin501,glm_2R_true.rolwin501,glm_3L_true.rolwin501,glm_3R_true.rolwin501,glm_X_true.rolwin501))


## How many sites remain after merging all filtered GLM results tables?
dim(glm_PAvSvSEvE_score.rolwin501)
#[1] 15197    20


### Random locus window scores ###
## Create indexes for each chromosome to slide over and calculate mean scores over 501-SNP 
## windows. If we don't do this per chromosome, the average will be calculated across 
## consecutive chromosomes.

## 2L
glm_2L_rand <- unique(glm_PAvSvSEvE_score_rand[glm_PAvSvSEvE_score_rand$CHROM=="2L",]) 
glm_2L_rand <- glm_2L_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2L_rand <- glm_2L_rand[is.finite(rowSums(glm_2L_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2L_rand$POS <- as.integer(glm_2L_rand$POS)
glm_2L_rand.rolwin501 <- glm_2L_rand[,c(1,2)]

for(i in 3:ncol(glm_2L_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2L_rand[is.finite(glm_2L_rand[,i])=="TRUE",i])+100
	glm_2L_rand[is.infinite(glm_2L_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2L_rand.rolwin501 <- cbind(glm_2L_rand.rolwin501,slide_mean(glm_2L_rand[,i], before=250, after=250, step = 100))
    colnames(glm_2L_rand.rolwin501)[i] <- paste(names(glm_2L_rand)[i], ".rolwin501", sep="")
}

## 2R
glm_2R_rand <- unique(glm_PAvSvSEvE_score_rand[glm_PAvSvSEvE_score_rand$CHROM=="2R",])
glm_2R_rand <- glm_2R_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2R_rand <- glm_2R_rand[is.finite(rowSums(glm_2R_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2R_rand$POS <- as.integer(glm_2R_rand$POS)
glm_2R_rand.rolwin501 <- glm_2R_rand[,c(1,2)]

for(i in 3:ncol(glm_2R_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2R_rand[is.finite(glm_2R_rand[,i])=="TRUE",i])+100
	glm_2R_rand[is.infinite(glm_2R_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2R_rand.rolwin501 <- cbind(glm_2R_rand.rolwin501,slide_mean(glm_2R_rand[,i], before=250, after=250, step = 100))
    colnames(glm_2R_rand.rolwin501)[i] <- paste(names(glm_2R_rand)[i], ".rolwin501", sep="")
}

## 3L
glm_3L_rand <- unique(glm_PAvSvSEvE_score_rand[glm_PAvSvSEvE_score_rand$CHROM=="3L",])
glm_3L_rand <- glm_3L_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3L_rand <- glm_3L_rand[is.finite(rowSums(glm_3L_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3L_rand$POS <- as.integer(glm_3L_rand$POS)
glm_3L_rand.rolwin501 <- glm_3L_rand[,c(1,2)]

for(i in 3:ncol(glm_3L_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3L_rand[is.finite(glm_3L_rand[,i])=="TRUE",i])+100
	glm_3L_rand[is.infinite(glm_3L_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3L_rand.rolwin501 <- cbind(glm_3L_rand.rolwin501,slide_mean(glm_3L_rand[,i], before=250, after=250, step = 100))
    colnames(glm_3L_rand.rolwin501)[i] <- paste(names(glm_3L_rand)[i], ".rolwin501", sep="")
}

## 3R
glm_3R_rand <- unique(glm_PAvSvSEvE_score_rand[glm_PAvSvSEvE_score_rand$CHROM=="3R",])
glm_3R_rand <- glm_3R_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3R_rand <- glm_3R_rand[is.finite(rowSums(glm_3R_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3R_rand$POS <- as.integer(glm_3R_rand$POS)
glm_3R_rand.rolwin501 <- glm_3R_rand[,c(1,2)]

for(i in 3:ncol(glm_3R_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3R_rand[is.finite(glm_3R_rand[,i])=="TRUE",i])+100
	glm_3R_rand[is.infinite(glm_3R_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3R_rand.rolwin501 <- cbind(glm_3R_rand.rolwin501,slide_mean(glm_3R_rand[,i], before=250, after=250, step = 100))
    colnames(glm_3R_rand.rolwin501)[i] <- paste(names(glm_3R_rand)[i], ".rolwin501", sep="")
}

## X
glm_X_rand <- unique(glm_PAvSvSEvE_score_rand[glm_PAvSvSEvE_score_rand$CHROM=="X",])
glm_X_rand <- glm_X_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_X_rand <- glm_X_rand[is.finite(rowSums(glm_X_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_X_rand$POS <- as.integer(glm_X_rand$POS) #reformat
glm_X_rand.rolwin501 <- glm_X_rand[,c(1,2)]

for(i in 3:ncol(glm_X_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_X_rand[is.finite(glm_X_rand[,i])=="TRUE",i])+100
	glm_X_rand[is.infinite(glm_X_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_X_rand.rolwin501 <- cbind(glm_X_rand.rolwin501,slide_mean(glm_X_rand[,i], before=250, after=250, step = 100))
    colnames(glm_X_rand.rolwin501)[i] <- paste(names(glm_X_rand)[i], ".rolwin501", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm_PAvSvSEvE_score_rand.rolwin501 <- na.omit(rbind(glm_2L_rand.rolwin501,glm_2R_rand.rolwin501,glm_3L_rand.rolwin501,glm_3R_rand.rolwin501,glm_X_rand.rolwin501))


## How many sites remain after merging all filtered GLM results tables?
dim(glm_PAvSvSEvE_score_rand.rolwin501)
#[1] 15197    20


### T-tests within windows ###
## Now we want to perform one-sided t-tests on scores within true locus windows and random 
## locus windows. Significant windows with true locus scores higher than random locus 
## scores will later be merged.

## Check all contrasts
for(i in 15:17) { #loop through score columns, only interested in contrasts with "E" pops
	## 2L
	glm_2L_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_2L_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_2L_ttest <- rbind(glm_2L_ttest,
    		cbind(glm_2L_true[j,c(1:2)],
    			#Set range
    			START=min(glm_2L_true[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),2]),
    			STOP=max(glm_2L_true[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),2]),
    			#Run test
    			t(t.test(glm_2L_true[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),i], glm_2L_rand[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_2L_true[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),i]),
    			med_score=median(glm_2L_true[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),i]),
    			max_score=max(glm_2L_true[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),i]),
    			rand_score=mean(glm_2L_rand[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),i]),
    			rand_med_score=median(glm_2L_true[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),i]),
    			rand_max_score=max(glm_2L_true[c(max(j-250,1):min(j+250,nrow(glm_2L_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_2L_ttest_",gsub(".score","",names(glm_2L_true)[i])), glm_2L_ttest)

	## 2R
	glm_2R_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_2R_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_2R_ttest <- rbind(glm_2R_ttest,
    		cbind(glm_2R_true[j,c(1:2)],
    			#Set range
    			START=min(glm_2R_true[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),2]),
    			STOP=max(glm_2R_true[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),2]),
    			#Run test
    			t(t.test(glm_2R_true[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),i], glm_2R_rand[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_2R_true[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),i]),
    			med_score=median(glm_2R_true[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),i]),
    			max_score=max(glm_2R_true[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),i]),
    			rand_score=mean(glm_2R_rand[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),i]),
    			rand_med_score=median(glm_2R_true[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),i]),
    			rand_max_score=max(glm_2R_true[c(max(j-250,1):min(j+250,nrow(glm_2R_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_2R_ttest_",gsub(".score","",names(glm_2R_true)[i])), glm_2R_ttest)

	## 3L
	glm_3L_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_3L_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_3L_ttest <- rbind(glm_3L_ttest,
    		cbind(glm_3L_true[j,c(1:2)],
    			#Set range
    			START=min(glm_3L_true[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),2]),
    			STOP=max(glm_3L_true[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),2]),
    			#Run test
    			t(t.test(glm_3L_true[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),i], glm_3L_rand[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_3L_true[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),i]),
    			med_score=median(glm_3L_true[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),i]),
    			max_score=max(glm_3L_true[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),i]),
    			rand_score=mean(glm_3L_rand[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),i]),
    			rand_med_score=median(glm_3L_true[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),i]),
    			rand_max_score=max(glm_3L_true[c(max(j-250,1):min(j+250,nrow(glm_3L_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_3L_ttest_",gsub(".score","",names(glm_3L_true)[i])), glm_3L_ttest)

	## 3R
	glm_3R_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_3R_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_3R_ttest <- rbind(glm_3R_ttest,
    		cbind(glm_3R_true[j,c(1:2)],
    			#Set range
    			START=min(glm_3R_true[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),2]),
    			STOP=max(glm_3R_true[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),2]),
    			#Run test
    			t(t.test(glm_3R_true[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),i], glm_3R_rand[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_3R_true[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),i]),
    			med_score=median(glm_3R_true[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),i]),
    			max_score=max(glm_3R_true[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),i]),
    			rand_score=mean(glm_3R_rand[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),i]),
    			rand_med_score=median(glm_3R_true[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),i]),
    			rand_max_score=max(glm_3R_true[c(max(j-250,1):min(j+250,nrow(glm_3R_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_3R_ttest_",gsub(".score","",names(glm_3R_true)[i])), glm_3R_ttest)

	## X
	glm_X_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_X_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_X_ttest <- rbind(glm_X_ttest,
    		cbind(glm_X_true[j,c(1:2)],
    			#Set range
    			START=min(glm_X_true[c(max(j-250,1):min(j+250,nrow(glm_X_true))),2]),
    			STOP=max(glm_X_true[c(max(j-250,1):min(j+250,nrow(glm_X_true))),2]),
    			#Run test
    			t(t.test(glm_X_true[c(max(j-250,1):min(j+250,nrow(glm_X_true))),i], glm_X_rand[c(max(j-250,1):min(j+250,nrow(glm_X_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_X_true[c(max(j-250,1):min(j+250,nrow(glm_X_true))),i]),
    			med_score=median(glm_X_true[c(max(j-250,1):min(j+250,nrow(glm_X_true))),i]),
    			max_score=max(glm_X_true[c(max(j-250,1):min(j+250,nrow(glm_X_true))),i]),
    			rand_score=mean(glm_X_rand[c(max(j-250,1):min(j+250,nrow(glm_X_true))),i]),
    			rand_med_score=median(glm_X_true[c(max(j-250,1):min(j+250,nrow(glm_X_true))),i]),
    			rand_max_score=max(glm_X_true[c(max(j-250,1):min(j+250,nrow(glm_X_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_X_ttest_",gsub(".score","",names(glm_X_true)[i])), glm_X_ttest)


	## Now, join stats from the the new rolling window tables for all chromosomes
    glm_PAvSvSEvE_score_ttest <- na.omit(rbind(glm_2L_ttest,glm_2R_ttest,glm_3L_ttest,glm_3R_ttest,glm_X_ttest))
    ## Reformat
    glm_PAvSvSEvE_score_ttest$statistic <- sapply(glm_PAvSvSEvE_score_ttest$statistic, toString)
    glm_PAvSvSEvE_score_ttest$parameter <- sapply(glm_PAvSvSEvE_score_ttest$parameter, toString)
    glm_PAvSvSEvE_score_ttest$p.value <- sapply(glm_PAvSvSEvE_score_ttest$p.value, toString)
    glm_PAvSvSEvE_score_ttest$stderr <- sapply(glm_PAvSvSEvE_score_ttest$stderr, toString)
    ## Correct t-test p-values
    glm_PAvSvSEvE_score_ttest$fdr <- p.adjust(glm_PAvSvSEvE_score_ttest$p.value, method = "fdr")    
    ## Name table for particular contrast 
    assign(paste0("glm_PAvSvSEvE_score_ttest_",
    	gsub(".score","",names(glm_X_true)[i])), glm_PAvSvSEvE_score_ttest)
    ## Save table for particular contrast 
    write.table(glm_PAvSvSEvE_score_ttest, file=paste0("glm_PAvSvSEvE_score_ttest_",
    	gsub(".score","",names(glm_X_true)[i]),".txt"),
    	sep = "\t", quote = FALSE, row.names = F)
    
    ## Filter windows for significant t-test results and minimum [max] score of 2.
    glm_PAvSvSEvE_score_ttest_filt <- glm_PAvSvSEvE_score_ttest[glm_PAvSvSEvE_score_ttest$fdr<0.05 & glm_PAvSvSEvE_score_ttest$max_score>=2,]
    
    ## We tried filtering by mean and median scores as well, but the EvSP had nothing.
    ## Using max scores within windows as filtering criteria worked for all contrasts.
    #glm_PAvSvSEvE_score_ttest_filt <- glm_PAvSvSEvE_score_ttest[glm_PAvSvSEvE_score_ttest$fdr<0.05 & glm_PAvSvSEvE_score_ttest$mean_score>=2,]
    
    ## Convert loci and window ranges to Genomic Ranges table
	gr <- GRanges(seqnames=glm_PAvSvSEvE_score_ttest_filt[,1],
		ranges=IRanges(as.integer(glm_PAvSvSEvE_score_ttest_filt[,3]), 	as.integer(glm_PAvSvSEvE_score_ttest_filt[,4])),
		strand="+",
		pos=glm_PAvSvSEvE_score_ttest_filt[,2])

	## Now merge all overlapping significant windows
	merged_gr <- reduce(gr)
    ## Name table for particular contrast 
	assign(paste0("glm_PAvSvSEvE_score_ttest_intervals_", 
		gsub(".score","",names(glm_X_true)[i])), merged_gr)
	## Save table for particular contrast 
	write.table(merged_gr, file=paste0("glm_PAvSvSEvE_score_ttest_intervals_",
		gsub(".score","",names(glm_X_true)[i]),".txt"),
		sep = "\t", quote = FALSE, row.names = F)
}

## How many merged clusters in each contrast?
dim(as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1))
#[1] 129   5
dim(as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1))
#[1] 35  5
dim(as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1))
#[1] 290   5


### Plot clusters per contrast ###

### E vs SE contrast
## Set each cluster as an alternating color
glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1)))

## Convert to normal data frame
glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1.df <- as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1)

## Use locus info to name each cluster
glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1.df$clust <- paste(glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1.df$seqnames,glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1.df$start,glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1.df$end,sep="_")

## Create a bed file using original SNP loci with "dummy" stop site
glm_PAvSvSEvE_bed <- as.data.frame(cbind(CHROM=as.character(contrast.PAvSvSEvE.table.logp$CHROM), POS=as.integer(contrast.PAvSvSEvE.table.logp$POS), STOP=as.integer(contrast.PAvSvSEvE.table.logp$POS + 1)))

## Sort bed file
glm_PAvSvSEvE_bed <- glm_PAvSvSEvE_bed[order(glm_PAvSvSEvE_bed[,1], as.numeric(glm_PAvSvSEvE_bed[,2])),]

## Create genomic ranges object with new bed file
glm_PAvSvSEvE_bed <- GRanges(seqnames=contrast.PAvSvSEvE.table.logp[,1],
		ranges=IRanges(as.integer(contrast.PAvSvSEvE.table.logp[,2]), 	as.integer(contrast.PAvSvSEvE.table.logp[,2]+1)),
		strand="+",
		pos=contrast.PAvSvSEvE.table.logp[,2])

## Merge original loci with intervals so we can assign cluster colors for plotting
glm_PAvSvSEvE_EvSE.T1_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvSEvE_bed, glm_PAvSvSEvE_score_ttest_intervals_EvSE.T1))

## Reformat merged table to retain only relevant columns
glm_PAvSvSEvE_EvSE.T1_overlap <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvSEvE_EvSE.T1_overlap[,1]), 
	POS=as.integer(glm_PAvSvSEvE_EvSE.T1_overlap[,2]),
		paste(glm_PAvSvSEvE_EvSE.T1_overlap[,8],
		glm_PAvSvSEvE_EvSE.T1_overlap[,9],
		glm_PAvSvSEvE_EvSE.T1_overlap[,10],sep="_"),
		color=glm_PAvSvSEvE_EvSE.T1_overlap[,13]))

## Rename cluster column header
names(glm_PAvSvSEvE_EvSE.T1_overlap)[3] = "clust"

## Finally, merge -log10(p) values for making manhattans
glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1 <- merge(contrast.PAvSvSEvE.table.logp, glm_PAvSvSEvE_EvSE.T1_overlap, by=c("CHROM","POS"))

## Save table for future use
write.table(glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1, file="glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1.txt", sep = "\t", quote = FALSE, row.names = F)

## Reload table if desired
#glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1.txt", header=TRUE)

# This manhattan plot shows the contrast between E and SE samples in TPT 1
manh.EvSE.T1.clust <- ggplot(contrast.PAvSvSEvE.table.logp, aes(POS, EvSE.T1.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Highlight clusters significant in SE vs E GLM contrast
	geom_point(data=glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1,aes(POS, EvSE.T1.logp), color = glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1$color, size = 0.1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(contrast.PAvSvSEvE.table.logp$EvSE.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SE (TPT1)")


## E vs SP
## Set each cluster as an alternating color
glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1)))

## Convert to normal data frame
glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1.df <- as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1)

## Use locus info to name each cluster
glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1.df$clust <- paste(glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1.df$seqnames,glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1.df$start,glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1.df$end,sep="_")

## Create a bed file using original SNP loci with "dummy" stop site
glm_PAvSvSEvE_bed <- as.data.frame(cbind(CHROM=as.character(contrast.PAvSvSEvE.table.logp$CHROM), POS=as.integer(contrast.PAvSvSEvE.table.logp$POS), STOP=as.integer(contrast.PAvSvSEvE.table.logp$POS + 1)))

## Sort bed file
glm_PAvSvSEvE_bed <- glm_PAvSvSEvE_bed[order(glm_PAvSvSEvE_bed[,1], as.numeric(glm_PAvSvSEvE_bed[,2])),]

## Create genomic ranges object with new bed file
glm_PAvSvSEvE_bed <- GRanges(seqnames=contrast.PAvSvSEvE.table.logp[,1],
		ranges=IRanges(as.integer(contrast.PAvSvSEvE.table.logp[,2]), 	as.integer(contrast.PAvSvSEvE.table.logp[,2]+1)),
		strand="+",
		pos=contrast.PAvSvSEvE.table.logp[,2])

## Merge original loci with intervals so we can assign cluster colors for plotting
glm_PAvSvSEvE_EvSP.T1_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvSEvE_bed, glm_PAvSvSEvE_score_ttest_intervals_EvSP.T1))

## Reformat merged table to retain only relevant columns
glm_PAvSvSEvE_EvSP.T1_overlap <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvSEvE_EvSP.T1_overlap[,1]), 
	POS=as.integer(glm_PAvSvSEvE_EvSP.T1_overlap[,2]),
		paste(glm_PAvSvSEvE_EvSP.T1_overlap[,8],
		glm_PAvSvSEvE_EvSP.T1_overlap[,9],
		glm_PAvSvSEvE_EvSP.T1_overlap[,10],sep="_"),
		color=glm_PAvSvSEvE_EvSP.T1_overlap[,13]))

## Rename cluster column header
names(glm_PAvSvSEvE_EvSP.T1_overlap)[3] = "clust"

## Finally, merge -log10(p) values for making manhattans
glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1 <- merge(contrast.PAvSvSEvE.table.logp, glm_PAvSvSEvE_EvSP.T1_overlap, by=c("CHROM","POS"))

## Save table for future use
write.table(glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1, file="glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1.txt", sep = "\t", quote = FALSE, row.names = F)

## Reload table if desired
#glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1.txt", header=TRUE)

# This manhattan plot shows the contrast between E and SP samples in TPT 1
manh.EvSP.T1.clust <- ggplot(contrast.PAvSvSEvE.table.logp, aes(POS, EvSP.T1.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Highlight clusters significant in SP vs E GLM contrast
	geom_point(data=glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1,aes(POS, EvSP.T1.logp), color = glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1$color, size = 0.1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(contrast.PAvSvSEvE.table.logp$EvSP.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SP (TPT1)")


## E vs PA
## Set each cluster as an alternating color
glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1)))

## Convert to normal data frame
glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1.df <- as.data.frame(glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1)

## Use locus info to name each cluster
glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1.df$clust <- paste(glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1.df$seqnames,glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1.df$start,glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1.df$end,sep="_")

## Create a bed file using original SNP loci with "dummy" stop site
glm_PAvSvSEvE_bed <- as.data.frame(cbind(CHROM=as.character(contrast.PAvSvSEvE.table.logp$CHROM), POS=as.integer(contrast.PAvSvSEvE.table.logp$POS), STOP=as.integer(contrast.PAvSvSEvE.table.logp$POS + 1)))

## Sort bed file
glm_PAvSvSEvE_bed <- glm_PAvSvSEvE_bed[order(glm_PAvSvSEvE_bed[,1], as.numeric(glm_PAvSvSEvE_bed[,2])),]

## Create genomic ranges object with new bed file
glm_PAvSvSEvE_bed <- GRanges(seqnames=contrast.PAvSvSEvE.table.logp[,1],
		ranges=IRanges(as.integer(contrast.PAvSvSEvE.table.logp[,2]), 	as.integer(contrast.PAvSvSEvE.table.logp[,2]+1)),
		strand="+",
		pos=contrast.PAvSvSEvE.table.logp[,2])

## Merge original loci with intervals so we can assign cluster colors for plotting
glm_PAvSvSEvE_EvPA.T1_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvSEvE_bed, glm_PAvSvSEvE_score_ttest_intervals_EvPA.T1))

## Reformat merged table to retain only relevant columns
glm_PAvSvSEvE_EvPA.T1_overlap <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvSEvE_EvPA.T1_overlap[,1]), 
	POS=as.integer(glm_PAvSvSEvE_EvPA.T1_overlap[,2]),
		paste(glm_PAvSvSEvE_EvPA.T1_overlap[,8],
		glm_PAvSvSEvE_EvPA.T1_overlap[,9],
		glm_PAvSvSEvE_EvPA.T1_overlap[,10],sep="_"),
		color=glm_PAvSvSEvE_EvPA.T1_overlap[,13]))

## Rename cluster column header
names(glm_PAvSvSEvE_EvPA.T1_overlap)[3] = "clust"

## Finally, merge -log10(p) values for making manhattans
glm_PAvSvSEvE_score_ttest_clusters_EvPA.T1 <- merge(contrast.PAvSvSEvE.table.logp, glm_PAvSvSEvE_EvPA.T1_overlap, by=c("CHROM","POS"))

## Save table for future use
write.table(glm_PAvSvSEvE_score_ttest_clusters_EvPA.T1, file="glm_PAvSvSEvE_score_ttest_clusters_EvPA.T1.txt", sep = "\t", quote = FALSE, row.names = F)

## Reload table if desired
#glm_PAvSvSEvE_score_ttest_clusters_EvPA.T1 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/glm_PAvSvSEvE_score_ttest_clusters_EvPA.T1.txt", header=TRUE)

# This manhattan plot shows the contrast between E and PA samples in TPT 1
manh.PAvE.T1.clust <- ggplot(contrast.PAvSvSEvE.table.logp, aes(POS, EvPA.T1.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Highlight clusters significant in PA vs E GLM contrast
	geom_point(data=glm_PAvSvSEvE_score_ttest_clusters_EvPA.T1,aes(POS, EvPA.T1.logp), color = glm_PAvSvSEvE_score_ttest_clusters_EvPA.T1$color, size = 0.1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(contrast.PAvSvSEvE.table.logp$EvPA.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA (TPT1)")


## Plot the three relevant pairwise TPT1 treatment comparisons together
pdf(file = "rudflies_2023_redo.EvSE.EvSP.PAvE.T1.rolwin501clust.glm.manh.pdf", width=7.5, height=6)
	ggarrange(manh.EvSE.T1.clust, manh.EvSP.T1.clust, manh.PAvE.T1.clust,
              ncol = 1, nrow = 3)
dev.off()



##########################################
### TPT4 S-cage leave-one-out analysis ###
##########################################

## These GLM analyses were performed while iteratively dropping each S population and 
## running TPT1 filtered pairwise treatment contrasts. Due to the low N for both SE 
## (extinct) and SP (persistent) populations, we can assess whether any indivual S sample 
## had an outsized impact on results, compared to the GLM results utilizing all samples.

## Load all GLM results generated via the separate R script:
## "glm.rudflies2023.PAvSvE.TPT4.leave1out.sh"

# Dropped cage 3 (SP)
glm_T4_no3 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvE.wLociGLMcontrast.TPT4.LOO.no3.txt", header=TRUE)
names(glm_T4_no3) <- c("CHROM","POS","EvPA.T4_no3","EvSP.T4_no3","PAvSP.T4_no3")

# Dropped cage 7 (SP)
glm_T4_no7 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvE.wLociGLMcontrast.TPT4.LOO.no7.txt", header=TRUE)
names(glm_T4_no7) <- c("CHROM","POS","EvPA.T4_no7","EvSP.T4_no7","PAvSP.T4_no7")

# Dropped cage 15 (SP)
glm_T4_no15 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvE.wLociGLMcontrast.TPT4.LOO.no15.txt", header=TRUE)
names(glm_T4_no15) <- c("CHROM","POS","EvPA.T4_no15","EvSP.T4_no15","PAvSP.T4_no15")

# Dropped cage 33 (SP)
glm_T4_no33 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvE.wLociGLMcontrast.TPT4.LOO.no33.txt", header=TRUE)
names(glm_T4_no33) <- c("CHROM","POS","EvPA.T4_no33","EvSP.T4_no33","PAvSP.T4_no33")

# Dropped cage 37 (SP)
glm_T4_no37 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_PAvSvE.wLociGLMcontrast.TPT4.LOO.no37.txt", header=TRUE)
names(glm_T4_no37) <- c("CHROM","POS","EvPA.T4_no37","EvSP.T4_no37","PAvSP.T4_no37")

# No dropped cages for reference
glm_PAvSvE_TPT4 <- merge(merge(cbind(contrastout.SvE.table[,c(1:2)],EvSP.T4=contrastout.SvE.table$SvE.T4),cbind(contrastout.PAvE.table[,c(1:2)],EvPA.T4=contrastout.PAvE.table$PAvE.T4),by=c("CHROM","POS")),cbind(contrast.PAvS[,c(1:2)],PAvSP.T4=contrast.PAvS$PAvS.T4),by=c("CHROM","POS"))

## Make a master file of all results by merging each table by locus
glm_PAvSvE_TPT4 <- merge(glm_PAvSvE_TPT4, glm_T4_no3, by=c("CHROM","POS"))
glm_PAvSvE_TPT4 <- merge(glm_PAvSvE_TPT4, glm_T4_no7, by=c("CHROM","POS"))
glm_PAvSvE_TPT4 <- merge(glm_PAvSvE_TPT4, glm_T4_no15, by=c("CHROM","POS"))
glm_PAvSvE_TPT4 <- merge(glm_PAvSvE_TPT4, glm_T4_no33, by=c("CHROM","POS"))
glm_PAvSvE_TPT4 <- merge(glm_PAvSvE_TPT4, glm_T4_no37, by=c("CHROM","POS"))

# filter out chromosome 4
glm_PAvSvE_TPT4 <- glm_PAvSvE_TPT4[glm_PAvSvE_TPT4$CHROM != "4",]

dim(glm_PAvSvE_TPT4)
#[1] 1425215      20

## P-value correction
x=ncol(glm_PAvSvE_TPT4) #use number of columns
for(i in c(3:x)) { #start after loci columns and append FDR cols at the end of the table
  glm_PAvSvE_TPT4[,i] <- as.numeric(glm_PAvSvE_TPT4[,i]) #format p-vals to numeric
  glm_PAvSvE_TPT4 <- cbind(glm_PAvSvE_TPT4, p.adjust(glm_PAvSvE_TPT4[,i], method = "fdr"))
  colnames(glm_PAvSvE_TPT4)[i+(x-2)] <- paste(names(glm_PAvSvE_TPT4)[i], ".fdr", sep="")
}

## Calculate -log10 for p-values
y=ncol(glm_PAvSvE_TPT4) #use number of columns
for(i in 3:x) { #start after loci columns and append logp cols at the end of the table
  #find minimum non-zero p-value and divide by 2 to reassign to zero values
  min.temp <- min(glm_PAvSvE_TPT4[glm_PAvSvE_TPT4[,i]!=0,i])/2 
  col.temp <- glm_PAvSvE_TPT4[,i]
  col.temp[col.temp==0] <- min.temp #reassign zero p-values for plotting only
  #perform -log10(p) calculations
  glm_PAvSvE_TPT4 <- cbind(glm_PAvSvE_TPT4, as.numeric(-log10(glm_PAvSvE_TPT4[,i])))
  colnames(glm_PAvSvE_TPT4)[i+(y-2)] <- paste(names(glm_PAvSvE_TPT4)[i], ".logp", sep="")
}


#################################################################
### Check correlations of full glm results vs. leave-one-out  ###
### results for the "EvSP" contrast.						  ###
#################################################################

## All samples vs. no cage 3 (SP)
cor.test(as.numeric(glm_PAvSvE_TPT4$EvSP.T4),
		 as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no3), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvE_TPT4$EvSP.T4) and as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no3)
#t = 1137.4, df = 1425213, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.6889275 0.6906487
#sample estimates:
#      cor 
#0.6897891 


## All samples vs. no cage 7 (SP)
cor.test(as.numeric(glm_PAvSvE_TPT4$EvSP.T4), as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no7), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvE_TPT4$EvSP.T4) and as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no7)
#t = 1317.3, df = 1425213, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7402303 0.7417110
#sample estimates:
#      cor 
#0.7409716 


## All samples vs. no cage 15 (SP)
cor.test(as.numeric(glm_PAvSvE_TPT4$EvSP.T4), as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no15), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvE_TPT4$EvSP.T4) and as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no15)
#t = 1005.4, df = 1425213, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.6432157 0.6451367
#sample estimates:
#      cor 
#0.6441772 


## All samples vs. no cage 33 (SP)
cor.test(as.numeric(glm_PAvSvE_TPT4$EvSP.T4), as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no33), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvE_TPT4$EvSP.T4) and as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no33)
#t = 1493.6, df = 1425213, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7805036 0.7817835
#sample estimates:
#      cor 
#0.7811444 


## All samples vs. no cage 37 (SP)
cor.test(as.numeric(glm_PAvSvE_TPT4$EvSP.T4), as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no37), method="pearson")
#	Pearson's product-moment correlation
#
#data:  as.numeric(glm_PAvSvE_TPT4$EvSP.T4) and as.numeric(glm_PAvSvE_TPT4$EvSP.T4_no37)
#t = 781.62, df = 1425213, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.5466094 0.5489077
#sample estimates:
#      cor 
#0.5477596


########################################################################
### Mean AF differences	for leave-one-out and leave-one-in contrasts ###
########################################################################

## Let's calculate mean frequencies with and without each SP cage, then calculate pairwise 
## differences from other conditions. The resulting values will reflect both the magnitude 
## and direction of AF change. Later, we'll use these values look at concordance of 
## left-out and left-in populations, as well as assign positive or negative signs to 
## -log10p values. 

## Make list of all SP cages in T4
SP_loo_cages <- haf.meta.T4filt[haf.meta.T4filt$treat.fix=="S",]$cage
## Cycle through all SP T4 cages
for(i in SP_loo_cages) { 
	## Calculate pairwise differences of mean frequencies leaving out one SE cage
	assign(paste("freq_diff_t4_no",i,sep=""),
		cbind(EvSP.T4.diff=rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="E"]) - rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="S" & haf.meta.T4filt$cage!=i]),
			EvPA.T4.diff=rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="E"]) - rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="PA"]),
			SPvPA.T4.diff=rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="S" & haf.meta.T4filt$cage!=i]) - rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="PA"])))
	
	## For use in bedtools, we'll create a bed file from AF differences
	tempdf <- get(paste("freq_diff_t4_no",i,sep="")) #grab above data frame
	tempdf2 <- cbind(haf.sites.T4filt,STOP=haf.sites.T4filt$POS+1, tempdf) #grab loci
	tempdf2 <- tempdf2[order(tempdf2[,1], tempdf2[,2]), ] #sort
	assign(paste("freq_diff_t4_no",i,"_bed",sep=""),tempdf2) 

	## Calculate pairwise differences of mean frequencies leaving in only one SP cage
	assign(paste("freq_diff_t4_only",i,sep=""),
		cbind(EvSP.T4.diff=rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="E"]) - haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="S" & haf.meta.T4filt$cage==i],
			EvPA.T4.diff=rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="E"]) - rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="PA"]),
			SPvPA.T4.diff=haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="S" & haf.meta.T4filt$cage==i] - rowMeans(haf.freq.T4filt[,haf.meta.T4filt$treat.fix=="PA"])))

	## For use in bedtools, we'll create a bed file from AF differences
	tempdf <- get(paste("freq_diff_t4_only",i,sep="")) #grab above data frame
	tempdf2 <- cbind(haf.sites.T4filt,STOP=haf.sites.T4filt$POS+1, tempdf) #grab loci
	tempdf2 <- tempdf2[order(tempdf2[,1], tempdf2[,2]), ] #sort
	assign(paste("freq_diff_t4_only",i,"_bed",sep=""),tempdf2) 
}


###############################################
### Try to find clusters for data reduction ###
###############################################

## Select locus and FDR cols only
glm_PAvSvE_TPT4.fdr <- cbind(glm_PAvSvE_TPT4[,c(1:2)],select(glm_PAvSvE_TPT4,contains("fdr")))

## Select locus and -log10(p) cols only
glm_PAvSvE_TPT4.logp <- cbind(glm_PAvSvE_TPT4[,c(1:2)],select(glm_PAvSvE_TPT4,contains("logp")))

## Join FDR and delta frequency cols that will be used to calculate scores for clustering
glm_PAvSvE_TPT4_score <- merge(glm_PAvSvE_TPT4.fdr[,c(1:5)], cbind(haf.sites.T4filt, abs(freq_diff_t4)), by=c("CHROM","POS"))
## Sort by locus
glm_PAvSvE_TPT4_score <- glm_PAvSvE_TPT4_score[order(glm_PAvSvE_TPT4_score[,1], glm_PAvSvE_TPT4_score[,2]), ]

## Assign scores based on a combination of significance level and mean frequency diffs
## These criteria were borrowed and modified from Rudman et. al, 2022:
## "Direct observation of adaptive tracking on ecological time scales in Drosophila"
for(f in c(3:5)) {
	glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score[,f] > 0.2,f+6] <- 0
	glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score[,f] < 0.2 | (glm_PAvSvE_TPT4_score[,f] > 0.2 & glm_PAvSvE_TPT4_score[,f+3] > 0.02),f+6] <- 1
	glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score[,f] < 0.05 & glm_PAvSvE_TPT4_score[,f+3] > 0.02,f+6] <- 2
	glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score[,f] < 0.01 & glm_PAvSvE_TPT4_score[,f+3] > 0.02,f+6] <- 3
	glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score[,f] < 0.001 & glm_PAvSvE_TPT4_score[,f+3] > 0.02,f+6] <- 4
	glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score[,f] < 0.0001 & glm_PAvSvE_TPT4_score[,f+3] > 0.02,f+6] <- 5
	glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score[,f] < 0.00001 & glm_PAvSvE_TPT4_score[,f+6] > 0.02,f+6] <- 6
	names(glm_PAvSvE_TPT4_score)[f+6] <- paste(names(glm_PAvSvE_TPT4_score[f]),".score",sep="")
	names(glm_PAvSvE_TPT4_score)[f+6] <- gsub(".fdr","",names(glm_PAvSvE_TPT4_score)[f+6])
}

## Main idea is testing whether window-based scores are significantly higher than randomly 
## assigned scores. Here, let's shuffle scores per chromosome to control for observed 
## variable signals of selective sweeps across chromosomes.
set.seed(42) 
glm_PAvSvE_TPT4_score_rand <- c() #initialize data frame

## Loop through all chromosomes
for(i in unique(glm_PAvSvE_TPT4_score$CHROM)) {
	shuffle_idx <- c()
	tempCHROM <- glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score$CHROM==i,]
	shuffle_idx <- append(shuffle_idx, sample(1:nrow(tempCHROM)))
	tempCHROM[, c("CHROM", "POS")] <- tempCHROM[shuffle_idx, c("CHROM", "POS")]
	glm_PAvSvE_TPT4_score_rand <- rbind(glm_PAvSvE_TPT4_score_rand,tempCHROM)
}

## Sort based on randomized loci
glm_PAvSvE_TPT4_score_rand <- glm_PAvSvE_TPT4_score_rand[order(glm_PAvSvE_TPT4_score_rand[,1], glm_PAvSvE_TPT4_score_rand[,2]), ]


### True loci window scores ###
## Create indexes for each chromosome to slide over and calculate mean scores over 501-SNP 
## windows. If we don't do this per chromosome, the average will be calculated across 
## consecutive chromosomes.

## 2L
glm_T4_2L_true <- unique(glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score$CHROM=="2L",]) 
glm_T4_2L_true <- glm_T4_2L_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_2L_true <- glm_T4_2L_true[is.finite(rowSums(glm_T4_2L_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_2L_true$POS <- as.integer(glm_T4_2L_true$POS)
glm_T4_2L_true.rolwin501 <- glm_T4_2L_true[,c(1,2)]

for(i in 3:ncol(glm_T4_2L_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_2L_true[is.finite(glm_T4_2L_true[,i])=="TRUE",i])+100
	glm_T4_2L_true[is.infinite(glm_T4_2L_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_2L_true.rolwin501 <- cbind(glm_T4_2L_true.rolwin501,slide_mean(glm_T4_2L_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_2L_true.rolwin501)[i] <- paste(names(glm_T4_2L_true)[i], ".rolwin501", sep="")
}

## 2R
glm_T4_2R_true <- unique(glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score$CHROM=="2R",])
glm_T4_2R_true <- glm_T4_2R_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_2R_true <- glm_T4_2R_true[is.finite(rowSums(glm_T4_2R_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_2R_true$POS <- as.integer(glm_T4_2R_true$POS)
glm_T4_2R_true.rolwin501 <- glm_T4_2R_true[,c(1,2)]

for(i in 3:ncol(glm_T4_2R_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_2R_true[is.finite(glm_T4_2R_true[,i])=="TRUE",i])+100
	glm_T4_2R_true[is.infinite(glm_T4_2R_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_2R_true.rolwin501 <- cbind(glm_T4_2R_true.rolwin501,slide_mean(glm_T4_2R_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_2R_true.rolwin501)[i] <- paste(names(glm_T4_2R_true)[i], ".rolwin501", sep="")
}

## 3L
glm_T4_3L_true <- unique(glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score$CHROM=="3L",])
glm_T4_3L_true <- glm_T4_3L_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_3L_true <- glm_T4_3L_true[is.finite(rowSums(glm_T4_3L_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_3L_true$POS <- as.integer(glm_T4_3L_true$POS)
glm_T4_3L_true.rolwin501 <- glm_T4_3L_true[,c(1,2)]

for(i in 3:ncol(glm_T4_3L_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_3L_true[is.finite(glm_T4_3L_true[,i])=="TRUE",i])+100
	glm_T4_3L_true[is.infinite(glm_T4_3L_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_3L_true.rolwin501 <- cbind(glm_T4_3L_true.rolwin501,slide_mean(glm_T4_3L_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_3L_true.rolwin501)[i] <- paste(names(glm_T4_3L_true)[i], ".rolwin501", sep="")
}

## 3R
glm_T4_3R_true <- unique(glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score$CHROM=="3R",])
glm_T4_3R_true <- glm_T4_3R_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_3R_true <- glm_T4_3R_true[is.finite(rowSums(glm_T4_3R_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_3R_true$POS <- as.integer(glm_T4_3R_true$POS)
glm_T4_3R_true.rolwin501 <- glm_T4_3R_true[,c(1,2)]

for(i in 3:ncol(glm_T4_3R_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_3R_true[is.finite(glm_T4_3R_true[,i])=="TRUE",i])+100
	glm_T4_3R_true[is.infinite(glm_T4_3R_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_3R_true.rolwin501 <- cbind(glm_T4_3R_true.rolwin501,slide_mean(glm_T4_3R_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_3R_true.rolwin501)[i] <- paste(names(glm_T4_3R_true)[i], ".rolwin501", sep="")
}

## X
glm_T4_X_true <- unique(glm_PAvSvE_TPT4_score[glm_PAvSvE_TPT4_score$CHROM=="X",])
glm_T4_X_true <- glm_T4_X_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_X_true <- glm_T4_X_true[is.finite(rowSums(glm_T4_X_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_X_true$POS <- as.integer(glm_T4_X_true$POS) #reformat
glm_T4_X_true.rolwin501 <- glm_T4_X_true[,c(1,2)]

for(i in 3:ncol(glm_T4_X_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_X_true[is.finite(glm_T4_X_true[,i])=="TRUE",i])+100
	glm_T4_X_true[is.infinite(glm_T4_X_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_X_true.rolwin501 <- cbind(glm_T4_X_true.rolwin501,slide_mean(glm_T4_X_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_X_true.rolwin501)[i] <- paste(names(glm_T4_X_true)[i], ".rolwin501", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm_PAvSvE_TPT4_score.rolwin501 <- na.omit(rbind(glm_T4_2L_true.rolwin501,glm_T4_2R_true.rolwin501,glm_T4_3L_true.rolwin501,glm_T4_3R_true.rolwin501,glm_T4_X_true.rolwin501))


## How many sites remain after merging all filtered GLM results tables?
dim(glm_PAvSvE_TPT4_score.rolwin501)
#[1] 14253    11


### Random locus window scores ###
## Create indexes for each chromosome to slide over and calculate mean scores over 501-SNP 
## windows. If we don't do this per chromosome, the average will be calculated across 
## consecutive chromosomes.

## 2L
glm_T4_2L_rand <- unique(glm_PAvSvE_TPT4_score_rand[glm_PAvSvE_TPT4_score_rand$CHROM=="2L",]) 
glm_T4_2L_rand <- glm_T4_2L_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_2L_rand <- glm_T4_2L_rand[is.finite(rowSums(glm_T4_2L_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_2L_rand$POS <- as.integer(glm_T4_2L_rand$POS)
glm_T4_2L_rand.rolwin501 <- glm_T4_2L_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_2L_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_2L_rand[is.finite(glm_T4_2L_rand[,i])=="TRUE",i])+100
	glm_T4_2L_rand[is.infinite(glm_T4_2L_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_2L_rand.rolwin501 <- cbind(glm_T4_2L_rand.rolwin501,slide_mean(glm_T4_2L_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_2L_rand.rolwin501)[i] <- paste(names(glm_T4_2L_rand)[i], ".rolwin501", sep="")
}

## 2R
glm_T4_2R_rand <- unique(glm_PAvSvE_TPT4_score_rand[glm_PAvSvE_TPT4_score_rand$CHROM=="2R",])
glm_T4_2R_rand <- glm_T4_2R_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_2R_rand <- glm_T4_2R_rand[is.finite(rowSums(glm_T4_2R_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_2R_rand$POS <- as.integer(glm_T4_2R_rand$POS)
glm_T4_2R_rand.rolwin501 <- glm_T4_2R_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_2R_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_2R_rand[is.finite(glm_T4_2R_rand[,i])=="TRUE",i])+100
	glm_T4_2R_rand[is.infinite(glm_T4_2R_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_2R_rand.rolwin501 <- cbind(glm_T4_2R_rand.rolwin501,slide_mean(glm_T4_2R_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_2R_rand.rolwin501)[i] <- paste(names(glm_T4_2R_rand)[i], ".rolwin501", sep="")
}

## 3L
glm_T4_3L_rand <- unique(glm_PAvSvE_TPT4_score_rand[glm_PAvSvE_TPT4_score_rand$CHROM=="3L",])
glm_T4_3L_rand <- glm_T4_3L_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_3L_rand <- glm_T4_3L_rand[is.finite(rowSums(glm_T4_3L_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_3L_rand$POS <- as.integer(glm_T4_3L_rand$POS)
glm_T4_3L_rand.rolwin501 <- glm_T4_3L_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_3L_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_3L_rand[is.finite(glm_T4_3L_rand[,i])=="TRUE",i])+100
	glm_T4_3L_rand[is.infinite(glm_T4_3L_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_3L_rand.rolwin501 <- cbind(glm_T4_3L_rand.rolwin501,slide_mean(glm_T4_3L_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_3L_rand.rolwin501)[i] <- paste(names(glm_T4_3L_rand)[i], ".rolwin501", sep="")
}

## 3R
glm_T4_3R_rand <- unique(glm_PAvSvE_TPT4_score_rand[glm_PAvSvE_TPT4_score_rand$CHROM=="3R",])
glm_T4_3R_rand <- glm_T4_3R_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_3R_rand <- glm_T4_3R_rand[is.finite(rowSums(glm_T4_3R_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_3R_rand$POS <- as.integer(glm_T4_3R_rand$POS)
glm_T4_3R_rand.rolwin501 <- glm_T4_3R_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_3R_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_3R_rand[is.finite(glm_T4_3R_rand[,i])=="TRUE",i])+100
	glm_T4_3R_rand[is.infinite(glm_T4_3R_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_3R_rand.rolwin501 <- cbind(glm_T4_3R_rand.rolwin501,slide_mean(glm_T4_3R_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_3R_rand.rolwin501)[i] <- paste(names(glm_T4_3R_rand)[i], ".rolwin501", sep="")
}

## X
glm_T4_X_rand <- unique(glm_PAvSvE_TPT4_score_rand[glm_PAvSvE_TPT4_score_rand$CHROM=="X",])
glm_T4_X_rand <- glm_T4_X_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_X_rand <- glm_T4_X_rand[is.finite(rowSums(glm_T4_X_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_X_rand$POS <- as.integer(glm_T4_X_rand$POS) #reformat
glm_T4_X_rand.rolwin501 <- glm_T4_X_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_X_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_X_rand[is.finite(glm_T4_X_rand[,i])=="TRUE",i])+100
	glm_T4_X_rand[is.infinite(glm_T4_X_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_X_rand.rolwin501 <- cbind(glm_T4_X_rand.rolwin501,slide_mean(glm_T4_X_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_X_rand.rolwin501)[i] <- paste(names(glm_T4_X_rand)[i], ".rolwin501", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm_PAvSvE_TPT4_score_rand.rolwin501 <- na.omit(rbind(glm_T4_2L_rand.rolwin501,glm_T4_2R_rand.rolwin501,glm_T4_3L_rand.rolwin501,glm_T4_3R_rand.rolwin501,glm_T4_X_rand.rolwin501))

## How many sites remain after merging all filtered GLM results tables?
dim(glm_PAvSvE_TPT4_score_rand.rolwin501)
#[1] 14253    11



### T-tests within windows ###
## Now we want to perform one-sided t-tests on scores within true locus windows and random 
## locus windows. Significant windows with true locus scores higher than random locus 
## scores will later be merged.

## Check all contrasts
for(i in 9:11) { #loop through score columns
	## 2L
	glm_T4_2L_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_2L_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_2L_ttest <- rbind(glm_T4_2L_ttest,
    		cbind(glm_T4_2L_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_2L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),2]),
    			STOP=max(glm_T4_2L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),2]),
    			#Run test
    			t(t.test(glm_T4_2L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),i], glm_T4_2L_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_2L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),i]),
    			med_score=median(glm_T4_2L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),i]),
    			max_score=max(glm_T4_2L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),i]),
    			rand_score=mean(glm_T4_2L_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),i]),
    			rand_med_score=median(glm_T4_2L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),i]),
    			rand_max_score=max(glm_T4_2L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_2L_ttest_",gsub(".score","",names(glm_T4_2L_true)[i])), glm_T4_2L_ttest)

	## 2R
	glm_T4_2R_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_2R_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_2R_ttest <- rbind(glm_T4_2R_ttest,
    		cbind(glm_T4_2R_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_2R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),2]),
    			STOP=max(glm_T4_2R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),2]),
    			#Run test
    			t(t.test(glm_T4_2R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),i], glm_T4_2R_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_2R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),i]),
    			med_score=median(glm_T4_2R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),i]),
    			max_score=max(glm_T4_2R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),i]),
    			rand_score=mean(glm_T4_2R_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),i]),
    			rand_med_score=median(glm_T4_2R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),i]),
    			rand_max_score=max(glm_T4_2R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_2R_ttest_",gsub(".score","",names(glm_T4_2R_true)[i])), glm_T4_2R_ttest)

	## 3L
	glm_T4_3L_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_3L_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_3L_ttest <- rbind(glm_T4_3L_ttest,
    		cbind(glm_T4_3L_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_3L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),2]),
    			STOP=max(glm_T4_3L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),2]),
    			#Run test
    			t(t.test(glm_T4_3L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),i], glm_T4_3L_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_3L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),i]),
    			med_score=median(glm_T4_3L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),i]),
    			max_score=max(glm_T4_3L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),i]),
    			rand_score=mean(glm_T4_3L_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),i]),
    			rand_med_score=median(glm_T4_3L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),i]),
    			rand_max_score=max(glm_T4_3L_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_3L_ttest_",gsub(".score","",names(glm_T4_3L_true)[i])), glm_T4_3L_ttest)

	## 3R
	glm_T4_3R_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_3R_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_3R_ttest <- rbind(glm_T4_3R_ttest,
    		cbind(glm_T4_3R_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_3R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),2]),
    			STOP=max(glm_T4_3R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),2]),
    			#Run test
    			t(t.test(glm_T4_3R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),i], glm_T4_3R_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_3R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),i]),
    			med_score=median(glm_T4_3R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),i]),
    			max_score=max(glm_T4_3R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),i]),
    			rand_score=mean(glm_T4_3R_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),i]),
    			rand_med_score=median(glm_T4_3R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),i]),
    			rand_max_score=max(glm_T4_3R_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_3R_ttest_",gsub(".score","",names(glm_T4_3R_true)[i])), glm_T4_3R_ttest)

	## X
	glm_T4_X_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_X_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_X_ttest <- rbind(glm_T4_X_ttest,
    		cbind(glm_T4_X_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_X_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),2]),
    			STOP=max(glm_T4_X_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),2]),
    			#Run test
    			t(t.test(glm_T4_X_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),i], glm_T4_X_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_X_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),i]),
    			med_score=median(glm_T4_X_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),i]),
    			max_score=max(glm_T4_X_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),i]),
    			rand_score=mean(glm_T4_X_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),i]),
    			rand_med_score=median(glm_T4_X_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),i]),
    			rand_max_score=max(glm_T4_X_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_X_ttest_",gsub(".score","",names(glm_T4_X_true)[i])), glm_T4_X_ttest)


	## Now, join stats from the the new rolling window tables for all chromosomes
    glm_PAvSvE_TPT4_score_ttest <- na.omit(rbind(glm_T4_2L_ttest,glm_T4_2R_ttest,glm_T4_3L_ttest,glm_T4_3R_ttest,glm_T4_X_ttest))
    ## Reformat
    glm_PAvSvE_TPT4_score_ttest$statistic <- sapply(glm_PAvSvE_TPT4_score_ttest$statistic, toString)
    glm_PAvSvE_TPT4_score_ttest$parameter <- sapply(glm_PAvSvE_TPT4_score_ttest$parameter, toString)
    glm_PAvSvE_TPT4_score_ttest$p.value <- sapply(glm_PAvSvE_TPT4_score_ttest$p.value, toString)
    glm_PAvSvE_TPT4_score_ttest$stderr <- sapply(glm_PAvSvE_TPT4_score_ttest$stderr, toString)
    ## Correct t-test p-values
    glm_PAvSvE_TPT4_score_ttest$fdr <- p.adjust(glm_PAvSvE_TPT4_score_ttest$p.value, method = "fdr")    
    ## Name table for particular contrast 
    assign(paste0("glm_PAvSvE_TPT4_score_ttest_",
    	gsub(".score","",names(glm_T4_X_true)[i])), glm_PAvSvE_TPT4_score_ttest)
    ## Save table for particular contrast 
    write.table(glm_PAvSvE_TPT4_score_ttest, file=paste0("glm_PAvSvE_TPT4_score_ttest_",
    	gsub(".score","",names(glm_T4_X_true)[i]),".txt"),
    	sep = "\t", quote = FALSE, row.names = F)
    
    ## Filter windows for significant t-test results and minimum [max] score of 2.
    glm_PAvSvE_TPT4_score_ttest_filt <- glm_PAvSvE_TPT4_score_ttest[glm_PAvSvE_TPT4_score_ttest$fdr<0.05 & glm_PAvSvE_TPT4_score_ttest$max_score>=2,]
    
    ## We tried filtering by mean and median scores as well, but the EvSP had nothing.
    ## Using max scores within windows as filtering criteria worked for all contrasts.
    #glm_PAvSvE_TPT4_score_ttest_filt <- glm_PAvSvE_TPT4_score_ttest[glm_PAvSvE_TPT4_score_ttest$fdr<0.05 & glm_PAvSvE_TPT4_score_ttest$mean_score>=2,]
    
    ## Convert loci and window ranges to Genomic Ranges table
	gr <- GRanges(seqnames=glm_PAvSvE_TPT4_score_ttest_filt[,1],
		ranges=IRanges(as.integer(glm_PAvSvE_TPT4_score_ttest_filt[,3]), 	as.integer(glm_PAvSvE_TPT4_score_ttest_filt[,4])),
		strand="+",
		pos=glm_PAvSvE_TPT4_score_ttest_filt[,2])

	## Now merge all overlapping significant windows
	merged_gr <- reduce(gr)
    ## Name table for particular contrast 
	assign(paste0("glm_PAvSvE_TPT4_score_ttest_intervals_", 
		gsub(".score","",names(glm_T4_X_true)[i])), merged_gr)
	## Save table for particular contrast 
	write.table(merged_gr, file=paste0("glm_PAvSvE_TPT4_score_ttest_intervals_",
		gsub(".score","",names(glm_T4_X_true)[i]),".txt"),
		sep = "\t", quote = FALSE, row.names = F)
}

## How many merged clusters in each contrast?
dim(as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4))
#[1] 372   5
dim(as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4))
#[1] 255   5
dim(as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4))
#[1] 291   5


### Plot clusters per contrast ###

## E vs SP
## Set each cluster as an alternating color
glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4)))

## Convert to normal data frame
glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4.df <- as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4)

## Use locus info to name each cluster
glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4.df$clust <- paste(glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4.df$seqnames,glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4.df$start,glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4.df$end,sep="_")

## Create a bed file using original SNP loci with "dummy" stop site
glm_PAvSvE_TPT4_bed <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvE_TPT4.logp$CHROM), POS=as.integer(glm_PAvSvE_TPT4.logp$POS), STOP=as.integer(glm_PAvSvE_TPT4.logp$POS + 1)))

## Sort bed file
glm_PAvSvE_TPT4_bed <- glm_PAvSvE_TPT4_bed[order(glm_PAvSvE_TPT4_bed[,1], as.numeric(glm_PAvSvE_TPT4_bed[,2])),]

## Create genomic ranges object with new bed file
glm_PAvSvE_TPT4_bed <- GRanges(seqnames=glm_PAvSvE_TPT4.logp[,1],
		ranges=IRanges(as.integer(glm_PAvSvE_TPT4.logp[,2]), 	as.integer(glm_PAvSvE_TPT4.logp[,2]+1)),
		strand="+",
		pos=glm_PAvSvE_TPT4.logp[,2])

## Merge original loci with intervals so we can assign cluster colors for plotting
glm_PAvSvE_TPT4_EvSP.T4_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvE_TPT4_bed, glm_PAvSvE_TPT4_score_ttest_intervals_EvSP.T4))

## Reformat merged table to retain only relevant columns
glm_PAvSvE_TPT4_EvSP.T4_overlap <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvE_TPT4_EvSP.T4_overlap[,1]), 
	POS=as.integer(glm_PAvSvE_TPT4_EvSP.T4_overlap[,2]),
		paste(glm_PAvSvE_TPT4_EvSP.T4_overlap[,8],
		glm_PAvSvE_TPT4_EvSP.T4_overlap[,9],
		glm_PAvSvE_TPT4_EvSP.T4_overlap[,10],sep="_"),
		color=glm_PAvSvE_TPT4_EvSP.T4_overlap[,13]))

## Rename cluster column header
names(glm_PAvSvE_TPT4_EvSP.T4_overlap)[3] = "clust"

## Finally, merge -log10(p) values for making manhattans
glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4 <- merge(glm_PAvSvE_TPT4.logp, glm_PAvSvE_TPT4_EvSP.T4_overlap, by=c("CHROM","POS"))

## Save table for future use
write.table(glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4, file="glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4.txt", sep = "\t", quote = FALSE, row.names = F)

## Reload table if desired
#glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4.txt", header=TRUE)

# This manhattan plot shows the contrast between E and SP samples in TPT 1
manh.EvSP.T4.clust <- ggplot(glm_PAvSvE_TPT4.logp, aes(POS, EvSP.T4.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Highlight clusters significant in SP vs E GLM contrast
	geom_point(data=glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4,aes(POS, EvSP.T4.logp), color = glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4$color, size = 0.1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvE_TPT4.logp$EvSP.T4.logp))+0.5)) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SP (TPT4)")


## E vs PA
## Set each cluster as an alternating color
glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4)))

## Convert to normal data frame
glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4.df <- as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4)

## Use locus info to name each cluster
glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4.df$clust <- paste(glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4.df$seqnames,glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4.df$start,glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4.df$end,sep="_")

## Create a bed file using original SNP loci with "dummy" stop site
glm_PAvSvE_TPT4_bed <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvE_TPT4.logp$CHROM), POS=as.integer(glm_PAvSvE_TPT4.logp$POS), STOP=as.integer(glm_PAvSvE_TPT4.logp$POS + 1)))

## Sort bed file
glm_PAvSvE_TPT4_bed <- glm_PAvSvE_TPT4_bed[order(glm_PAvSvE_TPT4_bed[,1], as.numeric(glm_PAvSvE_TPT4_bed[,2])),]

## Create genomic ranges object with new bed file
glm_PAvSvE_TPT4_bed <- GRanges(seqnames=glm_PAvSvE_TPT4.logp[,1],
		ranges=IRanges(as.integer(glm_PAvSvE_TPT4.logp[,2]), 	as.integer(glm_PAvSvE_TPT4.logp[,2]+1)),
		strand="+",
		pos=glm_PAvSvE_TPT4.logp[,2])

## Merge original loci with intervals so we can assign cluster colors for plotting
glm_PAvSvE_TPT4_EvPA.T4_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvE_TPT4_bed, glm_PAvSvE_TPT4_score_ttest_intervals_EvPA.T4))

## Reformat merged table to retain only relevant columns
glm_PAvSvE_TPT4_EvPA.T4_overlap <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvE_TPT4_EvPA.T4_overlap[,1]), 
	POS=as.integer(glm_PAvSvE_TPT4_EvPA.T4_overlap[,2]),
		paste(glm_PAvSvE_TPT4_EvPA.T4_overlap[,8],
		glm_PAvSvE_TPT4_EvPA.T4_overlap[,9],
		glm_PAvSvE_TPT4_EvPA.T4_overlap[,10],sep="_"),
		color=glm_PAvSvE_TPT4_EvPA.T4_overlap[,13]))

## Rename cluster column header
names(glm_PAvSvE_TPT4_EvPA.T4_overlap)[3] = "clust"

## Finally, merge -log10(p) values for making manhattans
glm_PAvSvE_TPT4_score_ttest_clusters_EvPA.T4 <- merge(glm_PAvSvE_TPT4.logp, glm_PAvSvE_TPT4_EvPA.T4_overlap, by=c("CHROM","POS"))

## Save table for future use
write.table(glm_PAvSvE_TPT4_score_ttest_clusters_EvPA.T4, file="glm_PAvSvE_TPT4_score_ttest_clusters_EvPA.T4.txt", sep = "\t", quote = FALSE, row.names = F)

## Reload table if desired
#glm_PAvSvE_TPT4_score_ttest_clusters_EvPA.T4 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/glm_PAvSvE_TPT4_score_ttest_clusters_EvPA.T4.txt", header=TRUE)

# This manhattan plot shows the contrast between E and PA samples in TPT 1
manh.EvPA.T4.clust <- ggplot(glm_PAvSvE_TPT4.logp, aes(POS, EvPA.T4.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Highlight clusters significant in PA vs E GLM contrast
	geom_point(data=glm_PAvSvE_TPT4_score_ttest_clusters_EvPA.T4,aes(POS, EvPA.T4.logp), color = glm_PAvSvE_TPT4_score_ttest_clusters_EvPA.T4$color, size = 0.1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvE_TPT4.logp$EvPA.T4.logp))+0.5)) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA (TPT4)")
	
	
## PA vs SP
## Set each cluster as an alternating color
glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4)))

## Convert to normal data frame
glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4.df <- as.data.frame(glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4)

## Use locus info to name each cluster
glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4.df$clust <- paste(glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4.df$seqnames,glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4.df$start,glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4.df$end,sep="_")

## Create a bed file using original SNP loci with "dummy" stop site
glm_PAvSvE_TPT4_bed <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvE_TPT4.logp$CHROM), POS=as.integer(glm_PAvSvE_TPT4.logp$POS), STOP=as.integer(glm_PAvSvE_TPT4.logp$POS + 1)))

## Sort bed file
glm_PAvSvE_TPT4_bed <- glm_PAvSvE_TPT4_bed[order(glm_PAvSvE_TPT4_bed[,1], as.numeric(glm_PAvSvE_TPT4_bed[,2])),]

## Create genomic ranges object with new bed file
glm_PAvSvE_TPT4_bed <- GRanges(seqnames=glm_PAvSvE_TPT4.logp[,1],
		ranges=IRanges(as.integer(glm_PAvSvE_TPT4.logp[,2]), 	as.integer(glm_PAvSvE_TPT4.logp[,2]+1)),
		strand="+",
		pos=glm_PAvSvE_TPT4.logp[,2])

## Merge original loci with intervals so we can assign cluster colors for plotting
glm_PAvSvE_TPT4_PAvSP.T4_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvE_TPT4_bed, glm_PAvSvE_TPT4_score_ttest_intervals_PAvSP.T4))

## Reformat merged table to retain only relevant columns
glm_PAvSvE_TPT4_PAvSP.T4_overlap <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvE_TPT4_PAvSP.T4_overlap[,1]), 
	POS=as.integer(glm_PAvSvE_TPT4_PAvSP.T4_overlap[,2]),
		paste(glm_PAvSvE_TPT4_PAvSP.T4_overlap[,8],
		glm_PAvSvE_TPT4_PAvSP.T4_overlap[,9],
		glm_PAvSvE_TPT4_PAvSP.T4_overlap[,10],sep="_"),
		color=glm_PAvSvE_TPT4_PAvSP.T4_overlap[,13]))

## Rename cluster column header
names(glm_PAvSvE_TPT4_PAvSP.T4_overlap)[3] = "clust"

## Finally, merge -log10(p) values for making manhattans
glm_PAvSvE_TPT4_score_ttest_clusters_PAvSP.T4 <- merge(glm_PAvSvE_TPT4.logp, glm_PAvSvE_TPT4_PAvSP.T4_overlap, by=c("CHROM","POS"))

## Save table for future use
write.table(glm_PAvSvE_TPT4_score_ttest_clusters_PAvSP.T4, file="glm_PAvSvE_TPT4_score_ttest_clusters_PAvSP.T4.txt", sep = "\t", quote = FALSE, row.names = F)

## Reload table if desired
#glm_PAvSvE_TPT4_score_ttest_clusters_PAvSP.T4 <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/glm_PAvSvE_TPT4_score_ttest_clusters_PAvSP.T4.txt", header=TRUE)

# This manhattan plot shows the contrast between PA and SP samples in TPT 1
manh.PAvSP.T4.clust <- ggplot(glm_PAvSvE_TPT4.logp, aes(POS, PAvSP.T4.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Highlight clusters significant in PA and SP GLM contrast
	geom_point(data=glm_PAvSvE_TPT4_score_ttest_clusters_PAvSP.T4,aes(POS, PAvSP.T4.logp), color = glm_PAvSvE_TPT4_score_ttest_clusters_PAvSP.T4$color, size = 0.1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm_PAvSvE_TPT4.logp$PAvSP.T4.logp))+0.5)) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("PA vs SP (TPT4)")


## Plot two most relevant pairwise TPT4 treatment comparisons together
pdf(file = "rudflies_2023_redo.EvSP.EvPA.T4.rolwin501clust.glm.manh.pdf", width=7.5, height=4)
	ggarrange(manh.EvSP.T4.clust, manh.EvPA.T4.clust,
              ncol = 1, nrow = 2)
dev.off()


##############################################################################
### Test for elevated allele frequency differences in top outliers in each ###
### contrast of interest (EvSE.T1, EvSP.T1, EvSPT4) within clusters from   ###
### each contrast of interest. There will be nine total comparisons (3X3). ###
###                                                                        ###
### There will be two versions of this reciprocal parallelism analysis.    ###
### Here, we will look for parallelism of the top outliers per cluster     ###
### based on the lowest all-sample GLM FDR values per respective contrast. ###
### Later, we will choose top outliers per leave-one-out GLMs results.     ###
##############################################################################

##############################################################
### "EvSE T1" AF difference per "EvSE T1" cluster top SNPs ###
##############################################################
f=4 #pick E vs SE field for analysis

## Make list of all SE cages in T1
SE_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SE",]$cage

for(i in SE_loo_cages) { #cycle through all SE T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SET1clust_EvSET1_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSE.T1 FDR value from full-sample GLM
		fdr_temp <- cbind(contrast.PAvSvSEvE.table[,c(1:2)], fdr=p.adjust(contrast.PAvSvSEvE.table[,4], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1[,c(1,2,9)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SET1clust_EvSET1_topsig_no",i,sep=""), rbind(get(paste("ttest_SET1clust_EvSET1_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SET1clust_EvSET1_topsig <- rbind(cbind("no11","SE",ttest_SET1clust_EvSET1_topsig_no11),cbind("no21","SE",ttest_SET1clust_EvSET1_topsig_no21),cbind("no27","SE",ttest_SET1clust_EvSET1_topsig_no27),cbind("no41","SE",ttest_SET1clust_EvSET1_topsig_no41),cbind("no45","SE",ttest_SET1clust_EvSET1_topsig_no45))

write.table(ttest_SET1clust_EvSET1_topsig, file="rudflies_2023_redo.ttest_SET1clust_EvSET1_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSE T1" AF difference per "EvSP T1" cluster top SNPs ###
##############################################################
f=4 #pick E vs SE field for analysis

## Make list of all SE cages in T1
SE_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SE",]$cage

for(i in SE_loo_cages) { #cycle through all SE T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SPT1clust_EvSET1_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSP.T1 FDR value from full-sample GLM
		fdr_temp <- cbind(contrast.PAvSvSEvE.table[,c(1:2)], fdr=p.adjust(contrast.PAvSvSEvE.table[,5], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1[,c(1,2,9)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SPT1clust_EvSET1_topsig_no",i,sep=""), rbind(get(paste("ttest_SPT1clust_EvSET1_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SPT1clust_EvSET1_topsig <- rbind(cbind("no11","SE",ttest_SPT1clust_EvSET1_topsig_no11),cbind("no21","SE",ttest_SPT1clust_EvSET1_topsig_no21),cbind("no27","SE",ttest_SPT1clust_EvSET1_topsig_no27),cbind("no41","SE",ttest_SPT1clust_EvSET1_topsig_no41),cbind("no45","SE",ttest_SPT1clust_EvSET1_topsig_no45))

write.table(ttest_SPT1clust_EvSET1_topsig, file="rudflies_2023_redo.ttest_SPT1clust_EvSET1_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSE T1" AF difference per "EvSP T4" cluster top SNPs ###
##############################################################
f=4 #pick E vs SE field for analysis

## Make list of all SE cages in T1
SE_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SE",]$cage

for(i in SE_loo_cages) { #cycle through all SE T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SPT4clust_EvSET1_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSP.T4 FDR value from full-sample GLM
		fdr_temp <- cbind(glm_PAvSvE_TPT4[,c(1:2)], fdr=p.adjust(glm_PAvSvE_TPT4[,3], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4[,c(1,2,21)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SPT4clust_EvSET1_topsig_no",i,sep=""), rbind(get(paste("ttest_SPT4clust_EvSET1_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SPT4clust_EvSET1_topsig <- rbind(cbind("no11","SE",ttest_SPT4clust_EvSET1_topsig_no11),cbind("no21","SE",ttest_SPT4clust_EvSET1_topsig_no21),cbind("no27","SE",ttest_SPT4clust_EvSET1_topsig_no27),cbind("no41","SE",ttest_SPT4clust_EvSET1_topsig_no41),cbind("no45","SE",ttest_SPT4clust_EvSET1_topsig_no45))

write.table(ttest_SPT4clust_EvSET1_topsig, file="rudflies_2023_redo.ttest_SPT4clust_EvSET1_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSP T1" AF difference per "EvSE T1" cluster top SNPs ###
##############################################################
f=5 #pick E vs SP field for analysis

## Make list of all SP cages in T1
SP_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SP",]$cage

for(i in SP_loo_cages) { #cycle through all SP T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SET1clust_EvSPT1_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSE.T1 FDR value from full-sample GLM
		fdr_temp <- cbind(contrast.PAvSvSEvE.table[,c(1:2)], fdr=p.adjust(contrast.PAvSvSEvE.table[,4], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1[,c(1,2,9)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SET1clust_EvSPT1_topsig_no",i,sep=""), rbind(get(paste("ttest_SET1clust_EvSPT1_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SET1clust_EvSPT1_topsig <- rbind(cbind("no3","SP",ttest_SET1clust_EvSPT1_topsig_no3),cbind("no7","SP",ttest_SET1clust_EvSPT1_topsig_no7),cbind("no15","SP",ttest_SET1clust_EvSPT1_topsig_no15),cbind("no33","SP",ttest_SET1clust_EvSPT1_topsig_no33),cbind("no37","SP",ttest_SET1clust_EvSPT1_topsig_no37))

write.table(ttest_SET1clust_EvSPT1_topsig, file="rudflies_2023_redo.ttest_SET1clust_EvSPT1_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSP T1" AF difference per "EvSP T1" cluster top SNPs ###
##############################################################
f=5 #pick E vs SP field for analysis

## Make list of all SP cages in T1
SP_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SP",]$cage

for(i in SP_loo_cages) { #cycle through all SP T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SPT1clust_EvSPT1_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSP.T1 FDR value from full-sample GLM
		fdr_temp <- cbind(contrast.PAvSvSEvE.table[,c(1:2)], fdr=p.adjust(contrast.PAvSvSEvE.table[,5], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1[,c(1,2,9)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SPT1clust_EvSPT1_topsig_no",i,sep=""), rbind(get(paste("ttest_SPT1clust_EvSPT1_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SPT1clust_EvSPT1_topsig <- rbind(cbind("no3","SP",ttest_SPT1clust_EvSPT1_topsig_no3),cbind("no7","SP",ttest_SPT1clust_EvSPT1_topsig_no7),cbind("no15","SP",ttest_SPT1clust_EvSPT1_topsig_no15),cbind("no33","SP",ttest_SPT1clust_EvSPT1_topsig_no33),cbind("no37","SP",ttest_SPT1clust_EvSPT1_topsig_no37))

write.table(ttest_SPT1clust_EvSPT1_topsig, file="rudflies_2023_redo.ttest_SPT1clust_EvSPT1_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSP T1" AF difference per "EvSP T4" cluster top SNPs ###
##############################################################
f=5 #pick E vs SP field for analysis

## Make list of all SP cages in T1
SP_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SP",]$cage

for(i in SP_loo_cages) { #cycle through all SP T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SPT4clust_EvSPT1_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSP.T4 FDR value from full-sample GLM
		fdr_temp <- cbind(glm_PAvSvE_TPT4[,c(1:2)], fdr=p.adjust(glm_PAvSvE_TPT4[,3], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4[,c(1,2,21)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SPT4clust_EvSPT1_topsig_no",i,sep=""), rbind(get(paste("ttest_SPT4clust_EvSPT1_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SPT4clust_EvSPT1_topsig <- rbind(cbind("no3","SP",ttest_SPT4clust_EvSPT1_topsig_no3),cbind("no7","SP",ttest_SPT4clust_EvSPT1_topsig_no7),cbind("no15","SP",ttest_SPT4clust_EvSPT1_topsig_no15),cbind("no33","SP",ttest_SPT4clust_EvSPT1_topsig_no33),cbind("no37","SP",ttest_SPT4clust_EvSPT1_topsig_no37))

write.table(ttest_SPT4clust_EvSPT1_topsig, file="rudflies_2023_redo.ttest_SPT4clust_EvSPT1_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSP T4" AF difference per "EvSE T1" cluster top SNPs ###
##############################################################
f=4 #pick E vs SP field for analysis

## Make list of all SP cages in T1
SP_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SP",]$cage

for(i in SP_loo_cages) { #cycle through all SP T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SET1clust_EvSPT4_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_t4_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_t4_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSE.T1 FDR value from full-sample GLM
		fdr_temp <- cbind(contrast.PAvSvSEvE.table[,c(1:2)], fdr=p.adjust(contrast.PAvSvSEvE.table[,4], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvSEvE_score_ttest_clusters_EvSE.T1[,c(1,2,9)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SET1clust_EvSPT4_topsig_no",i,sep=""), rbind(get(paste("ttest_SET1clust_EvSPT4_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_t4_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SET1clust_EvSPT4_topsig <- rbind(cbind("no3","SP",ttest_SET1clust_EvSPT4_topsig_no3),cbind("no7","SP",ttest_SET1clust_EvSPT4_topsig_no7),cbind("no15","SP",ttest_SET1clust_EvSPT4_topsig_no15),cbind("no33","SP",ttest_SET1clust_EvSPT4_topsig_no33),cbind("no37","SP",ttest_SET1clust_EvSPT4_topsig_no37))

write.table(ttest_SET1clust_EvSPT4_topsig, file="rudflies_2023_redo.ttest_SET1clust_EvSPT4_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSP T4" AF difference per "EvSP T1" cluster top SNPs ###
##############################################################
f=4 #pick E vs SP field for analysis

## Make list of all SP cages in T1
SP_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SP",]$cage

for(i in SP_loo_cages) { #cycle through all SP T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SPT1clust_EvSPT4_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_t4_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_t4_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSP.T1 FDR value from full-sample GLM
		fdr_temp <- cbind(contrast.PAvSvSEvE.table[,c(1:2)], fdr=p.adjust(contrast.PAvSvSEvE.table[,5], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvSEvE_score_ttest_clusters_EvSP.T1[,c(1,2,9)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SPT1clust_EvSPT4_topsig_no",i,sep=""), rbind(get(paste("ttest_SPT1clust_EvSPT4_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_t4_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SPT1clust_EvSPT4_topsig <- rbind(cbind("no3","SP",ttest_SPT1clust_EvSPT4_topsig_no3),cbind("no7","SP",ttest_SPT1clust_EvSPT4_topsig_no7),cbind("no15","SP",ttest_SPT1clust_EvSPT4_topsig_no15),cbind("no33","SP",ttest_SPT1clust_EvSPT4_topsig_no33),cbind("no37","SP",ttest_SPT1clust_EvSPT4_topsig_no37))

write.table(ttest_SPT1clust_EvSPT4_topsig, file="rudflies_2023_redo.ttest_SPT1clust_EvSPT4_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSP T4" AF difference per "EvSP T4" cluster top SNPs ###
##############################################################
f=4 #pick E vs SP field for analysis

## Make list of all SP cages in T1
SP_loo_cages <- haf.meta.T1filt[haf.meta.T1filt$condition=="SP",]$cage

for(i in SP_loo_cages) { #cycle through all SP T4 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_SPT4clust_EvSPT4_topsig_no",i,sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_t4_no",i,"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_t4_only",i,"_bed",sep=""))[,f])
		## Retrieve EvSP.T4 FDR value from full-sample GLM
		fdr_temp <- cbind(glm_PAvSvE_TPT4[,c(1:2)], fdr=p.adjust(glm_PAvSvE_TPT4[,3], method = "fdr"))
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, glm_PAvSvE_TPT4_score_ttest_clusters_EvSP.T4[,c(1,2,21)], by=c("CHROM","POS"))
		#vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
	
		## Select 1 random SNP per cluster (optional but not used)
		#vep_sig <- c()
		#for(c in unique(vep_temp2$clust)) { 
			#vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr < 0.05 & vep_temp2$clust == c,]),1))
		#}
		
		## Select 1k random SNPs across genome (optional but not used)
		#vep_sig <- sample_n(unique(vep_temp[vep_temp$fdr<0.05,]),min(nrow(unique(vep_temp[vep_temp$fdr<0.05,])),1000))
		
		## Select top 1K significant SNPs in order of FDR value (optional but not used)
		#vep_sig <- head(unique(vep_temp[order(vep_temp$fdr),]),1000)
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_SPT4clust_EvSPT4_topsig_no",i,sep=""), rbind(get(paste("ttest_SPT4clust_EvSPT4_topsig_no",i,sep="")),cbind(contrast=names(get(paste("freq_diff_t4_no",i,"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_SPT4clust_EvSPT4_topsig <- rbind(cbind("no3","SP",ttest_SPT4clust_EvSPT4_topsig_no3),cbind("no7","SP",ttest_SPT4clust_EvSPT4_topsig_no7),cbind("no15","SP",ttest_SPT4clust_EvSPT4_topsig_no15),cbind("no33","SP",ttest_SPT4clust_EvSPT4_topsig_no33),cbind("no37","SP",ttest_SPT4clust_EvSPT4_topsig_no37))

write.table(ttest_SPT4clust_EvSPT4_topsig, file="rudflies_2023_redo.ttest_SPT4clust_EvSPT4_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)




######################################################################################
### Leave-one-out parallelism analysis. Above, we looked at parallelism within     ###
### top sites within clusters so we'd have a uniform SNP set to compare different  ###
### treatments. Here, we will select top sites from unique clusters constructed    ###
### using each set of leave-one-out GLM results. This will provide more insight    ###
### into parallelism within treatments rather than across treatments. To do so,    ###
### we will cycle through each set of LOO GLM results, build clusters, and select  ###
### top sites for each SE.T1, SP.T1, and SP.T4 sample prior to running parallelism ###  
### tests.                                                                         ###
######################################################################################

### Start with EvSE.T1 and EvSP.T1 leave-one-out GLM results ###

###############################################
### Try to find clusters for data reduction ###
###############################################

### Create a table with leave-one-out p, fdr, and logp tables for EvSE.T1 & EvSP.T1
## Create lists of EvSE.T1 LOO header patterns
SE_loo_cages2 <- unlist(lapply(SE_loo_cages, function(x) paste0("EvSE.T1_no", x)))
SE_loo_cages2 <- paste(SE_loo_cages2, collapse = "|")

## Create lists of EvSP.T1 LOO header patterns
SP_loo_cages2 <- unlist(lapply(SP_loo_cages, function(x) paste0("EvSP.T1_no", x)))
SP_loo_cages2 <- paste(SP_loo_cages2, collapse = "|")

## Select columns that match any LOO cage/contrast in the list
glm_PAvSvSEvE_LOO <- glm_PAvSvSEvE_all %>% select(matches(SE_loo_cages2) | matches(SP_loo_cages2))

## Remove extra unnecessary matches
glm_PAvSvSEvE_LOO <- select(glm_PAvSvSEvE_LOO,-contains("SEvSP"))

## Add locus columns
glm_PAvSvSEvE_LOO <- cbind(glm_PAvSvSEvE_all[,c(1:2)],glm_PAvSvSEvE_LOO)

## Store duplicate table with corrected p-values for later
glm_PAvSvSEvE_LOO_fdr <- select(glm_PAvSvSEvE_LOO,contains("CHROM") | contains("POS") | contains("fdr"))

## Store duplicate table with -log10(P) values for later
glm_PAvSvSEvE_LOO_logp <- select(glm_PAvSvSEvE_LOO,contains("CHROM") | contains("POS") | contains("logp"))

## Remove unnecessary columns for clustering
glm_PAvSvSEvE_LOO <- select(glm_PAvSvSEvE_LOO,-contains("fdr"))
glm_PAvSvSEvE_LOO <- select(glm_PAvSvSEvE_LOO,-contains("logp"))


### Create a table with leave-one-out frequency difference tables for EvSE.T1 & EvSP.T1
## Create lists of SE LOO frequency difference tables
SE_loo_cages3 <- unlist(lapply(SE_loo_cages, function(x) paste0("freq_diff_no", x)))

## Create lists of SP LOO frequency difference tables
SP_loo_cages3 <- unlist(lapply(SP_loo_cages, function(x) paste0("freq_diff_no", x)))

## Grab locus columns to intialize master LOO frequency table
freq_diff_allLOO <- freq_diff_no3_bed[,c(1:2)]
## Cycle through all SE LOO tables and grab column 1 (EvSE.T1)
for(f in SE_loo_cages3) {
	freq_diff_allLOO <- cbind(freq_diff_allLOO,get(f)[,1] )
}
## Name columns
names(freq_diff_allLOO)[c(3:7)] <- SE_loo_cages3
names(freq_diff_allLOO) <- gsub("freq_diff","EvSE.T1",names(freq_diff_allLOO))

## Cycle through all SP LOO tables and grab column 2 (EvSP.T1)
for(f in SP_loo_cages3) {
	freq_diff_allLOO <- cbind(freq_diff_allLOO,get(f)[,2] )
}
## Name columns
names(freq_diff_allLOO)[c(8:12)] <- SP_loo_cages3
names(freq_diff_allLOO) <- gsub("freq_diff","EvSP.T1",names(freq_diff_allLOO))


## Join FDR and delta frequency cols that will be used to calculate scores for clustering
glm_PAvSvSEvE_LOO_score <- merge(glm_PAvSvSEvE_LOO_fdr, freq_diff_allLOO, by=c("CHROM","POS"))
## Sort by locus
glm_PAvSvSEvE_LOO_score <- glm_PAvSvSEvE_LOO_score[order(glm_PAvSvSEvE_LOO_score[,1], glm_PAvSvSEvE_LOO_score[,2]), ]

## Assign scores based on a combination of significance level and mean frequency diffs
## These criteria were borrowed and modified from Rudman et. al, 2022:
## "Direct observation of adaptive tracking on ecological time scales in Drosophila"
for(f in c(3:12)) {
	glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score[,f] > 0.2,f+20] <- 0
	glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score[,f] < 0.2 | (glm_PAvSvSEvE_LOO_score[,f] > 0.2 & glm_PAvSvSEvE_LOO_score[,f+10] > 0.02),f+20] <- 1
	glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score[,f] < 0.05 & glm_PAvSvSEvE_LOO_score[,f+10] > 0.02,f+20] <- 2
	glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score[,f] < 0.01 & glm_PAvSvSEvE_LOO_score[,f+10] > 0.02,f+20] <- 3
	glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score[,f] < 0.001 & glm_PAvSvSEvE_LOO_score[,f+10] > 0.02,f+20] <- 4
	glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score[,f] < 0.0001 & glm_PAvSvSEvE_LOO_score[,f+10] > 0.02,f+20] <- 5
	glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score[,f] < 0.00001 & glm_PAvSvSEvE_LOO_score[,f+10] > 0.02,f+20] <- 6
	names(glm_PAvSvSEvE_LOO_score)[f+20] <- paste(names(glm_PAvSvSEvE_LOO_score[f]),".score",sep="")
}

## Main idea is testing whether window-based scores are significantly higher than randomly 
## assigned scores. Here, let's shuffle scores per chromosome to control for observed 
## variable signals of selective sweeps across chromosomes.
set.seed(42) 
glm_PAvSvSEvE_LOO_score_rand <- c() #initialize data frame

## Loop through all chromosomes
for(i in unique(glm_PAvSvSEvE_LOO_score$CHROM)) {
	shuffle_idx <- c()
	tempCHROM <- glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score$CHROM==i,]
	shuffle_idx <- append(shuffle_idx, sample(1:nrow(tempCHROM)))
	tempCHROM[, c("CHROM", "POS")] <- tempCHROM[shuffle_idx, c("CHROM", "POS")]
	glm_PAvSvSEvE_LOO_score_rand <- rbind(glm_PAvSvSEvE_LOO_score_rand,tempCHROM)
}

## Sort based on randomized loci
glm_PAvSvSEvE_LOO_score_rand <- glm_PAvSvSEvE_LOO_score_rand[order(glm_PAvSvSEvE_LOO_score_rand[,1], glm_PAvSvSEvE_LOO_score_rand[,2]), ]


### True loci window scores ###
## Create indexes for each chromosome to slide over and calculate mean scores over 501-SNP 
## windows. If we don't do this per chromosome, the average will be calculated across 
## consecutive chromosomes.

## 2L
glm_2L_LOO_true <- unique(glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score$CHROM=="2L",]) 
glm_2L_LOO_true <- glm_2L_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2L_LOO_true <- glm_2L_LOO_true[is.finite(rowSums(glm_2L_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2L_LOO_true$POS <- as.integer(glm_2L_LOO_true$POS)
glm_2L_LOO_true.rolwin501 <- glm_2L_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_2L_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2L_LOO_true[is.finite(glm_2L_LOO_true[,i])=="TRUE",i])+100
	glm_2L_LOO_true[is.infinite(glm_2L_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2L_LOO_true.rolwin501 <- cbind(glm_2L_LOO_true.rolwin501,slide_mean(glm_2L_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_2L_LOO_true.rolwin501)[i] <- paste(names(glm_2L_LOO_true)[i], ".rolwin501", sep="")
}

## 2R
glm_2R_LOO_true <- unique(glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score$CHROM=="2R",])
glm_2R_LOO_true <- glm_2R_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2R_LOO_true <- glm_2R_LOO_true[is.finite(rowSums(glm_2R_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2R_LOO_true$POS <- as.integer(glm_2R_LOO_true$POS)
glm_2R_LOO_true.rolwin501 <- glm_2R_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_2R_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2R_LOO_true[is.finite(glm_2R_LOO_true[,i])=="TRUE",i])+100
	glm_2R_LOO_true[is.infinite(glm_2R_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2R_LOO_true.rolwin501 <- cbind(glm_2R_LOO_true.rolwin501,slide_mean(glm_2R_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_2R_LOO_true.rolwin501)[i] <- paste(names(glm_2R_LOO_true)[i], ".rolwin501", sep="")
}

## 3L
glm_3L_LOO_true <- unique(glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score$CHROM=="3L",])
glm_3L_LOO_true <- glm_3L_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3L_LOO_true <- glm_3L_LOO_true[is.finite(rowSums(glm_3L_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3L_LOO_true$POS <- as.integer(glm_3L_LOO_true$POS)
glm_3L_LOO_true.rolwin501 <- glm_3L_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_3L_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3L_LOO_true[is.finite(glm_3L_LOO_true[,i])=="TRUE",i])+100
	glm_3L_LOO_true[is.infinite(glm_3L_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3L_LOO_true.rolwin501 <- cbind(glm_3L_LOO_true.rolwin501,slide_mean(glm_3L_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_3L_LOO_true.rolwin501)[i] <- paste(names(glm_3L_LOO_true)[i], ".rolwin501", sep="")
}

## 3R
glm_3R_LOO_true <- unique(glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score$CHROM=="3R",])
glm_3R_LOO_true <- glm_3R_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3R_LOO_true <- glm_3R_LOO_true[is.finite(rowSums(glm_3R_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3R_LOO_true$POS <- as.integer(glm_3R_LOO_true$POS)
glm_3R_LOO_true.rolwin501 <- glm_3R_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_3R_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3R_LOO_true[is.finite(glm_3R_LOO_true[,i])=="TRUE",i])+100
	glm_3R_LOO_true[is.infinite(glm_3R_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3R_LOO_true.rolwin501 <- cbind(glm_3R_LOO_true.rolwin501,slide_mean(glm_3R_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_3R_LOO_true.rolwin501)[i] <- paste(names(glm_3R_LOO_true)[i], ".rolwin501", sep="")
}

## X
glm_X_LOO_true <- unique(glm_PAvSvSEvE_LOO_score[glm_PAvSvSEvE_LOO_score$CHROM=="X",])
glm_X_LOO_true <- glm_X_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_X_LOO_true <- glm_X_LOO_true[is.finite(rowSums(glm_X_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_X_LOO_true$POS <- as.integer(glm_X_LOO_true$POS) #reformat
glm_X_LOO_true.rolwin501 <- glm_X_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_X_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_X_LOO_true[is.finite(glm_X_LOO_true[,i])=="TRUE",i])+100
	glm_X_LOO_true[is.infinite(glm_X_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_X_LOO_true.rolwin501 <- cbind(glm_X_LOO_true.rolwin501,slide_mean(glm_X_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_X_LOO_true.rolwin501)[i] <- paste(names(glm_X_LOO_true)[i], ".rolwin501", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm_PAvSvSEvE_LOO_score.rolwin501 <- na.omit(rbind(glm_2L_LOO_true.rolwin501,glm_2R_LOO_true.rolwin501,glm_3L_LOO_true.rolwin501,glm_3R_LOO_true.rolwin501,glm_X_LOO_true.rolwin501))


## How many sites remain after merging all filtered GLM results tables?
dim(glm_PAvSvSEvE_LOO_score.rolwin501)
#[1] 15874    32


### Random locus window scores ###
## Create indexes for each chromosome to slide over and calculate mean scores over 501-SNP 
## windows. If we don't do this per chromosome, the average will be calculated across 
## consecutive chromosomes.

## 2L
glm_2L_LOO_rand <- unique(glm_PAvSvSEvE_LOO_score_rand[glm_PAvSvSEvE_LOO_score_rand$CHROM=="2L",]) 
glm_2L_LOO_rand <- glm_2L_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2L_LOO_rand <- glm_2L_LOO_rand[is.finite(rowSums(glm_2L_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2L_LOO_rand$POS <- as.integer(glm_2L_LOO_rand$POS)
glm_2L_LOO_rand.rolwin501 <- glm_2L_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_2L_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2L_LOO_rand[is.finite(glm_2L_LOO_rand[,i])=="TRUE",i])+100
	glm_2L_LOO_rand[is.infinite(glm_2L_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2L_LOO_rand.rolwin501 <- cbind(glm_2L_LOO_rand.rolwin501,slide_mean(glm_2L_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_2L_LOO_rand.rolwin501)[i] <- paste(names(glm_2L_LOO_rand)[i], ".rolwin501", sep="")
}

## 2R
glm_2R_LOO_rand <- unique(glm_PAvSvSEvE_LOO_score_rand[glm_PAvSvSEvE_LOO_score_rand$CHROM=="2R",])
glm_2R_LOO_rand <- glm_2R_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_2R_LOO_rand <- glm_2R_LOO_rand[is.finite(rowSums(glm_2R_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_2R_LOO_rand$POS <- as.integer(glm_2R_LOO_rand$POS)
glm_2R_LOO_rand.rolwin501 <- glm_2R_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_2R_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_2R_LOO_rand[is.finite(glm_2R_LOO_rand[,i])=="TRUE",i])+100
	glm_2R_LOO_rand[is.infinite(glm_2R_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_2R_LOO_rand.rolwin501 <- cbind(glm_2R_LOO_rand.rolwin501,slide_mean(glm_2R_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_2R_LOO_rand.rolwin501)[i] <- paste(names(glm_2R_LOO_rand)[i], ".rolwin501", sep="")
}

## 3L
glm_3L_LOO_rand <- unique(glm_PAvSvSEvE_LOO_score_rand[glm_PAvSvSEvE_LOO_score_rand$CHROM=="3L",])
glm_3L_LOO_rand <- glm_3L_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3L_LOO_rand <- glm_3L_LOO_rand[is.finite(rowSums(glm_3L_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3L_LOO_rand$POS <- as.integer(glm_3L_LOO_rand$POS)
glm_3L_LOO_rand.rolwin501 <- glm_3L_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_3L_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3L_LOO_rand[is.finite(glm_3L_LOO_rand[,i])=="TRUE",i])+100
	glm_3L_LOO_rand[is.infinite(glm_3L_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3L_LOO_rand.rolwin501 <- cbind(glm_3L_LOO_rand.rolwin501,slide_mean(glm_3L_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_3L_LOO_rand.rolwin501)[i] <- paste(names(glm_3L_LOO_rand)[i], ".rolwin501", sep="")
}

## 3R
glm_3R_LOO_rand <- unique(glm_PAvSvSEvE_LOO_score_rand[glm_PAvSvSEvE_LOO_score_rand$CHROM=="3R",])
glm_3R_LOO_rand <- glm_3R_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_3R_LOO_rand <- glm_3R_LOO_rand[is.finite(rowSums(glm_3R_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_3R_LOO_rand$POS <- as.integer(glm_3R_LOO_rand$POS)
glm_3R_LOO_rand.rolwin501 <- glm_3R_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_3R_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_3R_LOO_rand[is.finite(glm_3R_LOO_rand[,i])=="TRUE",i])+100
	glm_3R_LOO_rand[is.infinite(glm_3R_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_3R_LOO_rand.rolwin501 <- cbind(glm_3R_LOO_rand.rolwin501,slide_mean(glm_3R_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_3R_LOO_rand.rolwin501)[i] <- paste(names(glm_3R_LOO_rand)[i], ".rolwin501", sep="")
}

## X
glm_X_LOO_rand <- unique(glm_PAvSvSEvE_LOO_score_rand[glm_PAvSvSEvE_LOO_score_rand$CHROM=="X",])
glm_X_LOO_rand <- glm_X_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_X_LOO_rand <- glm_X_LOO_rand[is.finite(rowSums(glm_X_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_X_LOO_rand$POS <- as.integer(glm_X_LOO_rand$POS) #reformat
glm_X_LOO_rand.rolwin501 <- glm_X_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_X_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_X_LOO_rand[is.finite(glm_X_LOO_rand[,i])=="TRUE",i])+100
	glm_X_LOO_rand[is.infinite(glm_X_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_X_LOO_rand.rolwin501 <- cbind(glm_X_LOO_rand.rolwin501,slide_mean(glm_X_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_X_LOO_rand.rolwin501)[i] <- paste(names(glm_X_LOO_rand)[i], ".rolwin501", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm_PAvSvSEvE_LOO_score_rand.rolwin501 <- na.omit(rbind(glm_2L_LOO_rand.rolwin501,glm_2R_LOO_rand.rolwin501,glm_3L_LOO_rand.rolwin501,glm_3R_LOO_rand.rolwin501,glm_X_LOO_rand.rolwin501))

## How many sites remain after merging all filtered GLM results tables?
dim(glm_PAvSvSEvE_LOO_score_rand.rolwin501)
#[1] 15874    32


### T-tests within windows ###
## Now we want to perform one-sided t-tests on scores within true locus windows and random 
## locus windows. Significant windows with true locus scores higher than random locus 
## scores will later be merged.

## Check all contrasts
for(i in 23:32) { #loop through score columns, "SE" and "SP" contrasts with "E" pops
	## 2L
	glm_2L_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_2L_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_2L_LOO_ttest <- rbind(glm_2L_LOO_ttest,
    		cbind(glm_2L_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),2]),
    			STOP=max(glm_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),i], glm_2L_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),i]),
    			med_score=median(glm_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),i]),
    			max_score=max(glm_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),i]),
    			rand_score=mean(glm_2L_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),i]),
    			rand_med_score=median(glm_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),i]),
    			rand_max_score=max(glm_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2L_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_2L_LOO_ttest_",gsub(".score","",names(glm_2L_LOO_true)[i])), glm_2L_LOO_ttest)

	## 2R
	glm_2R_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_2R_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_2R_LOO_ttest <- rbind(glm_2R_LOO_ttest,
    		cbind(glm_2R_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),2]),
    			STOP=max(glm_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),i], glm_2R_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),i]),
    			med_score=median(glm_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),i]),
    			max_score=max(glm_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),i]),
    			rand_score=mean(glm_2R_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),i]),
    			rand_med_score=median(glm_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),i]),
    			rand_max_score=max(glm_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_2R_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_2R_LOO_ttest_",gsub(".score","",names(glm_2R_LOO_true)[i])), glm_2R_LOO_ttest)

	## 3L
	glm_3L_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_3L_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_3L_LOO_ttest <- rbind(glm_3L_LOO_ttest,
    		cbind(glm_3L_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),2]),
    			STOP=max(glm_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),i], glm_3L_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),i]),
    			med_score=median(glm_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),i]),
    			max_score=max(glm_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),i]),
    			rand_score=mean(glm_3L_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),i]),
    			rand_med_score=median(glm_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),i]),
    			rand_max_score=max(glm_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3L_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_3L_LOO_ttest_",gsub(".score","",names(glm_3L_LOO_true)[i])), glm_3L_LOO_ttest)

	## 3R
	glm_3R_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_3R_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_3R_LOO_ttest <- rbind(glm_3R_LOO_ttest,
    		cbind(glm_3R_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),2]),
    			STOP=max(glm_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),i], glm_3R_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),i]),
    			med_score=median(glm_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),i]),
    			max_score=max(glm_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),i]),
    			rand_score=mean(glm_3R_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),i]),
    			rand_med_score=median(glm_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),i]),
    			rand_max_score=max(glm_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_3R_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_3R_LOO_ttest_",gsub(".score","",names(glm_3R_LOO_true)[i])), glm_3R_LOO_ttest)

	## X
	glm_X_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_X_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_X_LOO_ttest <- rbind(glm_X_LOO_ttest,
    		cbind(glm_X_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),2]),
    			STOP=max(glm_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),i], glm_X_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),i]),
    			med_score=median(glm_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),i]),
    			max_score=max(glm_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),i]),
    			rand_score=mean(glm_X_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),i]),
    			rand_med_score=median(glm_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),i]),
    			rand_max_score=max(glm_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_X_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_X_LOO_ttest_",gsub(".score","",names(glm_X_LOO_true)[i])), glm_X_LOO_ttest)


	## Now, join stats from the the new rolling window tables for all chromosomes
    glm_PAvSvSEvE_LOO_score_ttest <- na.omit(rbind(glm_2L_LOO_ttest,glm_2R_LOO_ttest,glm_3L_LOO_ttest,glm_3R_LOO_ttest,glm_X_LOO_ttest))
    ## Reformat
    glm_PAvSvSEvE_LOO_score_ttest$statistic <- sapply(glm_PAvSvSEvE_LOO_score_ttest$statistic, toString)
    glm_PAvSvSEvE_LOO_score_ttest$parameter <- sapply(glm_PAvSvSEvE_LOO_score_ttest$parameter, toString)
    glm_PAvSvSEvE_LOO_score_ttest$p.value <- sapply(glm_PAvSvSEvE_LOO_score_ttest$p.value, toString)
    glm_PAvSvSEvE_LOO_score_ttest$stderr <- sapply(glm_PAvSvSEvE_LOO_score_ttest$stderr, toString)
    ## Correct t-test p-values
    glm_PAvSvSEvE_LOO_score_ttest$fdr <- p.adjust(glm_PAvSvSEvE_LOO_score_ttest$p.value, method = "fdr")    
    ## Name table for particular contrast 
    assign(paste0("glm_PAvSvSEvE_LOO_score_ttest_",
    	gsub(".score","",names(glm_X_LOO_true)[i])), glm_PAvSvSEvE_LOO_score_ttest)
    ## Save table for particular contrast 
    write.table(glm_PAvSvSEvE_LOO_score_ttest, file=paste0("glm_PAvSvSEvE_LOO_score_ttest_",
    	gsub(".score","",names(glm_X_LOO_true)[i]),".txt"),
    	sep = "\t", quote = FALSE, row.names = F)
    
    ## Filter windows for significant t-test results and minimum [max] score of 2.
    glm_PAvSvSEvE_LOO_score_ttest_filt <- glm_PAvSvSEvE_LOO_score_ttest[glm_PAvSvSEvE_LOO_score_ttest$fdr<0.05 & glm_PAvSvSEvE_LOO_score_ttest$max_score>=2,]
    
    ## We tried filtering by mean and median scores as well, but the EvSP had nothing.
    ## Using max scores within windows as filtering criteria worked for all contrasts.
    #glm_PAvSvSEvE_LOO_score_ttest_filt <- glm_PAvSvSEvE_LOO_score_ttest[glm_PAvSvSEvE_LOO_score_ttest$fdr<0.05 & glm_PAvSvSEvE_LOO_score_ttest$mean_score>=2,]
    
    ## Convert loci and window ranges to Genomic Ranges table
	gr <- GRanges(seqnames=glm_PAvSvSEvE_LOO_score_ttest_filt[,1],
		ranges=IRanges(as.integer(glm_PAvSvSEvE_LOO_score_ttest_filt[,3]), 	as.integer(glm_PAvSvSEvE_LOO_score_ttest_filt[,4])),
		strand="+",
		pos=glm_PAvSvSEvE_LOO_score_ttest_filt[,2])

	## Now merge all overlapping significant windows
	merged_gr_LOO <- reduce(gr)
    ## Name table for particular contrast 
	assign(paste0("glm_PAvSvSEvE_LOO_score_ttest_intervals_", 
		gsub(".fdr.score","",names(glm_X_LOO_true)[i])), merged_gr_LOO)
	## Save table for particular contrast 
	write.table(merged_gr_LOO, file=paste0("glm_PAvSvSEvE_LOO_score_ttest_intervals_",
		gsub(".fdr.score","",names(glm_X_LOO_true)[i]),".txt"),
		sep = "\t", quote = FALSE, row.names = F)
}

## How many merged clusters in each contrast?
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSE.T1_no11))
#[1] 138   5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSE.T1_no21))
#[1] 100   5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSE.T1_no27))
#[1] 133   5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSE.T1_no41))
#[1] 108   5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSE.T1_no45))
#[1] 120   5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSP.T1_no3))
#[1] 20  5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSP.T1_no7))
#[1] 35  5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSP.T1_no15))
#[1] 27  5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSP.T1_no33))
#[1] 60  5
dim(as.data.frame(glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSP.T1_no37))
#[1] 55  5


### Plot EvSE.T1 clusters per contrast ###

## Check all contrasts
for(i in 1:5) { #loop through SE cages
	## Grab data frame for left out cage
	tempdf <- get(paste("glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSE.T1_no",SE_loo_cages[i],sep="")) 

	### E vs SE contrast
	## Set each cluster as an alternating color
	tempdf$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(tempdf)))

	## Convert to normal data frame
	temp.df <- as.data.frame(tempdf)

	## Use locus info to name each cluster
	temp.df$clust <- paste(temp.df$seqnames,temp.df$start,temp.df$end,sep="_")

	## Create a bed file using original SNP loci with "dummy" stop site
	glm_PAvSvSEvE_LOO_bed <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvSEvE_LOO$CHROM), POS=as.integer(glm_PAvSvSEvE_LOO$POS), STOP=as.integer(glm_PAvSvSEvE_LOO$POS + 1)))

	## Sort bed file
	glm_PAvSvSEvE_LOO_bed <- glm_PAvSvSEvE_LOO_bed[order(glm_PAvSvSEvE_LOO_bed[,1], as.numeric(glm_PAvSvSEvE_LOO_bed[,2])),]

	## Create genomic ranges object with new bed file
	glm_PAvSvSEvE_LOO_bed <- GRanges(seqnames=glm_PAvSvSEvE_LOO_bed[,1],
		ranges=IRanges(as.integer(glm_PAvSvSEvE_LOO_bed[,2]), 	as.integer(glm_PAvSvSEvE_LOO_bed[,2])+1),
		strand="+",
		pos=glm_PAvSvSEvE_LOO_bed[,2])

	## Merge original loci with intervals so we can assign cluster colors for plotting
	tempdf_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvSEvE_LOO_bed, tempdf))

	## Reformat merged table to retain only relevant columns
	tempdf_overlap <- as.data.frame(cbind(CHROM=as.character(tempdf_overlap[,1]), 
		POS=as.integer(tempdf_overlap[,2]),
		paste(tempdf_overlap[,8],
		tempdf_overlap[,9],
		tempdf_overlap[,10],sep="_"),
		color=tempdf_overlap[,13]))

	## Rename cluster column header
	names(tempdf_overlap)[3] = "clust"

	## Finally, merge -log10(p) values for making manhattans
	tempdf_clusters <- merge(glm_PAvSvSEvE_LOO_logp[,c(1,2,i+2)], tempdf_overlap, by=c("CHROM","POS"))

	## Save table for future use
	assign(paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSE.T1_no", SE_loo_cages[i], sep=""), tempdf_clusters)
	## Write table
	write.table(tempdf_clusters, file=paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSE.T1_no",SE_loo_cages[i],".txt",sep=""), sep = "\t", quote = FALSE, row.names = F)

	## Reload table if desired
	#assign(paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSE.T1_no", SE_loo_cages[i], sep=""), read.table(paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSE.T1_no", SE_loo_cages[i], ".txt",sep=""), header=TRUE))

	## Grab logp values for current LOO cage
	tempdf_plot <- glm_PAvSvSEvE_LOO_logp[,c(1,2,i+2)]
	names(tempdf_plot)[3] <- "EvSE.T1.logp"
	names(tempdf_clusters)[3] <- "EvSE.T1.logp"

	# This manhattan plot shows the contrast between E and SE samples in TPT 1
	assign(paste("manh.EvSE.T1.clust_no",SE_loo_cages[i],sep=""),
	ggplot(tempdf_plot, aes(POS, EvSE.T1.logp)) +
		geom_line(alpha = 1, colour = "#CCCCCC") +
		## Highlight clusters significant in SE vs E GLM contrast
		geom_point(data=tempdf_clusters,aes(POS, EvSE.T1.logp), color = tempdf_clusters$color, size = 0.1) +
		## General formatting commands
		facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
		scale_y_continuous(limits = c(0, max(na.omit(tempdf_plot$EvSE.T1.logp)))) +
		scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
		labs(col="candidate\ngene\n-log10(p)") +
		xlab("chromosome position") +
		ylab("-log10(p)") +
		theme_classic() +
		theme(legend.position = "none") +
		ggtitle(paste("E vs SE (TPT1): no cage ",SE_loo_cages[i],sep="")))

}

## Plot the five EvSE.T1 LOO clustered manhattans together
pdf(file = "rudflies_2023_redo.EvSE.T1.rolwin501clust.LOO.glm.manh.pdf", width=7.5, height=10)
	ggarrange(manh.EvSE.T1.clust_no11, 
		manh.EvSE.T1.clust_no21, 
		manh.EvSE.T1.clust_no27, 
		manh.EvSE.T1.clust_no41, 
		manh.EvSE.T1.clust_no45,
        ncol = 1, nrow = 5)
dev.off()



### Plot EvSP.T1 clusters per contrast ###

## Check all contrasts
for(i in 1:5) { #loop through SP cages
	## Grab data frame for left out cage
	tempdf <- get(paste("glm_PAvSvSEvE_LOO_score_ttest_intervals_EvSP.T1_no",SP_loo_cages[i],sep="")) 

	### E vs SP contrast
	## Set each cluster as an alternating color
	tempdf$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(tempdf)))

	## Convert to normal data frame
	temp.df <- as.data.frame(tempdf)

	## Use locus info to name each cluster
	temp.df$clust <- paste(temp.df$seqnames,temp.df$start,temp.df$end,sep="_")

	## Create a bed file using original SNP loci with "dummy" stop site
	glm_PAvSvSEvE_LOO_bed <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvSEvE_LOO$CHROM), POS=as.integer(glm_PAvSvSEvE_LOO$POS), STOP=as.integer(glm_PAvSvSEvE_LOO$POS + 1)))

	## Sort bed file
	glm_PAvSvSEvE_LOO_bed <- glm_PAvSvSEvE_LOO_bed[order(glm_PAvSvSEvE_LOO_bed[,1], as.numeric(glm_PAvSvSEvE_LOO_bed[,2])),]

	## Create genomic ranges object with new bed file
	glm_PAvSvSEvE_LOO_bed <- GRanges(seqnames=glm_PAvSvSEvE_LOO_bed[,1],
		ranges=IRanges(as.integer(glm_PAvSvSEvE_LOO_bed[,2]), 	as.integer(glm_PAvSvSEvE_LOO_bed[,2])+1),
		strand="+",
		pos=glm_PAvSvSEvE_LOO_bed[,2])

	## Merge original loci with intervals so we can assign cluster colors for plotting
	tempdf_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvSEvE_LOO_bed, tempdf))

	## Reformat merged table to retain only relevant columns
	tempdf_overlap <- as.data.frame(cbind(CHROM=as.character(tempdf_overlap[,1]), 
		POS=as.integer(tempdf_overlap[,2]),
		paste(tempdf_overlap[,8],
		tempdf_overlap[,9],
		tempdf_overlap[,10],sep="_"),
		color=tempdf_overlap[,13]))

	## Rename cluster column header
	names(tempdf_overlap)[3] = "clust"

	## Finally, merge -log10(p) values for making manhattans
	tempdf_clusters <- merge(glm_PAvSvSEvE_LOO_logp[,c(1,2,i+7)], tempdf_overlap, by=c("CHROM","POS"))

	## Save table for future use
	assign(paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSP.T1_no", SP_loo_cages[i], sep=""), tempdf_clusters)
	## Write table
	write.table(tempdf_clusters, file=paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSP.T1_no",SP_loo_cages[i],".txt",sep=""), sep = "\t", quote = FALSE, row.names = F)

	## Reload table if desired
	#assign(paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSP.T1_no", SP_loo_cages[i], sep=""), read.table(paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSP.T1_no", SP_loo_cages[i], ".txt",sep=""), header=TRUE))

	## Grab logp values for current LOO cage
	tempdf_plot <- glm_PAvSvSEvE_LOO_logp[,c(1,2,i+7)]
	names(tempdf_plot)[3] <- "EvSP.T1.logp"
	names(tempdf_clusters)[3] <- "EvSP.T1.logp"

	# This manhattan plot shows the contrast between E and SP samples in TPT 1
	assign(paste("manh.EvSP.T1.clust_no",SP_loo_cages[i],sep=""),
	ggplot(tempdf_plot, aes(POS, EvSP.T1.logp)) +
		geom_line(alpha = 1, colour = "#CCCCCC") +
		## Highlight clusters significant in SP vs E GLM contrast
		geom_point(data=tempdf_clusters,aes(POS, EvSP.T1.logp), color = tempdf_clusters$color, size = 0.1) +
		## General formatting commands
		facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
		scale_y_continuous(limits = c(0, max(na.omit(tempdf_plot$EvSP.T1.logp)))) +
		scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
		labs(col="candidate\ngene\n-log10(p)") +
		xlab("chromosome position") +
		ylab("-log10(p)") +
		theme_classic() +
		theme(legend.position = "none") +
		ggtitle(paste("E vs SP (TPT1): no cage ",SP_loo_cages[i],sep="")))

}

## Plot the five LOO pairwise TPT1 treatment comparisons together
pdf(file = "rudflies_2023_redo.EvSP.T1.rolwin501clust.LOO.glm.manh.pdf", width=7.5, height=10)
	ggarrange(manh.EvSP.T1.clust_no3, 
		manh.EvSP.T1.clust_no7, 
		manh.EvSP.T1.clust_no15, 
		manh.EvSP.T1.clust_no33, 
		manh.EvSP.T1.clust_no37,
        ncol = 1, nrow = 5)
dev.off()



### Next let's analyze EvSP.T4 leave-one-out GLM results ###

###############################################
### Try to find clusters for data reduction ###
###############################################

### Create a table with leave-one-out p, fdr, and logp tables for EvSP.T4
## Create lists of EvSP.T4 LOO header patterns
SP_loo_cages4 <- unlist(lapply(SP_loo_cages, function(x) paste0("EvSP.T4_no", x)))
SP_loo_cages4 <- paste(SP_loo_cages4, collapse = "|")

## Select columns that match any LOO cage/contrast in the list
glm_PAvSvE_TPT4_LOO <- glm_PAvSvE_TPT4 %>% select(matches(SP_loo_cages4))

## Remove extra unnecessary matches
glm_PAvSvE_TPT4_LOO <- select(glm_PAvSvE_TPT4_LOO,-contains("SEvSP"))

## Add locus columns
glm_PAvSvE_TPT4_LOO <- cbind(glm_PAvSvE_TPT4[,c(1:2)],glm_PAvSvE_TPT4_LOO)

## Store duplicate table with corrected p-values for later
glm_PAvSvE_TPT4_LOO_fdr <- select(glm_PAvSvE_TPT4_LOO,contains("CHROM") | contains("POS") | contains("fdr"))

## Store duplicate table with -log10(P) values for later
glm_PAvSvE_TPT4_LOO_logp <- select(glm_PAvSvE_TPT4_LOO,contains("CHROM") | contains("POS") | contains("logp"))

## Remove unnecessary columns for clustering
glm_PAvSvE_TPT4_LOO <- select(glm_PAvSvE_TPT4_LOO,-contains("fdr"))
glm_PAvSvE_TPT4_LOO <- select(glm_PAvSvE_TPT4_LOO,-contains("logp"))


### Create a table with leave-one-out frequency difference tables for EvSP.T4
## Create lists of SP LOO frequency difference tables
SP_loo_cages5 <- unlist(lapply(SP_loo_cages, function(x) paste0("freq_diff_t4_no", x)))

## Grab locus columns to intialize master LOO frequency table
freq_diff_t4_allLOO <- freq_diff_t4_no3_bed[,c(1:2)]
## Cycle through all SP LOO tables and grab column 1 (EvSP.T4)
for(f in SP_loo_cages5) {
	freq_diff_t4_allLOO <- cbind(freq_diff_t4_allLOO,get(f)[,1] )
}
## Name columns
names(freq_diff_t4_allLOO)[c(3:7)] <- SP_loo_cages5
names(freq_diff_t4_allLOO) <- gsub("freq_diff_t4","EvSP.T4",names(freq_diff_t4_allLOO))


## Join FDR and delta frequency cols that will be used to calculate scores for clustering
glm_PAvSvE_TPT4_LOO_score <- merge(glm_PAvSvE_TPT4_LOO_fdr, freq_diff_t4_allLOO, by=c("CHROM","POS"))
## Sort by locus
glm_PAvSvE_TPT4_LOO_score <- glm_PAvSvE_TPT4_LOO_score[order(glm_PAvSvE_TPT4_LOO_score[,1], glm_PAvSvE_TPT4_LOO_score[,2]), ]

## Assign scores based on a combination of significance level and mean frequency diffs
## These criteria were borrowed and modified from Rudman et. al, 2022:
## "Direct observation of adaptive tracking on ecological time scales in Drosophila"
for(f in c(3:7)) {
	glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score[,f] > 0.2,f+10] <- 0
	glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score[,f] < 0.2 | (glm_PAvSvE_TPT4_LOO_score[,f] > 0.2 & glm_PAvSvE_TPT4_LOO_score[,f+5] > 0.02),f+10] <- 1
	glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score[,f] < 0.05 & glm_PAvSvE_TPT4_LOO_score[,f+5] > 0.02,f+10] <- 2
	glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score[,f] < 0.01 & glm_PAvSvE_TPT4_LOO_score[,f+5] > 0.02,f+10] <- 3
	glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score[,f] < 0.001 & glm_PAvSvE_TPT4_LOO_score[,f+5] > 0.02,f+10] <- 4
	glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score[,f] < 0.0001 & glm_PAvSvE_TPT4_LOO_score[,f+5] > 0.02,f+10] <- 5
	glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score[,f] < 0.00001 & glm_PAvSvE_TPT4_LOO_score[,f+5] > 0.02,f+10] <- 6
	names(glm_PAvSvE_TPT4_LOO_score)[f+10] <- paste(names(glm_PAvSvE_TPT4_LOO_score[f]),".score",sep="")
	names(glm_PAvSvE_TPT4_LOO_score)[f+10] <- gsub(".fdr",".score",names(glm_PAvSvE_TPT4_LOO_score)[f])
}

## Main idea is testing whether window-based scores are significantly higher than randomly 
## assigned scores. Here, let's shuffle scores per chromosome to control for observed 
## variable signals of selective sweeps across chromosomes.
set.seed(42) 
glm_PAvSvE_TPT4_LOO_score_rand <- c() #initialize data frame

## Loop through all chromosomes
for(i in unique(glm_PAvSvE_TPT4_LOO_score$CHROM)) {
	shuffle_idx <- c()
	tempCHROM <- glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score$CHROM==i,]
	shuffle_idx <- append(shuffle_idx, sample(1:nrow(tempCHROM)))
	tempCHROM[, c("CHROM", "POS")] <- tempCHROM[shuffle_idx, c("CHROM", "POS")]
	glm_PAvSvE_TPT4_LOO_score_rand <- rbind(glm_PAvSvE_TPT4_LOO_score_rand,tempCHROM)
}

## Sort based on randomized loci
glm_PAvSvE_TPT4_LOO_score_rand <- glm_PAvSvE_TPT4_LOO_score_rand[order(glm_PAvSvE_TPT4_LOO_score_rand[,1], glm_PAvSvE_TPT4_LOO_score_rand[,2]), ]


### True loci window scores ###
## Create indexes for each chromosome to slide over and calculate mean scores over 501-SNP 
## windows. If we don't do this per chromosome, the average will be calculated across 
## consecutive chromosomes.

## 2L
glm_T4_2L_LOO_true <- unique(glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score$CHROM=="2L",]) 
glm_T4_2L_LOO_true <- glm_T4_2L_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_2L_LOO_true <- glm_T4_2L_LOO_true[is.finite(rowSums(glm_T4_2L_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_2L_LOO_true$POS <- as.integer(glm_T4_2L_LOO_true$POS)
glm_T4_2L_LOO_true.rolwin501 <- glm_T4_2L_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_T4_2L_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_2L_LOO_true[is.finite(glm_T4_2L_LOO_true[,i])=="TRUE",i])+100
	glm_T4_2L_LOO_true[is.infinite(glm_T4_2L_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_2L_LOO_true.rolwin501 <- cbind(glm_T4_2L_LOO_true.rolwin501,slide_mean(glm_T4_2L_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_2L_LOO_true.rolwin501)[i] <- paste(names(glm_T4_2L_LOO_true)[i], ".rolwin501", sep="")
}

## 2R
glm_T4_2R_LOO_true <- unique(glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score$CHROM=="2R",])
glm_T4_2R_LOO_true <- glm_T4_2R_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_2R_LOO_true <- glm_T4_2R_LOO_true[is.finite(rowSums(glm_T4_2R_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_2R_LOO_true$POS <- as.integer(glm_T4_2R_LOO_true$POS)
glm_T4_2R_LOO_true.rolwin501 <- glm_T4_2R_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_T4_2R_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_2R_LOO_true[is.finite(glm_T4_2R_LOO_true[,i])=="TRUE",i])+100
	glm_T4_2R_LOO_true[is.infinite(glm_T4_2R_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_2R_LOO_true.rolwin501 <- cbind(glm_T4_2R_LOO_true.rolwin501,slide_mean(glm_T4_2R_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_2R_LOO_true.rolwin501)[i] <- paste(names(glm_T4_2R_LOO_true)[i], ".rolwin501", sep="")
}

## 3L
glm_T4_3L_LOO_true <- unique(glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score$CHROM=="3L",])
glm_T4_3L_LOO_true <- glm_T4_3L_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_3L_LOO_true <- glm_T4_3L_LOO_true[is.finite(rowSums(glm_T4_3L_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_3L_LOO_true$POS <- as.integer(glm_T4_3L_LOO_true$POS)
glm_T4_3L_LOO_true.rolwin501 <- glm_T4_3L_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_T4_3L_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_3L_LOO_true[is.finite(glm_T4_3L_LOO_true[,i])=="TRUE",i])+100
	glm_T4_3L_LOO_true[is.infinite(glm_T4_3L_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_3L_LOO_true.rolwin501 <- cbind(glm_T4_3L_LOO_true.rolwin501,slide_mean(glm_T4_3L_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_3L_LOO_true.rolwin501)[i] <- paste(names(glm_T4_3L_LOO_true)[i], ".rolwin501", sep="")
}

## 3R
glm_T4_3R_LOO_true <- unique(glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score$CHROM=="3R",])
glm_T4_3R_LOO_true <- glm_T4_3R_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_3R_LOO_true <- glm_T4_3R_LOO_true[is.finite(rowSums(glm_T4_3R_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_3R_LOO_true$POS <- as.integer(glm_T4_3R_LOO_true$POS)
glm_T4_3R_LOO_true.rolwin501 <- glm_T4_3R_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_T4_3R_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_3R_LOO_true[is.finite(glm_T4_3R_LOO_true[,i])=="TRUE",i])+100
	glm_T4_3R_LOO_true[is.infinite(glm_T4_3R_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_3R_LOO_true.rolwin501 <- cbind(glm_T4_3R_LOO_true.rolwin501,slide_mean(glm_T4_3R_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_3R_LOO_true.rolwin501)[i] <- paste(names(glm_T4_3R_LOO_true)[i], ".rolwin501", sep="")
}

## X
glm_T4_X_LOO_true <- unique(glm_PAvSvE_TPT4_LOO_score[glm_PAvSvE_TPT4_LOO_score$CHROM=="X",])
glm_T4_X_LOO_true <- glm_T4_X_LOO_true %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_X_LOO_true <- glm_T4_X_LOO_true[is.finite(rowSums(glm_T4_X_LOO_true[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_X_LOO_true$POS <- as.integer(glm_T4_X_LOO_true$POS) #reformat
glm_T4_X_LOO_true.rolwin501 <- glm_T4_X_LOO_true[,c(1,2)]

for(i in 3:ncol(glm_T4_X_LOO_true)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_X_LOO_true[is.finite(glm_T4_X_LOO_true[,i])=="TRUE",i])+100
	glm_T4_X_LOO_true[is.infinite(glm_T4_X_LOO_true[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_X_LOO_true.rolwin501 <- cbind(glm_T4_X_LOO_true.rolwin501,slide_mean(glm_T4_X_LOO_true[,i], before=250, after=250, step = 100))
    colnames(glm_T4_X_LOO_true.rolwin501)[i] <- paste(names(glm_T4_X_LOO_true)[i], ".rolwin501", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm_PAvSvE_TPT4_LOO_score.rolwin501 <- na.omit(rbind(glm_T4_2L_LOO_true.rolwin501,glm_T4_2R_LOO_true.rolwin501,glm_T4_3L_LOO_true.rolwin501,glm_T4_3R_LOO_true.rolwin501,glm_T4_X_LOO_true.rolwin501))


## How many sites remain after merging all filtered GLM results tables?
dim(glm_PAvSvE_TPT4_LOO_score.rolwin501)
#[1] 14253    17


### Random locus window scores ###
## Create indexes for each chromosome to slide over and calculate mean scores over 501-SNP 
## windows. If we don't do this per chromosome, the average will be calculated across 
## consecutive chromosomes.

## 2L
glm_T4_2L_LOO_rand <- unique(glm_PAvSvE_TPT4_LOO_score_rand[glm_PAvSvE_TPT4_LOO_score_rand$CHROM=="2L",]) 
glm_T4_2L_LOO_rand <- glm_T4_2L_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_2L_LOO_rand <- glm_T4_2L_LOO_rand[is.finite(rowSums(glm_T4_2L_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_2L_LOO_rand$POS <- as.integer(glm_T4_2L_LOO_rand$POS)
glm_T4_2L_LOO_rand.rolwin501 <- glm_T4_2L_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_2L_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_2L_LOO_rand[is.finite(glm_T4_2L_LOO_rand[,i])=="TRUE",i])+100
	glm_T4_2L_LOO_rand[is.infinite(glm_T4_2L_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_2L_LOO_rand.rolwin501 <- cbind(glm_T4_2L_LOO_rand.rolwin501,slide_mean(glm_T4_2L_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_2L_LOO_rand.rolwin501)[i] <- paste(names(glm_T4_2L_LOO_rand)[i], ".rolwin501", sep="")
}

## 2R
glm_T4_2R_LOO_rand <- unique(glm_PAvSvE_TPT4_LOO_score_rand[glm_PAvSvE_TPT4_LOO_score_rand$CHROM=="2R",])
glm_T4_2R_LOO_rand <- glm_T4_2R_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_2R_LOO_rand <- glm_T4_2R_LOO_rand[is.finite(rowSums(glm_T4_2R_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_2R_LOO_rand$POS <- as.integer(glm_T4_2R_LOO_rand$POS)
glm_T4_2R_LOO_rand.rolwin501 <- glm_T4_2R_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_2R_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_2R_LOO_rand[is.finite(glm_T4_2R_LOO_rand[,i])=="TRUE",i])+100
	glm_T4_2R_LOO_rand[is.infinite(glm_T4_2R_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_2R_LOO_rand.rolwin501 <- cbind(glm_T4_2R_LOO_rand.rolwin501,slide_mean(glm_T4_2R_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_2R_LOO_rand.rolwin501)[i] <- paste(names(glm_T4_2R_LOO_rand)[i], ".rolwin501", sep="")
}

## 3L
glm_T4_3L_LOO_rand <- unique(glm_PAvSvE_TPT4_LOO_score_rand[glm_PAvSvE_TPT4_LOO_score_rand$CHROM=="3L",])
glm_T4_3L_LOO_rand <- glm_T4_3L_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_3L_LOO_rand <- glm_T4_3L_LOO_rand[is.finite(rowSums(glm_T4_3L_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_3L_LOO_rand$POS <- as.integer(glm_T4_3L_LOO_rand$POS)
glm_T4_3L_LOO_rand.rolwin501 <- glm_T4_3L_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_3L_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_3L_LOO_rand[is.finite(glm_T4_3L_LOO_rand[,i])=="TRUE",i])+100
	glm_T4_3L_LOO_rand[is.infinite(glm_T4_3L_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_3L_LOO_rand.rolwin501 <- cbind(glm_T4_3L_LOO_rand.rolwin501,slide_mean(glm_T4_3L_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_3L_LOO_rand.rolwin501)[i] <- paste(names(glm_T4_3L_LOO_rand)[i], ".rolwin501", sep="")
}

## 3R
glm_T4_3R_LOO_rand <- unique(glm_PAvSvE_TPT4_LOO_score_rand[glm_PAvSvE_TPT4_LOO_score_rand$CHROM=="3R",])
glm_T4_3R_LOO_rand <- glm_T4_3R_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_3R_LOO_rand <- glm_T4_3R_LOO_rand[is.finite(rowSums(glm_T4_3R_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_3R_LOO_rand$POS <- as.integer(glm_T4_3R_LOO_rand$POS)
glm_T4_3R_LOO_rand.rolwin501 <- glm_T4_3R_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_3R_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_3R_LOO_rand[is.finite(glm_T4_3R_LOO_rand[,i])=="TRUE",i])+100
	glm_T4_3R_LOO_rand[is.infinite(glm_T4_3R_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_3R_LOO_rand.rolwin501 <- cbind(glm_T4_3R_LOO_rand.rolwin501,slide_mean(glm_T4_3R_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_3R_LOO_rand.rolwin501)[i] <- paste(names(glm_T4_3R_LOO_rand)[i], ".rolwin501", sep="")
}

## X
glm_T4_X_LOO_rand <- unique(glm_PAvSvE_TPT4_LOO_score_rand[glm_PAvSvE_TPT4_LOO_score_rand$CHROM=="X",])
glm_T4_X_LOO_rand <- glm_T4_X_LOO_rand %>%
  arrange(POS)
## remove rows with infinite values that will mess up calculations
#glm_T4_X_LOO_rand <- glm_T4_X_LOO_rand[is.finite(rowSums(glm_T4_X_LOO_rand[,-c(1:2)])),]
## create new empty file for filtered rows
glm_T4_X_LOO_rand$POS <- as.integer(glm_T4_X_LOO_rand$POS) #reformat
glm_T4_X_LOO_rand.rolwin501 <- glm_T4_X_LOO_rand[,c(1,2)]

for(i in 3:ncol(glm_T4_X_LOO_rand)) { #loop through all non-positional columns
	# infinite values not allowed, so create a ceiling of the max finite value plus 100
	# then assign to infinite values.
	max.temp <- max(glm_T4_X_LOO_rand[is.finite(glm_T4_X_LOO_rand[,i])=="TRUE",i])+100
	glm_T4_X_LOO_rand[is.infinite(glm_T4_X_LOO_rand[,i])=="TRUE",i] <- max.temp
	# now calculate rolling means
	glm_T4_X_LOO_rand.rolwin501 <- cbind(glm_T4_X_LOO_rand.rolwin501,slide_mean(glm_T4_X_LOO_rand[,i], before=250, after=250, step = 100))
    colnames(glm_T4_X_LOO_rand.rolwin501)[i] <- paste(names(glm_T4_X_LOO_rand)[i], ".rolwin501", sep="")
}

## Now, join the the new rolling window tables for all chromosomes
glm_PAvSvE_TPT4_LOO_score_rand.rolwin501 <- na.omit(rbind(glm_T4_2L_LOO_rand.rolwin501,glm_T4_2R_LOO_rand.rolwin501,glm_T4_3L_LOO_rand.rolwin501,glm_T4_3R_LOO_rand.rolwin501,glm_T4_X_LOO_rand.rolwin501))


## How many sites remain after merging all filtered GLM results tables?
dim(glm_PAvSvE_TPT4_LOO_score_rand.rolwin501)
#[1] 14253    17


### T-tests within windows ###
## Now we want to perform one-sided t-tests on scores within true locus windows and random 
## locus windows. Significant windows with true locus scores higher than random locus 
## scores will later be merged.

## Check all contrasts
for(i in 13:17) { #loop through score columns, "SP" contrasts with "E" pops
	## 2L
	glm_T4_2L_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_2L_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_2L_LOO_ttest <- rbind(glm_T4_2L_LOO_ttest,
    		cbind(glm_T4_2L_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),2]),
    			STOP=max(glm_T4_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_T4_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),i], glm_T4_2L_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),i]),
    			med_score=median(glm_T4_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),i]),
    			max_score=max(glm_T4_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),i]),
    			rand_score=mean(glm_T4_2L_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),i]),
    			rand_med_score=median(glm_T4_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),i]),
    			rand_max_score=max(glm_T4_2L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2L_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_2L_LOO_ttest_",gsub(".score","",names(glm_T4_2L_LOO_true)[i])), glm_T4_2L_LOO_ttest)

	## 2R
	glm_T4_2R_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_2R_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_2R_LOO_ttest <- rbind(glm_T4_2R_LOO_ttest,
    		cbind(glm_T4_2R_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),2]),
    			STOP=max(glm_T4_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_T4_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),i], glm_T4_2R_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),i]),
    			med_score=median(glm_T4_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),i]),
    			max_score=max(glm_T4_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),i]),
    			rand_score=mean(glm_T4_2R_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),i]),
    			rand_med_score=median(glm_T4_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),i]),
    			rand_max_score=max(glm_T4_2R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_2R_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_2R_LOO_ttest_",gsub(".score","",names(glm_T4_2R_LOO_true)[i])), glm_T4_2R_LOO_ttest)

	## 3L
	glm_T4_3L_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_3L_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_3L_LOO_ttest <- rbind(glm_T4_3L_LOO_ttest,
    		cbind(glm_T4_3L_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),2]),
    			STOP=max(glm_T4_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_T4_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),i], glm_T4_3L_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),i]),
    			med_score=median(glm_T4_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),i]),
    			max_score=max(glm_T4_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),i]),
    			rand_score=mean(glm_T4_3L_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),i]),
    			rand_med_score=median(glm_T4_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),i]),
    			rand_max_score=max(glm_T4_3L_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3L_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_3L_LOO_ttest_",gsub(".score","",names(glm_T4_3L_LOO_true)[i])), glm_T4_3L_LOO_ttest)

	## 3R
	glm_T4_3R_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_3R_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_3R_LOO_ttest <- rbind(glm_T4_3R_LOO_ttest,
    		cbind(glm_T4_3R_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),2]),
    			STOP=max(glm_T4_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_T4_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),i], glm_T4_3R_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),i]),
    			med_score=median(glm_T4_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),i]),
    			max_score=max(glm_T4_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),i]),
    			rand_score=mean(glm_T4_3R_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),i]),
    			rand_med_score=median(glm_T4_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),i]),
    			rand_max_score=max(glm_T4_3R_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_3R_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_3R_LOO_ttest_",gsub(".score","",names(glm_T4_3R_LOO_true)[i])), glm_T4_3R_LOO_ttest)

	## X
	glm_T4_X_LOO_ttest <- c()

    for(j in seq(from = 1, to = nrow(glm_T4_X_LOO_true), by = 100)) { #use 100-SNP step size
    	tryCatch({
    	glm_T4_X_LOO_ttest <- rbind(glm_T4_X_LOO_ttest,
    		cbind(glm_T4_X_LOO_true[j,c(1:2)],
    			#Set range
    			START=min(glm_T4_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),2]),
    			STOP=max(glm_T4_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),2]),
    			#Run test
    			t(t.test(glm_T4_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),i], glm_T4_X_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),i], alternative = "greater")[c(1,2,3,7)]),
    			#Save various score summary stats for filtering
    			mean_score=mean(glm_T4_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),i]),
    			med_score=median(glm_T4_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),i]),
    			max_score=max(glm_T4_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),i]),
    			rand_score=mean(glm_T4_X_LOO_rand[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),i]),
    			rand_med_score=median(glm_T4_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),i]),
    			rand_max_score=max(glm_T4_X_LOO_true[c(max(j-250,1):min(j+250,nrow(glm_T4_X_LOO_true))),i])))
    	}, error=function(e){})
    }
    assign(paste0("glm_T4_X_LOO_ttest_",gsub(".score","",names(glm_T4_X_LOO_true)[i])), glm_T4_X_LOO_ttest)


	## Now, join stats from the the new rolling window tables for all chromosomes
    glm_PAvSvE_TPT4_LOO_score_ttest <- na.omit(rbind(glm_T4_2L_LOO_ttest,glm_T4_2R_LOO_ttest,glm_T4_3L_LOO_ttest,glm_T4_3R_LOO_ttest,glm_T4_X_LOO_ttest))
    ## Reformat
    glm_PAvSvE_TPT4_LOO_score_ttest$statistic <- sapply(glm_PAvSvE_TPT4_LOO_score_ttest$statistic, toString)
    glm_PAvSvE_TPT4_LOO_score_ttest$parameter <- sapply(glm_PAvSvE_TPT4_LOO_score_ttest$parameter, toString)
    glm_PAvSvE_TPT4_LOO_score_ttest$p.value <- sapply(glm_PAvSvE_TPT4_LOO_score_ttest$p.value, toString)
    glm_PAvSvE_TPT4_LOO_score_ttest$stderr <- sapply(glm_PAvSvE_TPT4_LOO_score_ttest$stderr, toString)
    ## Correct t-test p-values
    glm_PAvSvE_TPT4_LOO_score_ttest$fdr <- p.adjust(glm_PAvSvE_TPT4_LOO_score_ttest$p.value, method = "fdr")    
    ## Name table for particular contrast 
    assign(paste0("glm_PAvSvE_TPT4_LOO_score_ttest_",
    	gsub(".score","",names(glm_T4_X_LOO_true)[i])), glm_PAvSvE_TPT4_LOO_score_ttest)
    ## Save table for particular contrast 
    write.table(glm_PAvSvE_TPT4_LOO_score_ttest, file=paste0("glm_PAvSvE_TPT4_LOO_score_ttest_",
    	gsub(".score","",names(glm_T4_X_LOO_true)[i]),".txt"),
    	sep = "\t", quote = FALSE, row.names = F)
    
    ## Filter windows for significant t-test results and minimum [max] score of 2.
    glm_PAvSvE_TPT4_LOO_score_ttest_filt <- glm_PAvSvE_TPT4_LOO_score_ttest[glm_PAvSvE_TPT4_LOO_score_ttest$fdr<0.05 & glm_PAvSvE_TPT4_LOO_score_ttest$max_score>=2,]
    
    ## We tried filtering by mean and median scores as well, but the EvSP had nothing.
    ## Using max scores within windows as filtering criteria worked for all contrasts.
    #glm_PAvSvE_TPT4_LOO_score_ttest_filt <- glm_PAvSvE_TPT4_LOO_score_ttest[glm_PAvSvE_TPT4_LOO_score_ttest$fdr<0.05 & glm_PAvSvE_TPT4_LOO_score_ttest$mean_score>=2,]
    
    ## Convert loci and window ranges to Genomic Ranges table
	gr <- GRanges(seqnames=glm_PAvSvE_TPT4_LOO_score_ttest_filt[,1],
		ranges=IRanges(as.integer(glm_PAvSvE_TPT4_LOO_score_ttest_filt[,3]), 	as.integer(glm_PAvSvE_TPT4_LOO_score_ttest_filt[,4])),
		strand="+",
		pos=glm_PAvSvE_TPT4_LOO_score_ttest_filt[,2])

	## Now merge all overlapping significant windows
	merged_gr_T4_LOO <- reduce(gr)
    ## Name table for particular contrast 
	assign(paste0("glm_PAvSvE_TPT4_LOO_score_ttest_intervals_", 
		gsub(".score","",names(glm_T4_X_LOO_true)[i])), merged_gr_T4_LOO)
	## Save table for particular contrast 
	write.table(merged_gr_T4_LOO, file=paste0("glm_PAvSvE_TPT4_LOO_score_ttest_intervals_",
		gsub(".score","",names(glm_T4_X_LOO_true)[i]),".txt"),
		sep = "\t", quote = FALSE, row.names = F)
}

## How many merged clusters in each contrast?
dim(as.data.frame(glm_PAvSvE_TPT4_LOO_score_ttest_intervals_EvSP.T4_no3))
#[1] 360  5
dim(as.data.frame(glm_PAvSvE_TPT4_LOO_score_ttest_intervals_EvSP.T4_no7))
#[1] 374  5
dim(as.data.frame(glm_PAvSvE_TPT4_LOO_score_ttest_intervals_EvSP.T4_no15))
#[1] 325  5
dim(as.data.frame(glm_PAvSvE_TPT4_LOO_score_ttest_intervals_EvSP.T4_no33))
#[1] 363  5
dim(as.data.frame(glm_PAvSvE_TPT4_LOO_score_ttest_intervals_EvSP.T4_no37))
#[1] 283  5


### Plot clusters per EvSP.T4 contrast ###

## Check all contrasts
for(i in 1:5) { #loop through SP cages
	## Grab data frame for left out cage
	tempdf <- get(paste("glm_PAvSvE_TPT4_LOO_score_ttest_intervals_EvSP.T4_no",SP_loo_cages[i],sep="")) 

	### E vs SP contrast
	## Set each cluster as an alternating color
	tempdf$color <- rep(c("blue4", "orange3"),length.out = nrow(as.data.frame(tempdf)))

	## Convert to normal data frame
	temp.df <- as.data.frame(tempdf)

	## Use locus info to name each cluster
	temp.df$clust <- paste(temp.df$seqnames,temp.df$start,temp.df$end,sep="_")

	## Create a bed file using original SNP loci with "dummy" stop site
	glm_PAvSvE_TPT4_LOO_bed <- as.data.frame(cbind(CHROM=as.character(glm_PAvSvE_TPT4_LOO$CHROM), POS=as.integer(glm_PAvSvE_TPT4_LOO$POS), STOP=as.integer(glm_PAvSvE_TPT4_LOO$POS + 1)))

	## Sort bed file
	glm_PAvSvE_TPT4_LOO_bed <- glm_PAvSvE_TPT4_LOO_bed[order(glm_PAvSvE_TPT4_LOO_bed[,1], as.numeric(glm_PAvSvE_TPT4_LOO_bed[,2])),]

	## Create genomic ranges object with new bed file
	glm_PAvSvE_TPT4_LOO_bed <- GRanges(seqnames=glm_PAvSvE_TPT4_LOO_bed[,1],
		ranges=IRanges(as.integer(glm_PAvSvE_TPT4_LOO_bed[,2]), 	as.integer(glm_PAvSvE_TPT4_LOO_bed[,2])+1),
		strand="+",
		pos=glm_PAvSvE_TPT4_LOO_bed[,2])

	## Merge original loci with intervals so we can assign cluster colors for plotting
	tempdf_overlap <- as.data.frame(mergeByOverlaps(glm_PAvSvE_TPT4_LOO_bed, tempdf))

	## Reformat merged table to retain only relevant columns
	tempdf_overlap <- as.data.frame(cbind(CHROM=as.character(tempdf_overlap[,1]), 
		POS=as.integer(tempdf_overlap[,2]),
		paste(tempdf_overlap[,8],
		tempdf_overlap[,9],
		tempdf_overlap[,10],sep="_"),
		color=tempdf_overlap[,13]))

	## Rename cluster column header
	names(tempdf_overlap)[3] = "clust"

	## Finally, merge -log10(p) values for making manhattans
	tempdf_clusters <- merge(glm_PAvSvE_TPT4_LOO_logp[,c(1,2,i+2)], tempdf_overlap, by=c("CHROM","POS"))

	## Save table for future use
	assign(paste("glm_PAvSvE_TPT4_LOO_score_ttest_clusters_EvSP.T4_no", SP_loo_cages[i], sep=""), tempdf_clusters)
	## Write table
	write.table(tempdf_clusters, file=paste("glm_PAvSvE_TPT4_LOO_score_ttest_clusters_EvSP.T4_no",SP_loo_cages[i],".txt",sep=""), sep = "\t", quote = FALSE, row.names = F)

	## Reload table if desired
	#assign(paste("glm_PAvSvE_TPT4_LOO_score_ttest_clusters_EvSP.T4_no", SP_loo_cages[i], sep=""), read.table(paste("glm_PAvSvE_TPT4_LOO_score_ttest_clusters_EvSP.T4_no", SP_loo_cages[i], ".txt",sep=""), header=TRUE))

	## Grab logp values for current LOO cage
	tempdf_plot <- glm_PAvSvE_TPT4_LOO_logp[,c(1,2,i+2)]
	names(tempdf_plot)[3] <- "EvSP.T4.logp"
	names(tempdf_clusters)[3] <- "EvSP.T4.logp"

	# This manhattan plot shows the contrast between E and SP samples in TPT 4
	assign(paste("manh.EvSP.T4.clust_no",SP_loo_cages[i],sep=""),
	ggplot(tempdf_plot, aes(POS, EvSP.T4.logp)) +
		geom_line(alpha = 1, colour = "#CCCCCC") +
		## Highlight clusters significant in SP vs E GLM contrast
		geom_point(data=tempdf_clusters,aes(POS, EvSP.T4.logp), color = tempdf_clusters$color, size = 0.1) +
		## General formatting commands
		facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
		scale_y_continuous(limits = c(0, max(na.omit(tempdf_plot[!grepl("Inf", tempdf_plot$EvSP.T4.logp), ]$EvSP.T4.logp)))) +
		scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
		labs(col="candidate\ngene\n-log10(p)") +
		xlab("chromosome position") +
		ylab("-log10(p)") +
		theme_classic() +
		theme(legend.position = "none") +
		ggtitle(paste("E vs SP (TPT4): no cage ",SP_loo_cages[i],sep="")))

}

## Plot the five LOO pairwise TPT4 treatment comparisons together
pdf(file = "rudflies_2023_redo.EvSP.T4.rolwin501clust.LOO.glm.manh.pdf", width=7.5, height=10)
	ggarrange(manh.EvSP.T4.clust_no3, 
		manh.EvSP.T4.clust_no7, 
		manh.EvSP.T4.clust_no15, 
		manh.EvSP.T4.clust_no33, 
		manh.EvSP.T4.clust_no37,
        ncol = 1, nrow = 5)
dev.off()




##############################################################################
### Test for elevated allele frequency differences of top outliers in each ###
### leave-one-out cluster set of interest: five samples from each of three ###
### contrasts (EvSE.T1, EvSP.T1, EvSP.T4). There will be nine total        ###
### comparisons (3 X 5).                                                   ###
##############################################################################

##############################################################
### "EvSE T1" AF difference per "EvSE T1" cluster top SNPs ###
##############################################################
f=4 #pick E vs SE field for analysis

for(i in 1:5) { #cycle through all SE T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_LOO_SET1clust_EvSET1_topsig_no",SE_loo_cages[i],sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_no",SE_loo_cages[i],"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_only",SE_loo_cages[i],"_bed",sep=""))[,f])
		## Retrieve EvSE.T1 FDR value from LOO GLM
		fdr_temp <- cbind(glm_PAvSvSEvE_LOO_fdr[,c(1:2)], fdr=glm_PAvSvSEvE_LOO_fdr[,i+2])
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, get(paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSE.T1_no", SE_loo_cages[i], sep=""))[,c(1,2,4)], by=c("CHROM","POS"))
		## Filter for significant sites
		vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]

		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_LOO_SET1clust_EvSET1_topsig_no",SE_loo_cages[i],sep=""), rbind(get(paste("ttest_LOO_SET1clust_EvSET1_topsig_no",SE_loo_cages[i],sep="")),cbind(contrast=names(get(paste("freq_diff_no",SE_loo_cages[i],"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_LOO_SET1clust_EvSET1_topsig <- rbind(cbind("no11","SE",ttest_LOO_SET1clust_EvSET1_topsig_no11),
	cbind("no21","SE",ttest_LOO_SET1clust_EvSET1_topsig_no21),
	cbind("no27","SE",ttest_LOO_SET1clust_EvSET1_topsig_no27),
	cbind("no41","SE",ttest_LOO_SET1clust_EvSET1_topsig_no41),
	cbind("no45","SE",ttest_LOO_SET1clust_EvSET1_topsig_no45))

write.table(ttest_LOO_SET1clust_EvSET1_topsig, file="rudflies_2023_redo.ttest_LOO_SET1clust_EvSET1_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)



##############################################################
### "EvSP T1" AF difference per "EvSP T1" cluster top SNPs ###
##############################################################
f=5 #pick E vs SP field for analysis

for(i in 1:5) { #cycle through all SP T1 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_LOO_SPT1clust_EvSPT1_topsig_no",SP_loo_cages[i],sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_no",SP_loo_cages[i],"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_only",SP_loo_cages[i],"_bed",sep=""))[,f])
		## Retrieve EvSP.T1 FDR value from LOO GLM
		fdr_temp <- cbind(glm_PAvSvSEvE_LOO_fdr[,c(1:2)], fdr=glm_PAvSvSEvE_LOO_fdr[,i+7])
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, get(paste("glm_PAvSvSEvE_LOO_score_ttest_clusters_EvSP.T1_no", SP_loo_cages[i], sep=""))[,c(1,2,4)], by=c("CHROM","POS"))
		## Filter for significant sites
		vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]

		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_LOO_SPT1clust_EvSPT1_topsig_no",SP_loo_cages[i],sep=""), rbind(get(paste("ttest_LOO_SPT1clust_EvSPT1_topsig_no",SP_loo_cages[i],sep="")),cbind(contrast=names(get(paste("freq_diff_no",SP_loo_cages[i],"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_LOO_SPT1clust_EvSPT1_topsig <- rbind(cbind("no3","SP",ttest_LOO_SPT1clust_EvSPT1_topsig_no3),
	cbind("no7","SP",ttest_LOO_SPT1clust_EvSPT1_topsig_no7),
	cbind("no15","SP",ttest_LOO_SPT1clust_EvSPT1_topsig_no15),
	cbind("no33","SP",ttest_LOO_SPT1clust_EvSPT1_topsig_no33),
	cbind("no37","SP",ttest_LOO_SPT1clust_EvSPT1_topsig_no37))

write.table(ttest_LOO_SPT1clust_EvSPT1_topsig, file="rudflies_2023_redo.ttest_LOO_SPT1clust_EvSPT1_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)


##############################################################
### "EvSP T4" AF difference per "EvSP T4" cluster top SNPs ###
##############################################################
f=4 #pick E vs SP field for analysis

for(i in 1:5) { #cycle through all SP T4 cages
	### T-tests for all contrasts each leave one out GLM ###
	assign(paste("ttest_LOO_SPT4clust_EvSPT4_topsig_no",SP_loo_cages[i],sep=""),c())
	## This run 100 iterations of focal and matched SNP selection
	for(z in 1:100) { #set number of iterations
		## Build dataframe with delta AF from left-out and left in samples for comparison
		sign_temp <- cbind(get(paste("freq_diff_t4_no",SP_loo_cages[i],"_bed",sep=""))[,c(1,2,f)],AFdiff_LOO=get(paste("freq_diff_t4_only",SP_loo_cages[i],"_bed",sep=""))[,f])
		## Retrieve EvSP.T4 FDR value from full-sample GLM
		fdr_temp <- cbind(glm_PAvSvE_TPT4_LOO[,c(1:2)], fdr=glm_PAvSvE_TPT4_LOO[,i+2])
		## Join them together
		sign_temp <- merge(sign_temp, fdr_temp, by=c("CHROM","POS"))
		## Add VEP info for finding matches
		vep_temp <- merge(sign_temp, vep_priority, by=c("CHROM","POS"))
		## Add cluster info
		vep_temp2 <- merge(vep_temp, get(paste("glm_PAvSvE_TPT4_LOO_score_ttest_clusters_EvSP.T4_no", SP_loo_cages[i], sep=""))[,c(1,2,4)], by=c("CHROM","POS"))
		## Filter for significant sites
		vep_temp2  <- vep_temp2[vep_temp2$fdr<0.05,]
	
		## Select top significant SNP per cluster
		vep_sig <- c()
		for(c in unique(vep_temp2$clust)) { 
			vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c,]$fdr) & vep_temp2$clust == c,]),1))
		}
		
		## Select non-significant SNPs to use as potential background matched SNPs
		vep_nonsig <- unique(vep_temp[vep_temp$fdr>0.05,])
		vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
		vep_sig_list <- unique(vep_sig[,c(1:2)])
	
		###First, determine how many matches so we can exclude those without any
		bg.samp.counts <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in c(1:nrow(vep_sig))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- nrow(vep_nonsig %>%
		 	filter( 
		    	Consequence==type,
				CHROM==chrom,
				((E < found*1.25) & (E > found*.75)),
				abs(POS - pos) > 50000))
			#append to background list
			bg.samp.counts <- rbind(bg.samp.counts,bg.samp.temp)
	 	 }
		
		#Establish background non-candidate set
	  	bg.samp.list <- c()
	  	#This loop finds a random set of matched background genes
	  	for(j in as.character(rownames(vep_sig[bg.samp.counts[,1]>4,]))) {
	  		#Select matched lists based on following criteria
			type <- vep_sig[j,]$Consequence
			chrom <- vep_sig[j,]$CHROM
			found <- vep_sig[j,]$E
			pos <- vep_sig[j,]$POS
			#Pull one match per candidate gene
			bg.samp.temp <- sample_n(vep_nonsig %>%
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
		focal.afdiff <- vep_sig[,4] * (abs(vep_sig[,3])/vep_sig[,3])
		bg.afdiff <- bg.samp.list[,4] * (abs(bg.samp.list[,3])/bg.samp.list[,3])

		# Run the t-test
		assign(paste("ttest_LOO_SPT4clust_EvSPT4_topsig_no",SP_loo_cages[i],sep=""), rbind(get(paste("ttest_LOO_SPT4clust_EvSPT4_topsig_no",SP_loo_cages[i],sep="")),cbind(contrast=names(get(paste("freq_diff_t4_no",SP_loo_cages[i],"_bed",sep="")))[f],mean_focal=mean(focal.afdiff),mean_BG=mean(bg.afdiff),t(t.test(focal.afdiff, bg.afdiff, alternative = "greater")[c(1,2,3,7)]))))
	}
}

## Join all t-test results together and save
ttest_LOO_SPT4clust_EvSPT4_topsig <- rbind(cbind("no3","SP",ttest_LOO_SPT4clust_EvSPT4_topsig_no3),
	cbind("no7","SP",ttest_LOO_SPT4clust_EvSPT4_topsig_no7),
	cbind("no15","SP",ttest_LOO_SPT4clust_EvSPT4_topsig_no15),
	cbind("no33","SP",ttest_LOO_SPT4clust_EvSPT4_topsig_no33),
	cbind("no37","SP",ttest_LOO_SPT4clust_EvSPT4_topsig_no37))

write.table(ttest_LOO_SPT4clust_EvSPT4_topsig, file="rudflies_2023_redo.ttest_LOO_SPT4clust_EvSPT4_topsig_results.txt", sep = "\t", quote = FALSE, row.names = F)