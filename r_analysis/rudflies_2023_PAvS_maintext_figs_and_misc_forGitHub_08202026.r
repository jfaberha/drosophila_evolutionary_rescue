##########################################################################################
### This script is associated with a study on the evolution of spinosad-resistance in  ###
### field populations of Drosophila melanogaster run in 2023 by the Rudman lab (WSUV). ###
### This particular script contains primary plotting and statistical analysis for the  ###
### main manscript. Several supplementary plots and results were generated here too.   ###
##########################################################################################							

### In R ###
#configure r environment
setwd("/scratch/user/jfaberha/20260915_042547/admera/gp_analysis/rudflies_2023_redo/r/scripts/dryad_upload")
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
vep <- read.table("filtered-all.annot.vcf.FLYCADD.tsv", header=TRUE) 
## Sample metadata table
haf.meta <- read.table("rudflies_2023_meta.tsv", header=TRUE)
## Hafpipe imputed allele frequency table
haf.freq <- read.delim("rudflies_2023_hafpipe.csv", header=TRUE, sep = ",")

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

## Merge annotations for filtering and highlighting
glm.all.annot <- merge(glm.all, vep, by=c("CHROM","POS"))

## Load spinosad-resistance candidate gene list
spino.cand <- read.table("spino.cand.list.txt", header=FALSE)
names(spino.cand) <- "Gene"

## Now merge to find All SNPs in and around candidate genes
glm.all.annot.spino <- merge(spino.cand, glm.all.annot, by="Gene")

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


################
### Figure 2 ###
################

## Panel D
# This manhattan plot shows the contrast between S and E samples combined across all TPTs
manh.SvE.all.rolwin21 <- ggplot(na.omit(glm.all.rolwin21), aes(POS, SvE.logp.rolwin21)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though 
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.fdr.rolwin21 < 0.01,]$SvE.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in S vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.fdr.rolwin21 < 0.01,]), aes(POS, SvE.logp.rolwin21), color = "#BD2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=na.omit(glm.all.rolwin21.annot.spino),aes(POS, SvE.logp.rolwin21), color = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=na.omit(glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$SvE.fdr.rolwin21 < 0.01,]), aes(POS, SvE.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$SvE.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs S")
  
## Panel E
# This manhattan plot shows the contrast between PA and E samples combined across all TPTs
manh.PAvE.all.rolwin21 <- ggplot(na.omit(glm.all.rolwin21), aes(POS, PAvE.logp.rolwin21)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.fdr.rolwin21 < 0.01,]$PAvE.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in PA vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.fdr.rolwin21 < 0.01,]),aes(POS, PAvE.logp.rolwin21), color = "#3F0100", size = 1) +
	## Highlight points significant in S vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.fdr.rolwin21 < 0.01,]),aes(POS, PAvE.logp.rolwin21), color = "#BD2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, PAvE.logp.rolwin21), color = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=na.omit(glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$PAvE.fdr.rolwin21 < 0.01,]), aes(POS, PAvE.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$PAvE.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA")
  
## Save them together
pdf(file = "rudflies_2023_redo.SvE.all.PAvE.all.glm.manh_rolwin21.pdf", width=7.5, height=4)
	ggarrange(manh.SvE.all.rolwin21, manh.PAvE.all.rolwin21,
              ncol = 1, nrow = 2)
dev.off()



### Figure S3 ###
## Alternative version Figure 2 without rolling window
# This manhattan plot shows the contrast between S and E samples combined across all TPTs
manh.SvE.all <- ggplot(na.omit(glm.all), aes(POS, SvE.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though 
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$SvE.fdr < 0.01,]$SvE.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in S vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$SvE.fdr < 0.01,]), aes(POS, SvE.logp), color = "#BD2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=na.omit(glm.all.annot.spino),aes(POS, SvE.logp), color = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=na.omit(glm.all.annot.spino[glm.all.annot.spino$SvE.fdr < 0.01,]), aes(POS, SvE.logp), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$SvE.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs S")
  
## Panel E
# This manhattan plot shows the contrast between PA and E samples combined across all TPTs
manh.PAvE.all <- ggplot(na.omit(glm.all), aes(POS, PAvE.logp)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$PAvE.fdr < 0.01,]$PAvE.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in PA vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$PAvE.fdr < 0.01,]),aes(POS, PAvE.logp), color = "#3F0100", size = 1) +
	## Highlight points significant in S vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$SvE.fdr < 0.01,]),aes(POS, PAvE.logp), color = "#BD2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.annot.spino,aes(POS, PAvE.logp), color = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=na.omit(glm.all.annot.spino[glm.all.annot.spino$PAvE.fdr < 0.01,]), aes(POS, PAvE.logp), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$PAvE.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA")
  
## Save them together
pdf(file = "rudflies_2023_redo.SvE.all.PAvE.all.glm.manh_no_rolwin21.pdf", width=7.5, height=4)
	ggarrange(manh.SvE.all, manh.PAvE.all,
              ncol = 1, nrow = 2)
dev.off()


### Figure 2 AF change side panels ###

## General idea is that for spinosad resistance candidate SNPs we'll be comparing absolute 
## AF differences in "S vs E" and "PA vs E" contrasts and compare those AF differences to 
## those of matched background SNPs. We'll check whether they we see higher AF differences 
## for spinosad candidates than for matched SNPs. This would tell us whether candidate AFs
## are changing in spino-exposed populations more than expected by chance, suggesting 
## possible adaptative response. We'll do the same analysis for "S vs E" and "PA vs E" 
## glm contrast top outliers.

## We'll also check AF means from "S vs E" contrast significant SNPs in the "PA vs E" 
## contrast to see if we can detect some adaptive convergence.

## Start by calculating mean allele frequencies for each treat by TPT and across TPTs
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
  
## Merge AF means only with annotation info
freq_means_annot <- merge(freq_means, vep, by=c("CHROM","POS"))

## Make sure we don't have duplicate loci, generated by overlapping annotations
freq_means_annot_nr <- freq_means_annot[!duplicated(freq_means_annot[,c(1:2)]), ]

## Now filter for spino candidates only
freq_means_annot_nr_spino <- freq_means_annot_nr[freq_means_annot_nr$Gene %in% spino.cand$Gene,]

## Now filter for non-spino candidates as a background SNP list
freq_means_annot_nr_bg <- freq_means_annot_nr[!freq_means_annot_nr$Gene %in% spino.cand$Gene,]

## Merge raw GLM results table (pre-rolling window average) with annotation info 
glm.all.annot <- merge(glm.all, vep, by=c("CHROM","POS"))
## remove duplicate
glm.all.annot.nr <- glm.all.annot[!duplicated(glm.all.annot[,c(1:2)]), ]
## Since there are more significant "S vs E" SNPs than all N spinosad candidate SNPs, 
## make same-sized list of the top N SNPs from "S vs E" contrast.
glm.all.annot.nr.topEvS <- head(glm.all.annot.nr[order(glm.all.annot.nr$SvE.fdr), ],
	nrow(freq_means_annot_nr_spino))
## Make AF mean frequency table for top "S vs E" SNP list
freq_means_annot_nr_EvS <- merge(glm.all.annot.nr.topEvS[,c(1:2)], 
								freq_means_annot_nr, 
								by=c("CHROM","POS"))
## Make AF mean frequency table for "S vs E" background SNPs (not top outliers)
freq_means_annot_nr_EvS_bg <- freq_means_annot_nr[!freq_means_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvS$Uploaded_variation,]

## Since there are more significant "PA vs E" SNPs than all N spinosad candidate SNPs, 
## make same-sized list of the top N SNPs from "PA vs E" contrast.
glm.all.annot.nr.topEvPA <- head(glm.all.annot.nr[order(glm.all.annot.nr$PAvE.fdr), ],nrow(freq_means_annot_nr_spino))
## Make AF mean frequency table for top "PA vs E" SNP list
freq_means_annot_nr_EvPA <- merge(glm.all.annot.nr.topEvPA[,c(1:2)],
								freq_means_annot_nr, 
								by=c("CHROM","POS"))
## Make AF mean frequency table for "PA vs E" background SNPs (not top outliers)
freq_means_annot_nr_EvPA_bg <- freq_means_annot_nr[!freq_means_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvPA$Uploaded_variation,]

## The following loop will go through the list of spino candidates, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_results <- read.table("rudflies_2023_redo.freq_means_results.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_results <- c() #create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_annot_nr_spino))) { #number of sites in spino candidate genes
  tryCatch({
  	temp_chrom <- freq_means_annot_nr_spino[g,1] #find chrom
  	temp_e <- freq_means_annot_nr_spino[g,3] #find E mean frequency
  	temp_cons <- freq_means_annot_nr_spino[g,24] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_annot_nr_bg[freq_means_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ] 
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2S=abs(temp_bg$E.af-temp_bg$S.af), #raw diff for S vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2S.perc=abs((temp_bg$E.af-temp_bg$S.af)/temp_bg$E.af), #perc diff for S vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af) #find means for all bg SNP differences
  	names(temp_means) <- paste0(names(temp_means), ".mean") 
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_results <- rbind(freq_means_results,append(temp_means,temp_vars))
  }, error=function(e){})
}  

#save table so we don't need to run this loop every time
colnames(freq_means_results) <- paste0(colnames(freq_means_results), ".spino")
write.table(freq_means_results, file="rudflies_2023_redo.freq_means_results.txt",sep = "\t", quote = FALSE, row.names = F)


## The following loop will go through the list of "S vs E" top outliers, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_results_EvS <- read.table("rudflies_2023_redo.freq_means_results_EvS.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_results_EvS <- c() #create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_annot_nr_spino))) { #number of sites in spino candidate genes
  tryCatch({
  	temp_chrom <- freq_means_annot_nr_EvS[g,1] #find chrom
  	temp_e <- freq_means_annot_nr_EvS[g,3] #find E mean frequency
  	temp_cons <- freq_means_annot_nr_EvS[g,24] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_annot_nr_EvS_bg[freq_means_annot_nr_EvS_bg$CHROM %in% temp_chrom[1] &
  		freq_means_annot_nr_EvS_bg$Consequence %in% temp_cons &
  		freq_means_annot_nr_EvS_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_annot_nr_EvS_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ] 
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2S=abs(temp_bg$E.af-temp_bg$S.af), #raw diff for S vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2S.perc=abs((temp_bg$E.af-temp_bg$S.af)/temp_bg$E.af), #perc diff for S vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_results_EvS <- rbind(freq_means_results_EvS,append(temp_means,temp_vars))
  }, error=function(e){})
}  

#save table so we don't need to run this loop every time
colnames(freq_means_results_EvS) <- paste0(colnames(freq_means_results_EvS), ".EvS")
write.table(freq_means_results_EvS, file="rudflies_2023_redo.freq_means_results_EvS.txt",sep = "\t", quote = FALSE, row.names = F)


## The following loop will go through the list of "PA vs E" top outliers, find all 
## matching background background SNPs, sample 100 of them, then calculate absolute 
## frequency differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_results_EvPA <- read.table("rudflies_2023_redo.freq_means_results_EvPA.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_results_EvPA <- c() #create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_annot_nr_spino))) { #number of sites in spino candidate genes
  tryCatch({
  	temp_chrom <- freq_means_annot_nr_EvPA[g,1] #find chrom
  	temp_e <- freq_means_annot_nr_EvPA[g,3] #find E mean frequency
  	temp_cons <- freq_means_annot_nr_EvPA[g,24] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_annot_nr_EvPA_bg[freq_means_annot_nr_EvPA_bg$CHROM %in% temp_chrom[1] &
  		freq_means_annot_nr_EvPA_bg$Consequence %in% temp_cons &
  		freq_means_annot_nr_EvPA_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_annot_nr_EvPA_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),] 
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ]
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2S=abs(temp_bg$E.af-temp_bg$S.af),
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af),
  		bg.E2S.perc=abs((temp_bg$E.af-temp_bg$S.af)/temp_bg$E.af),
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_results_EvPA <- rbind(freq_means_results_EvPA,append(temp_means,temp_vars))
  }, error=function(e){})
}  

#save table so we don't need to run this loop every time
colnames(freq_means_results_EvPA) <- paste0(colnames(freq_means_results_EvPA), ".EvPA")
write.table(freq_means_results_EvPA, file="rudflies_2023_redo.freq_means_results_EvPA.txt",sep = "\t", quote = FALSE, row.names = F)


### Now that we've calculated background AF difference, time to calculate the same for 
### actual spino candidates and GLM outlier SNPs, then make a master table for each.

## First, make spino candidate table
freq_means_annot_nr_spino_afdiff <- cbind(freq_means_annot_nr_spino,
	cand.E2S.spino=abs(freq_means_annot_nr_spino$E.af-freq_means_annot_nr_spino$S.af),
	cand.E2PA.spino=abs(freq_means_annot_nr_spino$E.af-freq_means_annot_nr_spino$PA.af),
	cand.E2S.perc.spino=abs((freq_means_annot_nr_spino$E.af - freq_means_annot_nr_spino$S.af)/freq_means_annot_nr_spino$E.af),
	cand.E2PA.perc.spino=abs((freq_means_annot_nr_spino$E.af - freq_means_annot_nr_spino$PA.af)/freq_means_annot_nr_spino$E.af),
	freq_means_results)

## Second, make "E vs S" GLM outlier table
freq_means_annot_nr_EvS_afdiff <- cbind(freq_means_annot_nr_EvS,
	cand.E2S.EvS=abs(freq_means_annot_nr_EvS$E.af - freq_means_annot_nr_EvS$S.af),
	cand.E2PA.EvS=abs(freq_means_annot_nr_EvS$E.af - freq_means_annot_nr_EvS$PA.af), cand.E2S.perc.EvS=abs((freq_means_annot_nr_EvS$E.af - freq_means_annot_nr_EvS$S.af)/freq_means_annot_nr_EvS$E.af),
	cand.E2PA.perc.EvS=abs((freq_means_annot_nr_EvS$E.af - freq_means_annot_nr_EvS$PA.af)/freq_means_annot_nr_EvS$E.af),
	freq_means_results_EvS)

## Third, make "E vs PA" GLM outlier table
freq_means_annot_nr_EvPA_afdiff <- cbind(freq_means_annot_nr_EvPA,
	cand.E2S.EvPA=abs(freq_means_annot_nr_EvPA$E.af-freq_means_annot_nr_EvPA$S.af),
	cand.E2PA.EvPA=abs(freq_means_annot_nr_EvPA$E.af-freq_means_annot_nr_EvPA$PA.af),
	cand.E2S.perc.EvPA=abs((freq_means_annot_nr_EvPA$E.af - freq_means_annot_nr_EvPA$S.af)/freq_means_annot_nr_EvPA$E.af),
	cand.E2PA.perc.EvPA=abs((freq_means_annot_nr_EvPA$E.af - freq_means_annot_nr_EvPA$PA.af)/freq_means_annot_nr_EvPA$E.af),
	freq_means_results_EvPA)
	
## Third, make "E vs PA" GLM outlier table
freq_means_annot_nr_EvPA_afdiff <- cbind(freq_means_annot_nr_EvPA,
	cand.E2S.EvPA=abs(freq_means_annot_nr_EvPA$E.af-freq_means_annot_nr_EvPA$S.af),
	cand.E2PA.EvPA=abs(freq_means_annot_nr_EvPA$E.af-freq_means_annot_nr_EvPA$PA.af),
	cand.E2S.perc.EvPA=abs((freq_means_annot_nr_EvPA$E.af - freq_means_annot_nr_EvPA$S.af)/freq_means_annot_nr_EvPA$E.af),
	cand.E2PA.perc.EvPA=abs((freq_means_annot_nr_EvPA$E.af - freq_means_annot_nr_EvPA$PA.af)/freq_means_annot_nr_EvPA$E.af),
	freq_means_results_EvPA)


## Combine all relevant AF difference stats into a master table
freq_means_annot_nr_afdiff <- cbind(freq_means_annot_nr_spino_afdiff[,c(33:ncol(freq_means_annot_nr_spino_afdiff))], freq_means_annot_nr_EvS_afdiff[,c(33:ncol(freq_means_annot_nr_EvS_afdiff))], freq_means_annot_nr_EvPA_afdiff[,c(33:ncol(freq_means_annot_nr_EvPA_afdiff))])


## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs S differences only
freq_means_annot_nr_E2Sdiff <- data.frame(rbind(cbind("candidate","E2S","spino","raw",freq_means_annot_nr_afdiff$cand.E2S.spino),
	   cbind("candidate","E2S","EvS","raw", freq_means_annot_nr_afdiff$cand.E2S.EvS),
	   cbind("matched","E2S","spino","mean", freq_means_annot_nr_afdiff$bg.E2S.mean.spino),
	   cbind("matched","E2S","EvS","mean", freq_means_annot_nr_afdiff$bg.E2S.mean.EvS),
	   cbind("matched","E2S","spino","perc", freq_means_annot_nr_afdiff$bg.E2S.perc.mean.spino),
	   cbind("matched","E2S","EvS","perc", freq_means_annot_nr_afdiff$bg.E2S.perc.mean.EvS)))
	   
## Rename cols
colnames(freq_means_annot_nr_E2Sdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_annot_nr_E2Sdiff[,5] <- as.numeric(freq_means_annot_nr_E2Sdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs S" 
## outliers), and type (raw, mean, or perc)
freq_means_annot_nr_E2Sdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs S AF-change in E vs S outliers vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2Sdiff[freq_means_annot_nr_E2Sdiff$type!="perc" &
	freq_means_annot_nr_E2Sdiff$category=="EvS",],
	alternative = 'greater')

#T-test for differences in absolute E vs S AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2Sdiff[freq_means_annot_nr_E2Sdiff$type!="perc" &
	freq_means_annot_nr_E2Sdiff$category=="spino",],
	alternative = 'greater')

### Panel F ###
# make plot of percent AF differences by set and category
E2Sdiff <- ggplot(data = freq_means_annot_nr_E2Sdiff[freq_means_annot_nr_E2Sdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs S") + 
	theme_classic()

## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs PA differences only
freq_means_annot_nr_E2PAdiff <- data.frame(rbind(cbind("candidate","E2PA","spino","raw",freq_means_annot_nr_afdiff$cand.E2PA.spino),
	   cbind("candidate","E2PA","EvS","raw",freq_means_annot_nr_afdiff$cand.E2PA.EvS),
	   cbind("matched","E2PA","spino","mean",freq_means_annot_nr_afdiff$bg.E2PA.mean.spino),
	   cbind("matched","E2PA","EvS","mean",freq_means_annot_nr_afdiff$bg.E2PA.mean.EvS),
	   cbind("matched","E2PA","spino","perc",freq_means_annot_nr_afdiff$bg.E2PA.perc.mean.spino),
	   cbind("matched","E2PA","EvS","perc",freq_means_annot_nr_afdiff$bg.E2PA.perc.mean.EvS)))

## Rename cols
colnames(freq_means_annot_nr_E2PAdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_annot_nr_E2PAdiff[,5] <- as.numeric(freq_means_annot_nr_E2PAdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs S" 
## outliers), and type (raw, mean, or perc)
freq_means_annot_nr_E2PAdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs PA AF-change in E vs S outliers vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2PAdiff[freq_means_annot_nr_E2PAdiff$type!="perc" &
	freq_means_annot_nr_E2PAdiff$category=="EvS",],
	alternative = 'greater')

#T-test for differences in absolute E vs PA AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2PAdiff[freq_means_annot_nr_E2PAdiff$type!="perc" &
	freq_means_annot_nr_E2PAdiff$category=="spino",],
	alternative = 'greater')


### Panel G ###
# make plot of percent AF differences by set and category, specifically EvS outliers and 
# spino candidates within EvPA contrast
E2PAdiff <- ggplot(data = freq_means_annot_nr_E2PAdiff[freq_means_annot_nr_E2PAdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs PA") + 
	theme_classic()


## Save panel F & G plots together
pdf(file = "rudflies_2023_redo.freq_means_annot_nr_E2Sdiff.mean.multiboxplot.pdf", width=3, height=6)
	ggarrange(E2Sdiff, E2PAdiff,
              ncol = 1, nrow = 2)
dev.off()


### Are there differences in AF changes overall? ###
## Calculate all allele frequency differences, both raw and percent changes
freq_means_afdiff <- as.data.frame(cbind(E2S=abs(freq_means$E.af-freq_means$S.af),
	E2PA=abs(freq_means$E.af-freq_means$PA.af),
	E2S.perc=abs((freq_means$E.af - freq_means$S.af)/freq_means$E.af),
	E2PA.perc=abs((freq_means$E.af - freq_means$PA.af)/freq_means$E.af)))

## Reformat for running T-tests
freq_means_afdiff <- data.frame(rbind(cbind(contrast="E2S",
		af_diff=freq_means_afdiff$E2S,
		perc_diff=freq_means_afdiff$E2S.perc),
	cbind(contrast="E2PA",
		af_diff=freq_means_afdiff$E2PA,
		perc_diff=freq_means_afdiff$E2PA.perc)))
		
## Set AF differences to numeric format
freq_means_afdiff[,2] <- as.numeric(freq_means_afdiff[,2])
freq_means_afdiff[,3] <- as.numeric(freq_means_afdiff[,3])

## T-test for differences in all pairwise absolute AF-changes
t.test(af_diff ~ contrast, 
	data = freq_means_afdiff,
	alternative = 'two.sided')

## T-test for differences in all pairwise percent AF-changes
t.test(perc_diff ~ contrast, 
	data = freq_means_afdiff,
	alternative = 'two.sided')


### Alternate Panel G ###
# make plot of percent AF differences by set and category, specifically EvPA outliers and 
# spino candidates within EvPA contrast

## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs PA differences only
freq_means_annot_nr_E2PAdiff_EvPA <- data.frame(rbind(cbind("candidate","E2PA","spino","raw",freq_means_annot_nr_afdiff$cand.E2PA.spino),
	   cbind("candidate","E2PA","EvPA","raw", freq_means_annot_nr_afdiff$cand.E2PA.EvPA),
	   cbind("matched","E2PA","spino","mean", freq_means_annot_nr_afdiff$bg.E2PA.mean.spino),
	   cbind("matched","E2PA","EvPA","mean", freq_means_annot_nr_afdiff$bg.E2PA.mean.EvPA),
	   cbind("matched","E2PA","spino","perc", freq_means_annot_nr_afdiff$bg.E2PA.perc.mean.spino),
	   cbind("matched","E2PA","EvPA","perc", freq_means_annot_nr_afdiff$bg.E2PA.perc.mean.EvPA)))
	   
## Rename cols
colnames(freq_means_annot_nr_E2PAdiff_EvPA) <- c("set","contrast","category","type","af_diff")
## PAet AF differences to numeric format
freq_means_annot_nr_E2PAdiff_EvPA[,5] <- as.numeric(freq_means_annot_nr_E2PAdiff_EvPA[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs PA" 
## outliers), and type (raw, mean, or perc)
freq_means_annot_nr_E2PAdiff_EvPA %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs PA AF-change in E vs PA outliers vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2PAdiff_EvPA[freq_means_annot_nr_E2PAdiff_EvPA$type!="perc" &
	freq_means_annot_nr_E2PAdiff_EvPA$category=="EvPA",],
	alternative = 'greater')

#T-test for differences in absolute E vs PA AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_annot_nr_E2PAdiff_EvPA[freq_means_annot_nr_E2PAdiff_EvPA$type!="perc" &
	freq_means_annot_nr_E2PAdiff_EvPA$category=="spino",],
	alternative = 'greater')


# make plot of percent AF differences by set and category
E2PAdiff_EvPA <- ggplot(data = freq_means_annot_nr_E2PAdiff_EvPA[freq_means_annot_nr_E2PAdiff_EvPA$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#3F0100")) +
	ggtitle("AF differences: E vs PA") + 
	theme_classic()


## Save panel F & alternate G plots together
pdf(file = "rudflies_2023_redo.freq_means_annot_nr_E2PAdiff_EvPA.mean.multiboxplot.pdf", width=3, height=6)
	ggarrange(E2Sdiff, E2PAdiff_EvPA,
              ncol = 1, nrow = 2)
dev.off()


	
	
################
### Figure 3 ###
################

### These first versions use unsmoothed GLM p-values, -log10(p) values, and FDR values.
### They have been relegated to the supplement.

## Panel D
# This manhattan plot shows the contrast between SE and E samples combined in TPT 1
manh.EvSE.T1 <- ggplot(glm.all, aes(POS, EvSE.T1.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$EvSE.T1.fdr < 0.05,]$EvSE.T1.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr < 0.05,]),aes(POS, EvSE.T1.logp), color = "#BF2D23", size = 1) +
	## Highlight points significant in SP vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$EvSP.T1.fdr<0.05,]),aes(POS, EvSE.T1.logp), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E and SP vs E GLM contrasts
	geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05 & glm.all$EvSP.T1.fdr<0.05,]),aes(POS, EvSE.T1.logp), color = "yellow", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	#geom_point(data=glm.all.annot.spino,aes(POS, EvSE.T1.logp), colour = "#C45C5C", size = 1) +
## Highlight significant spinosad resistance candidate SNPs
	#geom_point(data=glm.all.annot.spino[glm.all.annot.spino$EvSE.T1.fdr < 0.05,],aes(POS, EvSE.T1.logp), color = "black", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$EvSE.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SE (TPT1)")

## Panel E
# This manhattan plot shows the contrast between SP and E samples combined in TPT 1
manh.EvSP.T1 <- ggplot(glm.all, aes(POS, EvSP.T1.logp)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$EvSP.T1.fdr < 0.05,]$EvSP.T1.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SP vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$EvSP.T1.fdr<0.05,]),aes(POS, EvSP.T1.logp), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05,]),aes(POS, EvSP.T1.logp), color = "#BF2D23", size = 1) +
	## Highlight points significant in SE vs E and SP vs E GLM contrasts
	geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05 & glm.all$EvSP.T1.fdr<0.05,]),aes(POS, EvSP.T1.logp), color = "yellow", size = 1) +
	## Highlight spinosad resistance candidate SNPs
#geom_point(data=glm.all.annot.spino,aes(POS, EvSP.T1.logp), colour = "#C45C5C", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	#geom_point(data=glm.all.annot.spino[glm.all.annot.spino$EvSP.T1.fdr < 0.05,],aes(POS, EvSP.T1.logp), color = "black", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$EvSP.T1.logp))+0.5)) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SP (TPT1)")

## Unused in figure 3
# This manhattan plot shows the contrast between SE and SP samples combined in TPT 1
manh.SEvSP.T1 <- ggplot(glm.all, aes(POS, SEvSP.T1.logp)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$SEvSP.T1.fdr < 0.05,]$SEvSP.T1.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr < 0.05,]),aes(POS, SEvSP.T1.logp), color = "#BF2D23", size = 1) +
	## Highlight points significant in SP vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$EvSP.T1.fdr<0.05,]),aes(POS, SEvSP.T1.logp), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E and SP vs E GLM contrasts
	geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05 & glm.all$EvSP.T1.fdr<0.05,]),aes(POS, SEvSP.T1.logp), color = "yellow", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	#geom_point(data=glm.all.annot.spino,aes(POS, SEvSP.T1.logp), colour = "#C45C5C", size = 1) +
	#geom_point(data=glm.all.annot.spino[glm.all.annot.spino$SEvSP.T1.fdr < 0.05,],aes(POS, SEvSP.T1.logp), color = "black", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$SEvSP.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("SE vs SP (TPT1)")

## Panel F
# This manhattan plot shows the contrast between PA and E samples combined in TPT 1
manh.PAvE.T1 <- ggplot(glm.all, aes(POS, PAvE.T1.logp)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$PAvE.T1.fdr < 0.05,]$PAvE.T1.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in PA vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$PAvE.T1.fdr < 0.05,]),aes(POS, PAvE.T1.logp), color = "#3F0000", size = 1) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr < 0.05 & glm.all$PAvE.T1.fdr < 0.05,]),aes(POS, PAvE.T1.logp), color = "#BF2D23", size = 1) +
	## Highlight points significant in SP vs E GLM contrast
	geom_point(data=na.omit(glm.all[glm.all$EvSP.T1.fdr<0.05,]),aes(POS, PAvE.T1.logp), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E and SP vs E GLM contrasts
	geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05 & glm.all$EvSP.T1.fdr<0.05,]),aes(POS, PAvE.T1.logp), color = "yellow", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	#geom_point(data=glm.all.annot.spino,aes(POS, PAvE.T1.logp), colour = "#C45C5C", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	#geom_point(data=glm.all.annot.spino[glm.all.annot.spino$PAvE.T1.fdr < 0.05,],aes(POS, PAvE.T1.logp), color = "black", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$PAvE.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA (TPT1)")

## Plot all four pairwise TPT1 treatment comparisons together (not used)
pdf(file = "rudflies_2023_redo.EvSE.EvSP.SEvSP.PAvE.T1.glm.manh.pdf", width=7.5, height=8)
	ggarrange(manh.EvSE.T1, manh.EvSP.T1, manh.SEvSP.T1, manh.PAvE.T1,
              ncol = 1, nrow = 4)
dev.off()

## Plot three most relevant pairwise TPT1 treatment comparisons (not SP vs SE) together
pdf(file = "rudflies_2023_redo.EvSE.EvSP.PAvE.T1.glm.manh.pdf", width=7.5, height=6)
	ggarrange(manh.EvSE.T1, manh.EvSP.T1, manh.PAvE.T1,
              ncol = 1, nrow = 3)
dev.off()

### Alternate plot (Figure S4) that shows significance of spino candidates instead ###
### of GLM outliers.															   ###

## Figure S4 Panel A
# This manhattan plot shows the contrast between SE and E samples combined in TPT 1
manh.EvSE.T1.spino <- ggplot(glm.all, aes(POS, EvSE.T1.logp)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$EvSE.T1.fdr < 0.05,]$EvSE.T1.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr < 0.05,]),aes(POS, EvSE.T1.logp), color = "#BF2D23", size = 1) +
	## Highlight points significant in SP vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$EvSP.T1.fdr<0.05,]),aes(POS, EvSE.T1.logp), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E and SP vs E GLM contrasts
	#geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05 & glm.all$EvSP.T1.fdr<0.05,]),aes(POS, EvSE.T1.logp), color = "yellow", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.annot.spino,aes(POS, EvSE.T1.logp), colour = "#809BCA", size = 1) +
## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=glm.all.annot.spino[glm.all.annot.spino$EvSE.T1.fdr < 0.05,],aes(POS, EvSE.T1.logp), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$EvSE.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SE (TPT1)")

## Figure S4 Panel B
# This manhattan plot shows the contrast between SP and E samples combined in TPT 1
manh.EvSP.T1.spino <- ggplot(glm.all, aes(POS, EvSP.T1.logp)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$EvSP.T1.fdr < 0.05,]$EvSP.T1.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SP vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$EvSP.T1.fdr<0.05,]),aes(POS, EvSP.T1.logp), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05,]),aes(POS, EvSP.T1.logp), color = "#BF2D23", size = 1) +
	## Highlight points significant in SE vs E and SP vs E GLM contrasts
	#geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05 & glm.all$EvSP.T1.fdr<0.05,]),aes(POS, EvSP.T1.logp), color = "yellow", size = 1) +
	## Highlight spinosad resistance candidate SNPs
geom_point(data=glm.all.annot.spino,aes(POS, EvSP.T1.logp), colour = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=glm.all.annot.spino[glm.all.annot.spino$EvSP.T1.fdr < 0.05,],aes(POS, EvSP.T1.logp), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$EvSP.T1.logp))+0.5)) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SP (TPT1)")

## Unused in Figure S4
# This manhattan plot shows the contrast between SE and SP samples combined in TPT 1
manh.SEvSP.T1.spino <- ggplot(glm.all, aes(POS, SEvSP.T1.logp)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$SEvSP.T1.fdr < 0.05,]$SEvSP.T1.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr < 0.05,]),aes(POS, SEvSP.T1.logp), color = "#BF2D23", size = 1) +
	## Highlight points significant in SP vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$EvSP.T1.fdr<0.05,]),aes(POS, SEvSP.T1.logp), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E and SP vs E GLM contrasts
	#geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05 & glm.all$EvSP.T1.fdr<0.05,]),aes(POS, SEvSP.T1.logp), color = "yellow", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.annot.spino,aes(POS, SEvSP.T1.logp), colour = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=glm.all.annot.spino[glm.all.annot.spino$SEvSP.T1.fdr < 0.05,],aes(POS, SEvSP.T1.logp), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$SEvSP.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("SE vs SP (TPT1)")

## Figure S4 Panel C
# This manhattan plot shows the contrast between PA and E samples combined in TPT 1
manh.PAvE.T1.spino <- ggplot(glm.all, aes(POS, PAvE.T1.logp)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	geom_hline(yintercept=min(na.omit(glm.all[glm.all$PAvE.T1.fdr < 0.05,]$PAvE.T1.logp)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in PA vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$PAvE.T1.fdr < 0.05,]),aes(POS, PAvE.T1.logp), color = "#3F0000", size = 1) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr < 0.05 & glm.all$PAvE.T1.fdr < 0.05,]),aes(POS, PAvE.T1.logp), color = "#BF2D23", size = 1) +
	## Highlight points significant in SP vs E GLM contrast
	#geom_point(data=na.omit(glm.all[glm.all$EvSP.T1.fdr<0.05,]),aes(POS, PAvE.T1.logp), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E and SP vs E GLM contrasts
	#geom_point(data=na.omit(glm.all[glm.all$EvSE.T1.fdr<0.05 & glm.all$EvSP.T1.fdr<0.05,]),aes(POS, PAvE.T1.logp), color = "yellow", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.annot.spino,aes(POS, PAvE.T1.logp), colour = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=glm.all.annot.spino[glm.all.annot.spino$PAvE.T1.fdr < 0.05,],aes(POS, PAvE.T1.logp), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all$PAvE.T1.logp)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA (TPT1)")

## Plot all four pairwise TPT1 treatment comparisons together (not used)
pdf(file = "rudflies_2023_redo.EvSE.EvSP.SEvSP.PAvE.T1.glm.manh.spino.pdf", width=7.5, height=8)
	ggarrange(manh.EvSE.T1.spino, manh.EvSP.T1.spino, manh.SEvSP.T1.spino, manh.PAvE.T1.spino,
              ncol = 1, nrow = 4)
dev.off()

## Plot three most relevant pairwise TPT1 treatment comparisons (not SP vs SE) together
pdf(file = "rudflies_2023_redo.EvSE.EvSP.PAvE.T1.glm.manh.spino.pdf", width=7.5, height=6)
	ggarrange(manh.EvSE.T1.spino, manh.EvSP.T1.spino, manh.PAvE.T1.spino,
              ncol = 1, nrow = 3)
dev.off()


### These second versions use unsmoothed GLM p-values, -log10(p) values, and FDR values.
### The plots showing GLM outliers are in the main MS as figure 3, while the spino
### candidate versions will be in the supplement as figure S4???.

## Panel D
# This manhattan plot shows the contrast between SE and E samples combined in TPT 1
manh.EvSE.T1.rolwin21 <- ggplot(glm.all.rolwin21, aes(POS, EvSE.T1.logp.rolwin21)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05,]$EvSE.T1.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05,]),aes(POS, EvSE.T1.logp.rolwin21), color = "#BF2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	#geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, EvSE.T1.logp.rolwin21), colour = "#809BCA", size = 1) +
## Highlight significant spinosad resistance candidate SNPs
	#geom_point(data=glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$EvSE.T1.fdr.rolwin21 < 0.05,],aes(POS, EvSE.T1.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$EvSE.T1.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SE (TPT1): 21-SNP sliding window")

## Panel E
# This manhattan plot shows the contrast between SP and E samples combined in TPT 1
manh.EvSP.T1.rolwin21 <- ggplot(glm.all.rolwin21, aes(POS, EvSP.T1.logp.rolwin21)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSP.T1.fdr.rolwin21 < 0.05,]$EvSP.T1.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSP.T1.fdr.rolwin21<0.05,]),aes(POS, EvSP.T1.logp.rolwin21), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E GLM contrast but not in SP vs. E
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21$EvSP.T1.fdr.rolwin21>0.05,]),aes(POS, EvSP.T1.logp.rolwin21), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21$EvSP.T1.fdr.rolwin21<0.05,]),aes(POS, EvSP.T1.logp.rolwin21), color = "#BF2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
#geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, EvSP.T1.logp.rolwin21), colour = "#809BCA", size = 1) +
## Highlight significant spinosad resistance candidate SNPs
#geom_point(data=glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$EvSP.T1.fdr.rolwin21 < 0.05,],aes(POS, EvSP.T1.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$EvSP.T1.logp.rolwin21))+0.5)) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SP (TPT1): 21-SNP sliding window")

## Unused in figure 3
# This manhattan plot shows the contrast between SE and SP samples combined in TPT 1
manh.SEvSP.T1.rolwin21 <- ggplot(glm.all.rolwin21, aes(POS, SEvSP.T1.logp.rolwin21)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SEvSP.T1.fdr.rolwin21 < 0.05,]$SEvSP.T1.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05,]),aes(POS, SEvSP.T1.logp.rolwin21), color = "#BF2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	#geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, SEvSP.T1.logp.rolwin21), colour = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	#geom_point(data=glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$SEvSP.T1.fdr.rolwin21 < 0.05,],aes(POS, SEvSP.T1.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$SEvSP.T1.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("SE vs SP (TPT1): 21-SNP sliding window")

## Panel F
# This manhattan plot shows the contrast between PA and E samples combined in TPT 1
manh.PAvE.T1.rolwin21 <- ggplot(glm.all.rolwin21, aes(POS, PAvE.T1.logp.rolwin21)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T1.fdr.rolwin21 < 0.05,]$PAvE.T1.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in PA vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T1.fdr.rolwin21 < 0.05,]),aes(POS, PAvE.T1.logp.rolwin21), color = "#3F0000", size = 1) +
	## Highlight points significant in SE vs E GLM contrast but not in PA vs. E
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05 & glm.all.rolwin21$PAvE.T1.fdr.rolwin21 > 0.05,]),aes(POS, PAvE.T1.logp.rolwin21), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E GLM contrast and in PA vs. E
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05 & glm.all.rolwin21$PAvE.T1.fdr.rolwin21 < 0.05,]),aes(POS, PAvE.T1.logp.rolwin21), color = "#BF2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	#geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, PAvE.T1.logp.rolwin21), colour = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	#geom_point(data=glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$PAvE.T1.fdr.rolwin21 < 0.05,],aes(POS, PAvE.T1.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$PAvE.T1.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA (TPT1): 21-SNP sliding window")

## Plot all four pairwise TPT1 treatment comparisons together (not used)
pdf(file = "rudflies_2023_redo.EvSE.EvSP.SEvSP.PAvE.T1.rolwin21.glm.manh.pdf", width=7.5, height=8)
	ggarrange(manh.EvSE.T1.rolwin21, manh.EvSP.T1.rolwin21, manh.SEvSP.T1.rolwin21, manh.PAvE.T1.rolwin21,
              ncol = 1, nrow = 4)
dev.off()

## Plot three most relevant pairwise TPT1 treatment comparisons (not SP vs SE) together
pdf(file = "rudflies_2023_redo.EvSE.EvSP.PAvE.T1.rolwin21.glm.manh.pdf", width=7.5, height=6)
	ggarrange(manh.EvSE.T1.rolwin21, manh.EvSP.T1.rolwin21, manh.PAvE.T1.rolwin21,
              ncol = 1, nrow = 3)
dev.off()


### Alternate plot (Figure S4) that shows significance of spino candidates instead ###
### of GLM outliers.															   ###

## Figure S4 Panel A
# This manhattan plot shows the contrast between SE and E samples combined in TPT 1
manh.EvSE.T1.rolwin21.spino <- ggplot(glm.all.rolwin21, aes(POS, EvSE.T1.logp.rolwin21)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05,]$EvSE.T1.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05,]),aes(POS, EvSE.T1.logp.rolwin21), color = "#BF2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, EvSE.T1.logp.rolwin21), colour = "#809BCA", size = 1) +
## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$EvSE.T1.fdr.rolwin21 < 0.05,],aes(POS, EvSE.T1.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$EvSE.T1.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SE (TPT1): 21-SNP sliding window")

## Figure S4 Panel B
# This manhattan plot shows the contrast between SP and E samples combined in TPT 1
manh.EvSP.T1.rolwin21.spino <- ggplot(glm.all.rolwin21, aes(POS, EvSP.T1.logp.rolwin21)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSP.T1.fdr.rolwin21 < 0.05,]$EvSP.T1.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSP.T1.fdr.rolwin21<0.05,]),aes(POS, EvSP.T1.logp.rolwin21), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E GLM contrast but not in SP vs. E
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21$EvSP.T1.fdr.rolwin21>0.05,]),aes(POS, EvSP.T1.logp.rolwin21), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21$EvSP.T1.fdr.rolwin21<0.05,]),aes(POS, EvSP.T1.logp.rolwin21), color = "#BF2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, EvSP.T1.logp.rolwin21), colour = "#809BCA", size = 1) +
## Highlight significant spinosad resistance candidate SNPs
geom_point(data=glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$EvSP.T1.fdr.rolwin21 < 0.05,],aes(POS, EvSP.T1.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$EvSP.T1.logp.rolwin21))+0.5)) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs SP (TPT1): 21-SNP sliding window")

## Unused in figure S4
# This manhattan plot shows the contrast between SE and SP samples combined in TPT 1
manh.SEvSP.T1.rolwin21.spino <- ggplot(glm.all.rolwin21, aes(POS, SEvSP.T1.logp.rolwin21)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SEvSP.T1.fdr.rolwin21 < 0.05,]$SEvSP.T1.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05,]),aes(POS, SEvSP.T1.logp.rolwin21), color = "#BF2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, SEvSP.T1.logp.rolwin21), colour = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$SEvSP.T1.fdr.rolwin21 < 0.05,],aes(POS, SEvSP.T1.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$SEvSP.T1.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("SE vs SP (TPT1): 21-SNP sliding window")

## Figure S4 Panel C
# This manhattan plot shows the contrast between PA and E samples combined in TPT 1
manh.PAvE.T1.rolwin21.spino <- ggplot(glm.all.rolwin21, aes(POS, PAvE.T1.logp.rolwin21)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T1.fdr.rolwin21 < 0.05,]$PAvE.T1.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in PA vs E GLM contrast
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T1.fdr.rolwin21 < 0.05,]),aes(POS, PAvE.T1.logp.rolwin21), color = "#3F0000", size = 1) +
	## Highlight points significant in SE vs E GLM contrast but not in PA vs. E
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05 & glm.all.rolwin21$PAvE.T1.fdr.rolwin21 > 0.05,]),aes(POS, PAvE.T1.logp.rolwin21), color = "#DE9691", size = 1) +
	## Highlight points significant in SE vs E GLM contrast and in PA vs. E
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05 & glm.all.rolwin21$PAvE.T1.fdr.rolwin21 < 0.05,]),aes(POS, PAvE.T1.logp.rolwin21), color = "#BF2D23", size = 1) +
	## Highlight spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin21.annot.spino,aes(POS, PAvE.T1.logp.rolwin21), colour = "#809BCA", size = 1) +
	## Highlight significant spinosad resistance candidate SNPs
	geom_point(data=glm.all.rolwin21.annot.spino[glm.all.rolwin21.annot.spino$PAvE.T1.fdr.rolwin21 < 0.05,],aes(POS, PAvE.T1.logp.rolwin21), color = "#025197", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$PAvE.T1.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA (TPT1): 21-SNP sliding window")

## Plot all four pairwise TPT1 treatment comparisons together (not used)
pdf(file = "rudflies_2023_redo.EvSE.EvSP.SEvSP.PAvE.T1.rolwin21.spino.glm.manh.pdf", width=7.5, height=8)
	ggarrange(manh.EvSE.T1.rolwin21.spino, manh.EvSP.T1.rolwin21.spino, manh.SEvSP.T1.rolwin21.spino, manh.PAvE.T1.rolwin21.spino,
              ncol = 1, nrow = 4)
dev.off()

## Plot three most relevant pairwise TPT1 treatment comparisons (not SP vs SE) together
pdf(file = "rudflies_2023_redo.EvSE.EvSP.PAvE.T1.rolwin21.spino.glm.manh.pdf", width=7.5, height=6)
	ggarrange(manh.EvSE.T1.rolwin21.spino, manh.EvSP.T1.rolwin21.spino, manh.PAvE.T1.rolwin21.spino,
              ncol = 1, nrow = 3)
dev.off()



### Plot allele frequencies for missense outlier SNP in both EvSE and EvSP (no rolwin)
### *there are no SNPs significant in both these contrasts after window-based smoothing
## 2L_7692533 
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="2L" & haf.sites.filt$POS=="7692533",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.2L_7692533.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("2L:7692533") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()


### AF change side panels ###

## This time, for spinosad resistance candidate SNPs, we'll be comparing absolute 
## AF differences in TPT1 pairwise treatment contrasts and compare those AF differences to 
## those of matched background SNPs. We'll check whether they we see higher AF differences 
## for spinosad candidates than for matched SNPs. This would tell us whether candidate AFs
## are changing in spino-exposed populations more than expected by chance, suggesting 
## possible adaptative response. We'll do the same analysis for "SE vs E" and "PA vs E" 
## glm contrast top outliers. 

## In TPT1, the S treatment was divided into SE (extinct cages) and SP (persistent cages)

## We'll also check AF means from "SE vs E" contrast significant SNPs in other contrasts 
## to see if we can detect some adaptive convergence. We would have done the same with 
## "SP vs E" outliers, but there were none.

## Start by calculating mean allele frequencies for each treatment at TPT1
freq_means_t1 <- cbind(haf.sites.T1filt,
	E.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "E" & haf.meta.T1filt$tpt=="1"]),
	S.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "S" & haf.meta.T1filt$tpt=="1"]),
	SE.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition == "SE" & haf.meta.T1filt$tpt=="1"]),
	SP.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition == "SP" & haf.meta.T1filt$tpt=="1"]),
	PA.af=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix == "PA"&haf.meta.T1filt$tpt=="1"]))


## Merge AF means only with annotation info
freq_means_t1_annot <- merge(freq_means_t1, vep, by=c("CHROM","POS"))

## Make sure we don't have duplicate loci, generated by overlapping annotations
freq_means_t1_annot_nr <- freq_means_t1_annot[!duplicated(freq_means_t1_annot[,c(1:2)]), ]

## Now filter for spino candidates only
freq_means_t1_annot_nr_spino <- freq_means_t1_annot_nr[freq_means_t1_annot_nr$Gene %in% spino.cand$Gene,]

## Now filter for non-spino candidates as a background SNP list
freq_means_t1_annot_nr_bg <- freq_means_t1_annot_nr[!freq_means_t1_annot_nr$Gene %in% spino.cand$Gene,]

## Since the number of sig "SE vs E" SNPs is different from all N spinosad candidate SNPs, 
## make same-sized list of the top N SNPs from "SE vs E" contrast.
glm.all.annot.nr.topEvSE <- head(glm.all.annot.nr[order(glm.all.annot.nr$EvSE.T1.fdr), ],
	nrow(freq_means_t1_annot_nr_spino))
## Make AF mean frequency table for top "SE vs E" SNP list
freq_means_t1_annot_nr_EvSE <- merge(glm.all.annot.nr.topEvSE[,c(1:2)],
								freq_means_t1_annot_nr,
								by=c("CHROM","POS"))
## Make AF mean frequency table for "SE vs E" background SNPs (not top outliers)
freq_means_t1_annot_nr_EvSE_bg <- freq_means_t1_annot_nr[!freq_means_t1_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvSE$Uploaded_variation,]

## Since the number of sig "PA vs E" SNPs is different from all N spinosad candidate SNPs,
## make same-sized list of the top N SNPs from "PA vs E" contrast.
glm.all.annot.nr.topEvPA <- head(glm.all.annot.nr[order(glm.all.annot.nr$PAvE.T1.fdr), ],nrow(freq_means_t1_annot_nr_spino))
## Make AF mean frequency table for top "PA vs E" SNP list
freq_means_t1_annot_nr_EvPA <- merge(glm.all.annot.nr.topEvPA[,c(1:2)],
								freq_means_t1_annot_nr,
								by=c("CHROM","POS"))
## Make AF mean frequency table for "PA vs E" background SNPs (not top outliers)
freq_means_t1_annot_nr_EvPA_bg <- freq_means_t1_annot_nr[!freq_means_t1_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvPA$Uploaded_variation,]

## Since there are no sig "SP vs E" SNPs, make a list of the top N SNPs from "SP vs E" 
## contrast to match all N spinosad candidate SNPs.
glm.all.annot.nr.topEvSP <- head(glm.all.annot.nr[order(glm.all.annot.nr$EvSP.T1.fdr), ],
	nrow(freq_means_t1_annot_nr_spino))
## Make AF mean frequency table for top "SP vs E" SNP list
freq_means_t1_annot_nr_EvSP <- merge(glm.all.annot.nr.topEvSP[,c(1:2)],
								freq_means_t1_annot_nr,
								by=c("CHROM","POS"))
## Make AF mean frequency table for "SP vs E" background SNPs (not top outliers)
freq_means_t1_annot_nr_EvSP_bg <- freq_means_t1_annot_nr[!freq_means_t1_annot_nr$Uploaded_variation %in% glm.all.annot.nr.topEvSP$Uploaded_variation,]


## Instead of running full loop, you can just reload the table from a previous run
freq_means_t1_results <- read.table("rudflies_2023_redo.freq_means_t1_results.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_t1_results <- c()#create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_t1_annot_nr_spino))) { #number of sites in spino candidate genes
  tryCatch({
  	temp_chrom <- freq_means_t1_annot_nr_spino[g,1] #find chrom
  	temp_e <- freq_means_t1_annot_nr_spino[g,3] #find E mean frequency
  	temp_cons <- freq_means_t1_annot_nr_spino[g,14] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_t1_annot_nr_bg[freq_means_t1_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_t1_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_t1_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_t1_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ] #randomly sample 100 SNPs
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2SE=abs(temp_bg$E.af-temp_bg$SE.af),  #raw diff for SE vs E
  		bg.E2SP=abs(temp_bg$E.af-temp_bg$SP.af),  #raw diff for SP vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2SE.perc=abs((temp_bg$E.af-temp_bg$SE.af)/temp_bg$E.af),#perc diff for SE vs E
  		bg.E2SP.perc=abs((temp_bg$E.af-temp_bg$SP.af)/temp_bg$E.af),#perc diff for SP vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af) #find means for all bg SNP differences
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_t1_results <- rbind(freq_means_t1_results,append(temp_means,temp_vars))
  }, error=function(e){})
}

#save table so we don't need to run this loop every time
colnames(freq_means_t1_results) <- paste0(colnames(freq_means_t1_results), ".spino")
write.table(freq_means_t1_results, file="rudflies_2023_redo.freq_means_t1_results.txt",sep = "\t", quote = FALSE, row.names = F)

## The following loop will go through list of "SE vs E" top outliers, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_t1_results_EvSE <- read.table("rudflies_2023_redo.freq_means_t1_results_EvSE.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_t1_results_EvSE <- c()#create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_t1_annot_nr_EvSE))) { #number of sites in EvSE candidate genes
  tryCatch({
  	temp_chrom <- freq_means_t1_annot_nr_EvSE[g,1] #find chrom
  	temp_e <- freq_means_t1_annot_nr_EvSE[g,3] #find E mean frequency
  	temp_cons <- freq_means_t1_annot_nr_EvSE[g,14] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_t1_annot_nr_bg[freq_means_t1_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_t1_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_t1_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_t1_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ]
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2SE=abs(temp_bg$E.af-temp_bg$SE.af), #raw diff for SE vs E
  		bg.E2SP=abs(temp_bg$E.af-temp_bg$SP.af), #raw diff for SP vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2SE.perc=abs((temp_bg$E.af-temp_bg$SE.af)/temp_bg$E.af),#perc diff for SE vs E
  		bg.E2SP.perc=abs((temp_bg$E.af-temp_bg$SP.af)/temp_bg$E.af),#perc diff for SP vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_t1_results_EvSE <- rbind(freq_means_t1_results_EvSE,append(temp_means,temp_vars))
  }, error=function(e){})
}

#save table so we don't need to run this loop every time
colnames(freq_means_t1_results_EvSE) <- paste0(colnames(freq_means_t1_results_EvSE), ".EvSE.T1")
write.table(freq_means_t1_results_EvSE, file="rudflies_2023_redo.freq_means_t1_results_EvSE.txt",sep = "\t", quote = FALSE, row.names = F)


## The following loop will go through list of "SP vs E" top outliers, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_t1_results_EvSP <- read.table("rudflies_2023_redo.freq_means_t1_results_EvSP.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_t1_results_EvSP <- c()#create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_t1_annot_nr_EvSP))) { #number of sites in EvSP candidate genes
  tryCatch({
  	temp_chrom <- freq_means_t1_annot_nr_EvSP[g,1] #find chrom
  	temp_e <- freq_means_t1_annot_nr_EvSP[g,3] #find E mean frequency
  	temp_cons <- freq_means_t1_annot_nr_EvSP[g,14] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_t1_annot_nr_bg[freq_means_t1_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_t1_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_t1_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_t1_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ]
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2SE=abs(temp_bg$E.af-temp_bg$SE.af), #raw diff for SE vs E
  		bg.E2SP=abs(temp_bg$E.af-temp_bg$SP.af), #raw diff for SP vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2SE.perc=abs((temp_bg$E.af-temp_bg$SE.af)/temp_bg$E.af),#perc diff for SE vs E
  		bg.E2SP.perc=abs((temp_bg$E.af-temp_bg$SP.af)/temp_bg$E.af),#perc diff for SP vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_t1_results_EvSP <- rbind(freq_means_t1_results_EvSP,append(temp_means,temp_vars))
  }, error=function(e){})
}

#save table so we don't need to run this loop every time
colnames(freq_means_t1_results_EvSP) <- paste0(colnames(freq_means_t1_results_EvSP), ".EvSP.T1")
write.table(freq_means_t1_results_EvSP, file="rudflies_2023_redo.freq_means_t1_results_EvSP.txt",sep = "\t", quote = FALSE, row.names = F)


## The following loop will go through list of "SE vs E" top outliers, find all matching
## background background SNPs, sample 100 of them, then calculate absolute frequency
## differences for treatment contrasts of interest.

## Instead of running full loop, you can just reload the table from a previous run
freq_means_t1_results_EvPA <- read.table("rudflies_2023_redo.freq_means_t1_results_EvPA.txt", header=TRUE)

### DON'T RUN CODE BLOCK IF YOU JUST LOADED TABLE FROM PREVIOUS RUN ###
freq_means_t1_results_EvPA <- c()#create empty object
## Loop through all cols pval values for candidate genes
for(g in c(1:nrow(freq_means_t1_annot_nr_EvPA))) { #number of sites in EvPA candidate genes
  tryCatch({
  	temp_chrom <- freq_means_t1_annot_nr_EvPA[g,1] #find chrom
  	temp_e <- freq_means_t1_annot_nr_EvPA[g,3] #find E mean frequency
  	temp_pa <- freq_means_t1_annot_nr_EvPA[g,7] #find PA mean frequency
  	temp_cons <- freq_means_t1_annot_nr_EvPA[g,14] #find allele consequence
  	#select bg SNPs by matching chromosome, consequence, and approximate size
  	temp_bg <- freq_means_t1_annot_nr_bg[freq_means_t1_annot_nr_bg$CHROM %in% temp_chrom[1] &
  		freq_means_t1_annot_nr_bg$Consequence %in% temp_cons &
  		freq_means_t1_annot_nr_bg$E.af>(temp_e - (0.25 * min(temp_e,1-temp_e))) &
  		freq_means_t1_annot_nr_bg$E.af<(temp_e + (0.25 * min(temp_e,1-temp_e))),]
  	#randomly sample 100 SNPs, or use all samples if less than 100 matches
  	temp_bg <- temp_bg[sample(nrow(temp_bg), min(c(nrow(temp_bg),100))), ]
  	#calculate bg allele differences
  	temp_af <- cbind(bg.E2SE=abs(temp_bg$E.af-temp_bg$SE.af), #raw diff for SE vs E
  		bg.E2SP=abs(temp_bg$E.af-temp_bg$SP.af), #raw diff for SP vs E
  		bg.E2PA=abs(temp_bg$E.af-temp_bg$PA.af), #raw diff for PA vs E
  		bg.E2SE.perc=abs((temp_bg$E.af-temp_bg$SE.af)/temp_bg$E.af),#perc diff for SE vs E
  		bg.E2SP.perc=abs((temp_bg$E.af-temp_bg$SP.af)/temp_bg$E.af),#perc diff for SP vs E
  		bg.E2PA.perc=abs((temp_bg$E.af-temp_bg$PA.af)/temp_bg$E.af))#perc diff for PA vs E
  	#calculate means for all bg allele differences
  	temp_means <- colMeans(temp_af)
  	names(temp_means) <- paste0(names(temp_means), ".mean")
  	#calculate variances for all bg allele differences
  	temp_vars <- colVars(temp_af)
  	names(temp_vars) <- paste0(names(temp_vars), ".var")
  	#add row to table
  	freq_means_t1_results_EvPA <- rbind(freq_means_t1_results_EvPA,append(temp_means,temp_vars))
  }, error=function(e){})
}

#save table so we don't need to run this loop every time
colnames(freq_means_t1_results_EvPA) <- paste0(colnames(freq_means_t1_results_EvPA), ".EvPA.T1")
write.table(freq_means_t1_results_EvPA, file="rudflies_2023_redo.freq_means_t1_results_EvPA.txt",sep = "\t", quote = FALSE, row.names = F)


### Now that we've calculated background AF difference, time to calculate the same for 
### actual spino candidates and GLM outlier SNPs, then make a master table for each.

## First, make spino candidate table
freq_means_t1_annot_nr_spino_afdiff <- cbind(freq_means_t1_annot_nr_spino, 
	cand.E2SE.spino=abs(freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$SE.af),
	cand.E2SP.spino=abs(freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$SP.af), cand.E2PA.spino = abs(freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$PA.af),
	cand.E2SE.perc.spino=abs((freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$SE.af)/freq_means_t1_annot_nr_spino$E.af),
	cand.E2SP.perc.spino=abs((freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$SP.af)/freq_means_t1_annot_nr_spino$E.af),
	cand.E2PA.perc.spino=abs((freq_means_t1_annot_nr_spino$E.af - freq_means_t1_annot_nr_spino$PA.af)/freq_means_t1_annot_nr_spino$E.af),
	freq_means_t1_results)

## Second, make "E vs SE" GLM outlier table
freq_means_t1_annot_nr_EvSE_afdiff <- cbind(freq_means_t1_annot_nr_EvSE,
	cand.E2SE.EvSE=abs(freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$SE.af),
	cand.E2SP.EvSE=abs(freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$SP.af),
	cand.E2PA.EvSE=abs(freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$PA.af),
	cand.E2SE.perc.EvSE=abs((freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$SE.af)/freq_means_t1_annot_nr_EvSE$E.af),
	cand.E2SP.perc.EvSE=abs((freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$SP.af)/freq_means_t1_annot_nr_EvSE$E.af),
	cand.E2PA.perc.EvSE=abs((freq_means_t1_annot_nr_EvSE$E.af - freq_means_t1_annot_nr_EvSE$PA.af)/freq_means_t1_annot_nr_EvSE$E.af), 
	freq_means_t1_results_EvSE)

## Third, make "E vs SP" GLM outlier table
freq_means_t1_annot_nr_EvSP_afdiff <- cbind(freq_means_t1_annot_nr_EvSP,
	cand.E2SE.EvSP=abs(freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$SE.af),
	cand.E2SP.EvSP=abs(freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$SP.af),
	cand.E2PA.EvSP=abs(freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$PA.af),
	cand.E2SE.perc.EvSP=abs((freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$SE.af)/freq_means_t1_annot_nr_EvSP$E.af),
	cand.E2SP.perc.EvSP=abs((freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$SP.af)/freq_means_t1_annot_nr_EvSP$E.af),
	cand.E2PA.perc.EvSP=abs((freq_means_t1_annot_nr_EvSP$E.af - freq_means_t1_annot_nr_EvSP$PA.af)/freq_means_t1_annot_nr_EvSP$E.af),
	freq_means_t1_results_EvSP)

## Fourth, make "E vs PA" GLM outlier table
freq_means_t1_annot_nr_EvPA_afdiff <- cbind(freq_means_t1_annot_nr_EvPA,
	cand.E2SE.EvPA=abs(freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$SE.af),
	cand.E2SP.EvPA=abs(freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$SP.af),
	cand.E2PA.EvPA=abs(freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$PA.af),
	cand.E2SE.perc.EvPA=abs((freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$SE.af)/freq_means_t1_annot_nr_EvPA$E.af),
	cand.E2SP.perc.EvPA=abs((freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$SP.af)/freq_means_t1_annot_nr_EvPA$E.af),
	cand.E2PA.perc.EvPA=abs((freq_means_t1_annot_nr_EvPA$E.af - freq_means_t1_annot_nr_EvPA$PA.af)/freq_means_t1_annot_nr_EvPA$E.af),
	freq_means_t1_results_EvPA)

## Combine all relevant AF difference stats into a master table
freq_means_t1_annot_nr_afdiff <- cbind(freq_means_t1_annot_nr_spino_afdiff[,c(22:ncol(freq_means_t1_annot_nr_spino_afdiff))], freq_means_t1_annot_nr_EvSE_afdiff[,c(22:ncol(freq_means_t1_annot_nr_EvSE_afdiff))], freq_means_t1_annot_nr_EvPA_afdiff[,c(22:ncol(freq_means_t1_annot_nr_EvPA_afdiff))])

## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs SE differences only
freq_means_t1_annot_nr_E2SEdiff <- data.frame(rbind(cbind("candidate","E2SE","spino","raw",freq_means_t1_annot_nr_afdiff$cand.E2SE.spino),
	   cbind("candidate","E2SE","EvSE","raw",freq_means_t1_annot_nr_afdiff$cand.E2SE.EvSE),
	   cbind("matched","E2SE","spino","mean",freq_means_t1_annot_nr_afdiff$bg.E2SE.mean.spino),
	   cbind("matched","E2SE","EvSE","mean",freq_means_t1_annot_nr_afdiff$bg.E2SE.mean.EvSE),
	   cbind("matched","E2SE","spino","perc",freq_means_t1_annot_nr_afdiff$bg.E2SE.perc.mean.spino),
	   cbind("matched","E2SE","EvSE","perc",freq_means_t1_annot_nr_afdiff$bg.E2SE.perc.mean.EvSE)))

## Rename cols	   
colnames(freq_means_t1_annot_nr_E2SEdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_t1_annot_nr_E2SEdiff[,5] <- as.numeric(freq_means_t1_annot_nr_E2SEdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs SE" 
## outliers), and type (raw, mean, or perc)
freq_means_t1_annot_nr_E2SEdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs SE AF-change in E vs SE candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2SEdiff[freq_means_t1_annot_nr_E2SEdiff$type!="perc" &
	freq_means_t1_annot_nr_E2SEdiff$category=="EvSE",],
	alternative = 'greater')

#T-test for differences in absolute E vs SE AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2SEdiff[freq_means_t1_annot_nr_E2SEdiff$type!="perc" &
	freq_means_t1_annot_nr_E2SEdiff$category=="spino",],
	alternative = 'greater')

### Figure 3 Panel G ###
# make plot of percent AF differences by set and category
E2SEdiff <- ggplot(data = freq_means_t1_annot_nr_E2SEdiff[freq_means_t1_annot_nr_E2SEdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs SE") + 
	theme_classic()

## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs SP differences only
freq_means_t1_annot_nr_E2SPdiff <- data.frame(rbind(cbind("candidate","E2SP","spino","raw", freq_means_t1_annot_nr_afdiff$cand.E2SP.spino),
	   cbind("candidate","E2SP","EvSE","raw", freq_means_t1_annot_nr_afdiff$cand.E2SP.EvSE),
	   cbind("matched","E2SP","spino","mean", freq_means_t1_annot_nr_afdiff$bg.E2SP.mean.spino),
	   cbind("matched","E2SP","EvSE","mean", freq_means_t1_annot_nr_afdiff$bg.E2SP.mean.EvSE),
	   cbind("matched","E2SP","spino","perc", freq_means_t1_annot_nr_afdiff$bg.E2SP.perc.mean.spino),
	   cbind("matched","E2SP","EvSE","perc", freq_means_t1_annot_nr_afdiff$bg.E2SP.perc.mean.EvSE)))
	   
## Rename cols
colnames(freq_means_t1_annot_nr_E2SPdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_t1_annot_nr_E2SPdiff[,5] <- as.numeric(freq_means_t1_annot_nr_E2SPdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs SE" 
## outliers), and type (raw, mean, or perc)
freq_means_t1_annot_nr_E2SPdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs SP AF-change in E vs SE candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2SPdiff[freq_means_t1_annot_nr_E2SPdiff$type!="perc" &
	freq_means_t1_annot_nr_E2SPdiff$category=="EvSE",],
	 alternative = 'greater')

#T-test for differences in absolute E vs SP AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2SPdiff[freq_means_t1_annot_nr_E2SPdiff$type!="perc" &
	freq_means_t1_annot_nr_E2SPdiff$category=="spino",],
	 alternative = 'greater')
	 
### Figure 3 Panel H ###
# make plot of percent AF differences by set and category
E2SPdiff <- ggplot(data = freq_means_t1_annot_nr_E2SPdiff[freq_means_t1_annot_nr_E2SPdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs SP") + 
	theme_classic()
	
## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs PA differences only
freq_means_t1_annot_nr_E2PAdiff <- data.frame(rbind(cbind("candidate","E2PA","spino","raw", freq_means_t1_annot_nr_afdiff$cand.E2PA.spino),
	cbind("candidate","E2PA","EvSE","raw", freq_means_t1_annot_nr_afdiff$cand.E2PA.EvSE),
	cbind("matched","E2PA","spino","mean", freq_means_t1_annot_nr_afdiff$bg.E2PA.mean.spino),
	cbind("matched","E2PA","EvSE","mean", freq_means_t1_annot_nr_afdiff$bg.E2PA.mean.EvSE),
	cbind("matched","E2PA","spino","perc", freq_means_t1_annot_nr_afdiff$bg.E2PA.perc.mean.spino),
	cbind("matched","E2PA","EvSE","perc", freq_means_t1_annot_nr_afdiff$bg.E2PA.perc.mean.EvSE)))
	   
## Rename cols
colnames(freq_means_t1_annot_nr_E2PAdiff) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_t1_annot_nr_E2PAdiff[,5] <- as.numeric(freq_means_t1_annot_nr_E2PAdiff[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs SE" 
## outliers), and type (raw, mean, or perc)
freq_means_t1_annot_nr_E2PAdiff %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs PA AF-change in E vs SE candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2PAdiff[freq_means_t1_annot_nr_E2PAdiff$type!="perc" &
	freq_means_t1_annot_nr_E2PAdiff$category=="EvSE",],
	 alternative = 'greater')

#T-test for differences in absolute E vs PA AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2PAdiff[freq_means_t1_annot_nr_E2PAdiff$type!="perc" &
	freq_means_t1_annot_nr_E2PAdiff$category=="spino",],
	 alternative = 'greater')

### Figure 3 panel I (unused) ###
# make plot of percent AF differences by set and category
E2PAdiff <- ggplot(data = freq_means_t1_annot_nr_E2PAdiff[freq_means_t1_annot_nr_E2PAdiff$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs PA") + 
	theme_classic()


## Save all plots together
pdf(file = "rudflies_2023_redo.freq_means_t1_annot_nr_E2SEdiff.mean.multiboxplot.pdf", width=3, height=6)
	ggarrange(E2SEdiff, E2SPdiff, E2PAdiff,
              ncol = 1, nrow = 3)
dev.off()


### Alternate version of figure 3 panel I with E2PA AF differentiation for EvPA outliers 

## Reorganize master table for plots: all AF diff stats in one column plus metadata cols
## E vs PA differences only
freq_means_t1_annot_nr_E2PAdiff_EvPA <- data.frame(rbind(cbind("candidate","E2PA","spino","raw", freq_means_t1_annot_nr_afdiff$cand.E2PA.spino),
	cbind("candidate","E2PA","EvPA","raw", freq_means_t1_annot_nr_afdiff$cand.E2PA.EvPA),
	cbind("matched","E2PA","spino","mean", freq_means_t1_annot_nr_afdiff$bg.E2PA.mean.spino),
	cbind("matched","E2PA","EvPA","mean", freq_means_t1_annot_nr_afdiff$bg.E2PA.mean.EvPA),
	cbind("matched","E2PA","spino","perc", freq_means_t1_annot_nr_afdiff$bg.E2PA.perc.mean.spino),
	cbind("matched","E2PA","EvPA","perc", freq_means_t1_annot_nr_afdiff$bg.E2PA.perc.mean.EvPA)))
	   
## Rename cols
colnames(freq_means_t1_annot_nr_E2PAdiff_EvPA) <- c("set","contrast","category","type","af_diff")
## Set AF differences to numeric format
freq_means_t1_annot_nr_E2PAdiff_EvPA[,5] <- as.numeric(freq_means_t1_annot_nr_E2PAdiff_EvPA[,5])

##Calculate means of AF differences by set (cand or matched), category (spino or "E vs SE" 
## outliers), and type (raw, mean, or perc)
freq_means_t1_annot_nr_E2PAdiff_EvPA %>%
  group_by(set, category, type) %>%
  summarise(avg = mean(af_diff), med = median(af_diff), stdev = sd(af_diff), .groups="keep")

#T-test for differences in absolute E vs PA AF-change in E vs PA candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2PAdiff_EvPA[freq_means_t1_annot_nr_E2PAdiff_EvPA$type!="perc" &
	freq_means_t1_annot_nr_E2PAdiff_EvPA$category=="EvPA",],
	 alternative = 'greater')

#T-test for differences in absolute E vs PA AF-change in spino candidates vs background
t.test(af_diff ~ set, 
	data = freq_means_t1_annot_nr_E2PAdiff_EvPA[freq_means_t1_annot_nr_E2PAdiff_EvPA$type!="perc" &
	freq_means_t1_annot_nr_E2PAdiff_EvPA$category=="spino",],
	 alternative = 'greater')


# make plot of percent AF differences by set and category
E2PAdiff_EvPA <- ggplot(data = freq_means_t1_annot_nr_E2PAdiff_EvPA[freq_means_t1_annot_nr_E2PAdiff_EvPA$type!="perc",], aes(x=set, y=af_diff)) + 
	geom_boxplot(aes(colour=category)) + 
	scale_color_manual(values=c("#39568A","#C45C5C")) +
	ggtitle("AF differences: E vs PA") + 
	theme_classic()


## Save all plots together
pdf(file = "rudflies_2023_redo.freq_means_t1_annot_nr_E2SEdiff_E2PAdiff.mean.multiboxplot.pdf", width=3, height=4)
	ggarrange(E2SEdiff, E2PAdiff_EvPA,
              ncol = 1, nrow = 2)
dev.off()



### Are there differences in AF changes overall? ###
## Calculate all allele frequency differences, both raw and percent changes
freq_means_t1_afdiff <- as.data.frame(cbind(E2S=abs(freq_means_t1$E.af-freq_means_t1$S.af),
	E2SE=abs(freq_means_t1$E.af-freq_means_t1$SE.af),
	E2SP=abs(freq_means_t1$E.af-freq_means_t1$SP.af),
	E2PA=abs(freq_means_t1$E.af-freq_means_t1$PA.af),
	E2S.perc=abs((freq_means_t1$E.af - freq_means_t1$S.af)/freq_means_t1$E.af),
	E2SE.perc=abs((freq_means_t1$E.af - freq_means_t1$SE.af)/freq_means_t1$E.af),
	E2SP.perc=abs((freq_means_t1$E.af - freq_means_t1$SP.af)/freq_means_t1$E.af),
	E2PA.perc=abs((freq_means_t1$E.af - freq_means_t1$PA.af)/freq_means_t1$E.af)))

## Reformat for running T-tests
freq_means_t1_afdiff <- data.frame(rbind(cbind(contrast="E2S",
		af_diff=freq_means_t1_afdiff$E2S,
		perc_diff=freq_means_t1_afdiff$E2S.perc),
	cbind(contrast="E2SE",
		af_diff=freq_means_t1_afdiff$E2SE,
		perc_diff=freq_means_t1_afdiff$E2SE.perc),
	cbind(contrast="E2SP",
		af_diff=freq_means_t1_afdiff$E2SP,
		perc_diff=freq_means_t1_afdiff$E2SP.perc),
	cbind(contrast="E2PA",
		af_diff=freq_means_t1_afdiff$E2PA,
		perc_diff=freq_means_t1_afdiff$E2PA.perc)))
		
## Set AF differences to numeric format
freq_means_t1_afdiff[,2] <- as.numeric(freq_means_t1_afdiff[,2])
freq_means_t1_afdiff[,3] <- as.numeric(freq_means_t1_afdiff[,3])

	   
## T-test for differences in all pairwise absolute AF-changes
pairwise.t.test(freq_means_t1_afdiff$af_diff, 
	freq_means_t1_afdiff$contrast,
	alternative = 'two.sided')	
#	Pairwise comparisons using t tests with pooled SD 
#
#data:  freq_means_t1_afdiff$af_diff and freq_means_t1_afdiff$contrast 
#
#     E2PA   E2S    E2SE  
#E2S  <2e-16 -      -     
#E2SE <2e-16 <2e-16 -     
#E2SP <2e-16 <2e-16 <2e-16

## T-test for differences in all pairwise percent AF-changes
pairwise.t.test(freq_means_t1_afdiff$perc_diff, 
	freq_means_t1_afdiff$contrast,
	alternative = 'two.sided')
#	Pairwise comparisons using t tests with pooled SD 
#
#data:  freq_means_t1_afdiff$perc_diff and freq_means_t1_afdiff$contrast 
#
#     E2PA    E2S     E2SE   
#E2S  < 2e-16 -       -      
#E2SE < 2e-16 < 2e-16 -      
#E2SP < 2e-16 < 2e-16 2e-16

## Find group means and standard deviations
freq_means_t1_afdiff %>%
  group_by(contrast) %>%
  summarise(raw_avg = mean(af_diff), raw_stdev = sd(af_diff), perc_avg = mean(perc_diff), perc_stdev = sd(perc_diff), .groups="keep")
## A tibble: 4 × 5
## Groups:   contrast [4]
#  contrast raw_avg raw_stdev perc_avg perc_stdev
#  <chr>      <dbl>     <dbl>    <dbl>      <dbl>
#1 E2PA      0.0521    0.0524    0.372      1.04 
#2 E2S       0.0204    0.0178    0.118      0.167
#3 E2SE      0.0279    0.0279    0.153      0.211
#4 E2SP      0.0238    0.0194    0.140      0.214




###########################################################
### Non-parametric SNP rank correlations and RRHO plots ###
###########################################################

## This will allow us to visualize convergence between treatments via non-parametic
## rank correlations of different sets of GLM p-values. More specifically, want to test
## for convergence of PA with either SE or SP, relative to control E populations, at TPT1.
## This will be done by taking -log10(p) values for all SNPs the "E vs PA", "E vs SE", 
## and "E vs SP" glm contrasts, making them positive or negative according to the mean 
## change in allele frequency relative to control E populations, then running pairwise
## correlations between contrasts to check for concordance in SNP ranks. 

## These ranks will be closer to 1 if the AF changes were highly significant and changing
## in the positive direction in spinosad-exposed populations relative to controls. 
## Alternatively, ranks will be closer to N (# of SNPs) if AF changes were highly 
## significant and changing in the negative direction. 

## The plots will show areas of density in the quadrants where the SNP ranks in pairwise 
## contrasts are most clustered. If the clusters appear in the lower left or upper right ## quadrants, the contrasts are positively correlated, which suggests convergence in the
## exposed populations. If the clusters appear in the upper left or lower right 
## quadrants, the contrasts are negatively correlated. This may be biologically 
## interesting, but the patterns are discordant. Finally, if there are little or no areas
## of density in the plot, the glm contrasts are uncorrelated, suggesting a lack of 
## parallel adaptation to spinosad, relative to control E populations.

## Overall, the RRHO plots are visualization methods that complement non-parametric
## Spearman Rank correlation analysis, which we will also do in this script.


#Try this color pallette for RRHO plitting
# Define a vector of colors (using names or hex codes)
rrho_cols <- c("#E8C863","#DD9551","#D24D4D","#C03268","#684187","#34215F") #light to dark

# Create a palette function that interpolates between these colors
# The resulting 'color_palette' is a function that takes an integer 'n'
# and returns 'n' interpolated colors.
rrho_palette <- colorRampPalette(rrho_cols)

# Generate 200 colors from the continuous palette
rrho_shades <- rrho_palette(200)


## Let's calculate mean frequencies, then calculate pairwise differences
## We'll use these values to assign positive or negative signs to -log10p values
## The resulting values will reflect both the magnitude and direction of AF change
freq_diff <- cbind(EvSE.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]),
	EvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]),
	EvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="E"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SEvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SPvPA.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="PA"]),
	SEvSP.T1.diff=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SE"]) - rowMeans(haf.freq.T1filt[,haf.meta.T1filt$condition=="SP"]))

## For use in bedtools, we'll create a bed file from AF differences
freq_diff_bed <- cbind(haf.sites.T1filt,STOP=haf.sites.T1filt$POS+1, freq_diff)
freq_diff_bed <- freq_diff_bed[order(freq_diff_bed[,1], freq_diff_bed[,2]), ]
## save
options(scipen=20)
write.table(freq_diff_bed, file="rudflies_2023_redo.freq_diff.T1.bed",sep = "\t", quote = FALSE, row.names = F)
options(scipen=0)

## Join TPT1 GLM results with new freq_diff_bed table
freq_diff_PAvSvSEvE <- merge(contrast.PAvSvSEvE.table, freq_diff_bed[,-3], by=c("CHROM","POS"))

## And sort by locus
freq_diff_PAvSvSEvE <- freq_diff_PAvSvSEvE[order(freq_diff_PAvSvSEvE[,1], freq_diff_PAvSvSEvE[,2]), ]

## Extra filter to ensure we have no rows where p=0 of our contrasts of interest
## These 0 values will result in infinite -log10(p) values, which RRHOs can't handle
freq_diff_PAvSvSEvE <- freq_diff_PAvSvSEvE[rowProds(as.matrix(freq_diff_PAvSvSEvE[,c(3:5,8)]))!=0,]

## Now calculate AF means across treatments for additional filtering of table
haf.freq.T1filt.af_mean <- cbind(haf.sites.T1filt,af_mean=rowMeans(haf.freq.T1filt[,haf.meta.T1filt$treat.fix=="E" | haf.meta.T1filt$treat.fix=="S" | haf.meta.T1filt$treat.fix=="PA"]))
## Merge new means with GLM/freq_diff table
freq_diff_PAvSvSEvE <- merge(freq_diff_PAvSvSEvE, haf.freq.T1filt.af_mean, by=c("CHROM","POS"))

## Filter if AF means are less than 0.15 or greater than 0.85.
## This preserves ranking patterns for biologically relevant loci while removing those 
## that appear to have major treatment differences due to REF or ALT alleles being rare.
freq_diff_PAvSvSEvE <- freq_diff_PAvSvSEvE[freq_diff_PAvSvSEvE$af_mean>0.15 & freq_diff_PAvSvSEvE$af_mean<0.85,]

## Let create locus names to cross reference rankings across contrasts
RRHO.SNPs <- paste(freq_diff_PAvSvSEvE$CHROM, freq_diff_PAvSvSEvE$POS, sep = "_")

### Create a table for each contrast to use in correlations and RRHO plotting
### Here, we finally calculate signed -log10 p-values

## EvSE
EvSE.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs,logp=-log10(freq_diff_PAvSvSEvE$EvSE.T1)*(abs(freq_diff_PAvSvSEvE$EvSE.T1.diff)/freq_diff_PAvSvSEvE$EvSE.T1.diff)))
# reformat
EvSE.T1.signed.logp.snp.all$logp <- as.numeric(EvSE.T1.signed.logp.snp.all$logp)

## EvSP
EvSP.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs,logp=-log10(freq_diff_PAvSvSEvE$EvSP.T1)*(abs(freq_diff_PAvSvSEvE$EvSP.T1.diff)/freq_diff_PAvSvSEvE$EvSP.T1.diff)))
# reformat
EvSP.T1.signed.logp.snp.all$logp <- as.numeric(EvSP.T1.signed.logp.snp.all$logp)

## SPvSE
SPvSE.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs,logp=-log10(freq_diff_PAvSvSEvE$SEvSP.T1)*-(abs(freq_diff_PAvSvSEvE$SEvSP.T1.diff)/freq_diff_PAvSvSEvE$SEvSP.T1.diff)))
# reformat
SPvSE.T1.signed.logp.snp.all$logp <- as.numeric(SPvSE.T1.signed.logp.snp.all$logp)

## EvPA
EvPA.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs,logp=-log10(freq_diff_PAvSvSEvE$EvPA.T1)*(abs(freq_diff_PAvSvSEvE$EvPA.T1.diff)/freq_diff_PAvSvSEvE$EvPA.T1.diff)))
# reformat
EvPA.T1.signed.logp.snp.all$logp <- as.numeric(EvPA.T1.signed.logp.snp.all$logp)


### Spearman rank correlations
## Non-parametric rank correlation between EvSE.T1 and EvPA.T1
cor.test(as.numeric(EvSE.T1.signed.logp.snp.all[,2]), as.numeric(EvPA.T1.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvSE.T1.signed.logp.snp.all[, 2]) and as.numeric(EvPA.T1.signed.logp.snp.all[, 2])
#S = 8.5486e+16, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.5642118 

## Non-parametric rank correlation between EvSP.T1 and EvPA.T1
cor.test(as.numeric(EvSP.T1.signed.logp.snp.all[,2]), as.numeric(EvPA.T1.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvSP.T1.signed.logp.snp.all[, 2]) and as.numeric(EvPA.T1.signed.logp.snp.all[, 2])
#S = 1.6789e+17, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1441586 

## Non-parametric rank correlation between SPvSE.T1 and EvPA.T1
cor.test(as.numeric(SPvSE.T1.signed.logp.snp.all[,2]), as.numeric(EvPA.T1.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(SPvSE.T1.signed.logp.snp.all[, 2]) and as.numeric(EvPA.T1.signed.logp.snp.all[, 2])
#S = 1.2255e+17, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.3752749 

## Non-parametric rank correlation between EvSE.T1 and EvSP.T1
cor.test(as.numeric(EvSE.T1.signed.logp.snp.all[,2]), as.numeric(EvSP.T1.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvSE.T1.signed.logp.snp.all[, 2]) and as.numeric(EvSP.T1.signed.logp.snp.all[, 2])
#S = 1.6217e+17, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1732958 


### RRHO plotting ###

## Since the RRHO package is often used for gene expression studies with thousands of 
## genes, processing for millions of SNPs takes excessive computation time. For 
## simplicity, we sampled 25K SNPs multiple times and checked whether each RRHO arrived at 
## a consensus plot. For the record, they did. 

## Randomly sample master table
freq_diff_PAvSvSEvE_samp <- freq_diff_PAvSvSEvE[sample(nrow(freq_diff_PAvSvSEvE),25000), ]

## Let create locus names to cross reference rankings across contrasts
RRHO.SNPs.samp <- paste(freq_diff_PAvSvSEvE_samp$CHROM, freq_diff_PAvSvSEvE_samp$POS, sep = "_")

### Create a table for each contrast to use in correlations and RRHO plotting
### Here, we finally calculate signed -log10 p-values

## EvSE
EvSE.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.samp,logp=-log10(freq_diff_PAvSvSEvE_samp$EvSE.T1)*(abs(freq_diff_PAvSvSEvE_samp$EvSE.T1.diff)/freq_diff_PAvSvSEvE_samp$EvSE.T1.diff)))
# reformat
EvSE.T1.signed.logp.snp$logp <- as.numeric(EvSE.T1.signed.logp.snp$logp)

## EvSP
EvSP.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.samp,logp=-log10(freq_diff_PAvSvSEvE_samp$EvSP.T1)*(abs(freq_diff_PAvSvSEvE_samp$EvSP.T1.diff)/freq_diff_PAvSvSEvE_samp$EvSP.T1.diff)))
# reformat
EvSP.T1.signed.logp.snp$logp <- as.numeric(EvSP.T1.signed.logp.snp$logp)

## SPvSE
SPvSE.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.samp,logp=-log10(freq_diff_PAvSvSEvE_samp$SEvSP.T1)*-(abs(freq_diff_PAvSvSEvE_samp$SEvSP.T1.diff)/freq_diff_PAvSvSEvE_samp$SEvSP.T1.diff)))
# reformat
SPvSE.T1.signed.logp.snp$logp <- as.numeric(SPvSE.T1.signed.logp.snp$logp)

## EvPA
EvPA.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.samp,logp=-log10(freq_diff_PAvSvSEvE_samp$EvPA.T1)*(abs(freq_diff_PAvSvSEvE_samp$EvPA.T1.diff)/freq_diff_PAvSvSEvE_samp$EvPA.T1.diff)))
# reformat
EvPA.T1.signed.logp.snp$logp <- as.numeric(EvPA.T1.signed.logp.snp$logp)


### Finally, we can make RRHOs ###
## On first pass, we found the highest concordance when comparing "EvSE" and "EvPA",
## therefore we create this RRHO object first and normalize the heatmap color scale 
## maximums in all other RRHOs realative to this one. We also set the color gradient 
## to a custom one (rrho_shades).

### Figure 3 panel J
## Create RRHO object for comparison of "EvSE" and "EvPA" GLM contrasts
rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp <- RRHO2_initialize(EvSE.T1.signed.logp.snp, EvPA.T1.signed.logp.snp, labels = c("EvSE.T1", "EvPA.T1"), log10.ind=TRUE)
#save plot
pdf(file = "rudflies_2023_redo.EvSE.T1.EvPA.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()

### Figure 3 panel K
## Create RRHO object for comparison of "EvSP" and "EvPA" GLM contrasts
rrho2.EvSP.T1.EvPA.T1.signed.logp.25Ksnp <- RRHO2_initialize(EvSP.T1.signed.logp.snp, EvPA.T1.signed.logp.snp, labels = c("EvSP.T1", "EvPA.T1"), log10.ind=TRUE)
#save plot
pdf(file = "rudflies_2023_redo.EvSP.T1.EvPA.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvSP.T1.EvPA.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades) #use max from Fig. 3J 
dev.off()

### Not used
## Create RRHO object for comparison of "SPvSE" and "EvPA" GLM contrasts
rrho2.SPvSE.T1.EvPA.T1.signed.logp.25Ksnp <- RRHO2_initialize(SPvSE.T1.signed.logp.snp, EvPA.T1.signed.logp.snp, labels = c("SPvSE.T1", "EvPA.T1"), log10.ind=TRUE)
#save plot
pdf(file = "rudflies_2023_redo.SPvSE.T1.EvPA.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.SPvSE.T1.EvPA.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades) #use max from Fig. 3J 
dev.off()

### Figure S11 panel B
## Create RRHO object for comparison of "EvSE" and "EvSP" GLM contrasts
rrho2.EvSE.T1.EvSP.T1.signed.logp.25Ksnp <- RRHO2_initialize(EvSE.T1.signed.logp.snp, EvSP.T1.signed.logp.snp, labels = c("EvSE.T1", "EvSP.T1"), log10.ind=TRUE)
#save plot
pdf(file = "rudflies_2023_redo.EvSE.T1.EvSP.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvSE.T1.EvSP.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvSE.T1.EvPA.T1.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades) #use max from Fig. 3J 
dev.off()


####################################
### Bonus time-series RRHO plots ###
####################################

## Here we want to look at convergence of exposed S and PA pops at each TPT to visualize 
## potential temporal convergence with pre-adapted populations.

## Create a table that has AF means across all E, S, and PA pops, basically all remaining 
## samples but redundant filter is there for safety.
haf.freq.filt.af_mean <- cbind(haf.sites.filt,af_mean=rowMeans(haf.freq.filt[,haf.meta.filt$treat.fix=="E" | haf.meta.filt$treat.fix=="S" | haf.meta.filt$treat.fix=="PA"]))

## Now, filter for E vs PA: all timepoints
haf.meta.filt.PA.E <- haf.meta[haf.meta$batch == "a" & (haf.meta$treat.fix == "E" | haf.meta$treat.fix == "PA"),]
## Secondary filter to remove any heat-experiment related samples
haf.meta.filt.PA.E <- haf.meta.filt.PA.E[which((haf.meta.filt.PA.E$samp %in% c(1:9,30:171))==TRUE),]
## Match the frequency table samples to the filtered metadata samples
haf.freq.filt.PA.E <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt.PA.E$samp)==TRUE)]
## For safety, match the metadata samples to the frequency table samples (redundant)
haf.meta.filt.PA.E <- haf.meta.filt.PA.E[which((haf.meta.filt.PA.E$samp %in% names(haf.freq.filt.PA.E))==TRUE),]

## Additional filter to remove extreme low variance loci
haf.freq.filt.PA.E.loci <- haf.freq[rowVars(as.matrix(haf.freq.filt.PA.E))>0.001,c(1:2)]
haf.freq.filt.PA.E <- haf.freq.filt.PA.E[rowVars(as.matrix(haf.freq.filt.PA.E))>0.001,]

## Tables need a bit of reformatting for glm to run
haf.meta.filt.PA.E$treat.fix <- as.factor(haf.meta.filt.PA.E$treat.fix)
haf.meta.filt.PA.E$tpt <- as.factor(haf.meta.filt.PA.E$tpt)

## Create table with E vs PA allele frequency differences for all time-points combined 
## and at each individual time-point.
EvPA_allTPT_diff <- cbind(haf.freq.filt.PA.E.loci,
	EvPA.diff=rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="E"]) - rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="PA"]), #all TPT
	EvPA.T1.diff=rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="E" & haf.meta.filt.PA.E$tpt=="1"]) - rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="PA" & haf.meta.filt.PA.E$tpt=="1"]), #TPT1
	EvPA.T2.diff=rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="E" & haf.meta.filt.PA.E$tpt=="2"]) - rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="PA" & haf.meta.filt.PA.E$tpt=="2"]), #TPT2
	EvPA.T3.diff=rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="E" & haf.meta.filt.PA.E$tpt=="3"]) - rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="PA" & haf.meta.filt.PA.E$tpt=="3"]), #TPT3
	EvPA.T4.diff=rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="E" & haf.meta.filt.PA.E$tpt=="4"]) - rowMeans(haf.freq.filt.PA.E[,haf.meta.filt.PA.E$treat.fix=="PA" & haf.meta.filt.PA.E$tpt=="4"])) #TPT4


## Now, filter for E vs S: all timepoints
haf.meta.filt.S.E <- haf.meta[haf.meta$batch == "a" & (haf.meta$treat.fix == "E" | haf.meta$treat.fix == "S"),]
## Secondary filter to remove any heat-experiment related samples
haf.meta.filt.S.E <- haf.meta.filt.S.E[which((haf.meta.filt.S.E$samp %in% c(1:9,30:171))==TRUE),]
## Match the frequency table samples to the filtered metadata samples
haf.freq.filt.S.E <- haf.freq[, which((names(haf.freq) %in% haf.meta.filt.S.E$samp)==TRUE)]
## For safety, match the metadata samples to the frequency table samples (redundant)
haf.meta.filt.S.E <- haf.meta.filt.S.E[which((haf.meta.filt.S.E$samp %in% names(haf.freq.filt.S.E))==TRUE),]

## Additional filter to remove extreme low variance loci
haf.freq.filt.S.E.loci <- haf.freq[rowVars(as.matrix(haf.freq.filt.S.E))>0.001,c(1:2)]
haf.freq.filt.S.E <- haf.freq.filt.S.E[rowVars(as.matrix(haf.freq.filt.S.E))>0.001,]

## Tables need a bit of reformatting for glm to run
haf.meta.filt.S.E$treat.fix <- as.factor(haf.meta.filt.S.E$treat.fix)
haf.meta.filt.S.E$tpt <- as.factor(haf.meta.filt.S.E$tpt)

## Create table with E vs S allele frequency differences for all time-points combined 
## and at each individual time-point.
EvS_allTPT_diff <- cbind(haf.freq.filt.S.E.loci,
	EvS.diff=rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="E"]) - rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="S"]), #all TPT
	EvS.T1.diff=rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="E" & haf.meta.filt.S.E$tpt=="1"]) - rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="S" & haf.meta.filt.S.E$tpt=="1"]), #TPT1
	EvS.T2.diff=rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="E" & haf.meta.filt.S.E$tpt=="2"]) - rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="S" & haf.meta.filt.S.E$tpt=="2"]), #TPT2
	EvS.T3.diff=rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="E" & haf.meta.filt.S.E$tpt=="3"]) - rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="S" & haf.meta.filt.S.E$tpt=="3"]), #TPT3
	EvS.T4.diff=rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="E" & haf.meta.filt.S.E$tpt=="4"]) - rowMeans(haf.freq.filt.S.E[,haf.meta.filt.S.E$treat.fix=="S" & haf.meta.filt.S.E$tpt=="4"])) #TPT4


### Now we grab p-values for EvS and EvPA contrasts and combine with AF difference tables
## E vs S contrasts
freq_diff_EvS <- merge(contrastout.SvE.table, EvS_allTPT_diff, by=c("CHROM","POS"))
## E vs PA contrasts
freq_diff_EvPA <- merge(contrastout.PAvE.table, EvPA_allTPT_diff, by=c("CHROM","POS"))
## Merge both of the above tables
freq_diff_EvS_EvPA <- na.omit(merge(freq_diff_EvS, freq_diff_EvPA, by=c("CHROM","POS")))
## Finally, merge with all-sample AF mean table for filtering
freq_diff_EvS_EvPA <- merge(freq_diff_EvS_EvPA, haf.freq.filt.af_mean, by=c("CHROM","POS"))

## Filter for intermediate mean allele frequencies since rare alleles may have inflated
## significance levels from minor AF changes without necessarily relating to major 
## phenotypic impacts.
freq_diff_EvS_EvPA <- freq_diff_EvS_EvPA[freq_diff_EvS_EvPA$af_mean>0.15 & freq_diff_EvS_EvPA$af_mean<0.85,]

## Create a list of loci for correlation tables
RRHO.SNPs.tpt <- paste(freq_diff_EvS_EvPA$CHROM, freq_diff_EvS_EvPA$POS, sep = "_")

### Now build sign-corrected -log10(p) value tables. Strategy is multiplying -log10(p) 
### values by positive or negative 1 depending on the direction of AF differences between 
### control and exposed populations. This will result in more accurate rankings based on 
### both significance levels and direction.

## EvS all TPT
EvS.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$SvE)*(abs(freq_diff_EvS_EvPA$EvS.diff)/freq_diff_EvS_EvPA$EvS.diff)))

## EvS.T1
EvS.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$SvE.T1)*(abs(freq_diff_EvS_EvPA$EvS.T1.diff)/freq_diff_EvS_EvPA$EvS.T1.diff)))

## EvS.T2
EvS.T2.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$SvE.T2)*(abs(freq_diff_EvS_EvPA$EvS.T2.diff)/freq_diff_EvS_EvPA$EvS.T2.diff)))

## EvS.T3
EvS.T3.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$SvE.T3)*(abs(freq_diff_EvS_EvPA$EvS.T3.diff)/freq_diff_EvS_EvPA$EvS.T3.diff)))

## EvS.T4
EvS.T4.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$SvE.T4)*(abs(freq_diff_EvS_EvPA$EvS.T4.diff)/freq_diff_EvS_EvPA$EvS.T4.diff)))

## EvPA all TPT
EvPA.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$PAvE)*(abs(freq_diff_EvS_EvPA$EvPA.diff)/freq_diff_EvS_EvPA$EvPA.diff)))

## EvPA.T1
EvPA.T1.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$PAvE.T1)*(abs(freq_diff_EvS_EvPA$EvPA.T1.diff)/freq_diff_EvS_EvPA$EvPA.T1.diff)))

## EvPA.T2
EvPA.T2.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$PAvE.T2)*(abs(freq_diff_EvS_EvPA$EvPA.T2.diff)/freq_diff_EvS_EvPA$EvPA.T2.diff)))

## EvPA.T3
EvPA.T3.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$PAvE.T3)*(abs(freq_diff_EvS_EvPA$EvPA.T3.diff)/freq_diff_EvS_EvPA$EvPA.T3.diff)))

## EvPA.T4
EvPA.T4.signed.logp.snp.all <- data.frame(cbind(snps=RRHO.SNPs.tpt,logp=-log10(freq_diff_EvS_EvPA$PAvE.T4)*(abs(freq_diff_EvS_EvPA$EvPA.T4.diff)/freq_diff_EvS_EvPA$EvPA.T4.diff)))


### Run Spearman Rank correlations between sign-corrected -log10(p) values using all SNPs
## Non-parametric rank correlation between EvS and EvPA all TPTs
cor.test(as.numeric(EvS.signed.logp.snp.all[,2]), as.numeric(EvPA.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvS.signed.logp.snp.all[, 2]) and as.numeric(EvPA.signed.logp.snp.all[, 2])
#S = 6.9406e+16, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#     rho 
#0.639967 

## Non-parametric rank correlation between EvS.T1 and EvPA.T1
cor.test(as.numeric(EvS.T1.signed.logp.snp.all[,2]), as.numeric(EvPA.T1.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvS.T1.signed.logp.snp.all[, 2]) and as.numeric(EvPA.T1.signed.logp.snp.all[, 2])
#S = 9.722e+16, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.4956788 

## Non-parametric rank correlation between EvS.T2 and EvPA.T2
cor.test(as.numeric(EvS.T2.signed.logp.snp.all[,2]), as.numeric(EvPA.T2.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvS.T2.signed.logp.snp.all[, 2]) and as.numeric(EvPA.T2.signed.logp.snp.all[, 2])
#S = 8.0201e+16, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.5839663 

## Non-parametric rank correlation between EvS.T3 and EvPA.T3
cor.test(as.numeric(EvS.T3.signed.logp.snp.all[,2]), as.numeric(EvPA.T3.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvS.T3.signed.logp.snp.all[, 2]) and as.numeric(EvPA.T3.signed.logp.snp.all[, 2])
#S = 6.2047e+16, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.6781388 

## Non-parametric rank correlation between EvS.T4 and EvPA.T4
cor.test(as.numeric(EvS.T4.signed.logp.snp.all[,2]), as.numeric(EvPA.T4.signed.logp.snp.all[,2]), method="spearman")
#	Spearman's rank correlation rho
#
#data:  as.numeric(EvS.T4.signed.logp.snp.all[, 2]) and as.numeric(EvPA.T4.signed.logp.snp.all[, 2])
#S = 1.1507e+17, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.4030992 


### For creating RRHOs, we'll sample 25K SNPs. Since they are most efficient for 
### transcriptomic datasets with several thousand genes, not 1M+, runtime is exponential 
### for our locus count. To account for variation in sampling, we ran this script at least 
### 10 times to ensure consistency of patterns.  

## Sample 25K random SNPs
freq_diff_EvS_EvPA_samp <- freq_diff_EvS_EvPA[sample(nrow(freq_diff_EvS_EvPA), 25000), ]

## Create a list of loci for RRHO tables
RRHO.SNPs.tpt.samp <- paste(freq_diff_EvS_EvPA_samp$CHROM, freq_diff_EvS_EvPA_samp$POS, sep = "_")

### Calculate sign-corrected -log10(p) values, same as above
## EvS all TPT
EvS.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=as.numeric(-log10(freq_diff_EvS_EvPA_samp$SvE)*(abs(freq_diff_EvS_EvPA_samp$EvS.diff)/freq_diff_EvS_EvPA_samp$EvS.diff))))
# reformat
EvS.signed.logp.snp$logp <- as.numeric(EvS.signed.logp.snp$logp)

## EvS.T1
EvS.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=as.numeric(-log10(freq_diff_EvS_EvPA_samp$SvE.T1)*(abs(freq_diff_EvS_EvPA_samp$EvS.T1.diff)/freq_diff_EvS_EvPA_samp$EvS.T1.diff))))
# reformat
EvS.T1.signed.logp.snp$logp <- as.numeric(EvS.T1.signed.logp.snp$logp)

## EvS.T2
EvS.T2.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=as.numeric(-log10(freq_diff_EvS_EvPA_samp$SvE.T2)*(abs(freq_diff_EvS_EvPA_samp$EvS.T2.diff)/freq_diff_EvS_EvPA_samp$EvS.T2.diff))))
# reformat
EvS.T2.signed.logp.snp$logp <- as.numeric(EvS.T2.signed.logp.snp$logp)

## EvS.T3
EvS.T3.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=as.numeric(-log10(freq_diff_EvS_EvPA_samp$SvE.T3)*(abs(freq_diff_EvS_EvPA_samp$EvS.T3.diff)/freq_diff_EvS_EvPA_samp$EvS.T3.diff))))
# reformat
EvS.T3.signed.logp.snp$logp <- as.numeric(EvS.T3.signed.logp.snp$logp)

## EvS.T4
EvS.T4.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=as.numeric(-log10(freq_diff_EvS_EvPA_samp$SvE.T4)*(abs(freq_diff_EvS_EvPA_samp$EvS.T4.diff)/freq_diff_EvS_EvPA_samp$EvS.T4.diff))))
# reformat
EvS.T4.signed.logp.snp$logp <- as.numeric(EvS.T4.signed.logp.snp$logp)

## EvPA all TPT
EvPA.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=as.numeric(-log10(freq_diff_EvS_EvPA_samp$PAvE)*(abs(freq_diff_EvS_EvPA_samp$EvPA.diff)/freq_diff_EvS_EvPA_samp$EvPA.diff))))
# reformat
EvPA.signed.logp.snp$logp <- as.numeric(EvPA.signed.logp.snp$logp)

## EvPA.T1
EvPA.T1.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=-log10(freq_diff_EvS_EvPA_samp$PAvE.T1)*(abs(freq_diff_EvS_EvPA_samp$EvPA.T1.diff)/freq_diff_EvS_EvPA_samp$EvPA.T1.diff)))
# reformat
EvPA.T1.signed.logp.snp$logp <- as.numeric(EvPA.T1.signed.logp.snp$logp)

## EvPA.T2
EvPA.T2.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=-log10(freq_diff_EvS_EvPA_samp$PAvE.T2)*(abs(freq_diff_EvS_EvPA_samp$EvPA.T2.diff)/freq_diff_EvS_EvPA_samp$EvPA.T2.diff)))
# reformat
EvPA.T2.signed.logp.snp$logp <- as.numeric(EvPA.T2.signed.logp.snp$logp)

## EvPA.T3
EvPA.T3.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=-log10(freq_diff_EvS_EvPA_samp$PAvE.T3)*(abs(freq_diff_EvS_EvPA_samp$EvPA.T3.diff)/freq_diff_EvS_EvPA_samp$EvPA.T3.diff)))
# reformat
EvPA.T3.signed.logp.snp$logp <- as.numeric(EvPA.T3.signed.logp.snp$logp)

## EvPA.T4
EvPA.T4.signed.logp.snp <- data.frame(cbind(snps=RRHO.SNPs.tpt.samp,logp=-log10(freq_diff_EvS_EvPA_samp$PAvE.T4)*(abs(freq_diff_EvS_EvPA_samp$EvPA.T4.diff)/freq_diff_EvS_EvPA_samp$EvPA.T4.diff)))
# reformat
EvPA.T4.signed.logp.snp$logp <- as.numeric(EvPA.T4.signed.logp.snp$logp)

### Build RRHO objects and plots to visualize genome-wide contrast of contrasts
## Start with TPT3 since Spearman Rank correlations revealed the strongest relationships.
## After building RRHO object, we'll use to max values from this plot to scale others. 
rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp <- RRHO2_initialize(EvS.T3.signed.logp.snp, EvPA.T3.signed.logp.snp, labels = c("EvS.T3", "EvPA.T3"), log10.ind=TRUE)

## Now plot
pdf(file = "rudflies_2023_redo.EvS.T3.EvPA.T3.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()

## Build RRHO object for TPT1
rrho2.EvS.T1.EvPA.T1.signed.logp.25Ksnp <- RRHO2_initialize(EvS.T1.signed.logp.snp, EvPA.T1.signed.logp.snp, labels = c("EvS.T1", "EvPA.T1"), log10.ind=TRUE)

## Now plot
pdf(file = "rudflies_2023_redo.EvS.T1.EvPA.T1.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvS.T1.EvPA.T1.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()

## Build RRHO object for TPT2
rrho2.EvS.T2.EvPA.T2.signed.logp.25Ksnp <- RRHO2_initialize(EvS.T2.signed.logp.snp, EvPA.T2.signed.logp.snp, labels = c("EvS.T2", "EvPA.T2"), log10.ind=TRUE)

## Now plot
pdf(file = "rudflies_2023_redo.EvS.T2.EvPA.T2.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvS.T2.EvPA.T2.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()

## Build RRHO object for TPT4
rrho2.EvS.T4.EvPA.T4.signed.logp.25Ksnp <- RRHO2_initialize(EvS.T4.signed.logp.snp, EvPA.T4.signed.logp.snp, labels = c("EvS.T4", "EvPA.T4"), log10.ind=TRUE)

## Now plot
pdf(file = "rudflies_2023_redo.EvS.T4.EvPA.T4.signed.logp.25Ksnp.RRHO2.pdf", width=6, height=5.5)
	RRHO2_heatmap(rrho2.EvS.T4.EvPA.T4.signed.logp.25Ksnp, maximum = max(na.omit(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp$hypermat[is.finite(rrho2.EvS.T3.EvPA.T3.signed.logp.25Ksnp$hypermat)])), colorGradient = rrho_shades)
dev.off()


################
### Figure 4 ###
################

## For improved visualization, we aimed to mask genomic regions in Manhattan plots that
## do not appear to be convergent selective targets across spinosad-exposed treatments 
## (S and PA), relative to E control populations. To do so, we leveraged our experimental
## design to mask -log10(p) values in target GLM contrasts ("S vs E" and "PA vs E") 
## by subracting maximum-normalized values from the "PA vs S" contrast.

## The concept relies on hypothesis that PA and S genomes show major differences 
## since PA populations were split off from the E control populations prior to S and
## were exposed to spinosad selection almost immediately. Clustering analyses
## and Manhattan plots all confirm that PA is more divergent from both E and S pops
## than E and S populations are from each other. 

## Due to strong selection from spinosad, PA and S populations show putative strong 
## sweeps in similar regions relative to control populations, however genomic draft 
## obscures the actual loci subject to selection. PA and S pops should be most similar 
## at adaptive loci though while drift and draft are likely to have more disparate 
## impacts in non-target loci, resulting in more variable signals.  

## Using this information, we subtract the maximum-normalized "PA vs S" -log10(p) values
## from "S vs E" and "PA vs E" values. In effect, convergent adaptive loci should have 
## low -log10(p) values in "PA vs S" and high -log10(p) values in contrasts that include E 
## controls, therefore those regions will not be masked by subtraction and will remain 
## peaks. Conversely, regions where -log10(p) values are moderate to high in 
## "PA vs S" are not likely subject to convergent parallel adaptation, and subtraction 
## of these values will mask these non-target regions in "S vs E" and "PA vs E" GLM 
## contrasts. 



### 3-way GLM comparison
## To look for top empirical candidate loci subject to convergent selection between
## PA and S populations, relative to E, we looked for parallel adaptation in SE and PA
## in TPT1 and SP and PA in TPT4, then overlap in both those treatment/time-point combos.

## Join all relevant GLM results into a single filtering table
glm.all.rolwin21.adaptive.fdr <- na.omit(cbind(glm.all.rolwin21[,c(1:2)], glm.all.rolwin21$PAvE.T1.fdr.rolwin21, glm.all.rolwin21$PAvE.T4.fdr.rolwin21, glm.all.rolwin21$SvE.T4.fdr.rolwin21, glm.all.rolwin21$EvSE.T1.fdr.rolwin21))

## For filtering, calculate row maximums FDR to ensure all contrasts meet threshold
glm.all.rolwin21.adaptive.fdr$fdr_max <- rowMaxs(as.matrix(glm.all.rolwin21.adaptive.fdr[,c(3:6)]))

## Sort by maximum FDR
glm.all.rolwin21.adaptive.fdr <- glm.all.rolwin21.adaptive.fdr[order(glm.all.rolwin21.adaptive.fdr[,7]),]

## Look for loci where all contrasts of interest are significant at FDR < 0.01, control 
## allele frequencies at TPT1 are between 0.1 and 0.9, and AF difference between PA and E
## at TPT1 is greater than 0.1.
dim(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.01,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")))

## Look for loci where all contrasts of interest are significant at FDR < 0.05, control 
## allele frequencies at TPT1 are between 0.1 and 0.9, and AF difference between PA and E
## at TPT1 is greater than 0.1.
dim(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.05,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")))

## Since there's minimal overlap in SE SP and PA at FDR<0.01 or FDR<0.05, we'll find 
## the top genes (FDR < 0.1) that trend towards significance at a relaxed threshold.
glm.all.rolwin21.adaptive.fdr10 <- merge(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep, by=c("CHROM","POS"))

## Save gene list for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.adaptive.fdr10$Gene), file="glm.all.rolwin21.adaptive.fdr10.txt", sep = "\t", quote = FALSE, row.names = F)

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.adaptive.fdr10[glm.all.rolwin21.adaptive.fdr10$FLYCADD > 0.6,]$Gene), file="glm.all.rolwin21.adaptive.fdr10.flycadd60.txt", sep = "\t", quote = FALSE, row.names = F)

## Exclude flanking SNPs
glm.all.rolwin21.adaptive.fdr10_noflank <- unique(glm.all.rolwin21.adaptive.fdr10[glm.all.rolwin21.adaptive.fdr10$Consequence != "upstream_gene_variant" & glm.all.rolwin21.adaptive.fdr10$Consequence != "downstream_gene_variant",]$Gene)

## Save gene list, excluding flanking SNPs, for external GO enrichment analysis in BiNGO
write.table(glm.all.rolwin21.adaptive.fdr10_noflank, file="glm.all.rolwin21.adaptive.fdr10_noflank.txt", sep = "\t", quote = FALSE, row.names = F)

## Save gene list for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.adaptive.fdr10$Gene), file="glm.all.rolwin21.adaptive.fdr10.txt", sep = "\t", quote = FALSE, row.names = F)

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.adaptive.fdr10[glm.all.rolwin21.adaptive.fdr10$FLYCADD > 0.6,]$Gene), file="glm.all.rolwin21.adaptive.fdr10.flycadd60.txt", sep = "\t", quote = FALSE, row.names = F)

### TPT1 SE and PA convergent candidates ###

## Retrieve "TPT1 convergent candidate loci" using the following criteria:
#	TPT1 SE vs E GLM contrast significant at FDR < 0.05 (1st target threshold)
#	TPT1 PA vs E GLM contrast significant at FDR < 0.01 (2nd target threshold)
#	TPT4 PA vs E GLM contrast significant at FDR < 0.01 (ensures PA & E remain divergent)
#	TPT1 E allele frequency < 0.9 (not rare alleles with possible inflated AF diff)
#	TPT1 E allele frequency > 0.1 (not rare alleles with possible inflated AF diff)
#	TPT1 absolute AF difference between PA and E > 0.1 (sizable AF diff between PA & E)
glm.all.rolwin21.TPT1.convergent_candidates <- unique(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep, by=c("CHROM","POS")))

## Save gene list for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.TPT1.convergent_candidates$Gene), file="glm.all.rolwin21.TPT1.convergent_candidate_genes.txt", sep = "\t", quote = FALSE, row.names = F)

## Retrieve "Functional TPT1 candidate loci": 
# 	Same filters as above plus FlyCADD score > 0.6
#	FlyCADD details available at https://github.com/JuliaBeets/FlyCADD
glm.all.rolwin21.TPT1.functional_convergent_candidates <- unique(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep[vep$FLYCADD > 0.6,], by=c("CHROM","POS")))

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.TPT1.functional_convergent_candidates$Gene), file="glm.all.rolwin21.TPT1.functional_convergent_candidate_genes.txt", sep = "\t", quote = FALSE, row.names = F)

## Save gene list, excluding upstream and downstream features, for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.TPT1.convergent_candidates[glm.all.rolwin21.TPT1.convergent_candidates$Consequence!="upstream_gene_variant" & glm.all.rolwin21.TPT1.convergent_candidates$Consequence!="downstream_gene_variant",]$Gene), file="glm.all.rolwin21.TPT1.convergent_candidate_genes_noIntergenic.txt", sep = "\t", quote = FALSE, row.names = F)

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.TPT1.functional_convergent_candidates[glm.all.rolwin21.TPT1.convergent_candidates$Consequence!="upstream_gene_variant" & glm.all.rolwin21.TPT1.convergent_candidates$Consequence!="downstream_gene_variant",]$Gene), file="glm.all.rolwin21.TPT1.functional_convergent_candidate_genes_noIntergenic.txt", sep = "\t", quote = FALSE, row.names = F)


### TPT4 SP and PA convergent candidates ###

## Retrieve "TPT4 convergent candidate loci" using the following criteria:
#	TPT4 SP vs E GLM contrast significant at FDR < 0.05 (1st target threshold)
#	TPT4 PA vs E GLM contrast significant at FDR < 0.01 (2nd target threshold)
#	TPT1 PA vs E GLM contrast significant at FDR < 0.01 (ensures PA & E start divergent)
#	TPT1 E allele frequency < 0.9 (not rare alleles with possible inflated AF diff)
#	TPT1 E allele frequency > 0.1 (not rare alleles with possible inflated AF diff)
#	TPT1 absolute AF difference between PA and E > 0.1 (sizable AF diff between PA & E)
glm.all.rolwin21.TPT4.convergent_candidates <- unique(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep, by=c("CHROM","POS")))

## Save gene list for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.TPT4.convergent_candidates$Gene), file="glm.all.rolwin21.TPT4.convergent_candidate_genes.txt", sep = "\t", quote = FALSE, row.names = F)

## Retrieve "Functional TPT4 candidate loci": 
# 	Same filters as above plus FlyCADD score > 0.6
#	FlyCADD details available at https://github.com/JuliaBeets/FlyCADD
glm.all.rolwin21.TPT4.functional_convergent_candidates <- unique(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by=c("CHROM","POS")), vep[vep$FLYCADD > 0.6,], by=c("CHROM","POS")))

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.TPT4.functional_convergent_candidates$Gene), file="glm.all.rolwin21.TPT4.functional_convergent_candidate_genes.txt", sep = "\t", quote = FALSE, row.names = F)

## Save gene list, excluding upstream and downstream features, for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.TPT4.convergent_candidates[glm.all.rolwin21.TPT4.convergent_candidates$Consequence!="upstream_gene_variant" & glm.all.rolwin21.TPT4.convergent_candidates$Consequence!="downstream_gene_variant",]$Gene), file="glm.all.rolwin21.TPT4.convergent_candidate_genes_noIntergenic.txt", sep = "\t", quote = FALSE, row.names = F)

## Save gene list w/ FlyCADD score > 0.6 for external GO enrichment analysis in BiNGO
write.table(unique(glm.all.rolwin21.TPT4.functional_convergent_candidates[glm.all.rolwin21.TPT4.convergent_candidates$Consequence!="upstream_gene_variant" & glm.all.rolwin21.TPT4.convergent_candidates$Consequence!="downstream_gene_variant",]$Gene), file="glm.all.rolwin21.TPT4.functional_convergent_candidate_genes_noIntergenic.txt", sep = "\t", quote = FALSE, row.names = F)

### Let's build residual Manhattan plots and overlay convergent candidates ###

### TPT1 convergent candidates ###

## Try a subtraction method to calculate residuals, in reference to PA vs E, to look 
## for loci where SE is approaching PA by T1.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SE contrasts from PA vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = TPT1 PA vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = TPT1 PA vs SE -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin21.PAvE.logdiffT1 <- as.data.frame(na.omit(cbind(CHROM=glm.all.rolwin21$CHROM, POS=glm.all.rolwin21$POS, PAvE.T1minusPAvSE.T1=(as.numeric(glm.all.rolwin21$PAvE.T1.logp.rolwin21)-(as.numeric(glm.all.rolwin21$PAvSE.T1.logp.rolwin21)/max(na.omit(as.numeric(glm.all.rolwin21$PAvSE.T1.logp.rolwin21)))*max(na.omit(as.numeric(glm.all.rolwin21$PAvE.T1.logp.rolwin21))))),PAvE.T1.fdr.rolwin21=glm.all.rolwin21$PAvE.T1.fdr.rolwin21,PAvE.T1.logp.rolwin21=glm.all.rolwin21$PAvE.T1.logp.rolwin21,PAvSE.T1.fdr.rolwin21=glm.all.rolwin21$PAvSE.T1.fdr.rolwin21,PAvSE.T1.logp.rolwin21=glm.all.rolwin21$PAvSE.T1.logp.rolwin21)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin21.PAvE.logdiffT1$residual <- as.numeric(glm.all.rolwin21.PAvE.logdiffT1[,3])
glm.all.rolwin21.PAvE.logdiffT1$z <- ave(as.numeric(glm.all.rolwin21.PAvE.logdiffT1[,3]), glm.all.rolwin21.PAvE.logdiffT1[,1], FUN=scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.PAvE.T1minusPAvSE.T1.rolwin21.zscore_hist.pdf")
	hist(glm.all.rolwin21.PAvE.logdiffT1$z)
dev.off()


## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin21.PAvE.logdiffT1[glm.all.rolwin21.PAvE.logdiffT1[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin21.PAvE.logdiffT1[,3] <- as.numeric(glm.all.rolwin21.PAvE.logdiffT1[,3])
glm.all.rolwin21.PAvE.logdiffT1[,2] <- as.integer(glm.all.rolwin21.PAvE.logdiffT1[,2])
glm.all.rolwin21.PAvE.logdiffT1[,8] <- as.numeric(glm.all.rolwin21.PAvE.logdiffT1[,8])

## Now build plot
manh.PAvE.T1minusPAvSE.T1.rolwin21 <- ggplot(glm.all.rolwin21.PAvE.logdiffT1, aes(POS, PAvE.T1minusPAvSE.T1)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T1.fdr.rolwin21<0.05,]$PAvE.T1.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T1.fdr.rolwin21<0.01,]$PAvE.T1.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T1minusPAvSE.T1)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T1minusPAvSE.T1)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T1minusPAvSE.T1)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T1minusPAvSE.T1)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T1minusPAvSE.T1)), color = "#4C4C4C", size = 1) +
  	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin21.PAvE.logdiffT1[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs PA minus PA vs SE (TPT1): 21-SNP sliding window")


## Try a subtraction method to calculate residuals, in reference to SE vs E, to look 
## for loci where SE is approaching PA by T1.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SE contrasts from SE vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = TPT1 SE vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = TPT1 PA vs SE -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin21.SEvE.logdiffT1 <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin21$CHROM, POS = glm.all.rolwin21$POS, EvSE.T1minusPAvSE.T1 = (as.numeric(glm.all.rolwin21$EvSE.T1.logp.rolwin21) - (as.numeric(glm.all.rolwin21$PAvSE.T1.logp.rolwin21)/max(na.omit(as.numeric(glm.all.rolwin21$PAvSE.T1.logp.rolwin21)))*max(na.omit(as.numeric(glm.all.rolwin21$EvSE.T1.logp.rolwin21))))),EvSE.T1.fdr.rolwin21 = glm.all.rolwin21$EvSE.T1.fdr.rolwin21,EvSE.T1.logp.rolwin21 = glm.all.rolwin21$EvSE.T1.logp.rolwin21,PAvSE.T1.fdr.rolwin21 = glm.all.rolwin21$PAvSE.T1.fdr.rolwin21,PAvSE.T1.logp.rolwin21 = glm.all.rolwin21$PAvSE.T1.logp.rolwin21)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin21.SEvE.logdiffT1$residual <- as.numeric(glm.all.rolwin21.SEvE.logdiffT1[,3])
glm.all.rolwin21.SEvE.logdiffT1$z <- ave(as.numeric(glm.all.rolwin21.SEvE.logdiffT1[,3]), glm.all.rolwin21.SEvE.logdiffT1[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.SEvE.T1minusPAvSE.T1.rolwin21.zscore_hist.pdf")
	hist(glm.all.rolwin21.SEvE.logdiffT1$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin21.SEvE.logdiffT1[glm.all.rolwin21.SEvE.logdiffT1[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin21.SEvE.logdiffT1[,3] <- as.numeric(glm.all.rolwin21.SEvE.logdiffT1[,3])
glm.all.rolwin21.SEvE.logdiffT1[,2] <- as.integer(glm.all.rolwin21.SEvE.logdiffT1[,2])
glm.all.rolwin21.SEvE.logdiffT1[,8] <- as.numeric(glm.all.rolwin21.SEvE.logdiffT1[,8])

## Now build plot
manh.EvSE.T1minusPAvSE.T1.rolwin21 <- ggplot(glm.all.rolwin21.SEvE.logdiffT1, aes(POS, EvSE.T1minusPAvSE.T1)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21<0.05,]$EvSE.T1.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21<0.01,]$EvSE.T1.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SEvE.logdiffT1, by = c("CHROM","POS")),aes(POS, as.numeric(EvSE.T1minusPAvSE.T1)), color = "#4C4C4C", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin21.SEvE.logdiffT1[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("E vs SE minus PA vs SE (TPT1): 21-SNP sliding window")


### TPT4 convergent candidates ###

## Try a subtraction method to calculate residuals, in reference to PA vs E, to look 
## for loci where SP is approaching PA by T4.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SP contrasts from PA vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = TPT4 PA vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = TPT4 PA vs SP -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin21.PAvE.logdiffT4 <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin21$CHROM, POS = glm.all.rolwin21$POS, PAvE.T4minusPAvS.T4 = (as.numeric(glm.all.rolwin21$PAvE.T4.logp.rolwin21) - (as.numeric(glm.all.rolwin21$PAvS.T4.logp.rolwin21)/max(na.omit(as.numeric(glm.all.rolwin21$PAvS.T4.logp.rolwin21)))*max(na.omit(as.numeric(glm.all.rolwin21$PAvE.T4.logp.rolwin21))))),PAvE.T4.fdr.rolwin21 = glm.all.rolwin21$PAvE.T4.fdr.rolwin21,PAvE.T4.logp.rolwin21 = glm.all.rolwin21$PAvE.T4.logp.rolwin21,PAvS.T4.fdr.rolwin21 = glm.all.rolwin21$PAvS.T4.fdr.rolwin21,PAvS.T4.logp.rolwin21 = glm.all.rolwin21$PAvS.T4.logp.rolwin21)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin21.PAvE.logdiffT4$residual <- as.numeric(glm.all.rolwin21.PAvE.logdiffT4[,3])
glm.all.rolwin21.PAvE.logdiffT4$z <- ave(as.numeric(glm.all.rolwin21.PAvE.logdiffT4[,3]), glm.all.rolwin21.PAvE.logdiffT4[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.PAvE.T4minusPAvS.T4.rolwin21.zscore_hist.pdf")
	hist(glm.all.rolwin21.PAvE.logdiffT4$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin21.PAvE.logdiffT4[glm.all.rolwin21.PAvE.logdiffT4[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin21.PAvE.logdiffT4[,3] <- as.numeric(glm.all.rolwin21.PAvE.logdiffT4[,3])
glm.all.rolwin21.PAvE.logdiffT4[,2] <- as.integer(glm.all.rolwin21.PAvE.logdiffT4[,2])
glm.all.rolwin21.PAvE.logdiffT4[,8] <- as.numeric(glm.all.rolwin21.PAvE.logdiffT4[,8])

## Now build plot
manh.PAvE.T4minusPAvS.T4.rolwin21 <- ggplot(glm.all.rolwin21.PAvE.logdiffT4, aes(POS, PAvE.T4minusPAvS.T4)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T4.fdr.rolwin21<0.05,]$PAvE.T4.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T4.fdr.rolwin21<0.01,]$PAvE.T4.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(PAvE.T4minusPAvS.T4)), color = "#4C4C4C", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin21.PAvE.logdiffT4[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("PA vs E minus PA vs S (TPT4): 21-SNP sliding window")


## Try a subtraction method to calculate residuals, in reference to SP vs E, to look 
## for loci where SP is approaching PA by T4.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SP contrasts from SP vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = TPT4 SP vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = TPT4 PA vs SP -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin21.SPvE.logdiffT4 <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin21$CHROM, POS = glm.all.rolwin21$POS, SvE.T4minusPAvS.T4 = (as.numeric(glm.all.rolwin21$SvE.T4.logp.rolwin21) - (as.numeric(glm.all.rolwin21$PAvS.T4.logp.rolwin21)/max(na.omit(as.numeric(glm.all.rolwin21$PAvS.T4.logp.rolwin21)))*max(na.omit(as.numeric(glm.all.rolwin21$SvE.T4.logp.rolwin21))))),SvE.T4.fdr.rolwin21 = glm.all.rolwin21$SvE.T4.fdr.rolwin21,SvE.T4.logp.rolwin21 = glm.all.rolwin21$SvE.T4.logp.rolwin21,PAvS.T4.fdr.rolwin21 = glm.all.rolwin21$PAvS.T4.fdr.rolwin21,PAvS.T4.logp.rolwin21 = glm.all.rolwin21$PAvS.T4.logp.rolwin21)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin21.SPvE.logdiffT4$residual <- as.numeric(glm.all.rolwin21.SPvE.logdiffT4[,3])
glm.all.rolwin21.SPvE.logdiffT4$z <- ave(as.numeric(glm.all.rolwin21.SPvE.logdiffT4[,3]), glm.all.rolwin21.SPvE.logdiffT4[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.SPvE.T4minusPAvS.T4.rolwin21.zscore_hist.pdf")
	hist(glm.all.rolwin21.SPvE.logdiffT4$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin21.SPvE.logdiffT4[glm.all.rolwin21.SPvE.logdiffT4[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin21.SPvE.logdiffT4[,3] <- as.numeric(glm.all.rolwin21.SPvE.logdiffT4[,3])
glm.all.rolwin21.SPvE.logdiffT4[,2] <- as.integer(glm.all.rolwin21.SPvE.logdiffT4[,2])
glm.all.rolwin21.SPvE.logdiffT4[,8] <- as.numeric(glm.all.rolwin21.SPvE.logdiffT4[,8])

## Now build plot
manh.SPvE.T4minusPAvS.T4.rolwin21 <- ggplot(glm.all.rolwin21.SPvE.logdiffT4, aes(POS, SvE.T4minusPAvS.T4)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.T4.fdr.rolwin21<0.05,]$SvE.T4.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.T4.fdr.rolwin21<0.01,]$SvE.T4.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SPvE.logdiffT4, by = c("CHROM","POS")),aes(POS, as.numeric(SvE.T4minusPAvS.T4)), color = "#4C4C4C", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin21.SPvE.logdiffT4[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("S vs E minus PA vs S (TPT4): 21-SNP sliding window")



# This manhattan plot shows the contrast between SP and E samples combined in TPT 4
manh.EvS.T4.rolwin21 <- ggplot(glm.all.rolwin21, aes(POS, SvE.T4.logp.rolwin21)) + 
	geom_line(alpha = 0.5, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.T4.fdr.rolwin21 < 0.05,]$SvE.T4.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.T4.fdr.rolwin21 < 0.05,]),aes(POS, SvE.T4.logp.rolwin21), color = "#39568A", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
	ggtitle("E vs S (TPT4): 21-SNP sliding window")
	
## Plot both TPT4 treatment comparisons together
pdf(file = "rudflies_2023_redo.EvS.T4.rolwin21.glm.manh.pdf", width=10, height=3)
	manh.EvS.T4.rolwin21
dev.off()

## Panel E
# This manhattan plot shows the contrast between PA and E samples combined in TPT 1
manh.PAvE.T4.rolwin21 <- ggplot(glm.all.rolwin21, aes(POS, PAvE.T4.logp.rolwin21)) + 
	geom_line(alpha = 0.5, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T4.fdr.rolwin21 < 0.01,]$PAvE.T4.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.T4.fdr.rolwin21 < 0.01,]),aes(POS, PAvE.T4.logp.rolwin21), color = "#39568A", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
	ggtitle("E vs PA (TPT4): 21-SNP sliding window")

## Plot both TPT4 treatment comparisons together
pdf(file = "rudflies_2023_redo.PAvE.T4.rolwin21.glm.manh.pdf", width=10, height=3)
	manh.PAvE.T4.rolwin21
dev.off()


# This manhattan plot shows the contrast between PA and SP samples combined in TPT 4
manh.PAvS.T4.rolwin21 <- ggplot(glm.all.rolwin21, aes(POS, PAvS.T4.logp.rolwin21)) + 
	geom_line(alpha = 0.5, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvS.T4.fdr.rolwin21 < 0.01,]$PAvS.T4.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in SE vs E GLM contrast
	#geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvS.T4.fdr.rolwin21 < 0.01,]),aes(POS, PAvS.T4.logp.rolwin21), color = "#39568A", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
	ggtitle("PA vs S (TPT4): 21-SNP sliding window")
	
## Plot both TPT4 treatment comparisons together
pdf(file = "rudflies_2023_redo.PAvS.T4.rolwin21.glm.manh.pdf", width=10, height=3)
	manh.PAvS.T4.rolwin21
dev.off()



## Plot both TPT4 treatment comparisons together
pdf(file = "rudflies_2023_redo.PAvE.T4.EvS.T4.rolwin21.glm.manh.pdf", width=5.625, height = 4)
	ggarrange(manh.PAvE.T4.rolwin21, manh.EvS.T4.rolwin21, ncol = 1, nrow = 2)
dev.off()

## Plot three TPT4 treatment comparisons together
pdf(file = "rudflies_2023_redo.PAvE.T4.EvS.T4.PAvS.T4.rolwin21.glm.manh.pdf", width=5.625, height = 6)
	ggarrange(manh.PAvE.T4.rolwin21, manh.EvS.T4.rolwin21, manh.PAvS.T4.rolwin21, ncol = 1, nrow = 3)
dev.off()


  
### All TPT convergent candidates ###

## Try a subtraction method to calculate residuals, in reference to PA vs E, to look 
## for loci where S is approaching PA over all timepoints.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SE contrasts from PA vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = All TPT PA vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = All TPT PA vs S -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin21.PAvE.logdiff <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin21$CHROM, POS = glm.all.rolwin21$POS, PAvEminusPAvS = (as.numeric(glm.all.rolwin21$PAvE.logp.rolwin21) - (as.numeric(glm.all.rolwin21$PAvS.logp.rolwin21)/max(na.omit(as.numeric(glm.all.rolwin21$PAvS.logp.rolwin21)))*max(na.omit(as.numeric(glm.all.rolwin21$PAvE.logp.rolwin21))))),PAvE.fdr.rolwin21 = glm.all.rolwin21$PAvE.fdr.rolwin21,PAvE.logp.rolwin21 = glm.all.rolwin21$PAvE.logp.rolwin21,PAvS.fdr.rolwin21 = glm.all.rolwin21$PAvS.fdr.rolwin21,PAvS.logp.rolwin21 = glm.all.rolwin21$PAvS.logp.rolwin21)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin21.PAvE.logdiff$residual <- as.numeric(glm.all.rolwin21.PAvE.logdiff[,3])
glm.all.rolwin21.PAvE.logdiff$z <- ave(as.numeric(glm.all.rolwin21.PAvE.logdiff[,3]), glm.all.rolwin21.PAvE.logdiff[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.PAvEminusPAvS.rolwin21.zscore_hist.pdf")
	hist(glm.all.rolwin21.PAvE.logdiff$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin21.PAvE.logdiff[glm.all.rolwin21.PAvE.logdiff[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin21.PAvE.logdiff[,3] <- as.numeric(glm.all.rolwin21.PAvE.logdiff[,3])
glm.all.rolwin21.PAvE.logdiff[,2] <- as.integer(glm.all.rolwin21.PAvE.logdiff[,2])
glm.all.rolwin21.PAvE.logdiff[,8] <- as.numeric(glm.all.rolwin21.PAvE.logdiff[,8])

## Now build plot
manh.PAvEminusPAvS.rolwin21 <- ggplot(glm.all.rolwin21.PAvE.logdiff, aes(POS, PAvEminusPAvS)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.fdr.rolwin21<0.05,]$PAvE.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.fdr.rolwin21<0.01,]$PAvE.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(PAvEminusPAvS)), color = "#4C4C4C", size = 1) +
  	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin21.PAvE.logdiff[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("PA vs E minus PA vs S: 21-SNP sliding window")


## Try a subtraction method to calculate residuals, in reference to S vs E, to look 
## for loci where S is approaching PA over all timepoints.

## This is a really long R command, but basically we're creating a table with locus info,
## the residual values, the input -log10(p) values, and the corresponding FDR values. 

## To calculate residuals, we subtract the maximum-normalized -log10(p) values for 
## PA vs SE contrasts from PA vs E -log10(p) values. We use the following formula: 
#	residual = An - (Bn / Bmax * Amax)
#		An = All TPT S vs E -log10(p) value for locus n
#		Amax = maximum A value across all loci
#		Bn = All TPT PA vs S -log10(p) value for locus n
#		Bmax = maximum A value across all loci
glm.all.rolwin21.SvE.logdiff <- as.data.frame(na.omit(cbind(CHROM = glm.all.rolwin21$CHROM, POS = glm.all.rolwin21$POS, SvEminusPAvS = (as.numeric(glm.all.rolwin21$SvE.logp.rolwin21) - (as.numeric(glm.all.rolwin21$PAvS.logp.rolwin21)/max(na.omit(as.numeric(glm.all.rolwin21$PAvS.logp.rolwin21)))*max(na.omit(as.numeric(glm.all.rolwin21$SvE.logp.rolwin21))))),SvE.fdr.rolwin21 = glm.all.rolwin21$SvE.fdr.rolwin21,SvE.logp.rolwin21 = glm.all.rolwin21$SvE.logp.rolwin21,PAvS.fdr.rolwin21 = glm.all.rolwin21$PAvS.fdr.rolwin21,PAvS.logp.rolwin21 = glm.all.rolwin21$PAvS.logp.rolwin21)))

## Create new residual and z-score columns for exploratory analyses
glm.all.rolwin21.SvE.logdiff$residual <- as.numeric(glm.all.rolwin21.SvE.logdiff[,3])
glm.all.rolwin21.SvE.logdiff$z <- ave(as.numeric(glm.all.rolwin21.SvE.logdiff[,3]), glm.all.rolwin21.SvE.logdiff[,1], FUN = scale)

## Exploratory plot to visualize spread of z-scores. 
pdf(file = "rudflies_2023_redo.SvEminusPAvS.rolwin21.zscore_hist.pdf")
	hist(glm.all.rolwin21.SvE.logdiff$z)
dev.off()

## According to the motivating hypothesis, PA and S populations are both more 
## differentiated from E than PA and S are from each other at convergent adaptive loci. 
## These loci of interest will have high positive residual scores. Some residual scores 
## end up negative though, and by definition these would show higher divergence between 
## S vs P than either have from E control pops. Therefore, these loci are considered 
## non-convergent (divergent) between PA and S, and the values are floored at 0 for 
## plotting purposes. 
glm.all.rolwin21.SvE.logdiff[glm.all.rolwin21.SvE.logdiff[,3] < 0,3] <- 0

## Format data types for plotting
glm.all.rolwin21.SvE.logdiff[,3] <- as.numeric(glm.all.rolwin21.SvE.logdiff[,3])
glm.all.rolwin21.SvE.logdiff[,2] <- as.integer(glm.all.rolwin21.SvE.logdiff[,2])
glm.all.rolwin21.SvE.logdiff[,8] <- as.numeric(glm.all.rolwin21.SvE.logdiff[,8])

## Now build plot
manh.SvEminusPAvS.rolwin21 <- ggplot(glm.all.rolwin21.SvE.logdiff, aes(POS, SvEminusPAvS)) + geom_line(alpha = 0.5, colour = "grey") +  
	## Add threshold that all sig. SNPs (FDR<0.05) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.fdr.rolwin21<0.05,]$SvE.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept = min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.fdr.rolwin21<0.01,]$SvE.logp.rolwin21)), color = "black",linetype = "dashed",linewidth = .25) +
	## Highlight TPT4 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#97AA71", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")),glm.all.rolwin21.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#C95B5B", size = 1) +
	## Highlight TPT1 convergent candidate loci, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#78B0BF", size = 1) +
	## Highlight TPT4 convergent candidate loci w/ FlyCADD > 0.6, filters described above
	geom_point(data = merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD > 0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#E8B53E", size = 1) +
	## Highlight top convergent candidate loci, filters described above
	geom_point(data = merge(merge(glm.all.rolwin21.adaptive.fdr[glm.all.rolwin21.adaptive.fdr$fdr_max < 0.1,c(1:2)], freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), glm.all.rolwin21.SvE.logdiff, by = c("CHROM","POS")),aes(POS, as.numeric(SvEminusPAvS)), color = "#4C4C4C", size = 1) +
	#visual options
  	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  	scale_y_continuous(limits = c(0, round(max(as.numeric(na.omit(glm.all.rolwin21.SvE.logdiff[,3])))))) +
  	scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
  	labs(col = "candidate\ngene\n-log10(p)") +
  	xlab("chromosome position") +
  	ylab("-log10(p)") +
  	theme_classic() +
  	theme(legend.position = "none") +
  	ggtitle("S vs E minus PA vs S: 21-SNP sliding window")



### Lastly, let's make some reference manhattan panels that showcase non-residual -log10p
### values for GLMs across all timepoints. These are modified from figure 2.

# This manhattan plot shows the contrast between S and E samples combined across all TPTs
manh.SvE.logp.rolwin21 <- ggplot(na.omit(glm.all.rolwin21), aes(POS, SvE.logp.rolwin21)) +
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though 
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.fdr.rolwin21 < 0.01,]$SvE.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in S vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.fdr.rolwin21 < 0.01,]), aes(POS, SvE.logp.rolwin21), color = "darkmagenta", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$SvE.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs S")
  
# This manhattan plot shows the contrast between PA and E samples combined across all TPTs
manh.PAvE.logp.rolwin21 <- ggplot(na.omit(glm.all.rolwin21), aes(POS, PAvE.logp.rolwin21)) + 
	geom_line(alpha = 1, colour = "#CCCCCC") +  
	## Add threshold that all sig. SNPs (FDR<0.01) fall above, not all SNPs over line significant though
	#geom_hline(yintercept=min(na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.fdr.rolwin21 < 0.01,]$PAvE.logp.rolwin21)),color="black",linetype="dashed",linewidth=.25) +
	## Highlight points significant in PA vs E GLM contrast
	geom_point(data=na.omit(glm.all.rolwin21[glm.all.rolwin21$PAvE.fdr.rolwin21 < 0.01,]),aes(POS, PAvE.logp.rolwin21), color = "darkmagenta", size = 1) +
	## General formatting commands
	facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_y_continuous(limits = c(0, max(na.omit(glm.all.rolwin21$PAvE.logp.rolwin21)))) +
	scale_x_continuous(breaks=c(0, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000),guide = guide_axis(angle = 45)) + 
	labs(col="candidate\ngene\n-log10(p)") +
	xlab("chromosome position") +
	ylab("-log10(p)") +
	theme_classic() +
	theme(legend.position = "none") +
	ggtitle("E vs PA")


### Plot all reciprocal residual plots, modified for Figure 4 and Figure S8 ###

## 	In the plots below the highlighted points are as follows:
##	blue - TPT1 SE and PA convergent adaptive loci
##	yellow - TPT1 SE and PA convergent adaptive loci w/ FlyCADD score > 0.06
##	green - TPT4 SP and PA convergent adaptive loci
##	brown - TPT4 SP and PA convergent adaptive loci w/ FlyCADD score > 0.06
##	black - TPT1 SE, TPT4 SP, and PA convergent adaptive loci
pdf(file = "rudflies_2023_redo.convergent_evo_residuals.rolwin21.glm.manh.pdf", width = 11.25, height = 8)
	ggarrange(manh.PAvE.logp.rolwin21,
			  manh.SvE.logp.rolwin21,
			  manh.PAvE.T1minusPAvSE.T1.rolwin21, 
			  manh.EvSE.T1minusPAvSE.T1.rolwin21, 
			  manh.PAvE.T4minusPAvS.T4.rolwin21, 
			  manh.SPvE.T4minusPAvS.T4.rolwin21, 
			  manh.PAvEminusPAvS.rolwin21, 
			  manh.SvEminusPAvS.rolwin21, 
              ncol = 2, nrow = 4)
dev.off()


### Allele frequency plots over time for top candidates ###
## Find top candidates through rankings of GLM significance and FlyCADD score

## TPT1 SE & PA convergent adaptive candidates
## Merge functional candidate list with residual score tables
TPT1_top_candidates <- merge(merge(glm.all.rolwin21.TPT1.functional_convergent_candidates, glm.all.rolwin21.SEvE.logdiffT1, by=c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiffT1, by=c("CHROM","POS"))

## Rank by FlyCADD scores
TPT1_top_candidates$flycadd_rank <- NA
TPT1_top_candidates$flycadd_rank <- rank(TPT1_top_candidates$FLYCADD)
## Normalize FlyCADD ranks
TPT1_top_candidates$flycadd_rank <- TPT1_top_candidates$flycadd_rank/max(TPT1_top_candidates$flycadd_rank)*100
## Rank by TPT1 convergent candidate SE residual scores
TPT1_top_candidates$SE_rank <- NA
TPT1_top_candidates$SE_rank <- rank(TPT1_top_candidates$EvSE.T1minusPAvSE.T1)
## Normalize TPT1 SE residual ranks
TPT1_top_candidates$SE_rank <- TPT1_top_candidates$SE_rank/max(TPT1_top_candidates$SE_rank)*100
## Rank by TPT1 convergent candidate PA residual scores
TPT1_top_candidates$PA_rank <- NA
TPT1_top_candidates$PA_rank <- rank(TPT1_top_candidates$PAvE.T1minusPAvSE.T1)
## Normalize TPT1 PA residual ranks
TPT1_top_candidates$PA_rank <- TPT1_top_candidates$PA_rank/max(TPT1_top_candidates$PA_rank)*100
## Find mean rank
TPT1_top_candidates$mean_rank1 <- NA
TPT1_top_candidates$mean_rank1 <- (TPT1_top_candidates$flycadd_rank + TPT1_top_candidates$SE_rank + TPT1_top_candidates$PA_rank)/3
## Find top ranked loci
TPT1_top_candidates[TPT1_top_candidates$mean_rank1==max(TPT1_top_candidates$mean_rank1),]

## Save table for TPT1_top_candidates
write.table(TPT1_top_candidates, file="rudflies_2023_redo.TPT1_top_candidates.txt", sep = "\t", quote = FALSE, row.names = F)


### Plot allele frequencies of top candidates based on mean rank for TPT1

## 3R_15541740 - top TPT1 candidate in SE, PA, and ranked flycadd scores
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="15541740",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3R_15541740.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:15541740") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

### Calculate mean ranks again, this time without FlyCADD scores
TPT1_top_candidates$mean_rank2 <- (TPT1_top_candidates$SE_rank + TPT1_top_candidates$PA_rank)/2
TPT1_top_candidates[TPT1_top_candidates$mean_rank2==max(TPT1_top_candidates$mean_rank2),]


### Plot allele frequencies of top candidates based on mean rank 2 for TPT1

## 3R_15319344 - top TPT1 candidate in SE and PA (plus minimum flycadd score of 0.6)
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="15319344",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3R_15319344.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:15319344") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

## TPT4 SP & PA convergent adaptive candidates
TPT4_top_candidates <- merge(merge(glm.all.rolwin21.TPT4.functional_convergent_candidates, glm.all.rolwin21.SPvE.logdiffT4, by=c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiffT4, by=c("CHROM","POS"))
## Rank by FlyCADD scores
TPT4_top_candidates$flycadd_rank <- NA
TPT4_top_candidates$flycadd_rank <- rank(TPT4_top_candidates$FLYCADD)
## Normalize FlyCADD ranks
TPT4_top_candidates$flycadd_rank <- TPT4_top_candidates$flycadd_rank/max(TPT4_top_candidates$flycadd_rank)*100
## Rank by TPT4 convergent candidate SP residual scores
TPT4_top_candidates$SP_rank <- NA
TPT4_top_candidates$SP_rank <- rank(TPT4_top_candidates$SvE.T4minusPAvS.T4)
## Normalize TPT4 SP residual ranks
TPT4_top_candidates$SP_rank <- TPT4_top_candidates$SP_rank/max(TPT4_top_candidates$SP_rank)*100
## Rank by TPT4 convergent candidate PA residual scores
TPT4_top_candidates$PA_rank <- NA
TPT4_top_candidates$PA_rank <- rank(TPT4_top_candidates$PAvE.T4minusPAvS.T4)
## Normalize TPT4 PA residual ranks
TPT4_top_candidates$PA_rank <- TPT4_top_candidates$PA_rank/max(TPT4_top_candidates$PA_rank)*100
## Find mean rank
TPT4_top_candidates$mean_rank1 <- NA
TPT4_top_candidates$mean_rank1 <- (TPT4_top_candidates$flycadd_rank + TPT4_top_candidates$SP_rank + TPT4_top_candidates$PA_rank)/3
## Find top ranked loci
TPT4_top_candidates[TPT4_top_candidates$mean_rank1==max(TPT4_top_candidates$mean_rank1),]
## Find top ranked loci on 3L
TPT4_top_candidates[TPT4_top_candidates$CHROM == "3L" & TPT4_top_candidates$mean_rank1==max(TPT4_top_candidates[TPT4_top_candidates$CHROM == "3L",]$mean_rank1),]
## Find top ranked loci on 3R
TPT4_top_candidates[TPT4_top_candidates$CHROM == "3R" & TPT4_top_candidates$mean_rank1==max(TPT4_top_candidates[TPT4_top_candidates$CHROM == "3R",]$mean_rank1),]

## Save table for TPT4_top_candidates
write.table(TPT4_top_candidates, file="rudflies_2023_redo.TPT4_top_candidates.txt", sep = "\t", quote = FALSE, row.names = F)


### Plot allele frequencies of top candidates based on mean rank for TPT4

## 3R_7856392 - top TPT4 candidate in SP, PA, and ranked flycadd scores
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="7856392",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3R_7856392.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:7856392") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

## 3L_23061477 - top TPT4 candidate in SP, PA, and ranked flycadd scores on 3l
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3L" & haf.sites.filt$POS=="23061477",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3L_23061477.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3L:23061477") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

## 3L_21382635 - top TPT4 candidate in SP, PA, and ranked flycadd scores with 5_prime_UTR_variant on 3l
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3L" & haf.sites.filt$POS=="21382635",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3L_21382635.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3L:21382635") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

### Calculate mean ranks again, this time without FlyCADD scores
TPT4_top_candidates$mean_rank2 <- (TPT4_top_candidates$SP_rank + TPT4_top_candidates$PA_rank)/2
## Find top ranked loci
TPT4_top_candidates[TPT4_top_candidates$mean_rank2==max(TPT4_top_candidates$mean_rank2),]
## Find top ranked loci on 3L
TPT4_top_candidates[TPT4_top_candidates$CHROM == "3L" & TPT4_top_candidates$mean_rank2==max(TPT4_top_candidates[TPT4_top_candidates$CHROM == "3L",]$mean_rank2),]
## Find top ranked loci on 3R
TPT4_top_candidates[TPT4_top_candidates$CHROM == "3R" & TPT4_top_candidates$mean_rank2==max(TPT4_top_candidates[TPT4_top_candidates$CHROM == "3R",]$mean_rank2),]


### Plot allele frequencies of top candidates based on mean rank 2 for TPT4

## 3R_7856044 - top TPT4 candidate in SP and PA (plus minimum flycadd score of 0.6)
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="7856044",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3R_7856044.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:7856044") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()


## 3R_13128437 - top TPT4 candidate in SP and PA *Missense
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="13128437",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3R_13128437.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:13128437") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()


## All TPT S & PA convergent adaptive candidates
AllTPT_top_candidates <- merge(merge(glm.all.rolwin21.adaptive.fdr10, glm.all.rolwin21.SvE.logdiff, by=c("CHROM","POS")), glm.all.rolwin21.PAvE.logdiff, by=c("CHROM","POS"))
## Rank by FlyCADD scores
AllTPT_top_candidates$flycadd_rank <- NA
AllTPT_top_candidates$flycadd_rank <- rank(AllTPT_top_candidates$FLYCADD)
## Normalize FlyCADD ranks
AllTPT_top_candidates$flycadd_rank <- AllTPT_top_candidates$flycadd_rank/max(AllTPT_top_candidates$flycadd_rank)*100
## Rank by All TPT convergent candidate S residual scores
AllTPT_top_candidates$S_rank <- NA
AllTPT_top_candidates$S_rank <- rank(AllTPT_top_candidates$SvEminusPAvS)
## Normalize All TPT S residual ranks
AllTPT_top_candidates$S_rank <- AllTPT_top_candidates$S_rank/max(AllTPT_top_candidates$S_rank)*100
## Rank by All TPT convergent candidate PA residual scores
AllTPT_top_candidates$PA_rank <- NA
AllTPT_top_candidates$PA_rank <- rank(AllTPT_top_candidates$PAvEminusPAvS)
## Normalize All TPT PA residual ranks
AllTPT_top_candidates$PA_rank <- AllTPT_top_candidates$PA_rank/max(AllTPT_top_candidates$PA_rank)*100
## Find mean rank
AllTPT_top_candidates$mean_rank1 <- NA
AllTPT_top_candidates$mean_rank1 <- (AllTPT_top_candidates$flycadd_rank + AllTPT_top_candidates$S_rank + AllTPT_top_candidates$PA_rank)/3
## Find top ranked loci
AllTPT_top_candidates[AllTPT_top_candidates$mean_rank1==max(AllTPT_top_candidates$mean_rank1),]

## Save table for AllTPT_top_candidates
write.table(AllTPT_top_candidates, file="rudflies_2023_redo.AllTPT_top_candidates.txt", sep = "\t", quote = FALSE, row.names = F)

## Save table for AllTPT_top_candidates
write.table(unique(AllTPT_top_candidates[AllTPT_top_candidates$Consequence!="upstream_gene_variant" & AllTPT_top_candidates$Consequence!="downstream_gene_variant",]$Gene), file="rudflies_2023_redo.AllTPT_top_candidates_noIntergenic_genes.txt", sep = "\t", quote = FALSE, row.names = F)


### Plot allele frequencies of top candidates based on mean rank for all-TPTs

## 3L_18866746 - top candidate in all TPT in S, PA, and ranked flycadd scores
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3L" & haf.sites.filt$POS=="18866746",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3L_18866746.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3L:18866746") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

AllTPT_top_candidates$mean_rank2 <- (AllTPT_top_candidates$S_rank + AllTPT_top_candidates$PA_rank)/2
AllTPT_top_candidates[AllTPT_top_candidates$mean_rank2==max(AllTPT_top_candidates$mean_rank2),]


### Plot allele frequencies of top candidates based on mean rank 2 for all-TPTs

## 3L_23040611 - top candidate in all TPT in S, PA, and ranked flycadd scores
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3L" & haf.sites.filt$POS=="23040611",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3L_23040611.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3L:23040611") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()

AllTPT_top_candidates$mean_rank2 <- (AllTPT_top_candidates$S_rank + AllTPT_top_candidates$PA_rank)/2
AllTPT_top_candidates[AllTPT_top_candidates$mean_rank2==max(AllTPT_top_candidates$mean_rank2),]

## 3R_13398258 - top candidate in all TPT in S and PA (but no flycadd score weighting)
tempplot <- cbind(haf.meta.filt,t(na.omit(haf.freq.filt[haf.sites.filt$CHROM=="3R" & haf.sites.filt$POS=="13398258",])))
colnames(tempplot)[13] <- "af"
tempplot$treat.fix <- as.character(tempplot$treat.fix)
tempplot[tempplot$tpt=="1" & (tempplot$cage=="11" | tempplot$cage=="21" | tempplot$cage=="27" | tempplot$cage=="41" | tempplot$cage=="45"),8] <- "SE"
tempplot[tempplot$treat.fix=="S",8] <- "SP"
tempplot$treat.fix <- as.factor(tempplot$treat.fix)

tempplot.mean <- tempplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_at(vars("af"), mean)
  
tempplot.mean <- tempplot.mean[tempplot.mean$tpt != 0,]

##strip chart
pdf(file = "rudflies_2023_redo.haf.3R_13398258.geneplot.PAvSvE.strip.pdf", width=4, height=4)
	ggplot(tempplot) + geom_point(position=position_dodge(width=0.5), aes(x=tpt, y=af, color=treat.fix),
              width=0.2, height=0.1, alpha=0.4, size=1, shape=19) + theme_classic() +  ggtitle("3R:13398258") + scale_color_manual(values = sample_cols) + geom_line(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix)) + geom_point(data=tempplot.mean,aes(x=tpt,y=af, group=treat.fix, color=treat.fix), size=2.5)
dev.off()



################
### Figure 4 ###
################

### Try plotting all convergent candidates on a single plot, per list ###

## Filter for founder samples for PA, S, and E. 
## This will be the 1st plotted time-point.
haf.meta.found <- haf.meta[haf.meta$batch == "F" & (haf.meta$treat.fix == "PA" | haf.meta$treat.fix == "E"),]
haf.meta.found <- haf.meta[(haf.meta$batch == "a" | haf.meta$batch == "F") & haf.meta$experiment == "spino" & (haf.meta$treat.fix == "PA" | haf.meta$treat.fix == "S" | haf.meta$treat.fix == "E"),]
## Find matching sample names in the frequency table, post metadata filtering
haf.freq.found <- haf.freq[, which((names(haf.freq) %in% haf.meta.found$samp)==TRUE)]
## And run the reciprocal filter
haf.meta.found <- haf.meta.found[which((haf.meta.found$samp %in% names(haf.freq.found))==TRUE),]

## additional filter to remove low variance loci
haf.sites.found <- na.omit(haf.freq[rowVars(as.matrix(haf.freq.found))>0.001,c(1:2)])
haf.freq.found <- na.omit(haf.freq.found[rowVars(as.matrix(haf.freq.found))>0.001,])

## tables need a bit of reformatting for some downstream plotting
haf.meta.found$treat.fix <- as.factor(haf.meta.found$treat.fix)
haf.meta.found$tpt <- as.factor(haf.meta.found$tpt)


##### T4 SP Outliers ######

## Create a table that contains metadata and population allele frequencies for all T4 
## convergent candidate loci.
multilocusplot <- cbind(haf.meta.found[,c(5,4,8)],
t(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), cbind(haf.sites.found,haf.freq.found), by = c("CHROM","POS"))[,-c(1,2)]))

## Relabel treatment column for differentiating SE and SP
multilocusplot$treat.fix <- as.character(multilocusplot$treat.fix)
## Label SE cages
multilocusplot[multilocusplot$tpt=="1" & (multilocusplot$cage=="11" | multilocusplot$cage=="21" | multilocusplot$cage=="27" | multilocusplot$cage=="41" | multilocusplot$cage=="45"),3] <- "SE" 
## Label remaining S cages as SP
multilocusplot[multilocusplot$treat.fix=="S",3] <- "SP"
multilocusplot$treat.fix <- as.factor(multilocusplot$treat.fix)
## Drop unneeded cage column
multilocusplot <- multilocusplot[,-1]

## Calculate the group mean allele frequencies
multilocusplot.mean <- multilocusplot %>%
  group_by(treat.fix, tpt) %>% 
  summarise_all(mean)

## Add extra E founder rows, as these will be plotted as the SE and SP founders as well.  
multilocusplot.mean <- rbind(multilocusplot.mean[1,],multilocusplot.mean[1,],multilocusplot.mean)

## Relabel S founders
multilocusplot.mean[1,1] <- "SE"
multilocusplot.mean[2,1] <- "SP"

## Reformat AF mean table and transpose
multilocusplot.mean <- data.frame(multilocusplot.mean)
multilocusplot.flip <- t(multilocusplot.mean[,c(3:ncol(multilocusplot.mean))])

## For visualization, we want to orient all loci so that the starting PA allele 
## frequencies are higher than starting E allele frequencies. Which allele is considered 
## the reference or alternate allele in this study is arbitrary and originally designated 
## by the reference genome we aligned to. For loci where PA is lower than E, we will 
## recalculate as "1 - AF". 
multilocusplot.flip[multilocusplot.flip[,7] > multilocusplot.flip[,17],]  <- 1 - multilocusplot.flip[multilocusplot.flip[,7] > multilocusplot.flip[,17],] 

## Now, we join the corrected mean allele frequencies to the metadata
multilocusplot.mean <- cbind(multilocusplot.mean[,c(1,2)], t(multilocusplot.flip - multilocusplot.flip[,3]))
## Reshape for plotting
multilocusplot.long <- reshape2::melt(multilocusplot.mean, c("treat.fix","tpt"))
multilocusplot.long[order(multilocusplot.long$treat.fix, multilocusplot.long$tpt),]

## One final round of reformating
multilocusplot.long$treat.fix <- as.factor(multilocusplot.long$treat.fix)
multilocusplot.long$tpt <- as.factor(multilocusplot.long$tpt)
multilocusplot.long$variable <- as.character(multilocusplot.long$variable)
multilocusplot.long$value <- as.numeric(multilocusplot.long$value)
multilocusplot.long$grouping <- paste(multilocusplot.long$treat.fix, multilocusplot.long$variable)


### Figure 4B ### 
pdf(file = "rudflies_2023_redo.haf.T4_SP_convergent_candidate_AFplot.pdf", width=5, height=4)
	ggplot(multilocusplot.long[multilocusplot.long$treat.fix!="SE",], aes(x=tpt,y=value, group=grouping, color=treat.fix)) + 
	geom_path(alpha=0.8, size=0.2) +
	theme_classic() +  
	ggtitle("convergent candidate loci") + 
	scale_color_manual(values = sample_cols[-3])
dev.off()


### Alternate version with additional filter of FlyCADD score > 0.6 ###
## Create a table that contains metadata and population allele frequencies for all T4 
## convergent candidate loci with predicted functional impact (high FlyCadd score).
multilocusplot2 <- cbind(haf.meta.found[,c(5,4,8)],
t(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$SvE.T4.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD>0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), cbind(haf.sites.found,haf.freq.found), by = c("CHROM","POS"))[,-c(1,2)]))

## Relabel treatment column for differentiating SE and SP
multilocusplot2$treat.fix <- as.character(multilocusplot2$treat.fix)
## Label SE cages
multilocusplot2[multilocusplot2$tpt=="1" & (multilocusplot2$cage=="11" | multilocusplot2$cage=="21" | multilocusplot2$cage=="27" | multilocusplot2$cage=="41" | multilocusplot2$cage=="45"),3] <- "SE"
## Label remaining S cages as SP
multilocusplot2[multilocusplot2$treat.fix=="S",3] <- "SP"
multilocusplot2$treat.fix <- as.factor(multilocusplot2$treat.fix)
## Drop unneeded cage column
multilocusplot2 <- multilocusplot2[,-1]

## Calculate the group mean allele frequencies
multilocusplot2.mean <- multilocusplot2 %>%
  group_by(treat.fix, tpt) %>% 
  summarise_all(mean)
 
## Add extra E founder rows, as these will be plotted as the SE and SP founders as well.   
multilocusplot2.mean <- rbind(multilocusplot2.mean[1,],multilocusplot2.mean[1,],multilocusplot2.mean)

## Relabel S founders
multilocusplot2.mean[1,1] <- "SE"
multilocusplot2.mean[2,1] <- "SP"

## Reformat AF mean table and transpose
multilocusplot2.mean <- data.frame(multilocusplot2.mean)
multilocusplot2.flip <- t(multilocusplot2.mean[,c(3:ncol(multilocusplot2.mean))])

## For visualization, we want to orient all loci so that the starting PA allele 
## frequencies are higher than starting E allele frequencies. Which allele is considered 
## the reference or alternate allele in this study is arbitrary and originally designated 
## by the reference genome we aligned to. For loci where PA is lower than E, we will 
## recalculate as "1 - AF". 
multilocusplot2.flip[multilocusplot2.flip[,7] > multilocusplot2.flip[,17],]  <- 1 - multilocusplot2.flip[multilocusplot2.flip[,7] > multilocusplot2.flip[,17],] 

## Now, we join the corrected mean allele frequencies to the metadata
multilocusplot2.mean <- cbind(multilocusplot2.mean[,c(1,2)], t(multilocusplot2.flip - multilocusplot2.flip[,3]))
multilocusplot2.long <- reshape2::melt(multilocusplot2.mean, c("treat.fix","tpt"))
multilocusplot2.long[order(multilocusplot2.long$treat.fix, multilocusplot2.long$tpt),]

## One final round of reformating
multilocusplot2.long$treat.fix <- as.factor(multilocusplot2.long$treat.fix)
multilocusplot2.long$tpt <- as.factor(multilocusplot2.long$tpt)
multilocusplot2.long$variable <- as.character(multilocusplot2.long$variable)
multilocusplot2.long$value <- as.numeric(multilocusplot2.long$value)
multilocusplot2.long$variable <- as.character(multilocusplot2.long$variable)
multilocusplot2.long$grouping <- paste(multilocusplot2.long$treat.fix, multilocusplot2.long$variable)
 
 
### Unused alternate version of Figure 4B ###
pdf(file = "rudflies_2023_redo.haf.T4_SP_convergent_candidate_AFplot.flycadd60.pdf", width=5, height=4)
	ggplot(multilocusplot2.long[multilocusplot2.long$treat.fix!="SE",], aes(x=tpt,y=value, group=grouping, color=treat.fix)) + 
	geom_path(alpha=0.8, size=0.2) +
	theme_classic() +  
	ggtitle("convergent candidate loci (FlyCADD > 0.6)") + 
	scale_color_manual(values = sample_cols[-3])
dev.off()



##### T1 SE Outliers #######

## Create a table that contains metadata and population allele frequencies for all T1 
## convergent candidate loci.
multilocusplot3 <- cbind(haf.meta.found[,c(5,4,8)],
t(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), cbind(haf.sites.found,haf.freq.found), by = c("CHROM","POS"))[,-c(1,2)]))

## Relabel treatment column for differentiating SE and SP
multilocusplot3$treat.fix <- as.character(multilocusplot3$treat.fix)
## Label SE cages
multilocusplot3[multilocusplot3$tpt=="1" & (multilocusplot3$cage=="11" | multilocusplot3$cage=="21" | multilocusplot3$cage=="27" | multilocusplot3$cage=="41" | multilocusplot3$cage=="45"),3] <- "SE" 
## Label remaining S cages as SP
multilocusplot3[multilocusplot3$treat.fix=="S",3] <- "SP"
multilocusplot3$treat.fix <- as.factor(multilocusplot3$treat.fix)
## Drop unneeded cage column
multilocusplot3 <- multilocusplot3[,-1]

## Calculate the group mean allele frequencies
multilocusplot3.mean <- multilocusplot3 %>%
  group_by(treat.fix, tpt) %>% 
  summarise_all(mean)

## Add extra E founder rows, as these will be plotted as the SE and SP founders as well.  
multilocusplot3.mean <- rbind(multilocusplot3.mean[1,],multilocusplot3.mean[1,],multilocusplot3.mean)

## Relabel S founders
multilocusplot3.mean[1,1] <- "SE"
multilocusplot3.mean[2,1] <- "SP"

## Reformat AF mean table and transpose
multilocusplot3.mean <- data.frame(multilocusplot3.mean)
multilocusplot3.flip <- t(multilocusplot3.mean[,c(3:ncol(multilocusplot3.mean))])

## For visualization, we want to orient all loci so that the starting PA allele 
## frequencies are higher than starting E allele frequencies. Which allele is considered 
## the reference or alternate allele in this study is arbitrary and originally designated 
## by the reference genome we aligned to. For loci where PA is lower than E, we will 
## recalculate as "1 - AF". 
multilocusplot3.flip[multilocusplot3.flip[,7] > multilocusplot3.flip[,17],]  <- 1 - multilocusplot3.flip[multilocusplot3.flip[,7] > multilocusplot3.flip[,17],] 

## Now, we join the corrected mean allele frequencies to the metadata
multilocusplot3.mean <- cbind(multilocusplot3.mean[,c(1,2)], t(multilocusplot3.flip - multilocusplot3.flip[,3]))
## Reshape for plotting
multilocusplot3.long <- reshape2::melt(multilocusplot3.mean, c("treat.fix","tpt"))
multilocusplot3.long[order(multilocusplot3.long$treat.fix, multilocusplot3.long$tpt),]

## One final round of reformating
multilocusplot3.long$treat.fix <- as.factor(multilocusplot3.long$treat.fix)
multilocusplot3.long$tpt <- as.factor(multilocusplot3.long$tpt)
multilocusplot3.long$variable <- as.character(multilocusplot3.long$variable)
multilocusplot3.long$value <- as.numeric(multilocusplot3.long$value)
multilocusplot3.long$grouping <- paste(multilocusplot3.long$treat.fix, multilocusplot3.long$variable)


### Figure 4B ### 
pdf(file = "rudflies_2023_redo.haf.T1_SE_convergent_candidate_AFplot.pdf", width=5, height=4)
	ggplot(multilocusplot3.long[multilocusplot3.long$treat.fix!="SE",], aes(x=tpt,y=value, group=grouping, color=treat.fix)) + 
	geom_path(alpha=0.8, size=0.2) +
	theme_classic() +  
	ggtitle("convergent candidate loci") + 
	scale_color_manual(values = sample_cols[-3])
dev.off()


### Alternate version with additional filter of FlyCADD score > 0.6 ###
## Create a table that contains metadata and population allele frequencies for all T1 
## convergent candidate loci with predicted functional impact (high FlyCadd score).
multilocusplot4 <- cbind(haf.meta.found[,c(5,4,8)],
t(merge(merge(unique(glm.all.rolwin21.annot[glm.all.rolwin21.annot$PAvE.T1.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$PAvE.T4.fdr.rolwin21<0.01 & glm.all.rolwin21.annot$EvSE.T1.fdr.rolwin21<0.05 & glm.all.rolwin21.annot$FLYCADD>0.6, c(1:2)]), freq_means_t1[freq_means_t1$E.af < 0.9 & freq_means_t1$E.af > 0.1 & abs(freq_means_t1$E.af - freq_means_t1$PA.af) > 0.1, c(1,2)], by = c("CHROM","POS")), cbind(haf.sites.found,haf.freq.found), by = c("CHROM","POS"))[,-c(1,2)]))

## Relabel treatment column for differentiating SE and SP
multilocusplot4$treat.fix <- as.character(multilocusplot4$treat.fix)
## Label SE cages
multilocusplot4[multilocusplot4$tpt=="1" & (multilocusplot4$cage=="11" | multilocusplot4$cage=="21" | multilocusplot4$cage=="27" | multilocusplot4$cage=="41" | multilocusplot4$cage=="45"),3] <- "SE"
## Label remaining S cages as SP
multilocusplot4[multilocusplot4$treat.fix=="S",3] <- "SP"
multilocusplot4$treat.fix <- as.factor(multilocusplot4$treat.fix)
## Drop unneeded cage column
multilocusplot4 <- multilocusplot4[,-1]

## Calculate the group mean allele frequencies
multilocusplot4.mean <- multilocusplot4 %>%
  group_by(treat.fix, tpt) %>% 
  summarise_all(mean)
 
## Add extra E founder rows, as these will be plotted as the SE and SP founders as well.   
multilocusplot4.mean <- rbind(multilocusplot4.mean[1,],multilocusplot4.mean[1,],multilocusplot4.mean)

## Relabel S founders
multilocusplot4.mean[1,1] <- "SE"
multilocusplot4.mean[2,1] <- "SP"

## Reformat AF mean table and transpose
multilocusplot4.mean <- data.frame(multilocusplot4.mean)
multilocusplot4.flip <- t(multilocusplot4.mean[,c(3:ncol(multilocusplot4.mean))])

## For visualization, we want to orient all loci so that the starting PA allele 
## frequencies are higher than starting E allele frequencies. Which allele is considered 
## the reference or alternate allele in this study is arbitrary and originally designated 
## by the reference genome we aligned to. For loci where PA is lower than E, we will 
## recalculate as "1 - AF". 
multilocusplot4.flip[multilocusplot4.flip[,7] > multilocusplot4.flip[,17],]  <- 1 - multilocusplot4.flip[multilocusplot4.flip[,7] > multilocusplot4.flip[,17],] 

## Now, we join the corrected mean allele frequencies to the metadata
multilocusplot4.mean <- cbind(multilocusplot4.mean[,c(1,2)], t(multilocusplot4.flip - multilocusplot4.flip[,3]))
multilocusplot4.long <- reshape2::melt(multilocusplot4.mean, c("treat.fix","tpt"))
multilocusplot4.long[order(multilocusplot4.long$treat.fix, multilocusplot4.long$tpt),]

## One final round of reformating
multilocusplot4.long$treat.fix <- as.factor(multilocusplot4.long$treat.fix)
multilocusplot4.long$tpt <- as.factor(multilocusplot4.long$tpt)
multilocusplot4.long$variable <- as.character(multilocusplot4.long$variable)
multilocusplot4.long$value <- as.numeric(multilocusplot4.long$value)
multilocusplot4.long$variable <- as.character(multilocusplot4.long$variable)
multilocusplot4.long$grouping <- paste(multilocusplot4.long$treat.fix, multilocusplot4.long$variable)
 
 
### Unused alternate version of Figure 4B ###
pdf(file = "rudflies_2023_redo.haf.T1_SE_convergent_candidate_AFplot.flycadd60.pdf", width=5, height=4)
	ggplot(multilocusplot4.long[multilocusplot4.long$treat.fix!="SE",], aes(x=tpt,y=value, group=grouping, color=treat.fix)) + 
	geom_path(alpha=0.8, size=0.2) +
	theme_classic() +  
	ggtitle("convergent candidate loci (FlyCADD > 0.6)") + 
	scale_color_manual(values = sample_cols[-3])
dev.off()


#######################################################################
############### MISCELLANEOUS SUPPLEMENTARY ANALYSES ##################
#######################################################################

##################
### Figure S11 ###
##################

### Checking whether the hafpipe frequencies match the dgrp2 founder file frequencies ###
## Load mean frequencies of the relevant DGRP2 samples from the DGRP Freeze.2 VCF file
dgrp2.freq <- read.delim("../../dgrp2_founder_list.frq", header=TRUE, sep = "\t")
## Load mean frequencies of experimental samples from VCF generated by BCFtools
vcf.freq <- read.delim("../calling/filtered-all.frq", header=TRUE, sep = "\t")
## Load mean frequencies of experimental sample from frequency table generated by Hafpipe
haf.freq.means <- cbind(haf.freq[,c(1:2)],rowMeans(haf.freq[,c(3:ncol(haf.freq))]))
## Turns out Hafpipe frequencies are generated for the ALT alleles frequencies instead of 
## REF alleles, so we'll calculate REF frequencies (1 - ALT) to positively correlate with 
## other datasets.
colnames(haf.freq.means) <- c("CHROM","POS","hafpipe.means.A")
haf.freq.means$hafpipe.means.R <- 1-haf.freq.means$hafpipe.means.A
## Merge all three tables by locus
all.freq <- merge(dgrp2.freq, vcf.freq, by=c("CHROM","POS"))
all.freq <- merge(all.freq, haf.freq.means, by=c("CHROM","POS"))

## Create a density plot comparing DGRP2 founder frequencies and BCFtools VCF frequencies
pdf(file = "rudflies_2023_redo.dgrp2_vs_vcf_freqs.density.pdf", width=6.5, height=5)
	ggplot(all.freq, 
	  aes(x=R.FREQ.x, y=R.FREQ.y)) + 
	  geom_bin2d(bins = 50) +
	  scale_fill_continuous(type = "viridis") +
	  xlim(0.05, 0.95) +
	  ylim(0.05, 0.95) +
	  xlab("dgrp2 AF") + 
	  ylab("bcftools AF") +
	  theme_classic() +
	  theme_update(axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.y = element_blank())
dev.off()

## Test correlation
cor.test(as.numeric(all.freq$R.FREQ.x), as.numeric(all.freq$R.FREQ.y), method="pearson")

## Create a density plot comparing DGRP2 founder frequencies and Hafpipe frequencies
pdf(file = "rudflies_2023_redo.dgrp2_vs_hafpipe_freqs.density.pdf", width=6.5, height=5)
	ggplot(all.freq, 
	  aes(x=R.FREQ.x, y=hafpipe.means.R)) + 
	  geom_bin2d(bins = 50) +
	  scale_fill_continuous(type = "viridis") +
	  xlim(0.05, 0.95) +
	  ylim(0.05, 0.95) +
	  xlab("dgrp2 AF") + 
	  ylab("hafpipe AF") +
	  theme_classic() +
	  theme_update(axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.y = element_blank())
dev.off()

## Test correlation
cor.test(as.numeric(all.freq$R.FREQ.x), as.numeric(all.freq$hafpipe.means.R), method="pearson")

## Create a density plot comparing BCFtools VCF frequencies and Hafpipe frequencies
pdf(file = "rudflies_2023_redo.vcf_vs_hafpipe_freqs.density.pdf", width=6.5, height=5)
	ggplot(all.freq, 
	  aes(x=R.FREQ.y, y=hafpipe.means.R)) + 
	  geom_bin2d(bins = 50) +
	  scale_fill_continuous(type = "viridis") +
	  xlim(0.05, 0.95) +
	  ylim(0.05, 0.95) +
	  xlab("bcftools AF") + 
	  ylab("hafpipe AF") +
	  theme_classic() +
	  theme_update(axis.ticks.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.y = element_blank())
dev.off()

## Test correlation 
cor.test(as.numeric(all.freq$R.FREQ.y), as.numeric(all.freq$hafpipe.means.R), method="pearson")

## Overall takeaway, hafpipe output is in reference to the alt-allele, not the ref
## Plus, hafpipe frequencies are much better correlated with the founder file than the bcftools output.


##################
### Figure S12 ###
##################

### General premise for this analysis is checking to see whether loci that showed signs of 
### adaptation in the persistent S cages may have showed early signs of adaptation in 
### extinct S cages.

## Use this color scheme for "condition" side colors
sidecols <- haf.meta.T1filt$condition 
sidecols <- gsub("PA","#495184",sidecols)
sidecols <- gsub("SP","#D9B851",sidecols)
sidecols <- gsub("SE","#848556",sidecols)
sidecols <- gsub("E","#D26183",sidecols)

### Heatmaps for TPT1 SE outliers SNPs among TPT1 AF frequencies ###
## Baseline to see what adaptive signals look like the the pops they were detected within

## Select SNPs from AF table significant (FDR<0.05) in GLM contrasts between E vs SE:TPT1
haf.freq.SE.fdr05 <- as.matrix(merge(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05,c(1:2)]),cbind(haf.sites.T1filt,haf.freq.T1filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
haf.freq.SE.fdr05.heat <- heatmap.2(haf.freq.SE.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")

pdf(file = "rudflies_2023_redo.haf.freq.SE.fdr05.TPT1.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SE.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.01) in GLM contrasts between E vs SE:TPT1
haf.freq.SE.fdr01 <- as.matrix(merge(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.01,c(1:2)]),cbind(haf.sites.T1filt,haf.freq.T1filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
haf.freq.SE.fdr01.heat <- heatmap.2(haf.freq.SE.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")

pdf(file = "rudflies_2023_redo.haf.freq.SE.fdr01.TPT1.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SE.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")
dev.off()


### Heatmaps for TPT4 SP outliers SNPs among TPT1 AF frequencies ###
## Are SE pops adapting early in some of the same regions SP eventually shows signal?

## Select SNPs from AF table significant (FDR<0.05) in GLM contrasts between E vs SP:TPT4
haf.freq.SP.fdr05 <- as.matrix(merge(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.T4.fdr.rolwin21 < 0.05,c(1:2)]),cbind(haf.sites.T1filt,haf.freq.T1filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
haf.freq.SP.fdr05.heat <- heatmap.2(haf.freq.SP.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")

pdf(file = "rudflies_2023_redo.haf.freq.SP.fdr05.TPT1.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SP.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.01) in GLM contrasts between E vs SP:TPT4
haf.freq.SP.fdr01 <- as.matrix(merge(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.T4.fdr.rolwin21 < 0.01,c(1:2)]),cbind(haf.sites.T1filt,haf.freq.T1filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
haf.freq.SP.fdr01.heat <- heatmap.2(haf.freq.SP.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")

pdf(file = "rudflies_2023_redo.haf.freq.SP.fdr01.TPT1.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SP.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols, labRow=FALSE, trace = "none")
dev.off()


### Now we're switching to outlier heatmaps among TPT4-sample frequencies ###

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

## Use this color scheme for "treatment" side colors
sidecols2 <- haf.meta.T4filt$treat.fix 
sidecols2 <- gsub("PA","#495184",sidecols2)
sidecols2 <- gsub("S","#D9B851",sidecols2)
sidecols2 <- gsub("E","#D26183",sidecols2)

### Heatmaps for TPT4 SP outlier SNPs among TPT1 AF frequencies ###
## Are SE pops adapting early in some of the same regions SP eventually shows adaptation?

## Select SNPs from AF table significant (FDR<0.05) in GLM contrasts between E vs SE:TPT1
haf.freq.SE.t4.fdr05 <- as.matrix(merge(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.05,c(1:2)]),cbind(haf.sites.T4filt,haf.freq.T4filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
haf.freq.SE.t4.fdr05.heat <- heatmap.2(haf.freq.SE.t4.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")

pdf(file = "rudflies_2023_redo.haf.freq.SE.fdr05.TPT4.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SE.t4.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.01) in GLM contrasts between E vs SE:TPT1
haf.freq.SE.t4.fdr01 <- as.matrix(merge(na.omit(glm.all.rolwin21[glm.all.rolwin21$EvSE.T1.fdr.rolwin21 < 0.01,c(1:2)]),cbind(haf.sites.T4filt,haf.freq.T4filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
haf.freq.SE.t4.fdr01.heat <- heatmap.2(haf.freq.SE.t4.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")

pdf(file = "rudflies_2023_redo.haf.freq.SE.fdr01.TPT4.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SE.t4.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")
dev.off()

### Heatmaps for TPT4 SP outlier SNPs among TPT4 AF frequencies ###
## Baseline to see what adaptive signals look like the the pops they were detected within

## Select SNPs from AF table significant (FDR<0.05) in GLM contrasts between E vs SP:TPT4
haf.freq.SP.t4.fdr05 <- as.matrix(merge(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.T4.fdr.rolwin21 < 0.05,c(1:2)]),cbind(haf.sites.T4filt,haf.freq.T4filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
haf.freq.SP.t4.fdr05.heat <- heatmap.2(haf.freq.SP.t4.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")

pdf(file = "rudflies_2023_redo.haf.freq.SP.fdr05.TPT4.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SP.t4.fdr05, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")
dev.off()

## Select SNPs from AF table significant (FDR<0.01) in GLM contrasts between E vs SP:TPT4
haf.freq.SP.t4.fdr01 <- as.matrix(merge(na.omit(glm.all.rolwin21[glm.all.rolwin21$SvE.T4.fdr.rolwin21 < 0.01,c(1:2)]),cbind(haf.sites.T4filt,haf.freq.T4filt), by=c("CHROM","POS"))[,-c(1:2)])

## Make the heatmap
haf.freq.SP.t4.fdr01.heat <- heatmap.2(haf.freq.SP.t4.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")

pdf(file = "rudflies_2023_redo.haf.freq.SP.fdr01.TPT4.heatmap.pdf", width=10, height=10)
	heatmap.2(haf.freq.SP.t4.fdr01, Rowv=FALSE, dendrogram="col", ColSideColors=sidecols2, labRow=FALSE, trace = "none")
dev.off()
