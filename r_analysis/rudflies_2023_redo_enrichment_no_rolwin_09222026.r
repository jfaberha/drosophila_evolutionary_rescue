##########################################################################################
### This script is associated with a study on the evolution of spinosad-resistance in  ###
### field populations of Drosophila melanogaster run in 2023 by the Rudman lab (WSUV). ###
### This R script is a shorter companion script to:                                    ###
### "rudflies_2023_redo_enrichment_and_cost_of_adaptation_09222026.r"                  ###
###                                                                                    ###
### The longer script listed above contains the primary analysis associated with the   ### 
### question; "Is there enrichment of an a priori list of insecticide resistance-      ###
### associated genes in the various lists of outlier genes in our study?" and it uses  ###
### rolling-window smoothed GLM FDR values to answer it. For supplementary analysis,   ###
### we also wanted to check if results would be any different using raw/unsmoothed     ###
### FDR values to determine outlier SNP sets, and those analyses are run using the     ###
### current script.                                                                    ###
##########################################################################################

### In R ###
#configure r environment
setwd("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/enrichment_cost_analysis")
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
vep <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/filtered-all.annot.vcf.FLYCADD.tsv", header=TRUE) 
## More detailed gene information from gff
snp.gff.overlap <- read.delim("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_hafpipe_loci_genic_overlap_info.bed", sep="\t", header=TRUE)
## Sample metadata table
haf.meta <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_meta.tsv", header=TRUE)
## Hafpipe imputed allele frequency table
haf.freq <- read.delim("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_hafpipe.csv", header=TRUE, sep = ",")

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

## Load spinosad-resistance candidate gene list
spino.cand <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/spino.cand.list.txt", header=FALSE)
names(spino.cand) <- "Gene"

## Now merge to find All SNPs in and around candidate genes
glm.all.annot.spino <- merge(spino.cand, glm.all.annot, by="Gene")


####################################
### Hypergeometric overlap tests ###
####################################

## Calculate minimum p-values and FDR values per annotated spino candidate genes. We will 
## assign the most significant FDR value overlapping the gene, plus its immediate upstream 
## and downstream regions, to represent the gene in comparison with background gene sets.

## Pull FDR columns only for spino candidates
glm.all.fdr.annot.spino <- unique(cbind(glm.all.annot.spino[,c(1:3)],glm.all.annot.spino$Consequence,select(glm.all.annot.spino,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate spino candidate minimum values
glm.all.fdr.annot.spino.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.annot.spino)) {
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.annot.spino$Gene,glm.all.fdr.annot.spino[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"
    #assign minimum FDR values
    glm.all.fdr.annot.spino.min <- cbind(glm.all.fdr.annot.spino.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.annot.spino.min) <- names(glm.all.fdr.annot.spino)[c(5:ncol(glm.all.fdr.annot.spino))]
## Add gene names to the table
glm.all.fdr.annot.spino.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.annot.spino.min)
        	
## Save minimum FDR results for downstream use
write.table(glm.all.fdr.annot.spino.min, file="rudflies_2023_redo.glm.all.minfdr.spino.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)
        	
## Pull FDR columns only for all genes/features
glm.all.fdr.annot <- unique(cbind(Gene=glm.all.annot$Gene,glm.all.annot[,c(1,2)],Consequence=glm.all.annot$Consequence,select(glm.all.annot,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate minimum values for all genes
glm.all.fdr.annot.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.annot)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.annot$Gene,glm.all.fdr.annot[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.annot.min <- cbind(glm.all.fdr.annot.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.annot.min) <- names(glm.all.fdr.annot)[c(5:ncol(glm.all.fdr.annot))]
## Add gene names to the table
glm.all.fdr.annot.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.annot.min)


### Redo above with missense SNPs in candidates only ###

## Merge to find missense SNPs in and around candidate genes
glm.all.annot.spino.missense <- merge(spino.cand, glm.all.annot[glm.all.annot$Consequence=="missense_variant",], by="Gene")

## Pull FDR columns only for spino candidate missense SNPs
glm.all.fdr.annot.spino.missense <- unique(cbind(glm.all.annot.spino.missense[,c(1:3)],glm.all.annot.spino.missense$Consequence,select(glm.all.annot.spino.missense,contains("fdr"))))

## Loop through each GLM-results FDR column to calculate spino candidate minimum values
glm.all.fdr.annot.spino.missense.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.annot.spino.missense)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.annot.spino.missense$Gene,glm.all.fdr.annot.spino.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.annot.spino.missense.min <- cbind(glm.all.fdr.annot.spino.missense.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.annot.spino.missense.min) <- names(glm.all.fdr.annot.spino.missense)[c(5:ncol(glm.all.fdr.annot.spino.missense))]
## Add gene names to the table
glm.all.fdr.annot.spino.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.annot.spino.missense.min)

## Save minimum FDR results from missense SNPs for downstream use
write.table(glm.all.fdr.annot.spino.missense.min, file="rudflies_2023_redo.glm.all.minfdr.spino.missense.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)


## Pull FDR columns only for all genes/features, filtered for missense SNPs only
glm.all.fdr.annot.missense <- glm.all.fdr.annot[glm.all.fdr.annot$Consequence == "missense_variant",]

## Loop through each GLM-results FDR column to calculate minimum values for all genes
glm.all.fdr.annot.missense.min <- c() # create an object
#loop through all cols fdr values for candidate genes
for(g in 5:ncol(glm.all.fdr.annot.missense)) {
  tryCatch({
  	#create 2-column table with gene name and one FDR column at a time
    temp.min <- as.data.frame(cbind(glm.all.fdr.annot.missense$Gene,glm.all.fdr.annot.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "fdr"    
    #assign minimum FDR values
    glm.all.fdr.annot.missense.min <- cbind(glm.all.fdr.annot.missense.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.fdr.annot.missense.min) <- names(glm.all.fdr.annot.missense)[c(5:ncol(glm.all.fdr.annot.missense))]
## Add gene names to the table
glm.all.fdr.annot.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(fdr),
        	list(min = min)) %>%
        	pull(Gene),glm.all.fdr.annot.missense.min)


## Now we've found minimum FDR values for candidate gene and background gene lists, 
## let's run hypergeometric tests on the overlap of candidates and GLM outliers.

## run hypergeometric overlap test per GLM contrast, 20-SNP rolling windows
## First with all gene-related/adjacent SNPS
## Significance threshold FDR<0.05
phyper_all_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.annot.spino.min)) {
  q <- length(unique(glm.all.fdr.annot.spino.min[glm.all.fdr.annot.spino.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.annot.spino.min[,1]))
  n <- length(unique(glm.all.fdr.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.annot.min[glm.all.fdr.annot.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_contrasts05 <- append(phyper_all_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)    
}

## Create results summary table
phyper_spinoCand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.annot.spino.min[,2:ncol(glm.all.fdr.annot.spino.min)]),phyper_pval=phyper_all_contrasts05,phyper_sig=phyper_all_contrasts05 < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_all_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.annot.spino.min)) {
  q <- length(unique(glm.all.fdr.annot.spino.min[glm.all.fdr.annot.spino.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.annot.spino.min[,1]))
  n <- length(unique(glm.all.fdr.annot.min[,1])) - m
  k <- length(unique(glm.all.fdr.annot.min[glm.all.fdr.annot.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_all_contrasts01 <- append(phyper_all_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_spinoCand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.annot.spino.min[,2:ncol(glm.all.fdr.annot.spino.min)]),phyper_pval=phyper_all_contrasts01,phyper_sig=phyper_all_contrasts01 < 0.05,q_all,m_all,k_all,n_all)


### Next with only missense SNPS ###

## Significance threshold FDR<0.05
phyper_missense_contrasts05 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.annot.spino.missense.min)) {
  q <- length(unique(glm.all.fdr.annot.spino.missense.min[glm.all.fdr.annot.spino.missense.min[,contrast] < 0.05,1]))
  m <- length(unique(glm.all.fdr.annot.spino.missense.min[,1]))
  n <- length(unique(glm.all.fdr.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.annot.missense.min[glm.all.fdr.annot.missense.min[,contrast] < 0.05,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_contrasts05 <- append(phyper_missense_contrasts05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_spinoCand_FDR05 <- cbind(contrast=colnames(glm.all.fdr.annot.spino.missense.min[,2:ncol(glm.all.fdr.annot.spino.missense.min)]),phyper_pval=phyper_missense_contrasts05,phyper_sig=phyper_missense_contrasts05 < 0.05,q_all,m_all,k_all,n_all)

## Significance threshold FDR<0.01
phyper_missense_contrasts01 <- c() #create empty object
q_all <- c() #candidate - outlier overlap gene count
m_all <- c() #candidate gene count
n_all <- c() #non-candidate gene count
k_all <- c() #outlier gene count
## Loop through all GLM columns
for(contrast in 2:ncol(glm.all.fdr.annot.spino.missense.min)) {
  q <- length(unique(glm.all.fdr.annot.spino.missense.min[glm.all.fdr.annot.spino.missense.min[,contrast] < 0.01,1]))
  m <- length(unique(glm.all.fdr.annot.spino.missense.min[,1]))
  n <- length(unique(glm.all.fdr.annot.missense.min[,1])) - m
  k <- length(unique(glm.all.fdr.annot.missense.min[glm.all.fdr.annot.missense.min[,contrast] < 0.01,1]))
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  phyper_missense_contrasts01 <- append(phyper_missense_contrasts01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  ## save results as we go
  q_all <- append(q_all,q)
  m_all <- append(m_all,m)
  n_all <- append(n_all,n)
  k_all <- append(k_all,k)  
}

## Create results summary table
phyper_missense_spinoCand_FDR01 <- cbind(contrast=colnames(glm.all.fdr.annot.spino.missense.min[,2:ncol(glm.all.fdr.annot.spino.missense.min)]),phyper_pval=phyper_missense_contrasts01,phyper_sig=phyper_missense_contrasts01 < 0.05,q_all,m_all,k_all,n_all)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.05)
write.table(phyper_spinoCand_FDR05, file="rudflies_2023_redo.phyper_allSNPs_spinoCand_FDR05.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.01)
write.table(phyper_spinoCand_FDR01, file="rudflies_2023_redo.phyper_allSNPs_spinoCand_FDR01.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.05), missense only
write.table(phyper_missense_spinoCand_FDR05, file="rudflies_2023_redo.phyper_missenseSNPs_spinoCand_FDR05.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)

## Save results table for spinosyn candidates vs GLM outliers (FDR < 0.01), missense only
write.table(phyper_missense_spinoCand_FDR01, file="rudflies_2023_redo.phyper_missenseSNPs_spinoCand_FDR01.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)

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
glm.all.pval.annot.spino <- unique(cbind(glm.all.annot.spino[,c(1:3)],Consequence=glm.all.annot.spino$Consequence,select(glm.all.annot.spino,contains("rolwin20"),-contains("fdr"),-contains("logp"))))

## Loop through each GLM-results pval column to calculate spino candidate minimum values
glm.all.pval.annot.spino.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.annot.spino)) {
  tryCatch({
  	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.annot.spino$Gene,glm.all.pval.annot.spino[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum p-values    
    glm.all.pval.annot.spino.min <- cbind(glm.all.pval.annot.spino.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.annot.spino.min) <- names(glm.all.pval.annot.spino)[c(5:ncol(glm.all.pval.annot.spino))]
## Add gene names to the table
glm.all.pval.annot.spino.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.annot.spino.min)
        	
## Save minimum pval results for downstream use
write.table(glm.all.pval.annot.spino.min, file="rudflies_2023_redo.glm.all.minpval.spino.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)
        	
        	
## Pull pval columns only for all genes/features
glm.all.pval.annot <- unique(cbind(Gene=glm.all.annot$Gene,glm.all.annot[,c(1,2)],Consequence=glm.all.annot$Consequence,select(glm.all.annot,contains("rolwin20"),-contains("fdr"),-contains("logp"))))
glm.all.pval.annot <- na.omit(glm.all.pval.annot)

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.annot.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.annot)) {
  tryCatch({
    #create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.annot$Gene,glm.all.pval.annot[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"  
    #assign minimum p-values  
    glm.all.pval.annot.min <- cbind(glm.all.pval.annot.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.annot.min) <- names(glm.all.pval.annot)[c(5:ncol(glm.all.pval.annot))]
## Add gene names to the table
glm.all.pval.annot.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.annot.min)

## Calculate medians instead for plotting
glm.all.pval.annot.med <- unique(glm.all.pval.annot[,1])
glm.all.pval.annot.med <- glm.all.pval.annot.med[order(glm.all.pval.annot.med) ]
#loop through all cols p-values for candidate genes
for(g in c(3,5:ncol(glm.all.pval.annot))) {
  tryCatch({
    temp.med <- as.data.frame(cbind(glm.all.pval.annot$Gene,glm.all.pval.annot[,g]))
    colnames(temp.med)[1] <- "Gene" 
    colnames(temp.med)[2] <- "pval" 
    #assign median p-values
    glm.all.pval.annot.med <- cbind(glm.all.pval.annot.med,
    	temp.med %>%
			group_by(Gene) %>% 
			summarise(median = median(as.numeric(pval), na.rm = TRUE)) %>%
        	pull(median))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.annot.med)[c(2:ncol(glm.all.pval.annot.med))] <- names(glm.all.pval.annot)[c(3,5:ncol(glm.all.pval.annot))]
colnames(glm.all.pval.annot.med)[1]="Gene"

## Add chromosome info to table for plotting
chrom.list <- unique(glm.all.pval.annot[,c(1:2)])
chrom.list <- chrom.list[order(chrom.list$Gene),]
chrom.list <- chrom.list[-c(1:5),]
glm.all.pval.annot.med <- merge(chrom.list, glm.all.pval.annot.med, by="Gene")
        	
        	
## Pull pval columns only for all genes/features, filtered for missense SNPs only
glm.all.pval.annot.spino.missense <- unique(cbind(glm.all.annot.spino.missense[,c(1:3)],glm.all.annot.spino.missense$Consequence,select(glm.all.annot.spino.missense,contains("rolwin20"),-contains("fdr"),-contains("logp"))))

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.annot.spino.missense.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.annot.spino.missense)) { #nrow(haf.freq.cand)
  tryCatch({
	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.annot.spino.missense$Gene,glm.all.pval.annot.spino.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"
    #assign minimum FDR values    
    glm.all.pval.annot.spino.missense.min <- cbind(glm.all.pval.annot.spino.missense.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.annot.spino.missense.min) <- names(glm.all.pval.annot.spino.missense)[c(5:ncol(glm.all.pval.annot.spino.missense))]
## Add gene names to the table
glm.all.pval.annot.spino.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.annot.spino.missense.min)


## Save minimum pval results from missense SNPs for downstream use
write.table(glm.all.pval.annot.spino.missense.min, file="rudflies_2023_redo.glm.all.minpval.spino.missense.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)


## Pull pval columns only for all genes/features, filtered for missense SNPs only
glm.all.pval.annot.missense <- glm.all.pval.annot[glm.all.pval.annot$Consequence == "missense_variant",]

## Loop through each GLM-results pval column to calculate minimum values for all genes
glm.all.pval.annot.missense.min <- c() # create an object
#loop through all cols p-values for candidate genes
for(g in 5:ncol(glm.all.pval.annot.missense)) { #nrow(haf.freq.cand)
  tryCatch({
  	#create 2-column table with gene name and one pval column at a time
    temp.min <- as.data.frame(cbind(glm.all.pval.annot.missense$Gene,glm.all.pval.annot.missense[,g]))
    #name columns
    colnames(temp.min)[1] <- "Gene" 
    colnames(temp.min)[2] <- "pval"   
    #assign minimum p-values 
    glm.all.pval.annot.missense.min <- cbind(glm.all.pval.annot.missense.min,
    	as.numeric(temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(min)))
  }, error=function(e){})
}
## Rename columns of new table
colnames(glm.all.pval.annot.missense.min) <- names(glm.all.pval.annot.missense)[c(5:ncol(glm.all.pval.annot.missense))]
## Add gene names to the table
glm.all.pval.annot.missense.min <- cbind(Gene=temp.min %>%
  			group_by(Gene) %>%
  			summarise_at(vars(pval),
        	list(min = min)) %>%
        	pull(Gene),glm.all.pval.annot.missense.min)
        	
        	


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
glm.all.pval.annot.spino.min.gff <- merge(snp.gff.overlap, glm.all.pval.annot.spino.min, by.x="ID", by.y="Gene")

##Find gene characteristics of background genes from GFF for matching
glm.all.pval.annot.min.gff <- merge(snp.gff.overlap, glm.all.pval.annot.min, by.x="ID", by.y="Gene")

## Pull spino candidate gene IDs
spino_IDs <- glm.all.pval.annot.spino.min.gff$ID

## Remove spino candidates from potential background gene list
glm.all.pval.annot.bg.min.gff <- glm.all.pval.annot.min.gff[!(glm.all.pval.annot.min.gff$ID %in% spino_IDs),]

## Create empty results table with 1K rows (for each iteration) and one column for 
## each GLM contrast.
wilcox.p <- data.frame(matrix(NA, nrow = 1000, ncol = ncol(glm.all.pval.annot.spino.min.gff)-ncol(snp.gff.overlap)))
## Modify column names to reflect GLM contrasts
names(wilcox.p) <- names(glm.all.pval.annot.spino.min.gff)[c(9:ncol(glm.all.pval.annot.spino.min.gff))]

#generate candidate gene filtering tables
for(k in 1:1000) { #set number of iterations
  #Establish background non-candidate set
  bg.samp.list <- c()
  #This loop finds a random set of matched background genes
  for(i in c(1:nrow(glm.all.pval.annot.spino.min.gff))) {
  	#Select matched lists based on following criteria
	length <- glm.all.pval.annot.spino.min.gff[i,]$LENGTH
	count <- glm.all.pval.annot.spino.min.gff[i,]$SNP_COUNT
	type <- glm.all.pval.annot.spino.min.gff[i,]$TYPE
	chrom <- glm.all.pval.annot.spino.min.gff[i,]$CHROM
	#Pull one match per candidate gene
	bg.samp.temp <- sample_n(glm.all.pval.annot.bg.min.gff %>%
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
  for(j in 1:(ncol(glm.all.pval.annot.spino.min.gff)-ncol(snp.gff.overlap))) { #set number of iterations to the number of contrasts
    #Establish Rank Sum test sets for p values
    #candidate set
    candidates <- glm.all.pval.annot.spino.min.gff[,8+j]  
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
    wilcox.p[k,j] <- wilcox.test(candidates.samp, bg.samp, alternative = "less")$p.value
}
}

## Save all results
# Wilcoxon rank sum p-values for all iterations
write.table(wilcox.p, file="rudflies_2023_redo.glm.all.wilcox.p.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)

# Median Wilcoxon rank sum p-values across iterations
write.table(colMedians(as.matrix(wilcox.p)), file="rudflies_2023_redo.glm.all.wilcox.p.medians.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)

# Mean Wilcoxon rank sum p-values across iterations
write.table(colMeans(as.matrix(wilcox.p)), file="rudflies_2023_redo.glm.all.wilcox.p.means.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)

## Calculate the percent of p-values significant at a raw p<0.05
colMedians(as.matrix(wilcox.p))
colMeans(as.matrix(wilcox.p))
wilcox.contrast <- c()
wilcox.contrast.p05 <- wilcox.p<0.05
for(l in 1:ncol(wilcox.p)) {
   wilcox.contrast <- append(wilcox.contrast,sum(wilcox.contrast.p05[,l], na.rm=TRUE)*0.1)
}
cbind(contrast=names(wilcox.p),perc_sig=wilcox.contrast)

# Save percent of significant Wilcoxon rank sum p-values across iterations
write.table(cbind(contrast=names(wilcox.p),perc_sig=wilcox.contrast), file="rudflies_2023_redo.glm.all.wilcox.p.percsig.norolwin.txt", quote = FALSE, sep = "\t", row.names = F)


