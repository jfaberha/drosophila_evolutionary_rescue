##########################################################################################
### This script is associated with a study on the evolution of spinosad-resistance in  ###
### field populations of Drosophila melanogaster run in 2023 by the Rudman lab (WSUV). ###
###                                                                                    ###
### This R script was run to summarize Watterson's Theta diversity metrics as          ###
### calculated in Grenedalf using the halfpipe-imputed frequency table. The overall    ###
### goal was to check if major allele frequency differences between treatments are     ###
### a consequence of losses in diversity from genetic drift or bottleneck effects due  ###
### to strong selection from spinosad exposure.                                        ###
###                                                                                    ###
### Since there is a range of effective coverages after hafpipe imputation, we will    ###
### analyze diversity at both the high end (3500) and low end (500) of the ranges, as  ###
### estimated here:                                                                    ###
### https://ec-calculator.shinyapps.io/shinyapp/                                       ###
##########################################################################################

### In R ###
#configure r environment
setwd("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/main_analysis")
library(ggplot2)
library(matrixStats)
library(dplyr)
library(ggpubr)
library(emmeans)

#set treatment group color scheme right away for plotting
#order: E, PA, SE, SP
sample_cols <- c("#D26183","#495184","#848556","#D9B851")

#load required input file
haf.meta <- read.table("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/rudflies_2023_meta.tsv", header=TRUE)

#filter for founder samples for PA, S, and E
haf.meta.filt <- haf.meta[(haf.meta$batch == "a" | haf.meta$batch == "F") & (haf.meta$treat.fix == "PA" | haf.meta$treat.fix == "S" | haf.meta$treat.fix == "E") & haf.meta$experiment == "spino" ,-12]

## A couple commands to label SE (extinct) and SP (persistent) in a new condition column
haf.meta.filt$condition <- gsub("S","SP",haf.meta.filt$treat.fix)
# List of SE cages to filter by
se_cages <- c("11","21","27","41","45")
# Relabel SE cages
haf.meta.filt[haf.meta.filt$cage %in% se_cages,]$condition = "SE"

#tables need a bit of reformatting for glm to run
haf.meta.filt$treat.fix <- as.factor(haf.meta.filt$treat.fix)
haf.meta.filt$tpt <- as.factor(haf.meta.filt$tpt)


##################################
### Start with coverage ~ 3500 ###
##################################

##############################################
### First, analyze diversity by chromosome ###
##############################################

#The following table was generated in grenedalf with the following command:
#grenedalf diversity --frequency-table-path ../hafpipe/all.csv.gz --window-type chromosomes --pool-sizes 240 --filter-sample-min-count 2 --frequency-table-int-factor 3500 --measure theta-watterson --allow-file-overwriting --out-dir ./
rudflies_2023_diversity <- read.csv("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/diversity_chrom_eff3500.csv", header=TRUE)

## Calculate relative theta
## First, pull theta_watterson_rel stats from table
rudflies_2023_theta_rel <- rudflies_2023_diversity %>%
  select(matches("theta_watterson_rel"))
## headers need some relabeling to match sample names.
names(rudflies_2023_theta_rel) <- gsub("X","",names(rudflies_2023_theta_rel))
names(rudflies_2023_theta_rel) <- gsub(".theta_watterson_rel","",names(rudflies_2023_theta_rel))
## Since there's one row for each chromsome, we'll set them as row names
row.names(rudflies_2023_theta_rel) <- rudflies_2023_diversity[,1]
## Transpose to have chromosomes as columns
rudflies_2023_theta_rel <- t(rudflies_2023_theta_rel)
## Include sample names for merging
rudflies_2023_theta_rel <- cbind(samp=row.names(rudflies_2023_theta_rel),rudflies_2023_theta_rel)
## Reformat sample names and merge stats table with metadata for filtering
rudflies_2023_theta_rel[,1] <- as.factor(rudflies_2023_theta_rel[,1])
rudflies_2023_theta_rel_meta <- merge(haf.meta.filt,rudflies_2023_theta_rel,by="samp")
## Reformat columns for dplyr stats
rudflies_2023_theta_rel_meta$condition <- as.factor(rudflies_2023_theta_rel_meta$condition)
rudflies_2023_theta_rel_meta$"2L" <- as.numeric(rudflies_2023_theta_rel_meta$"2L")
rudflies_2023_theta_rel_meta$"2R" <- as.numeric(rudflies_2023_theta_rel_meta$"2R")
rudflies_2023_theta_rel_meta$"3L" <- as.numeric(rudflies_2023_theta_rel_meta$"3L")
rudflies_2023_theta_rel_meta$"3R" <- as.numeric(rudflies_2023_theta_rel_meta$"3R")
rudflies_2023_theta_rel_meta$"X" <- as.numeric(rudflies_2023_theta_rel_meta$"X")

## Now, calculate mean theta by timepoint, condition, and chromosome
rudflies_2023_theta_rel_mean <- rudflies_2023_theta_rel_meta %>%
  group_by(condition, tpt) %>% 
  summarise_at(vars(c("2L","2R","3L","3R","X")), mean)

## Reformat tibble to dataframe
rudflies_2023_theta_rel_mean <- as.data.frame(rudflies_2023_theta_rel_mean)
## Rename chrosomes so they are plotted properly in ggplot
rudflies_2023_theta_rel_mean_v2 <- rbind(cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chr2L",theta=rudflies_2023_theta_rel_mean[,3]),
cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chr2R",theta=rudflies_2023_theta_rel_mean[,4]),
cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chr3L",theta=rudflies_2023_theta_rel_mean[,5]),
cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chr3R",theta=rudflies_2023_theta_rel_mean[,6]),
cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chrX",theta=rudflies_2023_theta_rel_mean[,7]))
## Reformat columns for ggplot
rudflies_2023_theta_rel_mean_v2$tpt <- as.factor(rudflies_2023_theta_rel_mean_v2$tpt)
rudflies_2023_theta_rel_mean_v2$condition <- as.factor(rudflies_2023_theta_rel_mean_v2$condition)

#lineplot (all conditions over time, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel_mean.line.pdf", width=6, height=12)
	ggplot(rudflies_2023_theta_rel_mean_v2,aes(x=tpt,y=theta , group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

#boxplot (all conditions over time, all chromosomes combined)
pdf(file = "rudflies_2023_redo_theta_rel_mean.box.pdf")
		ggplot(rudflies_2023_theta_rel_mean_v2,aes(x=tpt,y=theta, fill=condition)) +
		geom_boxplot() +
		theme_classic2() +
		scale_fill_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Mutate the clunky way so all theta stats are in a single column for alternate plotting
rudflies_2023_theta_rel_meta <- as.data.frame(rudflies_2023_theta_rel_meta)
rudflies_2023_theta_rel_meta_v2 <- rbind(cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chr2L",theta=rudflies_2023_theta_rel_meta[,13]),
cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chr2R",theta=rudflies_2023_theta_rel_meta[,14]),
cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chr3L",theta=rudflies_2023_theta_rel_meta[,15]),
cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chr3R",theta=rudflies_2023_theta_rel_meta[,16]),
cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chrX",theta=rudflies_2023_theta_rel_meta[,18]))
## Reformat columns for ggplot
rudflies_2023_theta_rel_meta_v2$tpt <- as.factor(rudflies_2023_theta_rel_meta_v2$tpt)
rudflies_2023_theta_rel_meta_v2$condition <- as.factor(rudflies_2023_theta_rel_meta_v2$condition)

#boxplot (all conditions, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel.box.pdf", width=6, height=12)
	ggplot(rudflies_2023_theta_rel_meta_v2,aes(x=tpt,y=theta , group=interaction(tpt, condition), fill=condition)) +
		geom_boxplot() +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2() +
		scale_fill_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

#boxplot 2 (only SE and SP, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel2.box.pdf", width=6, height=3)
	ggplot(rudflies_2023_theta_rel_meta_v2[rudflies_2023_theta_rel_meta_v2$tpt=="1" & (rudflies_2023_theta_rel_meta_v2$condition=="SP" | rudflies_2023_theta_rel_meta_v2$condition=="SE"),],aes(x=chrom,y=theta , group=interaction(chrom, condition), fill=condition)) +
		geom_boxplot() +
		theme_classic2() +
		scale_fill_manual(values = sample_cols[c(3,4)]) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Because 3L and 3R show the most extreme signals of adaptation

#line (chromosome 3L only)
pdf(file = "rudflies_2023_redo_theta_rel_mean_3L.line.pdf")
	ggplot(rudflies_2023_theta_rel_mean_v2[rudflies_2023_theta_rel_mean_v2$chrom=="chr3L",],aes(x=tpt,y=theta, group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

#line (chromosome 3R only)
pdf(file = "rudflies_2023_redo_theta_rel_mean_3R.line.pdf")
	ggplot(rudflies_2023_theta_rel_mean_v2[rudflies_2023_theta_rel_mean_v2$chrom=="chr3R",],aes(x=tpt,y=theta, group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Save summary table 
write.table(rudflies_2023_theta_rel_mean, file="rudflies_2023_theta_rel_mean_byChrom.txt",sep = "\t", quote = FALSE, row.names = F)
## to reload
#rudflies_2023_theta_rel_mean <- read.table("rudflies_2023_theta_rel_mean_byChrom.txt", header=TRUE)

########################################################
### Second, analyze diversity across the full genome ###
########################################################

#The following table was generated in grenedalf with the following command:
#grenedalf diversity --frequency-table-path ../hafpipe/all.csv.gz --window-type genome --pool-sizes 240 --filter-sample-min-count 2 --frequency-table-int-factor 3500 --measure theta-watterson --allow-file-overwriting --out-dir ./
rudflies_2023_diversity_genome <- read.csv("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/diversity_genome_eff3500.csv", header=TRUE)

## Calculate relative theta
## First, pull theta_watterson_rel stats from table
rudflies_2023_theta_rel_genome <- rudflies_2023_diversity_genome %>%
  select(matches("theta_watterson_rel"))
## headers need some relabeling to match sample names.  
names(rudflies_2023_theta_rel_genome) <- gsub("X","",names(rudflies_2023_theta_rel_genome))
names(rudflies_2023_theta_rel_genome) <- gsub(".theta_watterson_rel","",names(rudflies_2023_theta_rel_genome))
## Since there's one row for each chromsome, we'll set them as row names
row.names(rudflies_2023_theta_rel_genome) <- "genome"
## Transpose to have chromosomes as columns
rudflies_2023_theta_rel_genome <- t(rudflies_2023_theta_rel_genome)
## Include sample names for merging
rudflies_2023_theta_rel_genome <- cbind(samp=row.names(rudflies_2023_theta_rel_genome),rudflies_2023_theta_rel_genome)
## Reformat sample names and merge stats table with metadata for filtering
rudflies_2023_theta_rel_genome[,1] <- as.factor(rudflies_2023_theta_rel_genome[,1])
rudflies_2023_theta_rel_genome_meta <- merge(haf.meta.filt,rudflies_2023_theta_rel_genome,by="samp")
rudflies_2023_theta_rel_genome_meta$condition <- as.factor(rudflies_2023_theta_rel_genome_meta$condition)
rudflies_2023_theta_rel_genome_meta$"genome" <- as.numeric(rudflies_2023_theta_rel_genome_meta$"genome")
names(rudflies_2023_theta_rel_genome_meta)[ncol(rudflies_2023_theta_rel_genome_meta)] <- "theta"

## Now, calculate mean theta by timepoint and condition
rudflies_2023_theta_rel_genome_mean <- rudflies_2023_theta_rel_genome_meta %>%
  group_by(condition, tpt) %>% 
  summarise_at(vars(theta), mean)

## Reformat tibble to dataframe
rudflies_2023_theta_rel_genome_mean <- as.data.frame(rudflies_2023_theta_rel_genome_mean)
## Reformat columns for ggplot
rudflies_2023_theta_rel_genome_mean$tpt <- as.factor(rudflies_2023_theta_rel_genome_mean$tpt)
rudflies_2023_theta_rel_genome_mean$condition <- as.factor(rudflies_2023_theta_rel_genome_mean$condition)

#lineplot (all conditions over time)
pdf(file = "rudflies_2023_redo_theta_rel_genome_mean.line.pdf", width=6, height=3)
	ggplot(rudflies_2023_theta_rel_genome_mean,aes(x=tpt,y=theta , group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Add a chromosome column and label full genome theta as "all"
rudflies_2023_theta_rel_genome_mean <- cbind(rudflies_2023_theta_rel_genome_mean[,c(1:2)],chrom="all",theta=rudflies_2023_theta_rel_genome_mean[,3])

#lineplot (all conditions over time, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel_mean_byChrom.line.pdf", width=6, height=15)
	ggplot(rbind(rudflies_2023_theta_rel_genome_mean, rudflies_2023_theta_rel_mean_v2),aes(x=tpt,y=theta , group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2()  +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Reformat metadata columns for ggplot
rudflies_2023_theta_rel_genome_meta$tpt <- as.character(rudflies_2023_theta_rel_genome_meta$tpt)
rudflies_2023_theta_rel_genome_meta$condition <- as.factor(rudflies_2023_theta_rel_genome_meta$condition)

#boxplot (all conditions over time, full genome)
pdf(file = "rudflies_2023_redo_theta_rel_genome.box.pdf", width=6, height=3)
	ggplot(rudflies_2023_theta_rel_genome_meta,aes(x=tpt,y=theta , group=interaction(tpt, condition), fill=condition)) +
		geom_boxplot() +
		theme_classic2() +
		scale_fill_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Append previous chromosome-level theta table with new genome-level table so we can plot 
## them all together
rudflies_2023_theta_rel_meta_v3 <- rbind(rudflies_2023_theta_rel_meta_v2,
cbind(rudflies_2023_theta_rel_genome_meta[,c(4,12)],chrom="all",theta=rudflies_2023_theta_rel_genome_meta$theta))

## Reformat metadata columns for ggplot
rudflies_2023_theta_rel_meta_v3$tpt <- as.factor(rudflies_2023_theta_rel_meta_v3$tpt)
rudflies_2023_theta_rel_meta_v3$condition <- as.factor(rudflies_2023_theta_rel_meta_v3$condition)

#boxplot (all conditions over time, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel_byChrom.box.pdf", width=6, height=15)
	ggplot(rudflies_2023_theta_rel_meta_v3,aes(x=tpt,y=theta , group=interaction(tpt, condition), fill=condition)) +
		geom_boxplot() +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2() +
		scale_fill_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()


#lineplot (all conditions over time, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel_mean_byChrom.line.pdf", width=6, height=15)
	ggplot(rbind(rudflies_2023_theta_rel_genome_mean, rudflies_2023_theta_rel_mean_v2),aes(x=tpt,y=theta , group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Save summary table 
write.table(rudflies_2023_theta_rel_genome_mean, file="rudflies_2023_theta_rel_genome_mean.txt",sep = "\t", quote = FALSE, row.names = F)
## To reload
#rudflies_2023_theta_rel_genome_mean <- read.table("rudflies_2023_theta_rel_genome_mean.txt", header=TRUE)



### GLMs on Diversity Stats ###
glmres <- c() #initialize glm object
contrastout <- c() #initialize p-val summary table

## Format table for running glms with chromosome- and genome-level theta values
glmdf <- cbind(rudflies_2023_theta_rel_meta, all=rudflies_2023_theta_rel_genome_meta$theta)
## Save this table for possible future use
write.table(glmdf, file="rudflies_2023_theta_rel_eff3500_wMeta.txt",sep = "\t", quote = FALSE, row.names = F)
## Prior to running glms, omit founder samples since they are unbalanced with other TPTs 
glmdf <- glmdf[glmdf$tpt != "0",]

## Now run GLMs
#loop through all rows of allele frequency table
for(g in c(13:19)) {
	glmdf_temp <- glmdf
	glmdf_temp[,g] <- round(as.numeric(glmdf_temp[,g])*100000)
    colnames(glmdf_temp)[g] <- "count"
    glmres[[g]] <- glm(count~condition, data = glmdf_temp, family = quasipoisson, na.action = na.exclude)
    emtemp <- emmeans(glmres[[g]], pairwise ~ condition)$contrasts
    #summary(emtemp)$p.value
    contrastout <- rbind(contrastout,summary(emtemp)$p.value)
}

## Rename rows and columns in GLM p-value table
colnames(contrastout) <- c("EvPA", "EvSP", "EvSE", "PAvSP", "PAvSE", "SEvSP")
rownames(contrastout) <- colnames(glmdf[,c(13:19)])
contrastout
#         EvPA      EvSP      EvSE     PAvSP     PAvSE     SEvSP
#2L  0.9651466 0.9999547 0.9651397 0.9930501 0.9997732 0.9902793
#2R  0.3107884 0.9432666 0.9999412 0.9950388 0.4489914 0.9423222
#3L  0.9753962 0.2321694 0.9999581 0.3233863 0.9791207 0.2597329
#3R  0.9548143 0.8795311 0.9999889 0.9581894 0.9677062 0.8868580
#4   0.7319415 0.9636141 0.8332761 1.0000000 1.0000000 1.0000000
#X   0.9998873 0.6713480 0.9907927 0.6891574 0.9844407 0.6025968
#all 0.8838644 0.7806106 0.9999669 0.9304013 0.9437276 0.8171487

## Save GLM p-value table 
write.table(glmdf, file="rudflies_2023_theta_rel_eff3500_wMeta.txt",sep = "\t", quote = FALSE, row.names = F)


##################################
### Next, check coverage ~ 500 ###
##################################

##############################################
### First, analyze diversity by chromosome ###
##############################################

#The following table was generated in grenedalf with the following command:
#grenedalf diversity --frequency-table-path ../hafpipe/all.csv.gz --window-type chromosomes --pool-sizes 240 --filter-sample-min-count 2 --frequency-table-int-factor 500 --measure theta-watterson --allow-file-overwriting --out-dir ./
rudflies_2023_diversity <- read.csv("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/diversity_chrom_eff500.csv", header=TRUE)

## Calculate relative theta
## First, pull theta_watterson_rel stats from table
rudflies_2023_theta_rel <- rudflies_2023_diversity %>%
  select(matches("theta_watterson_rel"))
## headers need some relabeling to match sample names.
names(rudflies_2023_theta_rel) <- gsub("X","",names(rudflies_2023_theta_rel))
names(rudflies_2023_theta_rel) <- gsub(".theta_watterson_rel","",names(rudflies_2023_theta_rel))
## Since there's one row for each chromsome, we'll set them as row names
row.names(rudflies_2023_theta_rel) <- rudflies_2023_diversity[,1]
## Transpose to have chromosomes as columns
rudflies_2023_theta_rel <- t(rudflies_2023_theta_rel)
## Include sample names for merging
rudflies_2023_theta_rel <- cbind(samp=row.names(rudflies_2023_theta_rel),rudflies_2023_theta_rel)
## Reformat sample names and merge stats table with metadata for filtering
rudflies_2023_theta_rel[,1] <- as.factor(rudflies_2023_theta_rel[,1])
rudflies_2023_theta_rel_meta <- merge(haf.meta.filt,rudflies_2023_theta_rel,by="samp")
## Reformat columns for dplyr stats
rudflies_2023_theta_rel_meta$condition <- as.factor(rudflies_2023_theta_rel_meta$condition)
rudflies_2023_theta_rel_meta$"2L" <- as.numeric(rudflies_2023_theta_rel_meta$"2L")
rudflies_2023_theta_rel_meta$"2R" <- as.numeric(rudflies_2023_theta_rel_meta$"2R")
rudflies_2023_theta_rel_meta$"3L" <- as.numeric(rudflies_2023_theta_rel_meta$"3L")
rudflies_2023_theta_rel_meta$"3R" <- as.numeric(rudflies_2023_theta_rel_meta$"3R")
rudflies_2023_theta_rel_meta$"X" <- as.numeric(rudflies_2023_theta_rel_meta$"X")

## Now, calculate mean theta by timepoint, condition, and chromosome
rudflies_2023_theta_rel_mean <- rudflies_2023_theta_rel_meta %>%
  group_by(condition, tpt) %>% 
  summarise_at(vars(c("2L","2R","3L","3R","X")), mean)

## Reformat tibble to dataframe
rudflies_2023_theta_rel_mean <- as.data.frame(rudflies_2023_theta_rel_mean)
## Rename chrosomes so they are plotted properly in ggplot
rudflies_2023_theta_rel_mean_v2 <- rbind(cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chr2L",theta=rudflies_2023_theta_rel_mean[,3]),
cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chr2R",theta=rudflies_2023_theta_rel_mean[,4]),
cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chr3L",theta=rudflies_2023_theta_rel_mean[,5]),
cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chr3R",theta=rudflies_2023_theta_rel_mean[,6]),
cbind(rudflies_2023_theta_rel_mean[,c(1,2)],chrom="chrX",theta=rudflies_2023_theta_rel_mean[,7]))
## Reformat columns for ggplot
rudflies_2023_theta_rel_mean_v2$tpt <- as.factor(rudflies_2023_theta_rel_mean_v2$tpt)
rudflies_2023_theta_rel_mean_v2$condition <- as.factor(rudflies_2023_theta_rel_mean_v2$condition)

#lineplot (all conditions over time, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel_mean.line.pdf", width=6, height=12)
	ggplot(rudflies_2023_theta_rel_mean_v2,aes(x=tpt,y=theta , group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

#boxplot (all conditions over time, all chromosomes combined)
pdf(file = "rudflies_2023_redo_theta_rel_mean.box.pdf")
		ggplot(rudflies_2023_theta_rel_mean_v2,aes(x=tpt,y=theta, fill=condition)) +
		geom_boxplot() +
		theme_classic2() +
		scale_fill_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Mutate the clunky way so all theta stats are in a single column for alternate plotting
rudflies_2023_theta_rel_meta <- as.data.frame(rudflies_2023_theta_rel_meta)
rudflies_2023_theta_rel_meta_v2 <- rbind(cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chr2L",theta=rudflies_2023_theta_rel_meta[,13]),
cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chr2R",theta=rudflies_2023_theta_rel_meta[,14]),
cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chr3L",theta=rudflies_2023_theta_rel_meta[,15]),
cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chr3R",theta=rudflies_2023_theta_rel_meta[,16]),
cbind(rudflies_2023_theta_rel_meta[,c(4,12)],chrom="chrX",theta=rudflies_2023_theta_rel_meta[,18]))
## Reformat columns for ggplot
rudflies_2023_theta_rel_meta_v2$tpt <- as.factor(rudflies_2023_theta_rel_meta_v2$tpt)
rudflies_2023_theta_rel_meta_v2$condition <- as.factor(rudflies_2023_theta_rel_meta_v2$condition)

#boxplot (all conditions, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel.box.pdf", width=6, height=12)
	ggplot(rudflies_2023_theta_rel_meta_v2,aes(x=tpt,y=theta , group=interaction(tpt, condition), fill=condition)) +
		geom_boxplot() +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2() +
		scale_fill_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

#boxplot 2 (only SE and SP, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel2.box.pdf", width=6, height=3)
	ggplot(rudflies_2023_theta_rel_meta_v2[rudflies_2023_theta_rel_meta_v2$tpt=="1" & (rudflies_2023_theta_rel_meta_v2$condition=="SP" | rudflies_2023_theta_rel_meta_v2$condition=="SE"),],aes(x=chrom,y=theta , group=interaction(chrom, condition), fill=condition)) +
		geom_boxplot() +
		theme_classic2() +
		scale_fill_manual(values = sample_cols[c(3,4)]) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Because 3L and 3R show the most extreme signals of adaptation

#line (chromosome 3L only)
pdf(file = "rudflies_2023_redo_theta_rel_mean_3L.line.pdf")
	ggplot(rudflies_2023_theta_rel_mean_v2[rudflies_2023_theta_rel_mean_v2$chrom=="chr3L",],aes(x=tpt,y=theta, group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

#line (chromosome 3R only)
pdf(file = "rudflies_2023_redo_theta_rel_mean_3R.line.pdf")
	ggplot(rudflies_2023_theta_rel_mean_v2[rudflies_2023_theta_rel_mean_v2$chrom=="chr3R",],aes(x=tpt,y=theta, group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Save summary table 
write.table(rudflies_2023_theta_rel_mean, file="rudflies_2023_theta_rel_mean_byChrom.txt",sep = "\t", quote = FALSE, row.names = F)
## to reload
#rudflies_2023_theta_rel_mean <- read.table("rudflies_2023_theta_rel_mean_byChrom.txt", header=TRUE)

########################################################
### Second, analyze diversity across the full genome ###
########################################################

#The following table was generated in grenedalf with the following command:
#grenedalf diversity --frequency-table-path ../hafpipe/all.csv.gz --window-type genome --pool-sizes 240 --filter-sample-min-count 2 --frequency-table-int-factor 500 --measure theta-watterson --allow-file-overwriting --out-dir ./
rudflies_2023_diversity_genome <- read.csv("/data/lab/rudman/gp_analysis/rudflies_2023_redo/r/r_input/diversity_genome_eff500.csv", header=TRUE)

## Calculate relative theta
## First, pull theta_watterson_rel stats from table
rudflies_2023_theta_rel_genome <- rudflies_2023_diversity_genome %>%
  select(matches("theta_watterson_rel"))
## headers need some relabeling to match sample names.  
names(rudflies_2023_theta_rel_genome) <- gsub("X","",names(rudflies_2023_theta_rel_genome))
names(rudflies_2023_theta_rel_genome) <- gsub(".theta_watterson_rel","",names(rudflies_2023_theta_rel_genome))
## Since there's one row for each chromsome, we'll set them as row names
row.names(rudflies_2023_theta_rel_genome) <- "genome"
## Transpose to have chromosomes as columns
rudflies_2023_theta_rel_genome <- t(rudflies_2023_theta_rel_genome)
## Include sample names for merging
rudflies_2023_theta_rel_genome <- cbind(samp=row.names(rudflies_2023_theta_rel_genome),rudflies_2023_theta_rel_genome)
## Reformat sample names and merge stats table with metadata for filtering
rudflies_2023_theta_rel_genome[,1] <- as.factor(rudflies_2023_theta_rel_genome[,1])
rudflies_2023_theta_rel_genome_meta <- merge(haf.meta.filt,rudflies_2023_theta_rel_genome,by="samp")
rudflies_2023_theta_rel_genome_meta$condition <- as.factor(rudflies_2023_theta_rel_genome_meta$condition)
rudflies_2023_theta_rel_genome_meta$"genome" <- as.numeric(rudflies_2023_theta_rel_genome_meta$"genome")
names(rudflies_2023_theta_rel_genome_meta)[ncol(rudflies_2023_theta_rel_genome_meta)] <- "theta"

## Now, calculate mean theta by timepoint and condition
rudflies_2023_theta_rel_genome_mean <- rudflies_2023_theta_rel_genome_meta %>%
  group_by(condition, tpt) %>% 
  summarise_at(vars(theta), mean)

## Reformat tibble to dataframe
rudflies_2023_theta_rel_genome_mean <- as.data.frame(rudflies_2023_theta_rel_genome_mean)
## Reformat columns for ggplot
rudflies_2023_theta_rel_genome_mean$tpt <- as.factor(rudflies_2023_theta_rel_genome_mean$tpt)
rudflies_2023_theta_rel_genome_mean$condition <- as.factor(rudflies_2023_theta_rel_genome_mean$condition)

#lineplot (all conditions over time)
pdf(file = "rudflies_2023_redo_theta_rel_genome_mean.line.pdf", width=6, height=3)
	ggplot(rudflies_2023_theta_rel_genome_mean,aes(x=tpt,y=theta , group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Add a chromosome column and label full genome theta as "all"
rudflies_2023_theta_rel_genome_mean <- cbind(rudflies_2023_theta_rel_genome_mean[,c(1:2)],chrom="all",theta=rudflies_2023_theta_rel_genome_mean[,3])

#lineplot (all conditions over time, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel_mean_byChrom.line.pdf", width=6, height=15)
	ggplot(rbind(rudflies_2023_theta_rel_genome_mean, rudflies_2023_theta_rel_mean_v2),aes(x=tpt,y=theta , group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2()  +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Reformat metadata columns for ggplot
rudflies_2023_theta_rel_genome_meta$tpt <- as.character(rudflies_2023_theta_rel_genome_meta$tpt)
rudflies_2023_theta_rel_genome_meta$condition <- as.factor(rudflies_2023_theta_rel_genome_meta$condition)

#boxplot (all conditions over time, full genome)
pdf(file = "rudflies_2023_redo_theta_rel_genome.box.pdf", width=6, height=3)
	ggplot(rudflies_2023_theta_rel_genome_meta,aes(x=tpt,y=theta , group=interaction(tpt, condition), fill=condition)) +
		geom_boxplot() +
		theme_classic2() +
		scale_fill_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Append previous chromosome-level theta table with new genome-level table so we can plot 
## them all together
rudflies_2023_theta_rel_meta_v3 <- rbind(rudflies_2023_theta_rel_meta_v2,
cbind(rudflies_2023_theta_rel_genome_meta[,c(4,12)],chrom="all",theta=rudflies_2023_theta_rel_genome_meta$theta))

## Reformat metadata columns for ggplot
rudflies_2023_theta_rel_meta_v3$tpt <- as.factor(rudflies_2023_theta_rel_meta_v3$tpt)
rudflies_2023_theta_rel_meta_v3$condition <- as.factor(rudflies_2023_theta_rel_meta_v3$condition)

#boxplot (all conditions over time, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel_byChrom.box.pdf", width=6, height=15)
	ggplot(rudflies_2023_theta_rel_meta_v3,aes(x=tpt,y=theta , group=interaction(tpt, condition), fill=condition)) +
		geom_boxplot() +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2() +
		scale_fill_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()


#lineplot (all conditions over time, by chromosome)
pdf(file = "rudflies_2023_redo_theta_rel_mean_byChrom.line.pdf", width=6, height=15)
	ggplot(rbind(rudflies_2023_theta_rel_genome_mean, rudflies_2023_theta_rel_mean_v2),aes(x=tpt,y=theta , group=condition, color=condition)) +
		geom_line() +
		geom_point(size = 3) +
		facet_wrap(~chrom, ncol = 1) +
		theme_classic2() +
		scale_color_manual(values = sample_cols) +
		ggtitle("Watterson's theta (relative)")
dev.off()

## Save summary table 
write.table(rudflies_2023_theta_rel_genome_mean, file="rudflies_2023_theta_rel_genome_mean.txt",sep = "\t", quote = FALSE, row.names = F)
## To reload
#rudflies_2023_theta_rel_genome_mean <- read.table("rudflies_2023_theta_rel_genome_mean.txt", header=TRUE)



### GLMs on Diversity Stats ###
glmres <- c() #initialize glm object
contrastout <- c() #initialize p-val summary table

## Format table for running glms with chromosome- and genome-level theta values
glmdf <- cbind(rudflies_2023_theta_rel_meta, all=rudflies_2023_theta_rel_genome_meta$theta)
## Save this table for possible future use
write.table(glmdf, file="rudflies_2023_theta_rel_eff500_wMeta.txt",sep = "\t", quote = FALSE, row.names = F)
## Prior to running glms, omit founder samples since they are unbalanced with other TPTs 
glmdf <- glmdf[glmdf$tpt != "0",]

## Now run GLMs
#loop through all rows of allele frequency table
for(g in c(13:19)) {
	glmdf_temp <- glmdf
	glmdf_temp[,g] <- round(as.numeric(glmdf_temp[,g])*100000)
    colnames(glmdf_temp)[g] <- "count"
    glmres[[g]] <- glm(count~condition, data = glmdf_temp, family = quasipoisson, na.action = na.exclude)
    emtemp <- emmeans(glmres[[g]], pairwise ~ condition)$contrasts
    #summary(emtemp)$p.value
    contrastout <- rbind(contrastout,summary(emtemp)$p.value)
}

## Rename rows and columns in GLM p-value table
colnames(contrastout) <- c("EvPA", "EvSP", "EvSE", "PAvSP", "PAvSE", "SEvSP")
rownames(contrastout) <- colnames(glmdf[,c(13:19)])
contrastout
#         EvPA      EvSP      EvSE     PAvSP     PAvSE     SEvSP
#2L  0.9651466 0.9999547 0.9651397 0.9930501 0.9997732 0.9902793
#2R  0.3107884 0.9432666 0.9999412 0.9950388 0.4489914 0.9423222
#3L  0.9753962 0.2321694 0.9999581 0.3233863 0.9791207 0.2597329
#3R  0.9548143 0.8795311 0.9999889 0.9581894 0.9677062 0.8868580
#4   0.7319415 0.9636141 0.8332761 1.0000000 1.0000000 1.0000000
#X   0.9998873 0.6713480 0.9907927 0.6891574 0.9844407 0.6025968
#all 0.8838644 0.7806106 0.9999669 0.9304013 0.9437276 0.8171487

## Save GLM p-value table 
write.table(glmdf, file="rudflies_2023_theta_rel_eff500_wMeta.txt",sep = "\t", quote = FALSE, row.names = F)