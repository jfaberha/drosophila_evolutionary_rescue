########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 11_hafpipe_validation.R — accuracy of HAF-pipe allele frequencies.

###  Produces:
###  Fig. S14 / SI 1.7  for every SNP found in all three sources, the mean reference-allele frequency from (i) the DGRP2 Freeze 2 genotypes of the founder lines, (ii) BCFtools calls from the pooled samples, and (iii) HAF-pipe (mean over all samples; HAF-pipe reports the ALT allele, so 1 - mean is used); Pearson correlations with the DGRP2 values and 2-D density plots (A: BCFtools, B: HAF-pipe).

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
suppressMessages(library(ggplot2))
OUT <- module_dir("11_hafpipe_validation"); setwd(OUT)
sink("11_hafpipe_validation_report.txt", split = TRUE)
FRQ <- list(dgrp2 = "dgrp2_founder_list.frq", bcftools = "filtered-all.frq")

haf <- load_haf(exclude = TRUE)
obs <- c(n_hafpipe_variants = nrow(haf$freq))
hp <- data.frame(haf$freq[, 1:2], hafpipe.means.A = rowMeans(haf$freq[, -(1:2)]))
hp$hafpipe.means.R <- 1 - hp$hafpipe.means.A
msg("HAF-pipe: %s variants; mean over %d samples", format(nrow(hp), big.mark = ","), ncol(haf$freq) - 2)
rm(haf); invisible(gc())

rd <- function(f) { d <- read.delim(file.path(PATHS$raw_input, f), header = TRUE, sep = "\t"); d$CHROM <- as.character(d$CHROM); d }
dg <- rd(FRQ$dgrp2); vc <- rd(FRQ$bcftools)
stopifnot("R.FREQ" %in% names(dg), "R.FREQ" %in% names(vc))
all.freq <- merge(dg, vc, by = c("CHROM", "POS"))          # R.FREQ.x = DGRP2, R.FREQ.y = BCFtools
all.freq <- merge(all.freq, hp, by = c("CHROM", "POS"))
msg("SNPs in all three sources: %s", format(nrow(all.freq), big.mark = ","))

r_bcf <- cor.test(as.numeric(all.freq$R.FREQ.x), as.numeric(all.freq$R.FREQ.y), method = "pearson")
r_hp  <- cor.test(as.numeric(all.freq$R.FREQ.x), as.numeric(all.freq$hafpipe.means.R), method = "pearson")
r_bh  <- cor.test(as.numeric(all.freq$R.FREQ.y), as.numeric(all.freq$hafpipe.means.R), method = "pearson")
res <- data.frame(comparison = c("DGRP2 vs BCFtools (Fig. S14A)", "DGRP2 vs HAF-pipe (Fig. S14B)", "BCFtools vs HAF-pipe (not in the paper)"),
	n_SNPs = nrow(all.freq), pearson_R = c(r_bcf$estimate, r_hp$estimate, r_bh$estimate), P = c(r_bcf$p.value, r_hp$p.value, r_bh$p.value))
write.csv(res, "FigS14_correlations.csv", row.names = FALSE); print(res, row.names = FALSE)
obs["R_dgrp2_vs_bcftools"] <- unname(r_bcf$estimate); obs["R_dgrp2_vs_hafpipe"] <- unname(r_hp$estimate)

dens <- function(y, ylab, tag) ggplot(all.freq, aes(x = R.FREQ.x, y = .data[[y]])) + geom_bin2d(bins = 50) +
	scale_fill_continuous(type = "viridis") + xlim(0.05, 0.95) + ylim(0.05, 0.95) +
	xlab("dgrp2 AF") + ylab(ylab) + ggtitle(tag) + theme_classic() +
	theme(axis.ticks = element_blank(), axis.text = element_blank())
pA <- dens("R.FREQ.y", "bcftools AF", "A"); pB <- dens("hafpipe.means.R", "hafpipe AF", "B")
pdf("FigS14A_dgrp2_vs_bcftools_density.pdf", width = 6.5, height = 5); suppressWarnings(print(pA)); invisible(dev.off())
pdf("FigS14B_dgrp2_vs_hafpipe_density.pdf", width = 6.5, height = 5); suppressWarnings(print(pB)); invisible(dev.off())
suppressWarnings(save_panels("FigS14_density.pdf", list(pA, pB), ncol = 2, nrow = 1, width = 13, height = 5))

ok <- check_expected("11_hafpipe_validation", obs)
sink()
if (!ok) quit(status = 1)
