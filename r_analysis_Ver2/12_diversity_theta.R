########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 12_diversity_theta.R — genetic diversity (Watterson's θ) in extinct, rescued, pre-adapted and unexposed populations.

###  Produces:
###  Fig. S10  relative Watterson's θ per sample, per chromosome and genome-wide, at every sampling date, at 500x (left) and 3,500x (right) effective coverage.
###  Table S5 / SI 1.8.5, 2.2.1  August samples: quasi-Poisson GLM of θ (x 10^5, rounded) on population group, one model per chromosome and one genome-wide; pairwise contrasts with emmeans.

###  Input: grenedalf (theta-watterson) tables computed from the HAF-pipe allele-frequency table with pool size 240, --filter-sample-min-count 2 and --frequency-table-int-factor 500 or 3500 (commands in grenedalf_theta.sh), in PATHS$raw_input.

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
suppressMessages({ library(ggplot2); library(emmeans) })
OUT <- module_dir("12_diversity_theta"); setwd(OUT)
sink("12_diversity_theta_report.txt", split = TRUE)
COVERAGES <- c("500", "3500")
GROUPS <- c("E", "PA", "SE", "SP")

## ---- Samples: spinosad experiment, founders (batch F) and field samples (batch a) --------
meta <- read.table(file.path(PATHS$raw_input, FILES$meta), header = TRUE)
meta <- meta[meta$batch %in% c("a", "F") & meta$treat.fix %in% c("PA", "S", "E") & meta$experiment == "spino", ]
meta <- meta[!meta$samp %in% EXCLUDED_SAMPLES, ]
meta$condition <- ifelse(meta$treat.fix == "S", ifelse(as.character(meta$cage) %in% SE_CAGES, "SE", "SP"), as.character(meta$treat.fix))
cat("Samples per group and time point (sample", paste(EXCLUDED_SAMPLES, collapse = ","), "excluded):\n")
print(table(meta$condition, meta$tpt))

## ---- Relative Watterson's θ per sample: rows = samples, columns = chromosomes + genome ----
read_theta <- function(f, label = NULL) {
	d <- read.csv(file.path(PATHS$raw_input, f), header = TRUE, check.names = FALSE)
	rc <- grep("theta_watterson_rel$", names(d), value = TRUE)
	x <- t(as.matrix(d[, rc, drop = FALSE]))
	rownames(x) <- sub("^X", "", sub("\\.theta_watterson_rel$", "", rc))
	colnames(x) <- if (is.null(label)) as.character(d[[1]]) else label
	x
}
theta <- lapply(setNames(COVERAGES, COVERAGES), function(cv) {
	x <- cbind(read_theta(paste0("diversity_chrom_eff", cv, ".csv")), read_theta(paste0("diversity_genome_eff", cv, ".csv"), "genome"))
	d <- merge(meta[, c("samp", "cage", "tpt", "condition")], data.frame(samp = rownames(x), x, check.names = FALSE), by = "samp")
	d$condition <- factor(d$condition, levels = GROUPS)
	write.csv(d, paste0("theta_rel_eff", cv, ".csv"), row.names = FALSE)
	msg("Effective coverage %sx: %d samples with theta; regions: %s", cv, nrow(d), paste(colnames(x), collapse = " "))
	d
})
regions <- setdiff(names(theta[[1]]), c("samp", "cage", "tpt", "condition"))

## ---- GLMs on August samples: one per region; theta x 10^5 rounded as the quasi-Poisson response
glm_theta <- function(d) do.call(rbind, lapply(regions, function(r) {
	dd <- d[!is.na(d[[r]]), ]; dd$count <- round(as.numeric(dd[[r]]) * 1e5); dd$condition <- droplevels(dd$condition)
	ct <- as.data.frame(emmeans(glm(count ~ condition, data = dd, family = quasipoisson), pairwise ~ condition)$contrasts)
	data.frame(region = r, contrast = sub(" - ", "v", ct$contrast), P = ct$p.value)   # labels from emmeans (E, PA, SE, SP)
}))
res <- do.call(rbind, lapply(COVERAGES, function(cv) cbind(coverage = paste0(cv, "x"), glm_theta(theta[[cv]][theta[[cv]]$tpt == 1, ]))))
cat("\nAugust samples per group:\n"); print(table(theta[[1]]$condition[theta[[1]]$tpt == 1]))
write.csv(res, "theta_glm_contrasts_long.csv", row.names = FALSE)

## Table S5: published layout (rows 2L, 2R, 3L, 3R, X, all; columns EvPA, EvSP, EvSE, PAvSP, PAvSE, SEvSP)
S5_COLS <- c("EvPA", "EvSP", "EvSE", "PAvSP", "PAvSE", "SEvSP"); S5_ROWS <- c("2L", "2R", "3L", "3R", "X", "genome")
s5 <- do.call(rbind, lapply(paste0(COVERAGES, "x"), function(cv) {
	x <- res[res$coverage == cv, ]
	m <- t(sapply(S5_ROWS, function(r) setNames(x$P[x$region == r], x$contrast[x$region == r])[S5_COLS]))
	data.frame(coverage = cv, region = sub("genome", "all", S5_ROWS), signif(m, 4), check.names = FALSE)
}))
rownames(s5) <- NULL
write.csv(s5, "TableS5.csv", row.names = FALSE)
cat("\n===== Table S5: GLM P-values, Watterson's θ, August samples =====\n"); print(s5, row.names = FALSE)

## ---- Fig. S10: every sampling date; rows all, 2L, 2R, 3L, 3R, X; 500x (left), 3,500x (right)
DATES <- c(`0` = "Jul 16", `1` = "Aug 24", `2` = "Sep 30", `3` = "Oct 27", `4` = "Nov 24")
GROUP_LAB <- c(PA = "pre-adapted", E = "unexposed", SP = "rescue", SE = "extinct")
FIG_ROWS <- c("genome", "2L", "2R", "3L", "3R", "X")
long <- do.call(rbind, lapply(COVERAGES, function(cv) {
	d <- theta[[cv]]
	do.call(rbind, lapply(FIG_ROWS, function(r) data.frame(coverage = paste0(cv, "x"), region = r, date = DATES[as.character(d$tpt)],
		group = GROUP_LAB[as.character(d$condition)], theta = as.numeric(d[[r]]))))
}))
long$region <- factor(sub("genome", "all", long$region), levels = sub("genome", "all", FIG_ROWS))
long$date <- factor(long$date, levels = DATES); long$group <- factor(long$group, levels = GROUP_LAB)
cols <- setNames(SAMPLE_COLS[names(GROUP_LAB)], GROUP_LAB)
pl <- lapply(paste0(COVERAGES, "x"), function(cv) ggplot(long[long$coverage == cv, ], aes(x = date, y = theta, fill = group)) +
	geom_boxplot(outlier.size = 0.6, linewidth = 0.3, position = position_dodge2(preserve = "single")) +
	facet_wrap(~region, ncol = 1, scales = "fixed") + scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
	ylab(expression("Watterson's" ~ theta)) + xlab(NULL) + ggtitle(bquote(.(sub("3500", "3,500", sub("x$", "", cv))) %*% "")) +
	theme_classic() + theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold")))
save_panels("FigS10_theta.pdf", pl, ncol = 2, nrow = 1, width = 8.5, height = 13)
write.csv(long, "FigS10_data.csv", row.names = FALSE)

g <- function(cv, r, cc) res$P[res$coverage == cv & res$region == r & res$contrast == cc]
obs <- c(S5_genome_SEvSP_500x = g("500x", "genome", "SEvSP"), S5_genome_SEvSP_3500x = g("3500x", "genome", "SEvSP"),
	S5_genome_EvPA_500x = g("500x", "genome", "EvPA"), S5_genome_EvPA_3500x = g("3500x", "genome", "EvPA"),
	S5_min_SEvSP = min(res$P[res$contrast == "SEvSP"]), S5_lowest_P = min(res$P))
ok <- check_expected("12_diversity_theta", obs)
msg("Finished")
sink()
if (!ok) quit(status = 1)
