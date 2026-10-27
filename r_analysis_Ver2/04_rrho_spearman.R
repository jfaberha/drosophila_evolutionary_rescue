########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 04_rrho_spearman.R — genome-wide concordance of August contrasts.

###  Produces:
###  Main text Spearman ρ: extinct vs pre-adapted (0.56); rescued vs pre-adapted (0.14), each relative to unexposed
###  Fig. S10B   RRHO extinct vs rescued
###  Fig. S10C   RRHO extinct vs pre-adapted
###  Fig. S10D   RRHO rescued vs pre-adapted

### Additional notes: Signed -log10(P) = -log10(P) x sign of the unexposed-minus-treatment frequency difference, from the August four-group GLM; SNPs with P = 0 in any contrast removed; SNPs with mean August frequency outside 0.15-0.85 removed. Spearman uses all remaining SNPs; RRHO uses a random 25,000 (seeded). Colour scales of all panels share the maximum of the extinct-vs-pre-adapted plot, as in the original. Not carried over (unused in the paper): 'rescued vs extinct' vs pre-adapted comparison and the time-series RRHOs for September-November.

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
if (!requireNamespace("RRHO2", quietly = TRUE)) stop("Package RRHO2 is required.")
suppressMessages(library(RRHO2))
OUT <- module_dir("04_rrho_spearman"); setwd(OUT)
sink("04_rrho_spearman_report.txt", split = TRUE)

## August four-group GLM table (not the merged master table)
T1 <- read.table(file.path(PATHS$glm_tables, "rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt"), header = TRUE)
fm1 <- load_rds("freq_means_t1"); fall <- load_rds("freq_t1_allmean")
## unexposed-minus-treatment differences in mean August frequency
dif <- data.frame(CHROM = fm1$CHROM, POS = fm1$POS,
	EvSE.T1.diff = fm1$E.af - fm1$SE.af, EvSP.T1.diff = fm1$E.af - fm1$SP.af, EvPA.T1.diff = fm1$E.af - fm1$PA.af)
d <- merge(T1, dif, by = c("CHROM", "POS"))
d <- d[order(d$CHROM, d$POS), ]
d <- d[rowProds(as.matrix(d[, c("EvPA.T1", "EvSE.T1", "EvSP.T1", "SEvSP.T1")])) != 0, ]   # no P = 0
d <- merge(d, fall, by = c("CHROM", "POS"))
d <- d[d$af_mean > 0.15 & d$af_mean < 0.85, ]
msg("SNPs used: %s", format(nrow(d), big.mark = ","))

## signed -log10(P), written as (sign = |diff| / diff)
signed <- function(x, ct) data.frame(snps = paste(x$CHROM, x$POS, sep = "_"),
	logp = as.numeric(as.character(-log10(x[[ct]]) * (abs(x[[paste0(ct, ".diff")]]) / x[[paste0(ct, ".diff")]]))))

## ---- Spearman rank correlations, all SNPs ----
sp <- function(a, b, label) { ct <- cor.test(signed(d, a)$logp, signed(d, b)$logp, method = "spearman", exact = FALSE)
	data.frame(comparison = label, rho = unname(ct$estimate), S = unname(ct$statistic), P = ct$p.value, n = nrow(d)) }
rho <- rbind(sp("EvSE.T1", "EvPA.T1", "extinct vs pre-adapted (each vs unexposed)"),
	sp("EvSP.T1", "EvPA.T1", "rescued vs pre-adapted (each vs unexposed)"),
	sp("EvSE.T1", "EvSP.T1", "extinct vs rescued (each vs unexposed)"))
cat("\n===== Spearman rank correlations of signed -log10(P), August =====\n"); print(rho, row.names = FALSE)
write.csv(rho, "spearman_August.csv", row.names = FALSE)

## ---- RRHO on 25,000 random SNPs ----
rrho_shades <- colorRampPalette(c("#E8C863", "#DD9551", "#D24D4D", "#C03268", "#684187", "#34215F"))(200)
set.seed(SEED)
s <- d[sample(nrow(d), min(25000, nrow(d))), ]   # 25,000 on the real data
EvSE <- signed(s, "EvSE.T1"); EvSP <- signed(s, "EvSP.T1"); EvPA <- signed(s, "EvPA.T1")
r_SE_PA <- RRHO2_initialize(EvSE, EvPA, labels = c("EvSE.T1", "EvPA.T1"), log10.ind = TRUE)
mx <- max(na.omit(r_SE_PA$hypermat[is.finite(r_SE_PA$hypermat)]))   # common colour maximum
r_SP_PA <- RRHO2_initialize(EvSP, EvPA, labels = c("EvSP.T1", "EvPA.T1"), log10.ind = TRUE)
r_SE_SP <- RRHO2_initialize(EvSE, EvSP, labels = c("EvSE.T1", "EvSP.T1"), log10.ind = TRUE)
for (x in list(list(r_SE_SP, "FigS10B_RRHO_extinct_vs_rescued_August.pdf"),
	list(r_SE_PA, "FigS10C_RRHO_extinct_vs_preadapted_August.pdf"),
	list(r_SP_PA, "FigS10D_RRHO_rescued_vs_preadapted_August.pdf"))) {
	pdf(x[[2]], width = 6, height = 5.5); RRHO2_heatmap(x[[1]], maximum = mx, colorGradient = rrho_shades); dev.off()
}
msg("RRHO plots written (common colour maximum = %.1f)", mx)

pass <- check_expected("04_rrho_spearman", list(
	rho_extinct_vs_preadapted = round(rho$rho[1], 2), rho_rescued_vs_preadapted = round(rho$rho[2], 2)))
msg("Finished")
sink()
if (!pass) quit(status = 1)
