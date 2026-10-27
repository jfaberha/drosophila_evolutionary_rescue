########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 06_cage_scores_founder_axis.R — how far each population moved toward the pre-adapted genomic state.

###  Produces:
###  Fig. 3J / SI 2.2.2 (SI 1.8.7)  August cages scored on the pre-adapted axis: SNPs with pre-adapted vs unexposed GLM FDR < 0.01 in November; direction of pre-adapted allelic change from September-November means; score = mean signed difference from the August unexposed mean (unexposed cages vs the other unexposed cages). Extinct vs rescued: exact permutation over all 252 splits, Welch, Wilcoxon.
###  Fig. 4C (SI 1.8.8)  every field population, August-November, scored on the founder axis: SNPs differing by > 0.1 between pre-adapted and unexposed founders with complete separation of technical replicates; each time point scored against the unexposed mean at that time point. Rescued vs unexposed at each time point and extinct vs unexposed in August: exact permutation tests.

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
suppressMessages(library(ggplot2))
OUT <- module_dir("06_cage_scores_founder_axis"); setwd(OUT)
sink("06_cage_scores_founder_axis_report.txt", split = TRUE)

glm.all <- load_rds("glm_all")
haf <- load_haf(exclude = TRUE)
fs <- field_samples(haf)
cat("Field samples: frequency columns in metadata order:", fs$positional_ok, "\n")
meta <- fs$meta; meta$cage <- as.character(meta$cage); meta$tpt <- as.character(meta$tpt); meta$treat.fix <- as.character(meta$treat.fix)
Fq <- as.matrix(fs$freq); key_f <- paste(fs$sites$CHROM, fs$sites$POS)

## exact permutation test, difference in means between two groups of cages
perm_test <- function(a, b) {
	pool <- c(a, b); obs <- mean(a) - mean(b)
	d <- apply(combn(length(pool), length(a)), 2, function(i) mean(pool[i]) - mean(pool[-i]))
	c(diff = obs, P_two_sided = mean(abs(d) >= abs(obs) - 1e-12), n_larger = sum(d > obs + 1e-12), n_perm = length(d))
}

## ======================================================================================
## Fig. 3J: August cages on the November pre-adapted axis
## ======================================================================================
gi <- match(paste(glm.all$CHROM, glm.all$POS), key_f); ok <- !is.na(gi)
G <- glm.all[ok, ]; F1 <- Fq[gi[ok], ]
late <- meta$tpt %in% c("2", "3", "4")
sgn <- sign(rowMeans(F1[, meta$treat.fix == "PA" & late]) - rowMeans(F1[, meta$treat.fix == "E" & late]))
E1c <- which(meta$treat.fix == "E" & meta$tpt == "1"); PA1c <- which(meta$treat.fix == "PA" & meta$tpt == "1"); S1c <- which(meta$treat.fix == "S" & meta$tpt == "1")
sel <- G$PAvE.T4.fdr < 0.01 & sgn != 0; s <- sgn[sel]
E1 <- rowMeans(F1[sel, E1c])
score <- function(j, ref) mean(s * (F1[sel, j] - ref))
cs <- rbind(
	data.frame(cage = meta$cage[S1c], group = ifelse(meta$cage[S1c] %in% SE_CAGES, "extinct", ifelse(meta$cage[S1c] %in% SP_CAGES, "rescued", "other")), score = sapply(S1c, score, ref = E1)),
	data.frame(cage = meta$cage[E1c], group = "unexposed", score = sapply(E1c, function(j) score(j, rowMeans(F1[sel, setdiff(E1c, j)])))),
	data.frame(cage = meta$cage[PA1c], group = "pre-adapted", score = sapply(PA1c, score, ref = E1)))
ex <- cs$score[cs$group == "extinct"]; re <- cs$score[cs$group == "rescued"]; un <- cs$score[cs$group == "unexposed"]; pa <- cs$score[cs$group == "pre-adapted"]
pt <- perm_test(ex, re); wt <- t.test(ex, re); wx <- wilcox.test(ex, re, exact = TRUE)
fig3j <- data.frame(n_loci = sum(sel), extinct_mean = mean(ex), extinct_min = min(ex), extinct_max = max(ex),
	rescued_mean = mean(re), rescued_min = min(re), rescued_max = max(re), preadapted_mean = mean(pa), preadapted_min = min(pa), preadapted_max = max(pa),
	unexposed_min = min(un), unexposed_max = max(un), extinct_above_unexposed_range = sum(ex > max(un)), rescued_within_unexposed_range = sum(re >= min(un) & re <= max(un)),
	permutation_P = pt[["P_two_sided"]], permutation_n_larger = pt[["n_larger"]], welch_t = unname(wt$statistic), welch_P = wt$p.value, wilcoxon_P = wx$p.value)
cat("\n===== Fig. 3J / SI 2.2.2 =====\n"); print(t(round(fig3j, 4)))
write.csv(cs, "Fig3J_cage_scores.csv", row.names = FALSE); write.csv(fig3j, "Fig3J_tests.csv", row.names = FALSE)

cols <- c("unexposed" = "#79B998", "rescued" = "#EBBE45", "extinct" = "#1C1C1C", "pre-adapted" = "#5B397D")
d <- cs[cs$group != "other", ]; d$group <- factor(d$group, levels = names(cols))
mns <- aggregate(score ~ group, d, mean)
top <- max(d$score[d$group %in% c("rescued", "extinct")]); span <- diff(range(d$score))
set.seed(1)
p3j <- ggplot(d, aes(x = group, y = score)) +
	geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
	geom_crossbar(data = mns, aes(ymin = score, ymax = score), width = 0.5, linewidth = 0.4, colour = "grey20") +
	geom_jitter(aes(fill = group), shape = 21, size = 3, width = 0.12, height = 0, colour = "white", stroke = 0.4) +
	annotate("segment", x = 2, xend = 3, y = top + 0.08 * span, yend = top + 0.08 * span, linewidth = 0.4) +
	annotate("text", x = 2.5, y = top + 0.14 * span, label = paste0("P = ", formatC(pt[["P_two_sided"]], digits = 2, format = "fg")), size = 3.6) +
	scale_fill_manual(values = cols, guide = "none") +
	labs(x = NULL, y = expression(paste("mean allele frequency shift (", Delta, italic(p), ")"))) + theme_classic(base_size = 12)
ggsave("Fig3J_cage_scores.pdf", p3j, width = 4.5, height = 4.2)

## ======================================================================================
## Fig. 4C: field populations on the founder (laboratory) axis
## ======================================================================================
meta$cage <- as.character(meta$cage.fix)   # corrected labels (six October samples were transposed; SI 1.7)
meta$group <- ifelse(meta$treat.fix == "E", "unexposed", ifelse(meta$treat.fix == "PA", "pre-adapted",
	ifelse(meta$cage %in% SE_CAGES, "extinct", ifelse(meta$cage %in% SP_CAGES, "rescued", "other"))))
fmeta <- haf$meta[haf$meta$batch == "F" & haf$meta$treat.fix %in% c("PA", "E") & as.character(haf$meta$samp) %in% names(haf$freq), ]
cat("\nFounder technical replicates:\n"); print(table(fmeta$treat.fix))
FF <- as.matrix(haf$freq[match(key_f, paste(haf$freq$CHROM, haf$freq$POS)), as.character(fmeta$samp)])
fPA <- FF[, fmeta$treat.fix == "PA", drop = FALSE]; fE <- FF[, fmeta$treat.fix == "E", drop = FALSE]
dF <- rowMeans(fPA) - rowMeans(fE)
sepF <- (apply(fPA, 1, min) > apply(fE, 1, max)) | (apply(fPA, 1, max) < apply(fE, 1, min))
selF <- which(!is.na(dF) & sepF & abs(dF) > 0.1 & sign(dF) != 0); sF <- sign(dF)[selF]
msg("Founder axis: %s SNPs", format(length(selF), big.mark = ","))
fa <- c()
for (tp in c("1", "2", "3", "4")) {
	ec <- which(meta$treat.fix == "E" & meta$tpt == tp); Ebar <- rowMeans(Fq[selF, ec])
	for (j in which(meta$tpt == tp & meta$group != "other")) {
		ref <- if (meta$group[j] == "unexposed") rowMeans(Fq[selF, setdiff(ec, j), drop = FALSE]) else Ebar
		fa <- rbind(fa, data.frame(tpt = tp, cage = meta$cage[j], group = meta$group[j], score = mean(sF * (Fq[selF, j] - ref))))
	}
}
fa <- aggregate(score ~ tpt + cage + group, fa, mean)   # average technical duplicates, if any
tests <- c()
for (tp in c("1", "2", "3", "4")) {
	y <- fa[fa$tpt == tp, ]
	pt <- perm_test(y$score[y$group == "rescued"], y$score[y$group == "unexposed"])
	tests <- rbind(tests, data.frame(comparison = paste0("rescued vs unexposed, tpt ", tp), mean_a = mean(y$score[y$group == "rescued"]), P = pt[["P_two_sided"]]))
	if (tp == "1") {
		pt <- perm_test(y$score[y$group == "extinct"], y$score[y$group == "unexposed"])
		tests <- rbind(tests, data.frame(comparison = "extinct vs unexposed, tpt 1", mean_a = mean(y$score[y$group == "extinct"]), P = pt[["P_two_sided"]]))
	}
}
gm <- aggregate(score ~ tpt + group, fa, mean)
cat("\n===== Fig. 4C: mean score by time point and group =====\n"); print(reshape(gm, idvar = "tpt", timevar = "group", direction = "wide"), row.names = FALSE)
cat("\n===== Fig. 4C tests (exact permutation, two-sided) =====\n"); print(tests, row.names = FALSE)
pa_oct <- gm$score[gm$tpt == "3" & gm$group == "pre-adapted"]; re_oct <- gm$score[gm$tpt == "3" & gm$group == "rescued"]
cat(sprintf("Rescued October shift as %% of pre-adapted October shift: %.0f%%\n", 100 * re_oct / pa_oct))
write.csv(fa, "Fig4C_founder_axis_scores.csv", row.names = FALSE); write.csv(tests, "Fig4C_tests.csv", row.names = FALSE)

months <- c("1" = "Aug", "2" = "Sep", "3" = "Oct", "4" = "Nov")
x <- fa; x$x <- as.numeric(x$tpt) + c("unexposed" = -0.14, "rescued" = 0, "extinct" = 0.12, "pre-adapted" = 0.24)[x$group]
x$group <- factor(x$group, levels = names(cols))
bd <- aggregate(score ~ tpt, x[x$group == "unexposed", ], function(v) c(lo = min(v), hi = max(v)))
bd <- data.frame(x = as.numeric(bd$tpt), lo = bd$score[, "lo"], hi = bd$score[, "hi"])
mn <- aggregate(score ~ x + group, x, mean)
sg <- tests[grepl("^rescued", tests$comparison), ]; sg$x <- 1:4; sg$y <- max(x$score) + 0.015
sg$label <- ifelse(sg$P < 0.001, "***", ifelse(sg$P < 0.01, "**", ifelse(sg$P < 0.05, "*", "n.s.")))
p4c <- ggplot(x, aes(x = x, y = score)) +
	geom_ribbon(data = bd, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "#79B998", alpha = 0.25) +
	geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
	geom_line(data = x[x$group == "rescued", ], aes(group = cage), colour = "#EBBE45", alpha = 0.45, linewidth = 0.5) +
	geom_line(data = mn[mn$group %in% c("rescued", "pre-adapted"), ], aes(group = group, colour = group), linewidth = 1.3) +
	geom_point(aes(fill = group), shape = 21, size = 2.6, colour = "white", stroke = 0.4, position = position_jitter(width = 0.04, height = 0, seed = 1)) +
	geom_point(data = mn[mn$group %in% c("rescued", "extinct", "pre-adapted"), ], aes(fill = group), shape = 23, size = 3.8, colour = "black", stroke = 0.6) +
	geom_text(data = sg, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 4) +
	scale_fill_manual(values = cols, name = NULL) + scale_colour_manual(values = cols, guide = "none") +
	scale_x_continuous(breaks = 1:4, labels = months, limits = c(0.75, 4.35)) +
	labs(x = NULL, y = expression(paste("shift in the direction of pre-adapted allelic change (", Delta, italic(p), ")"))) +
	theme_classic(base_size = 12) + theme(legend.position = "top") +
	ggtitle(paste0("Founder axis (", format(length(selF), big.mark = ","), " SNPs)")) + theme(plot.title = element_text(size = 11, face = "bold"))
ggsave("Fig4C_founder_axis_trajectory.pdf", p4c, width = 5.2, height = 4.6)

g4 <- function(tp, grp) round(gm$score[gm$tpt == tp & gm$group == grp], 3)
tp4 <- function(lbl) tests$P[tests$comparison == lbl]
pass <- check_expected("06_cage_scores_founder_axis", list(
	Fig3J_n_snps = fig3j$n_loci, Fig3J_extinct_mean = round(fig3j$extinct_mean, 3), Fig3J_rescued_mean = round(fig3j$rescued_mean, 3),
	Fig3J_preadapted_mean = round(fig3j$preadapted_mean, 3), Fig3J_permutation_P = round(fig3j$permutation_P, 3),
	SI223_extinct_min = round(fig3j$extinct_min, 3), SI223_extinct_max = round(fig3j$extinct_max, 3),
	SI223_rescued_min = round(fig3j$rescued_min, 3), SI223_rescued_max = round(fig3j$rescued_max, 3),
	SI223_preadapted_min = round(fig3j$preadapted_min, 3), SI223_preadapted_max = round(fig3j$preadapted_max, 3),
	SI223_unexposed_min = round(fig3j$unexposed_min, 3), SI223_unexposed_max = round(fig3j$unexposed_max, 3),
	SI223_extinct_above_unexposed = fig3j$extinct_above_unexposed_range, SI223_rescued_within_unexposed = fig3j$rescued_within_unexposed_range,
	SI223_welch_t = round(fig3j$welch_t, 2), SI223_welch_P = round(fig3j$welch_P, 3), SI223_wilcoxon_P = round(fig3j$wilcoxon_P, 3),
	SI223_permutations_larger = fig3j$permutation_n_larger,
	Fig4C_n_snps = length(selF),
	Fig4C_rescued_Aug = g4("1", "rescued"), Fig4C_rescued_Sep = g4("2", "rescued"), Fig4C_rescued_Oct = g4("3", "rescued"), Fig4C_rescued_Nov = g4("4", "rescued"),
	Fig4C_P_Aug = round(tp4("rescued vs unexposed, tpt 1"), 2), Fig4C_P_Sep = round(tp4("rescued vs unexposed, tpt 2"), 3),
	Fig4C_P_Oct = round(tp4("rescued vs unexposed, tpt 3"), 3), Fig4C_P_Nov = round(tp4("rescued vs unexposed, tpt 4"), 3),
	Fig4C_extinct_Aug = g4("1", "extinct"), Fig4C_extinct_Aug_P = round(tp4("extinct vs unexposed, tpt 1"), 3)))
msg("Finished")
sink()
if (!pass) quit(status = 1)
