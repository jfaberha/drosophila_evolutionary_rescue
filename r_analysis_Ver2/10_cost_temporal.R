########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 10_cost_temporal.R — is resistance costly? Temporal allele-frequency change at loci that differ between pre-adapted and unexposed populations.

###  Outputs:
###  Tables S1 / S9  hypergeometric overlap of parallel temporal-change SNPs (8 of 9 unexposed cages / 9 of 10 pre-adapted cages, > 3%) with August 'pre-adapted vs unexposed' outliers, in each direction
###  Table S12 Spearman correlations of sign-corrected -log10(P), temporal contrasts vs August 'pre-adapted vs unexposed'
###  Table S13 top-10K outliers of four focal contrasts vs matched SNPs (signed -log10(P) in each temporal contrast; 100 iterations)
###  Fig. S12 direction of change of the top-10K SNPs in pre-adapted populations

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
suppressMessages({ library(dplyr); library(ggplot2) })
OUT <- module_dir("10_cost_temporal"); setwd(OUT)
sink("10_cost_temporal_report.txt", split = TRUE)
TC <- c("T1vT2", "T1vT3", "T1vT4", "T2vT3", "T2vT4", "T3vT4")
E_TEMPORAL_CAGES  <- c("1", "5", "9", "13", "23", "29", "35", "39", "43")             # nine unexposed cages used for the temporal tests
PA_TEMPORAL_CAGES <- c("2", "6", "10", "14", "20", "24", "31", "36", "40", "44")       # the ten pre-adapted cages

glm.all <- load_rds("glm_all")
haf <- load_haf(exclude = TRUE)
fs <- field_samples(haf); m <- fs$meta; fr <- fs$freq; sites <- fs$sites
cat("All-time-point samples: frequency columns in metadata order:", fs$positional_ok, "\n")
cat("Unexposed cages in the metadata:", paste(sort(as.integer(unique(m$cage.fix[m$treat.fix == "E"]))), collapse = ","),
	"| used for the temporal tests:", paste(E_TEMPORAL_CAGES, collapse = ","), "\n")
mn <- function(sel) rowMeans(fr[, sel, drop = FALSE])
tf <- as.character(m$treat.fix); tp <- as.character(m$tpt); cg <- as.character(m$cage.fix)

## August samples (four groups), as module 02
mt <- merge(haf$meta, COND_T1, by = c("cage", "tpt"))
mt <- unique(mt[mt$batch == "a" & mt$treat %in% c("S", "PA", "E"), ]); mt <- mt[mt$samp %in% KEEP_SAMP, ]
ft <- haf$freq[, names(haf$freq) %in% mt$samp]
mt <- mt[order(mt$samp), ]
cat("August samples: frequency columns in metadata order:", identical(as.character(mt$samp), names(ft)), "\n")
mt <- mt[match(names(ft), as.character(mt$samp)), ]
v1 <- rowVars(as.matrix(ft)) > 0.001; sites1 <- haf$freq[v1, 1:2]; ft <- ft[v1, ]
Fmean <- cbind(haf$freq[, 1:2], E = rowMeans(haf$freq[, 3:5]), PA = rowMeans(haf$freq[, 6:8]))   # founders: sample columns 3-5 (E) and 6-8 (PA)
rm(haf); invisible(gc())

##########################################################################################
## 1. Sign-corrected -log10(P) table                                                   ##
##########################################################################################
dconv <- cbind(sites,
	EvSP.T4.diff = mn(tf == "S" & tp == "4") - mn(tf == "E" & tp == "4"),
	PAvE.T1.diff = mn(tf == "PA" & tp == "1") - mn(tf == "E" & tp == "1"),
	PAvE.diff = mn(tf == "PA") - mn(tf == "E"))
for (g in c("PA", "E")) for (t in TC) {
	a <- substr(t, 2, 2); b <- substr(t, 5, 5)
	dconv[[paste0(g, ".", t, ".diff")]] <- mn(tf == g & tp == b) - mn(tf == g & tp == a)
}
d1 <- cbind(sites1,
	EvSE.T1.diff = rowMeans(ft[, mt$condition == "SE", drop = FALSE]) - rowMeans(ft[, mt$treat.fix == "E", drop = FALSE]),
	EvSP.T1.diff = rowMeans(ft[, mt$condition == "SP", drop = FALSE]) - rowMeans(ft[, mt$treat.fix == "E", drop = FALSE]))
dconv <- merge(d1, dconv, by = c("CHROM", "POS"))
dsign <- cbind(dconv[, 1:2], dconv[, -(1:2)] / abs(dconv[, -(1:2)]))
logp_cols <- c(EvSE.T1 = "EvSE.T1.logp", EvSP.T1 = "EvSP.T1.logp", EvSP.T4 = "SvE.T4.logp", PAvE.T1 = "PAvE.T1.logp", PAvE = "PAvE.logp",
	setNames(paste0("PA.", TC, ".logp"), paste0("PA.", TC)), setNames(paste0("E.", TC, ".logp"), paste0("E.", TC)))
stopifnot(identical(sub("\\.diff$", "", names(dsign)[-(1:2)]), names(logp_cols)))
lp <- glm.all[, c("CHROM", "POS", logp_cols)]; names(lp)[-(1:2)] <- paste0(names(logp_cols), ".logp")
tmp <- merge(dsign, lp, by = c("CHROM", "POS"))
slogp <- tmp[, 1:2]
for (i in seq_along(logp_cols)) slogp[[paste0(names(logp_cols)[i], ".logp.sign")]] <- tmp[[2 + i]] * tmp[[2 + length(logp_cols) + i]]
## August mean frequency over all four groups; keep SNPs with mean 0.15-0.85
afm <- cbind(sites1, af_mean = rowMeans(ft[, mt$treat.fix %in% c("E", "S", "PA"), drop = FALSE]))
dsign <- merge(dsign, afm, by = c("CHROM", "POS")); dsign <- dsign[order(dsign[, 1], dsign[, 2]), ]
dsign <- dsign[dsign$af_mean > 0.15 & dsign$af_mean < 0.85, ]
slogp <- merge(slogp, afm, by = c("CHROM", "POS"))
slogp <- slogp[slogp$af_mean > 0.15 & slogp$af_mean < 0.85, ]; slogp <- slogp[order(slogp[, 1], slogp[, 2]), ]
msg("Sign-corrected table: %s SNPs with August mean frequency 0.15-0.85", format(nrow(slogp), big.mark = ","))

##########################################################################################
## 2. Table S12: Spearman correlations                                                 ##
##########################################################################################
s12 <- do.call(rbind, lapply(TC, function(t) {
	r <- lapply(c("PA", "E"), function(g) {
		ct <- suppressWarnings(cor.test(slogp[[paste0(g, ".", t, ".logp.sign")]], slogp$PAvE.T1.logp.sign, method = "spearman", exact = FALSE))
		c(rho = unname(ct$estimate), p = ct$p.value) })
	data.frame(temporal_contrast = sub("v", " vs. ", t), PA_rho = r[[1]][["rho"]], PA_P = r[[1]][["p"]], E_rho = r[[2]][["rho"]], E_P = r[[2]][["p"]], n_SNPs = nrow(slogp))
}))
write.csv(s12, "TableS12_spearman_values.csv", row.names = FALSE)
f12 <- function(r, p) paste0(sprintf("%.4f", r), ifelse(r < 0 & p < 0.05, "*", ""))
tabs <- list(S12 = data.frame(`temporal contrast` = s12$temporal_contrast, `pre-adapted` = f12(s12$PA_rho, s12$PA_P),
	unexposed = f12(s12$E_rho, s12$E_P), check.names = FALSE))

##########################################################################################
## 3. Tables S1 and S9: parallel temporal change vs August pre-adapted outliers        ##
##    A SNP changes in parallel if its sign-corrected change is > 3% in all but one     ##
##    cage (second-lowest value > 0.03: direction of pre-adapted allelic change;        ##
##    second-highest < -0.03: direction of unexposed allelic change)                    ##
##########################################################################################
pa_sign <- (mn(tf == "PA" & tp == "1") - mn(tf == "E" & tp == "1")) / abs(mn(tf == "PA" & tp == "1") - mn(tf == "E" & tp == "1"))
one_af <- function(g, t, c) {      # one population's frequencies; NA if that sample does not exist
	j <- which(tf == g & tp == t & cg == c)
	if (length(j) == 0) return(rep(NA_real_, nrow(fr)))
	if (length(j) > 1) stop("More than one ", g, " sample for cage ", c, " time ", t)
	fr[, j]
}
second <- function(x, from_top) {  # second lowest / second highest of the non-missing values in each row
	apply(x, 1, function(v) { v <- sort(v); n <- length(v); if (n < 2) NA_real_ else if (from_top) v[n - 1] else v[2] })
}
nsel <- function(l) sum(l, na.rm = TRUE) + sum(is.na(l))   # nrow(df[l, ]) counts NA rows
fdrs <- glm.all[, c("CHROM", "POS", "PAvE.T1.fdr")]
temporal_test <- function(g, cages) {
	res <- list()
	for (t in TC) {
		a <- substr(t, 2, 2); b <- substr(t, 5, 5)
		x <- sapply(cages, function(c) (one_af(g, b, c) - one_af(g, a, c)) * pa_sign)
		d <- data.frame(sites, lo2 = second(x, FALSE), hi2 = second(x, TRUE))
		d <- merge(d, fdrs, by = c("CHROM", "POS"))
		d <- merge(d, glm.all[, c("CHROM", "POS")], by = c("CHROM", "POS"))
		sig <- d$PAvE.T1.fdr < 0.05
		for (dir in c("positive", "negative")) {
			sel <- if (dir == "positive") d$lo2 > 0.03 else d$hi2 < -0.03
			q <- nsel(sel & sig); mm <- nsel(sel); n <- nrow(d) - mm; k <- nsel(sig)
			res[[length(res) + 1]] <- data.frame(group = g, temporal = t, direction = dir, q = q, m = mm, n = n, k = k,
				P = phyper(q - 1, mm, n, k, lower.tail = FALSE), cages_with_data = sum(colSums(!is.na(x)) > 0))
		}
	}
	do.call(rbind, res)
}
fp <- function(p) ifelse(p == 0, "0", ifelse(p < 0.001, sub("e", "E", sprintf("%.2e", p)), sub("\\.?0+$", "", sprintf("%.8f", p))))
lab <- setNames(sub("v", " vs. ", TC), TC)
hyp_table <- function(r) {
	r <- r[order(match(r$direction, c("positive", "negative")), match(r$temporal, TC)), ]
	data.frame(contrast = lab[r$temporal], direction = ifelse(r$direction == "positive", "pre-adapted", "unexposed"), overlap = r$q,
		`temporal outliers` = r$m, `non-temporal outliers` = r$n, `PAvE.T1 outliers` = r$k, `hypergeometric P-value` = fp(r$P), check.names = FALSE, row.names = NULL)
}
msg("Temporal tests: unexposed (8 of 9 cages; 7 of 8 for September) and pre-adapted (9 of 10 cages) ...")
rE <- temporal_test("E", E_TEMPORAL_CAGES); rPA <- temporal_test("PA", PA_TEMPORAL_CAGES)
write.csv(rbind(rE, rPA), "TableS1_S9_values.csv", row.names = FALSE)
cat("Cages with data per unexposed contrast:", paste(paste(rE$temporal, rE$cages_with_data)[rE$direction == "positive"], collapse = "; "), "\n")
tabs$S1 <- hyp_table(rE); tabs$S9 <- hyp_table(rPA)

##########################################################################################
## 4. Table S13: top-10K SNPs vs matched SNPs (signed -log10(P) in temporal contrasts)  ##
##    and Fig. S12 direction counts                                                    ##
##########################################################################################
vep <- load_vep()
vep_priority <- vep_priority_table(glm.all, vep, Fmean); rm(vep); invisible(gc())
msg("SNPs with a consequence for matching: %s", format(nrow(vep_priority), big.mark = ","))
glm.temp <- glm.all[which(glm.all$POS %in% dsign$POS), ]      # (matched on POS only)
SETS <- list(EvSE.T1 = list(fdr = "EvSE.T1.fdr", sign = "EvSE.T1.logp.sign", dsign = "EvSE.T1.diff", seed = 1, label = "SE"),
	EvSP.T1 = list(fdr = "EvSP.T1.fdr", sign = "EvSP.T1.logp.sign", dsign = "EvSP.T1.diff", seed = 2, label = "SP"),
	EvSP.T4 = list(fdr = "SvE.T4.fdr",  sign = "EvSP.T4.logp.sign", dsign = "EvSP.T4.diff", seed = 3, label = "SP"),
	PAvE.T1 = list(fdr = "PAvE.T1.fdr", sign = "PAvE.T1.logp.sign", dsign = "PAvE.T1.diff", seed = 4, label = "PA"))
tcols <- c(paste0("PA.", TC), paste0("E.", TC))
top_snps <- function(fdr_col, n) head(unique(glm.temp[order(glm.temp[[fdr_col]]), c("CHROM", "POS", fdr_col)]), n)
slogp_nf <- slogp[, names(slogp) != "af_mean"]
vep_all <- merge(slogp, vep_priority, by = c("CHROM", "POS"))
cost <- list()
for (sn in names(SETS)) {
	s <- SETS[[sn]]
	set.seed(20260924 + s$seed)
	vep_sig <- merge(merge(top_snps(s$fdr, 10000)[, 1:2], slogp_nf, by = c("CHROM", "POS")), vep_priority, by = c("CHROM", "POS"))
	vep_nonsig <- anti_join(vep_all, vep_sig, by = c("CHROM", "POS"))
	idx <- split(seq_len(nrow(vep_nonsig)), paste(vep_nonsig$CHROM, vep_nonsig$Consequence, sep = "|"))
	pool <- lapply(seq_len(nrow(vep_sig)), function(j) {
		ii <- idx[[paste(vep_sig$CHROM[j], vep_sig$Consequence[j], sep = "|")]]; if (is.null(ii)) return(integer(0))
		found <- vep_sig$E[j]
		ii[which((vep_nonsig$E[ii] < found * 1.25) & (vep_nonsig$E[ii] > found * .75) & (abs(vep_nonsig$POS[ii] - vep_sig$POS[j]) > 50000))]
	})
	use <- which(lengths(pool) > 4)
	msg("%s: %d focal SNPs, %d with >= 5 possible matches", sn, nrow(vep_sig), length(use))
	sg <- function(df) abs(df[[s$sign]]) / df[[s$sign]]
	focal <- lapply(tcols, function(cc) vep_sig[[paste0(cc, ".logp.sign")]] * sg(vep_sig)); names(focal) <- tcols
	res <- vector("list", 100 * length(tcols)); r <- 0
	for (z in 1:100) {
		b <- vapply(pool[use], function(ii) ii[sample.int(length(ii), 1)], 1L)
		bs <- sg(vep_nonsig[b, ])
		for (cc in tcols) {
			fo <- focal[[cc]]; bg <- vep_nonsig[[paste0(cc, ".logp.sign")]][b] * bs
			tt <- t.test(fo, bg, alternative = "less")
			r <- r + 1; res[[r]] <- data.frame(SNP_set = sn, contrast = cc, iteration = z, mean_focal = mean(na.omit(fo)), mean_BG = mean(na.omit(bg)),
				t = unname(tt$statistic), df = unname(tt$parameter), p.value = tt$p.value, stderr = tt$stderr)
		}
	}
	cost[[sn]] <- do.call(rbind, res)
	write.table(cost[[sn]], paste0("ttest_", sn, "_10kSNPs_vs_matched_cost.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}
ord <- c(paste0("E.", TC), paste0("PA.", TC))
tabs$S13 <- do.call(rbind, lapply(c("PAvE.T1", "EvSE.T1", "EvSP.T1", "EvSP.T4"), function(sn) do.call(rbind, lapply(ord, function(ct) {
	x <- cost[[sn]][cost[[sn]]$contrast == ct, ]
	data.frame(`SNP set` = paste0(sn, "_10K"), treatment = ifelse(grepl("^E\\.", ct), "unexposed", "pre-adapted"),
		`temporal contrast` = sub("v", " vs. ", sub("^(E|PA)\\.", "", ct)),
		`focal SNPs: median signed -log10(P)` = sprintf("%.4f", median(x$mean_focal)),
		`matched SNPs: median signed -log10(P)` = sprintf("%.4f", median(x$mean_BG)),
		`median t-test statistic` = sprintf("%.4f", median(x$t)), `median t-test P-value` = sprintf("%.2E", median(x$p.value)),
		`perc. sig` = sprintf("%.0f", 100 * mean(x$p.value < 0.05)), check.names = FALSE)
}))))

names(tabs$S13) <- sub("-log10", paste0(intToUtf8(0x2212), "log10"), names(tabs$S13), fixed = TRUE)   # minus sign as in the SI table

## Fig. S12: direction of temporal change at the top 1K / 10K SNPs, re-oriented to the focal contrast
counts <- list()
for (sn in c("PAvE.T1", "EvSE.T1", "EvSP.T1", "EvSP.T4")) for (N in c(`1K` = 1000, `10K` = 10000)) {
	s <- SETS[[sn]]
	st <- merge(top_snps(s$fdr, N)[, 1:2], dsign[, names(dsign) != "af_mean"], by = c("CHROM", "POS"))
	for (cc in tcols) {
		rel <- st[[paste0(cc, ".diff")]] * st[[s$dsign]]
		counts[[length(counts) + 1]] <- data.frame(SNP_set = sn, top_n = if (N == 1000) "1K" else "10K", n_SNPs = nrow(st), contrast = cc,
			toward_unexposed = sum(rel == -1, na.rm = TRUE), toward_exposed = sum(rel == 1, na.rm = TRUE), exposed_group = s$label)
	}
}
counts <- do.call(rbind, counts)
write.table(counts, "temporal_signed_counts_top1K_top10K.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cnt <- counts[counts$top_n == "10K" & grepl("^PA\\.", counts$contrast), ]
sets <- c("PAvE.T1", "EvSE.T1", "EvSP.T1", "EvSP.T4")
titles <- c("pre-adapted vs. unexposed (Aug)", "extinct vs. unexposed (Aug)", "rescued vs. unexposed (Aug)", "rescued vs. unexposed (Nov)")
exposed_fill <- c("direction of pre-adapted allelic change", "direction of extinct allelic change", "direction of rescued allelic change", "direction of rescued allelic change")
intervals <- c(PA.T1vT2 = "Aug-Sep", PA.T1vT3 = "Aug-Oct", PA.T1vT4 = "Aug-Nov", PA.T2vT3 = "Sep-Oct", PA.T2vT4 = "Sep-Nov", PA.T3vT4 = "Oct-Nov")
long <- do.call(rbind, lapply(seq_along(sets), function(i) {
	x <- cnt[cnt$SNP_set == sets[i], ]
	rbind(data.frame(panel = titles[i], interval = intervals[x$contrast], n = x$toward_unexposed, direction = "direction of unexposed allelic change"),
		data.frame(panel = titles[i], interval = intervals[x$contrast], n = x$toward_exposed, direction = exposed_fill[i]))
}))
long$panel <- factor(long$panel, levels = titles); long$interval <- factor(long$interval, levels = intervals)
long$direction <- factor(long$direction, levels = c("direction of unexposed allelic change", "direction of pre-adapted allelic change", "direction of extinct allelic change", "direction of rescued allelic change"))
fills <- c("direction of unexposed allelic change" = "#79B998", "direction of pre-adapted allelic change" = "#5B397D", "direction of extinct allelic change" = "#1C1C1C", "direction of rescued allelic change" = "#EBBE45")
letters_df <- data.frame(panel = factor(titles, levels = titles), label = c("A", "B", "C", "D"))
p <- ggplot(long, aes(x = interval, y = n, fill = direction)) + geom_col(width = 0.7) +
	geom_hline(yintercept = 5000, linetype = "dashed", linewidth = 0.6) +
	geom_text(data = letters_df, aes(x = 0.6, y = 11600, label = label), inherit.aes = FALSE, hjust = 0, size = 7, fontface = "bold") +
	facet_wrap(~ panel, nrow = 1) + scale_fill_manual(values = fills, name = "LEGEND A-D") +
	scale_y_continuous(limits = c(0, 12000), breaks = seq(0, 10000, 2500), expand = c(0, 0)) +
	labs(x = "temporal contrast (pre-adapted populations)", y = "number of SNPs") +
	guides(fill = guide_legend(nrow = 1, title.position = "top")) + theme_classic(base_size = 13) +
	theme(legend.position = "top", legend.justification = "left", legend.title = element_text(face = "bold"),
		strip.background = element_blank(), strip.text = element_text(face = "bold", size = 13), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("FigS12.pdf", p, width = 16, height = 6)

##########################################################################################
## 5. Tables, comparison with the published tables, checks                             ##
##########################################################################################
obs <- c()
cat("\n===== Tables compared with data/published_tables =====\n")
for (n in names(tabs)) {
	write_csv_utf8(tabs[[n]], paste0("Table", n, ".csv"))
	pf <- file.path(RF_HOME, "data", "published_tables", paste0("Table", n, ".csv"))
	if (!file.exists(pf)) next
	pub <- read.csv(pf, check.names = FALSE, colClasses = "character", encoding = "UTF-8")
	now <- read.csv(paste0("Table", n, ".csv"), check.names = FALSE, colClasses = "character", encoding = "UTF-8")
	cat(sprintf("  Table %-4s %d rows\n", n, nrow(now)))
	same <- identical(dim(pub), dim(now)) && all(as.matrix(pub) == as.matrix(now))
	cat(sprintf("  Table %-4s %s\n", n, if (same) "identical to the published table" else "DIFFERS from the published table"))
	if (!same && identical(dim(pub), dim(now))) {
		w <- which(as.matrix(pub) != as.matrix(now), arr.ind = TRUE)
		for (k in seq_len(min(nrow(w), 15))) cat(sprintf("    row %d (%s) %s: published %s, now %s\n", w[k, 1], paste(pub[w[k, 1], 1:3], collapse = " / "), names(pub)[w[k, 2]], pub[w[k, 1], w[k, 2]], now[w[k, 1], w[k, 2]]))
	}
	obs[paste0("Table", n, "_matches_published")] <- as.numeric(same)
}
x <- rPA[rPA$temporal == "T2vT4", ]
obs["S9_T2vT4_towardE_m"] <- x$m[x$direction == "negative"]; obs["S9_T2vT4_towardE_pct"] <- round(100 * x$q[x$direction == "negative"] / x$m[x$direction == "negative"], 1)
obs["S9_T2vT4_towardPA_m"] <- x$m[x$direction == "positive"]; obs["S9_T2vT4_towardPA_pct"] <- round(100 * x$q[x$direction == "positive"] / x$m[x$direction == "positive"], 1)
obs["S1_T2vT4_towardE_P"] <- signif(rE$P[rE$temporal == "T2vT4" & rE$direction == "negative"], 2)
obs["S12_PA_rho_min"] <- round(min(s12$PA_rho), 2); obs["S12_PA_rho_max"] <- round(max(s12$PA_rho), 2)
obs["S12_E_rho_min_neg"] <- round(max(s12$E_rho[s12$E_rho < 0]), 2); obs["S12_E_rho_max_neg"] <- round(min(s12$E_rho), 2)
n10 <- unique(counts[counts$top_n == "10K", c("SNP_set", "n_SNPs")])
obs["FigS12_n_min"] <- min(n10$n_SNPs); obs["FigS12_n_max"] <- max(n10$n_SNPs)
ok <- check_expected("10_cost_temporal", obs)
sink()
if (!ok) quit(status = 1)
