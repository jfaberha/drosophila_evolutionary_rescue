########################################

### Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
### Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###*Corresponding authors:
#René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
#Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

### In this file: 01_glm/glm_models.R — every per-SNP GLM used in the paper.

### Every model is reproduced:
###   count = round(allele frequency x 100); quasi-Poisson GLM per SNP; P-values from emmeans; SNPs kept if allele-frequency variance across the model's samples is > 0.001.
### Additional notes:
###   - samples are paired with frequency columns by sample name, not by position;
###   - a failed fit gives an NA row instead of being silently dropped (the originals used tryCatch + rbind, which would shift every later P-value onto the wrong SNP if a fit failed);
###   - SNPs run in parallel and results are saved in chunks, so a rerun resumes.
###   - exception: the founders model (PAvE.F) uses the allele frequency itself as the response in a quasi-binomial GLM, on every SNP without missing values (no variance filter).

### Usage:
###   Rscript glm_models.R list                                 # list models
###   Rscript glm_models.R <model> validate [n_snps] [cores]    # check vs paper tables
###   Rscript glm_models.R <model> full [0] [cores]             # all SNPs
### 'validate' fits n_snps (default 3000) SNPs spread evenly along the genome and compares them, matched by CHROM/POS, with the table used in the paper.

########################################

suppressMessages(library(emmeans))
args <- commandArgs(TRUE)
f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(file.path(dirname(sub("^--file=", "", f_arg)), "..")) else getwd()
source(file.path(RF_HOME, "00_config.R"))

## ---- Sample selections ----
sel_groups <- function(groups, keep = TRUE) { force(groups); force(keep); function(meta) {
	m <- meta[meta$batch == "a" & meta$treat.fix %in% groups, ]
	if (keep) m <- m[m$samp %in% KEEP_SAMP, ]
	m
}}
## August four-group model (E, PA, SE, SP); optionally leaving one exposed cage out
sel_T1 <- function(drop_cage = NULL) { force(drop_cage); function(meta) {
	ct <- COND_T1; if (!is.null(drop_cage)) ct <- ct[ct$cage != drop_cage, ]
	m <- merge(meta, ct, by = c("cage", "tpt"))
	m <- unique(m[m$batch == "a" & m$treat %in% c("S", "PA", "E"), ])
	m <- m[m$samp %in% KEEP_SAMP, ]
	m[order(m$samp), ]
}}
## November three-group model (E, PA, S), leaving one rescued cage out
sel_T4 <- function(drop_cage) { force(drop_cage); function(meta) {
	m <- meta[meta$tpt == 4 & meta$batch == "a" & meta$treat.fix %in% c("S", "PA", "E"), ]
	m <- m[!(m$treat.fix == "S" & m$cage == drop_cage), ]
	m <- unique(m[m$rep != "dupl", ])
	m <- m[m$samp %in% KEEP_SAMP, ]
	m[order(m$samp), ]
}}

## Founders (16 July): pre-adapted vs unexposed founder samples
sel_founders <- function(meta) meta[as.character(meta$cage) == "F" & meta$treat.fix %in% c("E", "PA"), ]

## ---- Per-SNP statistics ----
st_pair_treat <- function(d) as.data.frame(emmeans(glm(count ~ treat.fix, data = d, family = quasipoisson, na.action = na.exclude), pairwise ~ treat.fix)$contrasts)$p.value
st_pair_tpt   <- function(d) as.data.frame(emmeans(glm(count ~ tpt, data = d, family = quasipoisson, na.action = na.exclude), pairwise ~ tpt)$contrasts)$p.value
st_treat_by_tpt <- function(d) joint_tests(glm(count ~ treat.fix * tpt, data = d, family = quasipoisson, na.action = na.exclude), by = "tpt")$p.value
st_int <- function(d) {
	emm <- emmeans(glm(count ~ treat.fix * tpt, data = d, family = quasipoisson, na.action = na.exclude), ~ tpt * treat.fix)
	as.data.frame(contrast(emm, interaction = c("pairwise", "consec"), by = NULL))$p.value
}
## PA vs S: treatment test at each time point (treatment x time-point model), then
## time-point pairs (additive treatment + time-point model)
st_PAvS_add <- function(d) {
	gi <- glm(count ~ treat.fix * tpt, data = d, family = quasipoisson, na.action = na.exclude)
	ga <- glm(count ~ treat.fix + tpt, data = d, family = quasipoisson, na.action = na.exclude)
	c(joint_tests(gi, by = "tpt")$p.value, as.data.frame(emmeans(ga, pairwise ~ tpt)$contrasts)$p.value)
}
## founders: allele frequency itself as the response, quasi-binomial GLM
st_founders <- function(d) as.data.frame(emmeans(glm(count ~ treat.fix, data = d, family = quasibinomial, na.action = na.exclude), pairwise ~ treat.fix)$contrasts)$p.value
st_condition <- function(d) as.data.frame(emmeans(glm(count ~ condition, data = d, family = quasipoisson, na.action = na.exclude), pairwise ~ condition)$contrasts)$p.value

TP <- c("T1vT2", "T1vT3", "T1vT4", "T2vT3", "T2vT4", "T3vT4")
out1 <- function(file, cols, afmean = FALSE, idx = seq_along(cols)) list(file = file, cols = cols, afmean = afmean, idx = idx)

## ---- Model registry --------------------------------------------------------------------
## Each model: sample selection, per-SNP statistic, and the table(s) it writes (file name as
## used by 02_build_glm_master.R; column headers are copied from the paper's table if present).
MODELS <- list(
	SvE_all  = list(sel = sel_groups(c("E", "S")),  stat = st_pair_treat, out = list(out1("rudflies_2023_SvE_treat.all.table.GLMcontrast.txt",  "contrastout.treat.S.E",  TRUE))),
	PAvE_all = list(sel = sel_groups(c("E", "PA")), stat = st_pair_treat, out = list(out1("rudflies_2023_PAvE_treat.all.table.GLMcontrast.txt", "contrastout.treat.PA.E", TRUE))),
	PAvS_all = list(sel = sel_groups(c("PA", "S")), stat = st_pair_treat, out = list(out1("rudflies_2023_PAvS_treat.all.table.GLMcontrast.txt", "contrastout.treat.PA.S", TRUE))),
	SvE_tpt  = list(sel = sel_groups(c("E", "S")),  stat = st_treat_by_tpt, out = list(out1("rudflies_2023_SvE_treat.GLMcontrast.table.txt",   paste0("SvE.T", 1:4)))),
	PAvE_tpt = list(sel = sel_groups(c("E", "PA")), stat = st_treat_by_tpt, out = list(out1("rudflies_2023_PAvE_treat.table.GLMcontrast2.txt", paste0("PAvE.T", 1:4), TRUE))),
	SvE_int  = list(sel = sel_groups(c("E", "S")),  stat = st_int, out = list(out1("rudflies_2023_SvE_int.GLMcontrast.table.txt", paste0("int.SvE.", TP)))),
	PAvS_add = list(sel = sel_groups(c("PA", "S"), keep = FALSE), stat = st_PAvS_add, out = list(
		out1("rudflies_2023_PAvS_treat.GLMcontrast.txt", paste0("PAvS.T", 1:4), idx = 1:4),
		out1("rudflies_2023_PAvS_tpt.GLMcontrast.txt", paste0("PA.S.", TP), idx = 5:10))),
	PAvS_int = list(sel = sel_groups(c("PA", "S"), keep = FALSE), stat = st_int, out = list(out1("rudflies_2023_PAvS_int.GLMcontrast.txt", paste0("int.PAvS.", TP)))),
	S_tpt    = list(sel = sel_groups("S",  keep = FALSE), stat = st_pair_tpt, out = list(out1("rudflies_2023_S_tpt.table.GLMcontrast.txt",  paste0("S.", TP),  TRUE))),
	PA_tpt   = list(sel = sel_groups("PA", keep = FALSE), stat = st_pair_tpt, out = list(out1("rudflies_2023_PA_tpt.table.GLMcontrast.txt", paste0("PA.", TP), TRUE))),
	E_tpt    = list(sel = sel_groups("E"), stat = st_pair_tpt, out = list(out1("rudflies_2023_E_tpt.table.GLMcontrast.txt", paste0("E.", TP), TRUE))),
	PAvSvSEvE_T1 = list(sel = sel_T1(), stat = st_condition, out = list(out1("rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt",
		c("EvPA.T1", "EvSE.T1", "EvSP.T1", "PAvSE.T1", "PAvSP.T1", "SEvSP.T1"))))
)
## Leave-one-out models (Tables S3, S4, S6)
for (cg in c(SE_CAGES, SP_CAGES)) MODELS[[paste0("T1_LOO_no", cg)]] <- list(sel = sel_T1(cg), stat = st_condition,
	out = list(out1(paste0("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no", cg, ".txt"), c("EvPA.T1", "EvSE.T1", "EvSP.T1", "PAvSE.T1", "PAvSP.T1", "SEvSP.T1"))))
for (cg in SP_CAGES) MODELS[[paste0("T4_LOO_no", cg)]] <- list(sel = sel_T4(cg), stat = st_pair_treat,
	out = list(out1(paste0("rudflies_2023_PAvSvE.wLociGLMcontrast.TPT4.LOO.no", cg, ".txt"), c("EvPA.T4", "EvSP.T4", "PAvSP.T4"))))
## PA vs E founders model ('PAvE.F'): allele frequencies (not counts), quasi-binomial; every SNP without missing values (no variance filter)
MODELS$PAvE_founders <- list(sel = sel_founders, stat = st_founders, response = "freq", varfilter = FALSE,
	out = list(out1("rudflies_2023_PAvE_treat.table.GLMcontrast.founders.txt", "PAvE.F", TRUE)))

if (length(args) == 0 || args[1] == "list") { cat(names(MODELS), sep = "\n"); quit(save = "no") }
model <- args[1]; mode <- args[2]
n_val <- if (length(args) >= 3 && as.integer(args[3]) > 0) as.integer(args[3]) else 3000L
cores <- if (length(args) >= 4) as.integer(args[4]) else as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (!model %in% names(MODELS)) stop("Unknown model '", model, "'. Run with 'list' to see the models.")
stopifnot(mode %in% c("validate", "full"))
M <- MODELS[[model]]
K <- max(unlist(lapply(M$out, `[[`, "idx")))
OUT <- module_dir("01_glm"); setwd(OUT)
msg("Model %s (%s) on %d cores", model, mode, cores)

## ---- Inputs and sample selection -------------------------------------------------------
haf <- load_haf(exclude = TRUE)
meta <- M$sel(haf$meta)
meta <- meta[as.character(meta$samp) %in% names(haf$freq), ]
cols <- match(as.character(meta$samp), names(haf$freq))          # pair by NAME
## Diagnostic: does the file order agree with name pairing?
file_order <- names(haf$freq)[sort(cols)]
cat("Samples:", nrow(meta), "| frequency columns in metadata order:", identical(file_order, as.character(meta$samp)), "\n")
print(table(group = if ("condition" %in% names(meta)) meta$condition else meta$treat.fix, tpt = meta$tpt))
meta$treat.fix <- factor(meta$treat.fix); meta$tpt <- factor(meta$tpt)
if ("condition" %in% names(meta)) meta$condition <- factor(meta$condition)
meta <- meta[, intersect(c("samp", "treat.fix", "tpt", "condition"), names(meta))]

AF <- as.matrix(haf$freq[, cols]); loci <- haf$freq[, 1:2]; rm(haf); invisible(gc())
v <- rowSums((AF - rowMeans(AF))^2) / (ncol(AF) - 1)
ok <- if (isFALSE(M$varfilter)) !is.na(v) else !is.na(v) & v > 0.001
AF <- AF[ok, , drop = FALSE]; loci <- loci[ok, ]
msg("SNPs %s: %s", if (isFALSE(M$varfilter)) "without missing values" else "passing the variance filter", format(nrow(AF), big.mark = ","))
resp <- if (is.null(M$response)) "count" else M$response   # count = round(frequency x 100); freq = frequency
if (mode == "validate") {
	keep <- unique(round(seq(1, nrow(AF), length.out = min(n_val, nrow(AF)))))   # spread along the genome
	AF <- AF[keep, , drop = FALSE]; loci <- loci[keep, ]
}

## ---- Fit (parallel, chunked, resumable) -------------------------------------------------
fit_block <- function(B) {
	d <- meta; out <- matrix(NA_real_, nrow(B), K)
	for (i in seq_len(nrow(B))) {
		d$count <- if (resp == "freq") B[i, ] else round(B[i, ] * 100)
		p <- tryCatch(stat(d), error = function(e) rep(NA_real_, K))
		if (length(p) == K) out[i, ] <- p
	}
	out
}
stat <- M$stat
chunk_dir <- paste0("chunks_", model, "_", mode)
if (mode == "validate") unlink(chunk_dir, recursive = TRUE)   # validation always refits; only full runs resume
dir.create(chunk_dir, showWarnings = FALSE)
CHUNK <- if (mode == "validate") nrow(AF) else 50000L
starts <- seq(1L, nrow(AF), by = CHUNK)
cl <- parallel::makePSOCKcluster(cores)
invisible(parallel::clusterEvalQ(cl, suppressMessages(library(emmeans))))
parallel::clusterExport(cl, c("meta", "stat", "K", "fit_block", "resp"))
t1 <- Sys.time()
for (st in starts) {
	f <- file.path(chunk_dir, sprintf("chunk_%09d.rds", st))
	if (file.exists(f)) next
	idx <- st:min(st + CHUNK - 1L, nrow(AF))
	blocks <- lapply(split(idx, ceiling(seq_along(idx) / 200L)), function(ix) AF[ix, , drop = FALSE])
	P <- do.call(rbind, parallel::parLapplyLB(cl, blocks, fit_block))
	stopifnot(nrow(P) == length(idx))
	saveRDS(P, f); rm(blocks, P); invisible(gc())
	msg("  SNPs %d-%d done", st, max(idx))
}
parallel::stopCluster(cl)
P <- do.call(rbind, lapply(file.path(chunk_dir, sprintf("chunk_%09d.rds", starts)), readRDS))
stopifnot(nrow(P) == nrow(AF))
## Bookkeeping self-check: refit 40 SNPs here, one at a time, and compare with the worker results
chk <- unique(round(seq(1, nrow(AF), length.out = 40)))
d <- max(abs(fit_block(AF[chk, , drop = FALSE]) - P[chk, , drop = FALSE]), na.rm = TRUE)
if (d > 1e-12) stop("Self-check failed: fitted results do not match a refit of the same SNPs (in a full run, delete ", chunk_dir, " if the samples or model changed).")
el <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
msg("Fitted %s SNPs in %.0f s; failed fits (NA rows): %d", format(nrow(P), big.mark = ","), el, sum(!complete.cases(P)))

## ---- Write / compare -----------------------------------------------------------------
read_tab <- function(f) as.data.frame(if (requireNamespace("data.table", quietly = TRUE)) data.table::fread(f) else read.delim(f))
all_ok <- TRUE; n_cmp <- 0
for (o in M$out) {
	tab <- data.frame(CHROM = loci$CHROM, POS = loci$POS)
	if (o$afmean) tab$af.mean <- rowMeans(AF)
	tab <- cbind(tab, setNames(as.data.frame(P[, o$idx, drop = FALSE]), o$cols))
	paper_f <- file.path(PATHS$glm_tables, o$file)
	paper <- if (file.exists(paper_f)) read_tab(paper_f) else NULL
	if (mode == "full") {
		## use the paper table's column headers so downstream code sees the same names
		if (!is.null(paper) && ncol(paper) == ncol(tab)) names(tab) <- names(paper)
		write.table(tab, o$file, sep = "\t", quote = FALSE, row.names = FALSE)
		msg("Wrote %s (%s SNPs)", file.path(OUT, o$file), format(nrow(tab), big.mark = ","))
	}
	if (is.null(paper)) { cat("No paper table to compare for", o$file, "\n"); next }
	n_cmp <- n_cmp + 1
	## compare the P-value columns (last length(cols) columns of the paper table), matched by locus
	pc <- (ncol(paper) - length(o$cols) + 1):ncol(paper)
	pp <- setNames(paper[, c(1, 2, pc)], c("CHROM", "POS", paste0("old_", o$cols)))
	m <- merge(tab, pp, by = c("CHROM", "POS"))
	cat(sprintf("\n===== %s | %s =====\n  paper table: %s SNPs, headers: %s\n  this run: %s SNPs; %s matched by CHROM/POS\n",
		model, o$file, format(nrow(paper), big.mark = ","), paste(names(paper)[pc], collapse = " "),
		format(nrow(tab), big.mark = ","), format(nrow(m), big.mark = ",")))
	if (mode == "full") cat(sprintf("  SNP sets identical: %s\n", nrow(m) == nrow(tab) && nrow(m) == nrow(paper)))
	for (cc in o$cols) {
		a <- m[[cc]]; b <- m[[paste0("old_", cc)]]
		dd <- max(abs(a - b), na.rm = TRUE); r <- suppressWarnings(cor(a, b, use = "complete.obs"))
		ok_c <- is.finite(dd) && dd < 1e-6   # same model and data; allows for emmeans version rounding
		all_ok <- all_ok && ok_c
		cat(sprintf("  %-5s %-16s max |diff| = %.2e  cor = %.8f  NA this/paper = %d/%d\n", if (ok_c) "PASS" else "FAIL", cc, dd, r, sum(is.na(a)), sum(is.na(b))))
	}
	if (nrow(m) < 0.99 * nrow(tab)) { all_ok <- FALSE; cat("  FAIL  fewer than 99% of SNPs matched by position\n") }
}
cat(sprintf("\n%s: %s\n", model, if (n_cmp == 0) "NOTHING COMPARED (no paper table found)" else if (all_ok) "ALL PASS — reproduces the table used in the paper" else "DIFFERENCES FOUND — see above"))
if (mode == "validate") msg("Projected full run: %.1f h on %d cores (%s SNPs)", el / nrow(P) * sum(ok) / 3600, cores, format(sum(ok), big.mark = ","))
msg("Finished")
