########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 00_config.R — shared settings and helpers for every module.

###  Every module starts with  source("00_config.R")  and nothing else is hard-coded. To run on another system, edit only the PATHS block (or set the RF_* environment variables before running).                                                         

########################################

## ---- Paths -----------------------------------------------------------------------------
## Defaults are folders inside this repository; set the RF_* environment variables to use other locations.
if (!exists("RF_HOME")) RF_HOME <- getwd()
PATHS <- list(
	## Inputs: HAF-pipe allele frequencies, sample metadata, VEP/FlyCADD annotation, spinosad candidate gene list, genic-overlap BED
	raw_input  = Sys.getenv("RF_RAW_INPUT",  file.path(RF_HOME, "data", "input")),
	## Per-SNP GLM result tables (written by 01_glm/glm_models.R)
	glm_tables = Sys.getenv("RF_GLM_TABLES", file.path(RF_HOME, "data", "glm_tables")),
	## Where every module writes its outputs (one sub-folder per module)
	out_root   = Sys.getenv("RF_OUT_ROOT",   file.path(RF_HOME, "output")),
	## GO ontology and FlyBase annotation used for Table S8 (go-basic.obo, fb.gaf.gz; see get_GO_files.sh)
	go_files   = Sys.getenv("RF_GO_FILES",   file.path(RF_HOME, "data", "GO"))
)
PATHS$rds <- file.path(PATHS$out_root, "rds")   # shared objects written by 02 and read by later modules

FILES <- list(
	hafpipe   = "rudflies_2023_hafpipe.csv",
	meta      = "rudflies_2023_meta.tsv",
	vep       = "filtered-all.annot.vcf.FLYCADD.tsv",
	spino     = "spino.cand.list.txt",
	genic_bed = "rudflies_2023_hafpipe_loci_genic_overlap_info.bed"
)

## ---- Samples ---------------------------------------------------------------------------
## Sample 90 = unexposed (E) cage 29, September 30 (tpt 2): failed low-diversity library, excluded from every analysis.
EXCLUDED_SAMPLES <- "90"
## GLM scripts restrict experimental samples to these sample numbers
KEEP_SAMP <- c(1:9, 30:171)
## August (tpt 1) outcome of each exposed cage, and the T1 grouping used in the four-group GLM
SE_CAGES <- c("11", "21", "27", "41", "45")   # extinct
SP_CAGES <- c("3", "7", "15", "33", "37")     # rescued (persisted)
PA_CAGES <- c("2", "6", "10", "14", "20", "24", "31", "36", "40", "44")
E_CAGES  <- c("1", "5", "9", "13", "19", "23", "29", "35", "39", "43")
COND_T1 <- data.frame(cage = c(SE_CAGES, SP_CAGES, PA_CAGES, E_CAGES), tpt = "1",
	condition = rep(c("SE", "SP", "PA", "E"), c(5, 5, 10, 10)), stringsAsFactors = FALSE)

## Plot colours (treatment order E, PA, SE, SP), as in the published figures
SAMPLE_COLS <- c(E = "#79B998", PA = "#5B397D", SE = "#1C1C1C", SP = "#EBBE45")
SEED <- 2023

## ---- Package helpers (base-R versions used if a package is missing) -------------------
if (requireNamespace("matrixStats", quietly = TRUE)) {
	rowVars <- matrixStats::rowVars; colVars <- matrixStats::colVars; rowMaxs <- matrixStats::rowMaxs; rowMins <- matrixStats::rowMins; rowProds <- matrixStats::rowProds
} else {
	rowVars <- function(x, ...) { x <- as.matrix(x); rowSums((x - rowMeans(x))^2) / (ncol(x) - 1) }
	colVars <- function(x, ...) { x <- as.matrix(x); colSums((t(t(x) - colMeans(x)))^2) / (nrow(x) - 1) }
	rowProds <- function(x, ...) apply(as.matrix(x), 1, prod)
	rowMaxs <- function(x, ...) apply(as.matrix(x), 1, max)
	rowMins <- function(x, ...) apply(as.matrix(x), 1, min)
}
if (requireNamespace("slider", quietly = TRUE)) {
	slide_mean <- slider::slide_mean
} else {
	## centred window, partial windows at the ends, NA if any value in the window is NA (as slider)
	slide_mean <- function(x, before = 0, after = 0) {
		n <- length(x); i <- seq_len(n)
		lo <- pmax(1, i - before); hi <- pmin(n, i + after)
		x0 <- ifelse(is.na(x), 0, x); cs <- c(0, cumsum(x0)); cn <- c(0, cumsum(is.na(x)))
		out <- (cs[hi + 1] - cs[lo]) / (hi - lo + 1)
		out[(cn[hi + 1] - cn[lo]) > 0] <- NA
		out
	}
}

## ---- Small helpers ---------------------------------------------------------------------
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")

## Output folder for a module, created if needed
module_dir <- function(module) { d <- file.path(PATHS$out_root, module); dir.create(d, FALSE, TRUE); d }

## Read the HAF-pipe table and metadata, clean the headers, and remove the excluded sample(s). Returns list(meta, freq); freq columns 1-2 are CHROM, POS.
load_haf <- function(exclude = TRUE) {
	meta <- read.table(file.path(PATHS$raw_input, FILES$meta), header = TRUE)
	f <- file.path(PATHS$raw_input, FILES$hafpipe)
	freq <- if (requireNamespace("data.table", quietly = TRUE)) {
		as.data.frame(data.table::fread(f, sep = ",", header = TRUE))
	} else read.delim(f, header = TRUE, sep = ",")
	names(freq) <- gsub("X", "", gsub("[.]af", "", names(freq)))
	names(freq)[1:2] <- c("CHROM", "POS")
	if (exclude) {
		x <- meta[meta$samp %in% EXCLUDED_SAMPLES, ]
		stopifnot(nrow(x) == 1, x$treat.fix == "E", x$tpt == 2, x$cage.fix == 29, "90" %in% names(freq))
		meta <- meta[!meta$samp %in% EXCLUDED_SAMPLES, ]
		freq <- freq[, !names(freq) %in% EXCLUDED_SAMPLES]
	}
	msg("HAF-pipe table: %s SNPs x %d samples; metadata rows: %d; excluded: %s",
		format(nrow(freq), big.mark = ","), ncol(freq) - 2, nrow(meta), if (exclude) paste(EXCLUDED_SAMPLES, collapse = ",") else "none")
	list(meta = meta, freq = freq)
}

## Field samples of the spinosad experiment (PA, S, E; batch a) with the all-sample variance filter (> 0.001).
## Frequency columns are put in metadata-row order (paired by sample name).
field_samples <- function(haf) {
	m <- haf$meta[haf$meta$batch == "a" & haf$meta$experiment == "spino" & haf$meta$treat.fix %in% c("PA", "S", "E"), ]
	fr <- haf$freq[, names(haf$freq) %in% m$samp]
	m <- m[m$samp %in% names(fr), ]
	positional_ok <- identical(as.character(m$samp), names(fr))   # are the frequency columns in metadata order?
	m <- m[match(names(fr), as.character(m$samp)), ]
	vf <- rowVars(as.matrix(fr)) > 0.001
	list(meta = m, sites = na.omit(haf$freq[vf, 1:2]), freq = na.omit(fr[vf, ]), positional_ok = positional_ok)
}

## VEP/FlyCADD annotation (several rows per SNP where features overlap) and candidate genes
load_vep <- function() {
	f <- file.path(PATHS$raw_input, FILES$vep)
	v <- if (requireNamespace("data.table", quietly = TRUE)) as.data.frame(data.table::fread(f, header = TRUE, colClasses = list(character = "CHROM"))) else read.table(f, header = TRUE)
	v$CHROM <- as.character(v$CHROM); v
}
load_spino <- function() { s <- read.table(file.path(PATHS$raw_input, FILES$spino), header = FALSE); names(s) <- "Gene"; s }
## Objects saved by 02_build_glm_master.R
load_rds <- function(name) {
	f <- file.path(PATHS$rds, paste0(name, ".rds"))
	if (!file.exists(f)) stop("Missing ", f, " - run 02_build_glm_master.R first.")
	readRDS(f)
}

## One VEP consequence per SNP, chosen by impact (HIGH > MODERATE > LOW > MODIFIER) and then the consequence's mean FlyCADD score, with E/PA founder means (vep_priority table; used to pick matched background SNPs). keys = CHROM/POS of the SNPs to keep (glm.all).
vep_priority_table <- function(keys, vep, Fmean) {
	annot <- merge(keys[, c("CHROM", "POS")], vep, by = c("CHROM", "POS"))
	vep_cons <- unique(cbind(annot[, c(1, 2)], Consequence = annot$Consequence))
	sc <- as.data.frame(cbind(Consequence = vep$Consequence, FLYCADD = vep$FLYCADD, Extra = vep$Extra))
	sc$Extra <- sub(".*=", "", sub(";.*", "", sc$Extra))
	sc$FLYCADD <- suppressWarnings(as.numeric(sc$FLYCADD))
	m <- as.data.frame(dplyr::summarise_all(dplyr::group_by(sc, Consequence, Extra), mean))
	m$Impact <- as.numeric(c(MODIFIER = 1, LOW = 2, MODERATE = 3, HIGH = 4)[m$Extra])
	m <- m[order(-m[, 4], -m[, 3]), ]
	out <- c(); u <- merge(vep_cons, Fmean, by = c("CHROM", "POS"))
	for (r in seq_len(nrow(m))) {
		tmp <- unique(u[u$Consequence == m[r, 1], ])
		out <- rbind(out, tmp)
		u <- dplyr::anti_join(u, tmp, by = c("CHROM", "POS"))
	}
	out
}

## Write a data frame as CSV in UTF-8 whatever the session locale (Table S13 headers contain the minus sign U+2212). Character columns are quoted, numbers are not, as write.csv does.
write_csv_utf8 <- function(df, file) {
	f <- function(x) if (is.numeric(x)) as.character(x) else paste0('"', gsub('"', '""', enc2utf8(as.character(x)), useBytes = TRUE), '"')
	lines <- c(paste0('"', gsub('"', '""', enc2utf8(names(df)), useBytes = TRUE), '"', collapse = ","), do.call(paste, c(unname(lapply(df, f)), sep = ",")))
	con <- file(file, "wb"); on.exit(close(con)); writeBin(charToRaw(paste0(paste(enc2utf8(lines), collapse = "\n"), "\n")), con)
}

## Arrange ggplots in a grid and save to PDF (ggpubr::ggarrange if installed)
save_panels <- function(file, plots, ncol = 1, nrow = length(plots), width = 7.5, height = 6) {
	pdf(file, width = width, height = height)
	if (requireNamespace("ggpubr", quietly = TRUE)) print(ggpubr::ggarrange(plotlist = plots, ncol = ncol, nrow = nrow)) else {
		grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
		for (k in seq_along(plots)) print(plots[[k]], vp = grid::viewport(layout.pos.row = (k - 1) %/% ncol + 1, layout.pos.col = (k - 1) %% ncol + 1))
	}
	invisible(dev.off())
}

## -log10(P) with P = 0 set to half the smallest non-zero P in that column (for plotting)
neglog10 <- function(p) { m <- min(p[p != 0]) / 2; p[p == 0] <- m; -log10(p) }

## Compare observed values with the verified values in expected_values.csv.
## Prints PASS/FAIL per value, writes <module>_checks.csv, and returns TRUE if all pass.
check_expected <- function(module, observed) {
	ev <- read.csv(file.path(RF_HOME, "expected_values.csv"), stringsAsFactors = FALSE, colClasses = c(expected = "character"))
	## significant digits of each value as written (for TEXT comparisons)
	sig_digits <- function(x) { m <- sub("[eE].*$", "", sub("^-", "", x)); m <- gsub("\\.", "", m); nchar(sub("^0+", "", m)) }
	ev$digits <- sig_digits(ev$expected); ev$expected <- as.numeric(ev$expected)
	ev <- ev[ev$module == module, ]
	if (!nrow(ev)) { cat("\n===== Checks (", module, "): no verified values listed in expected_values.csv =====\n", sep = ""); return(invisible(TRUE)) }
	res <- do.call(rbind, lapply(seq_len(nrow(ev)), function(i) {
		k <- ev$key[i]; o <- if (k %in% names(observed)) observed[[k]] else NA
		## tolerance NA = value not yet verified:
		## reported as TEXT (match) or TEXT? (differs) and never counted as a failure
		st <- if (is.na(ev$tolerance[i])) { if (!is.na(o) && signif(o, max(1, ev$digits[i])) == ev$expected[i]) "TEXT" else "TEXT?" } else
			if (!is.na(o) && abs(o - ev$expected[i]) <= ev$tolerance[i] + 1e-9) "PASS" else "FAIL"
		data.frame(key = k, expected = ev$expected[i], observed = o, status = st,
			reported_in = ev$reported_in[i], stringsAsFactors = FALSE)
	}))
	cat("\n===== Checks against verified values (", module, ") =====\n", sep = "")
	for (i in seq_len(nrow(res))) cat(sprintf("  %-5s %-45s expected %12s  observed %12s   [%s]\n", res$status[i], res$key[i],
		format(res$expected[i], big.mark = ","), format(res$observed[i], big.mark = ","), res$reported_in[i]))
	write.csv(res, file.path(module_dir(module), paste0(module, "_checks.csv")), row.names = FALSE)
	nt <- sum(res$status %in% c("PASS", "FAIL"))
	cat(sprintf("  %d of %d checks passed", sum(res$status == "PASS"), nt))
	if (nt < nrow(res)) cat(sprintf("; %d values compared with the current text only (TEXT = same, TEXT? = differs: update the text)", nrow(res) - nt))
	cat("\n")
	invisible(!any(res$status == "FAIL"))
}

