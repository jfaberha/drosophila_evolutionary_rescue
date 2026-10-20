########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 09_GO_enrichment.R.

###  Produces:
###  — Table S8: GO enrichment of the convergent candidate gene lists.

###  Test:
###  - hypergeometric test for over-representation
###  - reference set = whole annotation (all D. melanogaster genes with >= 1 GO term)
###  - genes annotated to a term also count for all its ancestors (is_a, part_of)
###  - Benjamini-Hochberg across all terms annotated in the test set; keep corrected P < 0.1
###  - test set = list genes with >= 1 GO annotation
###  Annotation: GO Consortium go-basic.obo and fb.gaf.gz in PATHS$go_files (download with get_GO_files.sh on the login node). 'NOT' annotations are excluded; all evidence codes are used. Table S8 was made with the files downloaded on 12 Oct 2026; newer files can change the table slightly.
###  Gene lists: written by 05_convergent_candidates.R
###  T1 SE   = August (extinct) convergent candidates
###  T4 SP   = November (rescued) convergent candidates
###  T1 & T4 = relaxed three-way set (all four contrasts FDR < 0.1)

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
OUT <- module_dir("09_GO_enrichment"); setwd(OUT)
sink("09_GO_enrichment_report.txt", split = TRUE)
FDR_CUT <- 0.1
IN05 <- file.path(PATHS$out_root, "05_convergent_candidates")
LISTS <- c(`T1 SE` = "glm.all.rolwin21.TPT1.convergent_candidate_genes.txt",
	`T4 SP` = "glm.all.rolwin21.TPT4.convergent_candidate_genes.txt",
	`T1 & T4` = "glm.all.rolwin21.adaptive.fdr10.txt")
gf <- file.path(PATHS$go_files, c("go-basic.obo", "fb.gaf.gz"))
if (!all(file.exists(gf))) stop("GO files not found in ", PATHS$go_files, " - run get_GO_files.sh on the login node first.")
cat("Annotation files:", paste(basename(gf), format(file.info(gf)$mtime, "%Y-%m-%d"), collapse = "; "), "\n")

## ---- Ontology ------------------------------------------------------------------------
obo <- readLines(gf[1])
stanza <- which(startsWith(obo, "["))
starts <- which(obo == "[Term]")
ends <- vapply(starts, function(s) { nx <- stanza[stanza > s]; if (length(nx)) nx[1] - 1 else length(obo) }, 1)
terms <- lapply(seq_along(starts), function(i) {
	b <- obo[starts[i]:ends[i]]
	g <- function(tag) sub(paste0("^", tag, ": "), "", b[startsWith(b, paste0(tag, ": "))])
	list(id = g("id")[1], name = g("name")[1], ns = g("namespace")[1], obsolete = any(b == "is_obsolete: true"), alt = g("alt_id"),
		parents = c(sub(" .*", "", g("is_a")), sub("^part_of (GO:[0-9]+).*", "\\1", grep("^part_of ", g("relationship"), value = TRUE))))
})
terms <- terms[!vapply(terms, `[[`, logical(1), "obsolete")]
ids <- vapply(terms, `[[`, "", "id")
tname <- setNames(vapply(terms, `[[`, "", "name"), ids)
tns <- setNames(vapply(terms, `[[`, "", "ns"), ids)
parents <- setNames(lapply(terms, `[[`, "parents"), ids)
alt <- unlist(lapply(terms, function(t) setNames(rep(t$id, length(t$alt)), t$alt)))
cat("GO terms (non-obsolete):", length(ids), "\n")
anc_env <- new.env(hash = TRUE)
ancestors <- function(id) {
	if (!is.null(r <- anc_env[[id]])) return(r)
	p <- parents[[id]]; p <- p[p %in% ids]
	r <- unique(c(id, unlist(lapply(p, ancestors))))
	assign(id, r, envir = anc_env); r
}

## ---- Annotation, propagated to ancestors ---------------------------------------------
gaf <- read.delim(gzfile(gf[2]), header = FALSE, comment.char = "!", quote = "", colClasses = "character")
gaf <- gaf[startsWith(gaf$V2, "FBgn") & !grepl("NOT", gaf$V4), c("V2", "V5")]
names(gaf) <- c("gene", "go")
gaf$go <- ifelse(gaf$go %in% names(alt), alt[gaf$go], gaf$go)
gaf <- unique(gaf[gaf$go %in% ids, ])
direct <- split(gaf$go, gaf$gene)
gene2go <- lapply(direct, function(g) unique(unlist(lapply(unique(g), ancestors))))
N <- length(gene2go)
pairs <- data.frame(gene = rep(names(gene2go), lengths(gene2go)), go = unlist(gene2go, use.names = FALSE))
K <- table(pairs$go)
cat("Reference set (whole annotation):", N, "genes\n")

## ---- Enrichment ----------------------------------------------------------------------
enrich <- function(genes, label) {
	genes <- unique(genes[genes %in% names(gene2go)]); n <- length(genes)
	p <- pairs[pairs$gene %in% genes, ]; x <- table(p$go); go <- names(x)
	pv <- phyper(as.numeric(x) - 1, as.numeric(K[go]), N - as.numeric(K[go]), n, lower.tail = FALSE)
	padj <- p.adjust(pv, "BH"); keep <- padj < FDR_CUT
	cat(sprintf("%-8s %d annotated genes in the test set; %d terms at corrected P < %s\n", label, n, sum(keep), FDR_CUT))
	if (!any(keep)) return(NULL)
	gl <- tapply(p$gene, p$go, function(v) paste(toupper(sort(v)), collapse = ", "))
	r <- data.frame(`reciprocal outlier list` = label, `GO-ID` = as.integer(sub("GO:", "", go[keep])),
		`P-value` = sprintf("%.2E", pv[keep]), `corr P-value` = sprintf("%.2E", padj[keep]),
		`test set` = paste0(as.numeric(x[keep]), "/", n), `background set` = paste0(as.numeric(K[go[keep]]), "/", N),
		description = tname[go[keep]], `genes in test set` = gl[go[keep]], ontology = tns[go[keep]],
		check.names = FALSE, stringsAsFactors = FALSE)
	r[order(padj[keep], pv[keep]), ]
}
obs <- c()
res <- do.call(rbind, lapply(names(LISTS), function(l) {
	f <- file.path(IN05, LISTS[[l]])
	if (!file.exists(f)) stop("Missing ", f, " - run 05_convergent_candidates.R first.")
	x <- readLines(f); g <- unique(regmatches(x, regexpr("FBgn[0-9]+", x)))
	cat(l, ": ", length(g), " FBgn IDs in the list\n", sep = "")
	r <- enrich(g, l); obs[paste0("TableS8_rows_", gsub("[^A-Za-z0-9]", "", l))] <<- if (is.null(r)) 0 else nrow(r); r
}))
if (is.null(res)) stop("No GO term passed corrected P < ", FDR_CUT, " in any list.")
write.csv(res, "TableS8_with_ontology.csv", row.names = FALSE)
s8 <- res[, names(res) != "ontology"]; rownames(s8) <- NULL
write.csv(s8, "TableS8.csv", row.names = FALSE)
print(s8[, 1:7], row.names = FALSE)

pf <- file.path(RF_HOME, "data", "published_tables", "TableS8.csv")
if (file.exists(pf)) {
	pub <- read.csv(pf, check.names = FALSE, colClasses = "character")
	now <- read.csv("TableS8.csv", check.names = FALSE, colClasses = "character")
	same <- identical(dim(pub), dim(now)) && all(as.matrix(pub) == as.matrix(now))
	cat("\nTable S8:", if (same) "identical to the published table" else "DIFFERS from the published table", "\n")
	if (!same) {
		kp <- paste(pub[[1]], pub[[2]]); kn <- paste(now[[1]], now[[2]])
		if (length(setdiff(kp, kn))) cat("  terms only in the published table:", paste(setdiff(kp, kn), collapse = "; "), "\n")
		if (length(setdiff(kn, kp))) cat("  terms only in the new table:", paste(setdiff(kn, kp), collapse = "; "), "\n")
	}
	obs["TableS8_matches_published"] <- as.numeric(same)
}
ok <- check_expected("09_GO_enrichment", obs)
sink()
if (!ok) quit(status = 1)
