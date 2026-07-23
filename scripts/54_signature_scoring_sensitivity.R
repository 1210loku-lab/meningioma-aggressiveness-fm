# Cross-platform scoring sensitivity for the grade-associated program.
# Tests common-platform genes, single-sample rank scoring, multiple program
# sizes, and unbalanced discovery-derived gene sets. No validation-cohort
# outcome is used to select genes, scoring method, or direction.

suppressPackageStartupMessages({
  library(GEOquery)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(illuminaHumanv4.db)
})

root <- normalizePath(".")
outdir <- file.path(root, "results", "audit_submission")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

program <- readRDS(file.path(root, "results", "deg",
                            "GSE136661_aggressiveness_program.rds"))
res <- program$res
res$ensembl <- rownames(res)
res$symbol <- unname(mapIds(org.Hs.eg.db, keys=res$ensembl,
                            column="SYMBOL", keytype="ENSEMBL",
                            multiVals="first"))

ordered_symbols <- function(direction, n=NULL) {
  if (direction == "up") {
    keep <- !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1
  } else {
    keep <- !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1
  }
  z <- unique(res$symbol[keep & !is.na(res$symbol) & res$symbol != ""])
  if (!is.null(n)) z <- head(z, n)
  z
}

balanced.sets <- lapply(c(50,100,200,300,500), function(total) {
  per.direction <- total/2
  list(up=ordered_symbols("up", per.direction),
       down=ordered_symbols("down", per.direction))
})
names(balanced.sets) <- paste0("balanced_",c(50,100,200,300,500))

# Exact historical 200-entry program: map the originally selected 100 up and
# 100 down Ensembl entries without backfilling entries that lack a symbol.
fixed.up <- unique(na.omit(unname(mapIds(org.Hs.eg.db, keys=program$up,
                                        column="SYMBOL", keytype="ENSEMBL",
                                        multiVals="first"))))
fixed.down <- unique(na.omit(unname(mapIds(org.Hs.eg.db, keys=program$dn,
                                          column="SYMBOL", keytype="ENSEMBL",
                                          multiVals="first"))))

eligible <- !is.na(res$padj) & res$padj < 0.05 &
  abs(res$log2FoldChange) > 1 & !is.na(res$symbol) & res$symbol != ""
unbalanced <- res[eligible,]
unbalanced <- unbalanced[!duplicated(unbalanced$symbol),]
top200.unbalanced <- head(unbalanced, 200)
gene.sets <- c(list(fixed_200_original=list(up=fixed.up,down=fixed.down)),
               balanced.sets, list(
  top200_unbalanced=list(
    up=top200.unbalanced$symbol[top200.unbalanced$log2FoldChange > 0],
    down=top200.unbalanced$symbol[top200.unbalanced$log2FoldChange < 0]
  ),
  threshold_all_unbalanced=list(
    up=unbalanced$symbol[unbalanced$log2FoldChange > 0],
    down=unbalanced$symbol[unbalanced$log2FoldChange < 0]
  )
))

load_gse16581 <- function() {
  es <- getGEO(filename=file.path(root,"data","raw",
                                 "GSE16581_series_matrix.txt.gz"),
               getGPL=FALSE)
  gpl <- getGEO(filename=file.path(root,"data","raw","GPL570.soft.gz"))
  x <- Biobase::exprs(es)
  pd <- Biobase::pData(es)
  if (max(x,na.rm=TRUE) > 100) x <- log2(x+1)
  tab <- GEOquery::Table(gpl)
  fd <- tab[match(rownames(x),tab$ID),,drop=FALSE]
  symbol.col <- grep("symbol",names(fd),ignore.case=TRUE,value=TRUE)[1]
  symbol <- as.character(fd[[symbol.col]])
  keep <- !is.na(symbol) & symbol != ""
  x <- x[keep,,drop=FALSE]; symbol <- symbol[keep]
  ord <- order(-rowMeans(x)); x <- x[ord,,drop=FALSE]; symbol <- symbol[ord]
  keep <- !duplicated(symbol); x <- x[keep,,drop=FALSE]
  rownames(x) <- symbol[keep]
  recurrence.frequency <- suppressWarnings(as.integer(pd[["recurrence_frequency:ch1"]]))
  event <- as.integer(recurrence.frequency > 0)
  grade <- suppressWarnings(as.integer(pd[["who grade:ch1"]]))
  list(expr=x, grade_ord=grade,
       outcomes=list(recurrence=list(keep=!is.na(event),event=event)))
}

load_gse74385 <- function() {
  mx <- read.delim(gzfile(file.path(root,"data","raw","GSE74385",
                                   "GSE74385_normalized.txt.gz")),
                   check.names=FALSE,stringsAsFactors=FALSE)
  rownames(mx) <- mx[,1]; mx <- mx[,-1]
  mx <- mx[,!grepl("Detection",colnames(mx)),drop=FALSE]
  symbol <- mapIds(illuminaHumanv4.db, keys=rownames(mx),
                   column="SYMBOL",keytype="PROBEID",multiVals="first")
  keep <- !is.na(symbol) & symbol != ""
  x <- as.matrix(sapply(mx[keep,,drop=FALSE],as.numeric))
  if (nrow(x) != sum(keep)) x <- as.matrix(mx[keep,,drop=FALSE])
  storage.mode(x) <- "numeric"
  symbol <- as.character(symbol[keep])
  ord <- order(-rowMeans(x)); x <- x[ord,,drop=FALSE]; symbol <- symbol[ord]
  keep2 <- !duplicated(symbol); x <- x[keep2,,drop=FALSE]
  rownames(x) <- symbol[keep2]
  meta <- read.csv(file.path(root,"results","deg",
                             "GSE74385_program_recurrence.csv"),
                   stringsAsFactors=FALSE)
  x <- x[,match(meta$title,colnames(x)),drop=FALSE]
  stopifnot(!anyNA(match(meta$title,colnames(mx))))
  recurrence <- as.integer(meta$outcome == "R")
  keep.recurrence <- meta$outcome %in% c("NR","R")
  composite <- as.integer(meta$outcome %in% c("R","M"))
  keep.composite <- meta$outcome %in% c("NR","R","M")
  list(expr=x, grade_ord=as.integer(meta$grade),
       outcomes=list(
         recurrence_only=list(keep=keep.recurrence,event=recurrence),
         recurrence_or_progression=list(keep=keep.composite,event=composite)
       ))
}

cohorts <- list(GSE16581=load_gse16581(), GSE74385=load_gse74385())
common.platform.genes <- Reduce(intersect,lapply(cohorts,function(z) rownames(z$expr)))

score_program <- function(expr, genes, method=c("z","rank"), common.only=FALSE) {
  method <- match.arg(method)
  if (common.only) {
    genes$up <- intersect(genes$up, common.platform.genes)
    genes$down <- intersect(genes$down, common.platform.genes)
  }
  up <- intersect(genes$up,rownames(expr))
  down <- intersect(genes$down,rownames(expr))
  stopifnot(length(up)>0,length(down)>0)
  if (method == "z") {
    transformed <- t(scale(t(expr)))
  } else {
    transformed <- apply(expr,2,rank,ties.method="average") / nrow(expr)
    rownames(transformed) <- rownames(expr)
    colnames(transformed) <- colnames(expr)
  }
  score <- colMeans(transformed[up,,drop=FALSE],na.rm=TRUE) -
    colMeans(transformed[down,,drop=FALSE],na.rm=TRUE)
  list(score=score,n_up=length(up),n_down=length(down))
}

auc_binary <- function(y, score) {
  pos <- score[y==1]; neg <- score[y==0]
  mean(outer(pos,neg,">")) + 0.5*mean(outer(pos,neg,"=="))
}

specifications <- list()
for (set.name in names(gene.sets)) {
  for (method in c("z","rank")) {
    specifications[[paste(set.name,method,"available",sep="__")]] <-
      list(set_name=set.name,method=method,common_only=FALSE)
  }
}
for (method in c("z","rank")) {
  specifications[[paste("fixed_200_original",method,"common",sep="__")]] <-
    list(set_name="fixed_200_original",method=method,common_only=TRUE)
}

rows <- list(); k <- 1L
for (cohort.name in names(cohorts)) {
  cohort <- cohorts[[cohort.name]]
  for (spec.name in names(specifications)) {
    spec <- specifications[[spec.name]]
    scored <- score_program(cohort$expr,gene.sets[[spec$set_name]],
                            method=spec$method,common.only=spec$common_only)
    grade.rho <- cor(scored$score,cohort$grade_ord,method="spearman")
    for (outcome.name in names(cohort$outcomes)) {
      outcome <- cohort$outcomes[[outcome.name]]
      keep <- outcome$keep & !is.na(outcome$event) & !is.na(scored$score)
      y <- outcome$event[keep]; s <- scored$score[keep]
      rows[[k]] <- data.frame(
        cohort=cohort.name,
        outcome=outcome.name,
        specification=spec.name,
        discovery_gene_set=spec$set_name,
        scoring_method=spec$method,
        common_platform_genes_only=spec$common_only,
        n=sum(keep),
        events=sum(y),
        n_up=scored$n_up,
        n_down=scored$n_down,
        grade_spearman=grade.rho,
        recurrence_auc=auc_binary(y,s),
        wilcoxon_p=wilcox.test(s[y==1],s[y==0],exact=FALSE)$p.value,
        stringsAsFactors=FALSE
      )
      k <- k+1L
    }
  }
}
results <- do.call(rbind,rows)

gene.set.summary <- do.call(rbind,lapply(names(gene.sets),function(nm) {
  g <- gene.sets[[nm]]
  data.frame(
    discovery_gene_set=nm,
    discovery_up=length(g$up),
    discovery_down=length(g$down),
    common_platform_up=length(intersect(g$up,common.platform.genes)),
    common_platform_down=length(intersect(g$down,common.platform.genes))
  )
}))

write.csv(results,
          file.path(outdir,"signature_scoring_sensitivity.csv"),
          row.names=FALSE)
write.csv(gene.set.summary,
          file.path(outdir,"signature_gene_set_sensitivity_summary.csv"),
          row.names=FALSE)

headline <- results[results$specification %in%
  c("fixed_200_original__z__available",
    "fixed_200_original__rank__available",
    "fixed_200_original__z__common",
    "fixed_200_original__rank__common"),]
lines <- c(
  "Cross-platform signature-scoring sensitivity",
  sprintf("Genes measurable on both arrays: %d",length(common.platform.genes)),
  apply(headline,1,function(x) sprintf(
    "%s %s %s: %s (up/down %s/%s), grade rho=%.3f, recurrence AUC=%.3f, p=%.3g",
    x[["cohort"]],x[["outcome"]],x[["specification"]],
    ifelse(x[["common_platform_genes_only"]]=="TRUE","common genes","platform-available genes"),
    x[["n_up"]],x[["n_down"]],as.numeric(x[["grade_spearman"]]),
    as.numeric(x[["recurrence_auc"]]),as.numeric(x[["wilcoxon_p"]]))),
  "Interpretation: rank scoring is computable for a single sample and common-gene scoring uses an identical gene set on both arrays."
)
writeLines(lines,file.path(outdir,"signature_scoring_sensitivity.txt"))
cat(paste(lines,collapse="\n"),"\n")
