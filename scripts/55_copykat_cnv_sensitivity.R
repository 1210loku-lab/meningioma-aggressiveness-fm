# CopyKAT sensitivity analysis for GSE206647 tumour/tumour-like cells.
#
# Usage:
#   Rscript scripts/55_copykat_cnv_sensitivity.R               # pilot sample
#   Rscript scripts/55_copykat_cnv_sensitivity.R --all         # all 16 tumours
#   Rscript scripts/55_copykat_cnv_sensitivity.R GSM7476446 ...
#
# Each patient is analysed separately with immune cells from the same sample as
# the diploid reference. To keep the analysis computationally bounded and
# patient-balanced, at most 1,500 tumour/tumour-like and 500 immune reference
# cells are sampled per patient. Raw counts are never modified.

suppressPackageStartupMessages({
  library(Matrix)
  library(copykat)
})

set.seed(20260717)
root <- normalizePath(".")
rawdir <- file.path(root,"data","raw","GSE206647_ex")
outroot <- file.path(root,"results","audit_submission","copykat_GSE206647")
dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
meta <- read.csv(file.path(root,"results","scrna","GSE206647_cellmeta.csv"),
                 stringsAsFactors=FALSE)

tumour.samples <- sort(unique(meta$gsm[meta$grade %in% c("I","II","III")]))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0L) {
  samples <- "GSM7476446"
} else if (identical(args,"--all")) {
  samples <- tumour.samples
} else {
  samples <- intersect(args,tumour.samples)
}
stopifnot(length(samples)>0L)

read_sample_counts <- function(gsm) {
  matrix.file <- list.files(rawdir,pattern=paste0("^",gsm,"_.*_matrix\\.mtx\\.gz$"),
                            full.names=TRUE)
  barcode.file <- list.files(rawdir,pattern=paste0("^",gsm,"_.*_barcodes\\.tsv\\.gz$"),
                             full.names=TRUE)
  feature.file <- list.files(rawdir,pattern=paste0("^",gsm,"_.*_features\\.tsv\\.gz$"),
                             full.names=TRUE)
  stopifnot(length(matrix.file)==1L,length(barcode.file)==1L,length(feature.file)==1L)
  x <- readMM(gzfile(matrix.file))
  barcodes <- read.delim(gzfile(barcode.file),header=FALSE,stringsAsFactors=FALSE)[,1]
  features <- read.delim(gzfile(feature.file),header=FALSE,stringsAsFactors=FALSE)
  gene <- features[,2]
  if (anyDuplicated(gene)) {
    unique.gene <- unique(gene)
    collapse <- sparseMatrix(i=match(gene,unique.gene),
                             j=seq_along(gene),x=1,
                             dims=c(length(unique.gene),length(gene)))
    x <- collapse %*% x
    gene <- unique.gene
  }
  rownames(x) <- gene
  colnames(x) <- paste0(gsm,"_",barcodes)
  x
}

run_one <- function(gsm) {
  message("CopyKAT sample: ",gsm)
  md.all <- meta[meta$gsm==gsm,,drop=FALSE]
  sample.dir <- file.path(outroot,gsm)
  dir.create(sample.dir,recursive=TRUE,showWarnings=FALSE)
  cached.file <- file.path(sample.dir,"copykat_prediction.csv")
  if (file.exists(cached.file)) {
    prediction <- read.csv(cached.file,stringsAsFactors=FALSE)
    pred.cell.col <- grep("cell",names(prediction),ignore.case=TRUE,value=TRUE)[1]
    pred.class.col <- grep("pred",names(prediction),ignore.case=TRUE,value=TRUE)[1]
    pred.cell <- as.character(prediction[[pred.cell.col]])
    pred.class <- tolower(as.character(prediction[[pred.class.col]]))
    cell.type <- setNames(md.all$celltype,md.all$cell)[pred.cell]
    tumour.pred <- pred.class[cell.type=="Meningioma"]
    ref.pred <- pred.class[cell.type%in%c("Myeloid","Tcell","Bplasma")]
    return(data.frame(
      gsm=gsm,status="ok_cached",grade=unique(md.all$grade),
      n_tumour_input=length(tumour.pred),n_reference_input=length(ref.pred),
      n_predicted=length(pred.class),
      n_aneuploid=sum(pred.class=="aneuploid",na.rm=TRUE),
      n_diploid=sum(pred.class=="diploid",na.rm=TRUE),
      tumour_aneuploid_fraction=mean(tumour.pred=="aneuploid",na.rm=TRUE),
      reference_diploid_fraction=mean(ref.pred=="diploid",na.rm=TRUE),
      message="cached compact prediction"
    ))
  }

  x <- read_sample_counts(gsm)
  md <- md.all[md.all$cell %in% colnames(x),,drop=FALSE]
  tumour <- md$cell[md$celltype=="Meningioma"]
  immune <- md$cell[md$celltype %in% c("Myeloid","Tcell","Bplasma")]
  stopifnot(length(tumour)>=100L,length(immune)>=100L)
  set.seed(20260717 + match(gsm,tumour.samples))
  tumour.use <- sample(tumour,min(1500L,length(tumour)))
  immune.use <- sample(immune,min(500L,length(immune)))
  cells.use <- c(tumour.use,immune.use)
  rawmat <- x[,cells.use,drop=FALSE]
  # copykat 1.2.5 calls data.frame methods internally and does not reliably
  # accept a dgCMatrix under current Matrix versions.
  rawmat <- as.matrix(rawmat)
  rm(x); gc()

  # CopyKAT writes very large gene-by-cell intermediates. Run it in an exact
  # task-specific temporary directory and retain only the compact prediction.
  work.dir <- tempfile(pattern=paste0("copykat_",gsm,"_"),tmpdir=tempdir())
  dir.create(work.dir,recursive=TRUE,showWarnings=FALSE)
  oldwd <- getwd(); setwd(work.dir)
  fit <- tryCatch(
    copykat(rawmat=rawmat,id.type="S",cell.line="no",ngene.chr=5,
            min.gene.per.cell=200,LOW.DR=0.05,UP.DR=0.1,win.size=25,
            norm.cell.names=immune.use,KS.cut=0.1,sam.name=gsm,
            distance="euclidean",test.emd="FALSE",output.seg=FALSE,
            plot.genes=FALSE,genome="hg20",n.cores=4),
    error=function(e)e
  )
  setwd(oldwd)
  if (inherits(fit,"error")) {
    writeLines(conditionMessage(fit),file.path(sample.dir,"ERROR.txt"))
    unlink(work.dir,recursive=TRUE,force=TRUE)
    return(data.frame(gsm=gsm,status="error",grade=unique(md$grade),
                      n_tumour_input=length(tumour.use),
                      n_reference_input=length(immune.use),
                      n_predicted=NA,n_aneuploid=NA,n_diploid=NA,
                      tumour_aneuploid_fraction=NA,
                      reference_diploid_fraction=NA,
                      message=conditionMessage(fit)))
  }

  prediction <- fit$prediction
  write.csv(prediction,file.path(sample.dir,"copykat_prediction.csv"),row.names=FALSE)
  unlink(work.dir,recursive=TRUE,force=TRUE)
  pred.cell.col <- grep("cell",names(prediction),ignore.case=TRUE,value=TRUE)[1]
  pred.class.col <- grep("pred",names(prediction),ignore.case=TRUE,value=TRUE)[1]
  pred.cell <- as.character(prediction[[pred.cell.col]])
  pred.class <- tolower(as.character(prediction[[pred.class.col]]))
  cell.type <- setNames(md$celltype,md$cell)[pred.cell]
  tumour.pred <- pred.class[cell.type=="Meningioma"]
  ref.pred <- pred.class[cell.type%in%c("Myeloid","Tcell","Bplasma")]

  data.frame(
    gsm=gsm,status="ok",grade=unique(md$grade),
    n_tumour_input=length(tumour.use),
    n_reference_input=length(immune.use),
    n_predicted=length(pred.class),
    n_aneuploid=sum(pred.class=="aneuploid",na.rm=TRUE),
    n_diploid=sum(pred.class=="diploid",na.rm=TRUE),
    tumour_aneuploid_fraction=mean(tumour.pred=="aneuploid",na.rm=TRUE),
    reference_diploid_fraction=mean(ref.pred=="diploid",na.rm=TRUE),
    message=""
  )
}

summary <- do.call(rbind,lapply(samples,run_one))
summary.file <- if (length(samples)==length(tumour.samples)) {
  file.path(outroot,"copykat_all_patient_summary.csv")
} else {
  file.path(outroot,paste0("copykat_summary_",paste(samples,collapse="_"),".csv"))
}
write.csv(summary,summary.file,row.names=FALSE)
if (length(samples)==length(tumour.samples)) {
  write.csv(summary,file.path(root,"results","audit_submission",
                              "copykat_all_patient_summary.csv"),row.names=FALSE)
}
print(summary)
cat("Summary:",summary.file,"\n")

# When the full tumour cohort is available, recompute the program on the exact
# CopyKAT-classified tumour cells. This is a sensitivity analysis: CopyKAT was
# run on a patient-balanced subsample, so it does not replace the full-cell
# tumour/tumour-like pseudobulk used in the primary analysis.
if (length(samples)==length(tumour.samples) && all(summary$status %in% c("ok","ok_cached"))) {
  program <- read.csv(file.path(root,"docs","Table_S2_aggressiveness_program_genes.csv"),
                      stringsAsFactors=FALSE)
  program <- program[nzchar(program$symbol) & !grepl("^ENSG",program$symbol),]
  up <- program$symbol[program$direction=="up_in_WHO_II_vs_I"]
  down <- program$symbol[program$direction=="down_in_WHO_II_vs_I"]

  broad.profiles <- list()
  cnv.profiles <- list()
  cell.summary <- list()
  aggregate_cpm <- function(x,cells) {
    counts <- Matrix::rowSums(x[,cells,drop=FALSE])
    log1p(counts/sum(counts)*1e6)
  }

  for (gsm in tumour.samples) {
    pred <- read.csv(file.path(outroot,gsm,"copykat_prediction.csv"),stringsAsFactors=FALSE)
    pred.cell <- as.character(pred[[grep("cell",names(pred),ignore.case=TRUE,value=TRUE)[1]]])
    pred.class <- tolower(as.character(pred[[grep("pred",names(pred),ignore.case=TRUE,value=TRUE)[1]]]))
    md <- meta[meta$gsm==gsm,,drop=FALSE]
    cell.type <- setNames(md$celltype,md$cell)[pred.cell]
    broad.cells <- pred.cell[cell.type=="Meningioma" & pred.class %in% c("aneuploid","diploid")]
    cnv.cells <- pred.cell[cell.type=="Meningioma" & pred.class=="aneuploid"]
    x <- read_sample_counts(gsm)
    broad.cells <- intersect(broad.cells,colnames(x))
    cnv.cells <- intersect(cnv.cells,colnames(x))
    stopifnot(length(broad.cells)>=10L)
    broad.profiles[[gsm]] <- aggregate_cpm(x,broad.cells)
    if (length(cnv.cells)>=10L) {
      cnv.profiles[[gsm]] <- aggregate_cpm(x,cnv.cells)
    }
    cell.summary[[gsm]] <- data.frame(gsm=gsm,grade=unique(md$grade),
                                      n_copykat_tumour=length(broad.cells),
                                      n_copykat_aneuploid_tumour=length(cnv.cells),
                                      aneuploid_fraction=length(cnv.cells)/length(broad.cells))
    rm(x); gc()
  }

  eligible <- names(cnv.profiles)
  stopifnot(length(eligible)>=10L)
  genes <- Reduce(intersect,c(lapply(broad.profiles[eligible],names),
                              lapply(cnv.profiles,names)))
  broad.mat <- do.call(cbind,lapply(broad.profiles[eligible],function(x)x[genes]))
  cnv.mat <- do.call(cbind,lapply(cnv.profiles,function(x)x[genes]))
  colnames(broad.mat) <- colnames(cnv.mat) <- eligible
  score_matrix <- function(m) {
    z <- t(scale(t(m)))
    up.use <- intersect(up,rownames(z)); down.use <- intersect(down,rownames(z))
    up.use <- up.use[apply(z[up.use,,drop=FALSE],1,function(v)all(is.finite(v)))]
    down.use <- down.use[apply(z[down.use,,drop=FALSE],1,function(v)all(is.finite(v)))]
    colMeans(z[up.use,,drop=FALSE])-colMeans(z[down.use,,drop=FALSE])
  }
  broad.score <- score_matrix(broad.mat)
  cnv.score <- score_matrix(cnv.mat)
  cell.summary <- do.call(rbind,cell.summary)
  cell.summary$grade_ord <- match(cell.summary$grade,c("I","II","III"))
  cell.summary$cnv_eligible <- cell.summary$gsm %in% eligible
  cell.summary$broad_balanced_pseudobulk_score <- broad.score[cell.summary$gsm]
  cell.summary$copykat_aneuploid_pseudobulk_score <- cnv.score[cell.summary$gsm]
  write.csv(cell.summary,file.path(outroot,"copykat_cnv_confirmed_pseudobulk_scores.csv"),row.names=FALSE)
  write.csv(cell.summary,file.path(root,"results","audit_submission",
                                   "copykat_cnv_confirmed_pseudobulk_scores.csv"),row.names=FALSE)

  analysis.summary <- cell.summary[cell.summary$cnv_eligible,,drop=FALSE]
  broad.cor <- suppressWarnings(cor.test(analysis.summary$broad_balanced_pseudobulk_score,
                                         analysis.summary$grade_ord,method="spearman",exact=FALSE))
  cnv.cor <- suppressWarnings(cor.test(analysis.summary$copykat_aneuploid_pseudobulk_score,
                                       analysis.summary$grade_ord,method="spearman",exact=FALSE))
  broad.aov <- summary(aov(broad_balanced_pseudobulk_score~factor(grade),data=analysis.summary))[[1]][1,"Pr(>F)"]
  cnv.aov <- summary(aov(copykat_aneuploid_pseudobulk_score~factor(grade),data=analysis.summary))[[1]][1,"Pr(>F)"]
  score.cor <- suppressWarnings(cor.test(analysis.summary$broad_balanced_pseudobulk_score,
                                         analysis.summary$copykat_aneuploid_pseudobulk_score,
                                         method="spearman",exact=FALSE))
  lines <- c(
    "CopyKAT patient-balanced CNV sensitivity analysis",
    sprintf("Patients classified: %d; patients with >=10 aneuploid tumour cells: %d; median tumour cells classified per patient: %.0f; median aneuploid tumour fraction: %.3f",
            nrow(cell.summary),nrow(analysis.summary),median(cell.summary$n_copykat_tumour),median(cell.summary$aneuploid_fraction)),
    sprintf("Balanced broad tumour/tumour-like pseudobulk: grade rho=%.3f, p=%.4g; ANOVA p=%.4g",
            unname(broad.cor$estimate),broad.cor$p.value,broad.aov),
    sprintf("CopyKAT-aneuploid tumour pseudobulk: grade rho=%.3f, p=%.4g; ANOVA p=%.4g",
            unname(cnv.cor$estimate),cnv.cor$p.value,cnv.aov),
    sprintf("Broad versus CopyKAT-aneuploid patient scores: rho=%.3f, p=%.4g",
            unname(score.cor$estimate),score.cor$p.value),
    "Interpretation: CopyKAT is a sensitivity classifier on a balanced cell subsample and does not establish malignancy for every cell."
  )
  writeLines(lines,file.path(outroot,"copykat_cnv_sensitivity_summary.txt"))
  writeLines(lines,file.path(root,"results","audit_submission",
                             "copykat_cnv_sensitivity_summary.txt"))
  cat(paste(lines,collapse="\n"),"\n")
}
