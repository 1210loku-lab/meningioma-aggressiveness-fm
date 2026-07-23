suppressMessages(library(GEOquery))
options(timeout=300)
for (acc in c("GSE183655","GSE206647")) {
  cat("\n=====", acc, "=====\n")
  r <- tryCatch({
    g <- getGEO(acc, GSEMatrix=TRUE, getGPL=FALSE, destdir="data/raw")
    es <- g[[1]]; pd <- Biobase::pData(es)
    cat("n samples:", nrow(pd), "| platform:", as.character(annotation(es)), "\n")
    chc <- grep("title|characteristics|source|grade|recur|library|strateg", names(pd), ignore.case=TRUE, value=TRUE)
    for(cc in chc){v<-unique(as.character(pd[[cc]])); if(length(v)<=14) cat("[",cc,"]:",paste(v,collapse=" | "),"\n") else cat("[",cc,"]:(",length(v),"uniq;",paste(head(v,5),collapse=" | "),")\n")}
    # list supplementary file names (size only, NO download of big tars)
    cat("--- suppl files (GEO):\n")
    sf <- tryCatch(getGEOSuppFiles(acc, fetch_files=FALSE), error=function(e) NULL)
    if(!is.null(sf)) print(sf[,c("fname","size")]) else cat("(suppl list unavailable via API)\n")
    TRUE
  }, error=function(e){cat("ERR:",conditionMessage(e),"\n");FALSE})
}
