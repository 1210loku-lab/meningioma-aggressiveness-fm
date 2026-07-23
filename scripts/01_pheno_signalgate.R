suppressMessages({library(GEOquery); library(limma)})
options(timeout=900); set.seed(42)
acc <- "GSE16581"
g <- getGEO(acc, GSEMatrix=TRUE, getGPL=FALSE, destdir="data/raw")
es <- g[[1]]; pd <- Biobase::pData(es); X <- Biobase::exprs(es)
cat("==",acc,"== dim:",nrow(X),"x",ncol(X)," platform:",annotation(es),"\n")
# show all characteristics columns + their unique values
chc <- grep("characteristics|title|grade|who|recur|surviv|status|outcome", names(pd), ignore.case=TRUE, value=TRUE)
for(cc in chc){ v<-unique(as.character(pd[[cc]])); if(length(v)<=20) cat("[",cc,"]:",paste(v,collapse=" | "),"\n") else cat("[",cc,"]:(",length(v),"uniq;",paste(head(v,4),collapse=" | "),")\n")}

signal_gate <- function(X, grp, name){
  ok<-!is.na(grp); X<-X[,ok]; grp<-factor(grp[ok])
  if(nlevels(grp)<2){cat("  ",name,": <2 levels, skip\n");return(invisible())}
  if(max(X,na.rm=TRUE)>100) X<-log2(X+1)
  fit<-eBayes(lmFit(X,model.matrix(~grp)))
  ns<-sum(topTable(fit,coef=2,number=Inf)$adj.P.Val<0.05,na.rm=TRUE)
  perm<-replicate(200,{pg<-sample(grp);sum(topTable(eBayes(lmFit(X,model.matrix(~pg))),coef=2,number=Inf)$adj.P.Val<0.05,na.rm=TRUE)})
  ep<-(sum(perm>=ns)+1)/201
  cat(sprintf("  [GATE %s] n=%s | real adjP<0.05=%d | perm mean=%.1f max=%d | emp p=%.4f -> %s\n",
      name, paste(table(grp),collapse="v"), ns, mean(perm), max(perm), ep,
      ifelse(ep<0.05 & ns>=20,"PASS","FAIL/WEAK")))
  invisible(list(ns=ns,ep=ep))
}
# build grade groups from whichever column holds WHO grade
gcol <- chc[grep("grade|who",chc,ignore.case=TRUE)][1]
if(!is.na(gcol)){
  gr<-as.character(pd[[gcol]]); cat("grade raw sample:",paste(head(unique(gr),8),collapse=" | "),"\n")
  gnum<-ifelse(grepl("1|I\\b|grade: 1|^1",gr),"GI",ifelse(grepl("2|3|II|III",gr),"GII_III",NA))
  cat("grade table:\n");print(table(gnum,useNA="always"))
  signal_gate(X,gnum,"GradeI_vs_higher")
} else cat("no grade column auto-found; inspect chc above\n")
