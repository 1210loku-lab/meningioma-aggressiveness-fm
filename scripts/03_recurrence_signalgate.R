suppressMessages({library(GEOquery); library(limma)}); options(timeout=600); set.seed(42)
g <- getGEO("GSE16581", GSEMatrix=TRUE, getGPL=FALSE, destdir="data/raw")
es <- g[[1]]; pd <- Biobase::pData(es); X <- Biobase::exprs(es)
rf <- suppressWarnings(as.integer(as.character(pd[["recurrence_frequency:ch1"]])))
grade <- as.character(pd[["who grade:ch1"]])
cat("recurrence_frequency table:\n"); print(table(rf, useNA="always"))
rec <- ifelse(is.na(rf), NA, ifelse(rf>0,"Recur","NoRecur"))
cat("\nRecur vs NoRecur (all grades):\n"); print(table(rec, useNA="always"))
cat("Recur status within WHO grade I:\n"); print(table(grade, rec, useNA="always"))
if(max(X,na.rm=TRUE)>100) X<-log2(X+1)

gate <- function(X,grp,name){
  ok<-!is.na(grp); Xi<-X[,ok]; grp<-factor(grp[ok])
  if(nlevels(grp)<2||min(table(grp))<3){cat(" [",name,"] insufficient (",paste(table(grp),collapse="v"),")\n");return()}
  ns<-sum(topTable(eBayes(lmFit(Xi,model.matrix(~grp))),coef=2,number=Inf)$adj.P.Val<0.05,na.rm=TRUE)
  nr<-sum(topTable(eBayes(lmFit(Xi,model.matrix(~grp))),coef=2,number=Inf)$P.Value<0.05,na.rm=TRUE)
  perm<-replicate(300,{pg<-sample(grp);sum(topTable(eBayes(lmFit(Xi,model.matrix(~pg))),coef=2,number=Inf)$adj.P.Val<0.05,na.rm=TRUE)})
  ep<-(sum(perm>=ns)+1)/301
  cat(sprintf(" [GATE %s] n=%s | adjP<0.05=%d rawP<0.05=%d | perm mean=%.1f max=%d | emp p=%.4f -> %s\n",
    name,paste(table(grp),collapse="v"),ns,nr,mean(perm),max(perm),ep,ifelse(ep<0.05&ns>=15,"PASS",ifelse(nr>2000,"WEAK(raw signal only)","FAIL"))))
}
cat("\n== Recurrence signal gates ==\n")
gate(X, rec, "Recur_vs_NoRecur_allGrades")
recI <- ifelse(grade=="1", rec, NA)   # within grade I (the hard clinical case)
gate(X, recI, "Recur_vs_NoRecur_gradeI_only")
