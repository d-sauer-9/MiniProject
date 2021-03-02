#running Sleuth for Mini Project
#just a general idea until I get kallisto working

#import created table from kallisto
library(sleuth)
library(dplyr)
stab<-read.table("kalTable.txt"),header=TRUE,stringsAsFactors=FALSE)
sleuthOb<-sleuth_prep(stab)

#trying to find most significant results
sleuth_table<-slueth_results(sluethOb,'reduced:full','lrt',show_all=FALSE)
sleuth_sig<-dplyr::filter(sleuth_table,qval<=.05)%>% dplyr::arrange(pval)
head(dplyr::select(sleuth_sig,target_id,test_stat,pval,qval))
write.table(sleuth_sig,file = "MiniProjectLog.log",quote=FALSE,row.names=FALSE)