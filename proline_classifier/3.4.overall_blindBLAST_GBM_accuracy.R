

blindBLAST=do.call(rbind,conf_tables_all_loops_blindBLAST)
gbm=do.call(rbind,conf_tables_all_loops_gbm)


a=sum(blindBLAST[as.character(blindBLAST$Var1)==as.character(blindBLAST$Var2),"Freq"])/sum(blindBLAST[,"Freq"])
b=sum(gbm[as.character(gbm$Var1)==as.character(gbm$Var2),"Freq"])/sum(gbm[,"Freq"])
print(c("From all available data, blindBLAST accuracy is ",a," and gbm accuracy is ",b))




common_n=intersect(names(conf_tables_all_loops_blindBLAST),names(conf_tables_all_loops_gbm))
conf_tables_all_loops_blindBLAST_c=conf_tables_all_loops_blindBLAST[common_n]
blindBLAST=do.call(rbind,conf_tables_all_loops_blindBLAST_c)
conf_tables_all_loops_gbm_c=conf_tables_all_loops_gbm[common_n]
gbm=do.call(rbind,conf_tables_all_loops_gbm_c)


a=sum(blindBLAST[as.character(blindBLAST$Var1)==as.character(blindBLAST$Var2),"Freq"])/sum(blindBLAST[,"Freq"])
b=sum(gbm[as.character(gbm$Var1)==as.character(gbm$Var2),"Freq"])/sum(gbm[,"Freq"])
print(c("From data available for common loop, blindBLAST accuracy is ",a," and gbm accuracy is ",b))



common_n=intersect(names(conf_tables_all_loops_blindBLAST),names(conf_tables_all_loops_gbm))
conf_tables_all_loops_blindBLAST_c=conf_tables_all_loops_blindBLAST[common_n]
blindBLAST=do.call(rbind,conf_tables_all_loops_blindBLAST_c)
blindBLAST=blindBLAST[!grepl("none",blindBLAST$Var1),]
conf_tables_all_loops_gbm_c=conf_tables_all_loops_gbm[common_n]
gbm=do.call(rbind,conf_tables_all_loops_gbm_c)
gbm=gbm[!grepl("none",gbm$Var1),]


a=sum(blindBLAST[as.character(blindBLAST$Var1)==as.character(blindBLAST$Var2),"Freq"])/sum(blindBLAST[,"Freq"])
b=sum(gbm[as.character(gbm$Var1)==as.character(gbm$Var2),"Freq"])/sum(gbm[,"Freq"])
print(c("From data available for common loop and do not consider cases when queries are in none cluster, blindBLAST accuracy is ",a," and gbm accuracy is ",b))
