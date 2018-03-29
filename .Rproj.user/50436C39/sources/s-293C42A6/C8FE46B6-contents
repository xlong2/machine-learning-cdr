
#ten_foldcv_blindblastlist
conf_tables_all_loops_blindBLAST=list()
conf_tables_all_loops_blindBLAST_diff=list()
new_accuracy_list=list()
for(loop in names(ten_foldcv_blindblastlist)){
  
  model_re=ten_foldcv_blindblastlist[[loop]]
  repeats=ceiling(length(model_re)%/%10)
  print(length(model_re))
  print(repeats)
  model_re_com=do.call(rbind,model_re)
  model_re_com[,1]=sapply(strsplit(model_re_com[,1],"\\."),"[[",1)
  model_re_com[,2]=sapply(strsplit(model_re_com[,2],"\\."),"[[",1)
  unique_cluster=gsub("\\*","none",unique(data_by_loop_type_list_unduplicated[[loop]][[1]]$cluster_type))
  which_1=which(!model_re_com[,1]%in%unique_cluster)
  model_re_com[which_1,1]=rep(paste(c(strsplit(loop,"_")[[1]],"none"),collapse="-"),length(which_1))
  which_2=which(!model_re_com[,2]%in%unique_cluster)
  model_re_com[which_2,2]=rep(paste(c(strsplit(loop,"_")[[1]],"none"),collapse="-"),length(which_2))
  
  
  conf_t=as.data.frame(table(model_re_com[,1],model_re_com[,2]))
  conf_t$Freq=conf_t$Freq/repeats
  
  new_accuracy_list[[loop]]=sum(conf_t[as.character(conf_t$Var1)==as.character(conf_t$Var2),"Freq"])/sum(conf_t[,"Freq"])
  print(loop)
  
  
  conf_tables_all_loops_blindBLAST[[loop]]=conf_t
  
  conf_t=conf_t[conf_t$Freq>0 & as.character(conf_t[,1])!=as.character(conf_t[,2]),]
  conf_tables_all_loops_blindBLAST_diff[[loop]]=conf_t
}
conf_tables_all_loops_blindBLAST_diff
save_file("conf_tables_all_loops_blindBLAST_diff")
save_file("conf_tables_all_loops_blindBLAST")