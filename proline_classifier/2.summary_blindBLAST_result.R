
#ten_foldcv_blindblastlist
conf_tables_all_loops_blindBLAST=list()
conf_tables_all_loops_blindBLAST_diff=list()
new_accuracy_list=list()




all_folds_sd_list=list()
for(loop in names(ten_foldcv_blindblastlist)){
  
  model_re=ten_foldcv_blindblastlist[[loop]]
  repeats=ceiling(length(model_re)%/%10)
  print(length(model_re))
  print(repeats)
  unique_cluster=gsub("\\*","none",unique(data_by_loop_type_list_unduplicated[[loop]][[1]]$cluster_type))
  
  #individual accuracy
  all_folds_acc=unlist(lapply(model_re,function(ind_re){
    get_conf_from_blindBLAST(ind_re,unique_cluster,1)
  }))
  all_folds_acc=all_folds_acc[!is.na(all_folds_acc)]
  all_folds_sd_list[[loop]]=sd(all_folds_acc)
  # total accu
  model_re_com=do.call(rbind,model_re)
  
  accuracy_result=get_conf_from_blindBLAST(model_re_com,unique_cluster,repeats)
  new_accuracy_list[[loop]]=accuracy_result
  print(loop)
  
  
  conf_tables_all_loops_blindBLAST[[loop]]=conf_t
  
  conf_t=conf_t[conf_t$Freq>0 & as.character(conf_t[,1])!=as.character(conf_t[,2]),]
  conf_tables_all_loops_blindBLAST_diff[[loop]]=conf_t
}
conf_tables_all_loops_blindBLAST_diff
save_file("all_folds_sd_list")

save_file("conf_tables_all_loops_blindBLAST_diff")
save_file("conf_tables_all_loops_blindBLAST")