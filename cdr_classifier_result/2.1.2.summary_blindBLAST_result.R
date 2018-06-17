# this script calculate the misclassification error counts for blindBLAST result.


current_d = getwd()
if (grepl("cdr_classifier_result", current_d)) {
  source("0.load_function_and_data.R")
  
}


#ten_foldcv_blindblastlist
conf_tables_all_loops_blindBLAST = list()
conf_tables_all_loops_blindBLAST_diff = list()
new_accuracy_list = list()
blindBLAST_errorcount_lists=list()

all_folds_sd_list = list()
for (loop in names(ten_foldcv_blindblastlist)) {
  print(loop)
  model_re = ten_foldcv_blindblastlist[[loop]]

  model_re_by_repeats=lapply(chunk2(1:length(model_re),3),function(x){do.call(rbind,model_re[x])})
  
  repeats = ceiling(length(model_re) %/% n_folds)
  print(length(model_re))
  print(repeats)
  unique_cluster = gsub("\\*",
                        "none",
                        unique(data_by_loop_type_list_unduplicated[[loop]][[1]]$cluster_type))
  re_by_repeats=lapply(model_re_by_repeats,function(x){
    get_conf_from_blindBLAST(x, unique_cluster, 1)
    conf_t_r = as.data.frame(table(sapply(strsplit(
      x[, 1], "\\."
    ), "[[", 1), sapply(strsplit(
      x[, 2], "\\."
    ), "[[", 1)) )
  })
  error_by_repeats=lapply(re_by_repeats,function(x){
    x=as.data.frame(x)
     error_c=x[as.character(x[,1])!=as.character(x[,2]),]
     sum(error_c$Freq)
    
    
  })
  
  blindBLAST_errorcount_lists[[loop]]=unlist(error_by_repeats)
  
  frame=re_by_repeats[[1]]
  for(ind in 2:(length(re_by_repeats))){
    print(re_by_repeats[[ind]])
    frame=as.data.frame(merge(frame,re_by_repeats[[ind]],by=c("Var1","Var2")))
  }
  
 this_loop_acc_by_repeats=lapply(3:dim(frame)[2],function(x){
    sum(frame[as.character(frame[,1])==as.character(frame[,2]),x])/sum(frame[,x])
  })
 accuracy_av=mean(unlist(this_loop_acc_by_repeats))
 sd_by_repeats=sd(unlist(this_loop_acc_by_repeats))
  
each_classification_sds=apply(frame,1,function(x){sd(unlist(x[3:5]))})
  sd_f=cbind(frame[,1:2],each_classification_sds)
  #individual accuracy
  all_folds_acc = unlist(lapply(model_re, function(ind_re) {
    get_conf_from_blindBLAST(ind_re, unique_cluster, 1)
  }))
  all_folds_acc = all_folds_acc[!is.na(all_folds_acc)]
  all_folds_sd_list[[loop]] = sd(all_folds_acc)
  # total accu

  model_re_com = do.call(rbind, model_re)
  
  accuracy_result = get_conf_from_blindBLAST(model_re_com, unique_cluster, repeats)
  new_accuracy_list[[loop]] = accuracy_result
  print(loop)
  
  conf_t = as.data.frame(table(sapply(strsplit(
    model_re_com[, 1], "\\."
  ), "[[", 1), sapply(strsplit(
    model_re_com[, 2], "\\."
  ), "[[", 1)) / repeats)
  conf_t=merge(conf_t,sd_f,by=c("Var1","Var2"))
  conf_tables_all_loops_blindBLAST[[loop]] = conf_t
  
  conf_t = conf_t[conf_t$Freq > 0 &
                    as.character(conf_t[, 1]) != as.character(conf_t[, 2]), ]
  conf_tables_all_loops_blindBLAST_diff[[loop]] = conf_t
}
conf_tables_all_loops_blindBLAST_diff

blindblast_by_loop=lapply(ten_foldcv_blindblastlist,function(x){if(n_repeats>1){split_c=chunk2(1:length(x),n_repeats)}else{split_c=list();split_c[[1]]=1:length(x)}; lapply(split_c,function(y){do.call(rbind,x[y])})})
blindblast_by_loop_bind=lapply(blindblast_by_loop,function(x){
  accuracies=lapply(x,function(y){
    all=y
    all[,1]=split_vector_and_replace(all[,1],"\\.",1,1,"")
    all[,2]=split_vector_and_replace(all[,2],"\\.",1,1,"")
    accu=dim(all[all[,1]==all[,2],])[1]/dim(all)[1]
  })
  accuracies=unlist(accuracies)
  
})

blindBLAST_mean_accu=lapply(blindblast_by_loop_bind,mean)
blindBLAST_accu_std = lapply(blindblast_by_loop_bind,sd)

blindblast_by_repeats=lapply(1:3,function(x){
  all_l_for_this_repeats=lapply(blindblast_by_loop,function(y){y[[x]]}); all=do.call(rbind, all_l_for_this_repeats)
all[,1]=split_vector_and_replace(all[,1],"\\.",1,1,"")
all[,2]=split_vector_and_replace(all[,2],"\\.",1,1,"")
accu=dim(all[all[,1]==all[,2],])[1]/dim(all)[1]
})
sd(unlist(blindblast_by_repeats))






save_file("all_folds_sd_list")
save_file("blindBLAST_mean_accu")
save_file("blindBLAST_accu_std")

save_file("conf_tables_all_loops_blindBLAST_diff")
save_file("conf_tables_all_loops_blindBLAST")
