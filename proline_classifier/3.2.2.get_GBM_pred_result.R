
all_models=list.files(pattern="*extra_test.rds",path="./proline_classifier/rmsd_cluster_hits_rmsd",full.names = FALSE)
lapply(all_models_list_by_loop,function(x){do.call(rbind,x)})
conf_tables_all_loops_gbm_diff=list()
conf_tables_all_loops_gbm=list()

gbm_folds_sd=list()

for(loop in names(all_models_list_by_loop)){
  paras=best_parameters_each_loop[[loop]]
  best_para=paste(c(paras[1,1],paras[1,2],paras[1,3],paras[1,4]),collapse="-")
  file_n=grep(loop,grep(best_para,all_models,value=TRUE),value=TRUE)
  if(length(file_n)>1){
    file_n=grep("xgboost",file_n,value=TRUE)
  }
  model_file=paste(c("./proline_classifier/rmsd_cluster_hits_rmsd/",file_n),collapse="")
  model=readRDS(model_file)
  gbm_folds_sd[[loop]]= model$results$AccuracySD
  
  model_re=model$pred
  
  blindBLAST_unique_clusters=gsub("\\*","none",unique(data_by_loop_type_list_unduplicated_for_blindBLAST[[loop]][[1]]$cluster_type))
  obs=factor(gsub("_","-",as.character(model_re$obs)),levels=blindBLAST_unique_clusters)
  pred=factor(gsub("_","-",as.character(model_re$pred)),levels=blindBLAST_unique_clusters)

  obs_which=which(!obs%in%blindBLAST_unique_clusters)
  pred_which=which(!pred%in%blindBLAST_unique_clusters)
  pred[pred_which]=rep(paste(c(strsplit(loop,"_")[[1]],"none"),collapse="-"),length(pred_which))
  
  obs[obs_which]=rep(paste(c(strsplit(loop,"_")[[1]],"none"),collapse="-"),length(obs_which))
  repeats=ceiling(as.numeric(gsub("Resample","",model_re[dim(model_re)[1],"Resample"]))%/%10)
  print(repeats)
  conf_t=as.data.frame(table(obs,pred))
  colnames(conf_t)[1:2]=c("Var1","Var2")
  
  conf_tables_all_loops_gbm[[loop]]=conf_t
  conf_t=conf_t[conf_t$Freq>0 & conf_t[,1]!=conf_t[,2],]
  conf_t$Freq=conf_t$Freq/repeats
  conf_t[,1]=gsub("_","-",conf_t[,1])
  conf_t[,2]=gsub("_","-",conf_t[,2])
  conf_tables_all_loops_gbm_diff[[loop]]=conf_t
}
save_file("conf_tables_all_loops_gbm_diff")
save_file("conf_tables_all_loops_gbm")
save_file("gbm_folds_sd")

#compare the conf_table from gbm to that of blindBLAST




