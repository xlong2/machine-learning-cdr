#This script read all the models trained by GBM and extracting the prediction results for all the leave out cases
# for all the prediction result, average the misclassification error count by the repeat number

current_d=getwd()
if(grepl("proline_classifier",current_d)){
  source("0.load_function_and_data.R")
  
}

all_models=list.files(pattern="*extra_test.rds",path="./rmsd_cluster_hits_rmsd",full.names = FALSE)
lapply(all_models_list_by_loop,function(x){do.call(rbind,x)})
conf_tables_all_loops_gbm_diff=list()
conf_tables_all_loops_gbm=list()
all_gbm_pred_result=list()
gbm_folds_sd=list()
gbm_errorcount_list=list()
gbm_sd_by_loops=list()
gbm_accuracy_by_loops=list()
gbm_errorcount_sd_list=list()
# 
for(loop in names(all_models_list_by_loop)){
  paras=best_parameters_each_loop[[loop]]
  paras=lapply(paras,as.character)
  best_para=paste(unlist(paras),collapse="-")
  file_n=grep(loop,grep(best_para,all_models,value=TRUE),value=TRUE)
  print(file_n)
  if(length(file_n)>3){
    file_n=file_n[1:3]
  }
  if(loop=="H2_10"){
    file_n=file_n[1]
  }
  #if(length(file_n)>1){
  #  file_n=grep("xgboost",file_n,value=TRUE)
  #}

  model_file=paste("./rmsd_cluster_hits_rmsd/",file_n,sep = "")
  model=readRDS(model_file[1])
    model_results_list=lapply(model_file,function(x){
    model=readRDS(x)
    model_result=model$pred
    re= model_result[,c("pred", "obs","Resample")]; #print(re)
    names(re)=c("pred", "obs","Resample")
    return(re)
    })
    
    if(length(model_results_list)==1){
    model_re=do.call(rbind, model_results_list)
    model_re_resample = split(model_re, model_re$Resample)
    resample_chunks=chunk2(1:length(model_re_resample),3)
    results_by_repeats=lapply(resample_chunks,function(x){
      result=do.call(rbind,model_re_resample[x])
      
    })
    
    }else{
  # process irregular cases, overwrite the repeats numbering 
    results_by_repeats=model_results_list
    
  }
      
      
      for(l in 1:length(results_by_repeats)){
        the_f = results_by_repeats[[l]]
        the_f$repeats  = rep(l,dim(the_f)[1])
        results_by_repeats[[l]]=the_f
      }
      model_re=do.call(rbind, results_by_repeats)
    
    # split the results into 3 repeats 

  gbm_folds_sd[[loop]]= model$results$AccuracySD
  


  blindBLAST_unique_clusters=gsub("\\*","none",unique(data_by_loop_type_list_unduplicated_for_blindBLAST[[loop]][[1]]$cluster_type))

 returned_l=convert_none(model_re,blindBLAST_unique_clusters)
 obs=returned_l[[1]]
 pred=returned_l[[2]]
  obs=factor(obs,levels=blindBLAST_unique_clusters)
  pred=factor(pred,levels=blindBLAST_unique_clusters)
  
  model_re[,1]=obs
  model_re[,2]=pred
  model_re_resplit=split(model_re,model_re$repeats)
  all_gbm_pred_result[[loop]]=model_re

  accuracies= lapply(model_re_resplit,function(x){
    dim(x[x[,1]==x[,2],])[1]/dim(x)[1]
  })
  
  errorcount_sd= lapply(model_re_resplit,function(x){
    returned_l=convert_none(x,blindBLAST_unique_clusters)
    obs=returned_l[[1]]
    pred=returned_l[[2]] 
    obs=factor(obs,levels=blindBLAST_unique_clusters)
    pred=factor(pred,levels=blindBLAST_unique_clusters)
    conf_t=as.data.frame(table(obs,pred))
    conf_t=conf_t[conf_t$Freq>0 & conf_t[,1]!=conf_t[,2],]
    conf_t[,1]=gsub("_","-",conf_t[,1])
    conf_t[,2]=gsub("_","-",conf_t[,2])
    sum(conf_t$Freq)
    })
  
  errorcount_by_mis_sd= lapply(model_re_resplit,function(x){
    returned_l=convert_none(x,blindBLAST_unique_clusters)
    obs=returned_l[[1]]
    pred=returned_l[[2]] 
    obs=factor(obs,levels=blindBLAST_unique_clusters)
    pred=factor(pred,levels=blindBLAST_unique_clusters)
    conf_t=as.data.frame(table(obs,pred))
    conf_t=conf_t[conf_t$Freq>0 & conf_t[,1]!=conf_t[,2],]
    conf_t[,1]=gsub("_","-",conf_t[,1])
    conf_t[,2]=gsub("_","-",conf_t[,2])
    return(conf_t)
  })
  errorcount_by_mis_fr=errorcount_by_mis_sd[[1]]
  for(x in 2:length(errorcount_by_mis_sd)){
    errorcount_by_mis_fr= merge(errorcount_by_mis_fr,errorcount_by_mis_sd[[x]],by=c("obs","pred"))
  }
  sds=apply(errorcount_by_mis_fr,1,function(x){
   sd(as.numeric( unlist(x[3:(3+length(errorcount_by_mis_sd)-1)])))
  })
  errorcount_by_mis_fr=cbind(errorcount_by_mis_fr[,1:2],sds)
  
  gbm_errorcount_list[[loop]]=unlist(errorcount_sd)
  gbm_errorcount_sd_list[[loop]]=sd(unlist(errorcount_sd))
  print(accuracies)
  gbm_sd_by_loops[[loop]]=sd(unlist(accuracies))
  gbm_accuracy_by_loops[[loop]]=mean(unlist(accuracies))
  conf_t=as.data.frame(table(obs,pred))
  colnames(conf_t)[1:2]=c("Var1","Var2")
  conf_t$Freq= conf_t$Freq/3   # average by the repeat number 
  
  conf_tables_all_loops_gbm[[loop]]=conf_t
  conf_t=conf_t[conf_t$Freq>0 & conf_t[,1]!=conf_t[,2],]
  colnames(errorcount_by_mis_fr)=c("Var1", "Var2" ,"sds")
  conf_t=merge(conf_t,errorcount_by_mis_fr)
  conf_t[,1]=gsub("_","-",conf_t[,1])
  conf_t[,2]=gsub("_","-",conf_t[,2])
  conf_tables_all_loops_gbm_diff[[loop]]=conf_t
}

save_file("gbm_errorcount_list")
#result_lists=calculate_accuracy_mean_std_by_repeats(all_gbm_pred_result)
save_file("conf_tables_all_loops_gbm_diff")
save_file("conf_tables_all_loops_gbm")
save_file("gbm_folds_sd")
save_file("gbm_sd_by_loops")
 save_file("gbm_accuracy_by_loops")
#compare the conf_table from gbm to that of blindBLAST
rbind_all_gbm_pred_result=do.call(rbind, all_gbm_pred_result)
gbm_pred_result_by_repeats=split(rbind_all_gbm_pred_result,rbind_all_gbm_pred_result$repeats)
gbm_acc_by_repeats=lapply(gbm_pred_result_by_repeats,function(x){
  dim(x[x[,"pred"]==x[,"obs"],])[1]/dim(x)[1]
})
sd(unlist(gbm_acc_by_repeats))
mean(unlist(gbm_acc_by_repeats))
#all_dist_matrix[[1]]["1a0q","5j1s"]

