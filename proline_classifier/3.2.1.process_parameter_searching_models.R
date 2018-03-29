# Purpose
#   This script harvest all the machine learning models for searching the grids, and find the best tuned grids out of all
#   the models. 
all_models=list.files(pattern="*extra_test.rds",path="./proline_classifier/rmsd_cluster_hits_rmsd",full.names = FALSE)

all_models_list_by_loop=list()
index_i=1
all_t=list()
for(each_file in all_models){
  each_file_r=paste(c("./proline_classifier/rmsd_cluster_hits_rmsd/",each_file),collapse="")
  the_model=readRDS(each_file_r)
  
  index_i=index_i+1
print(each_file)
each_file=gsub("xgboost-north","gbm_test-north",each_file)
info=strsplit(each_file,"\\/")[[1]][length(strsplit(each_file,"\\/")[[1]])]
  loop_type=paste(strsplit(info, "_")[[1]][1:2],collapse="_")
  all_t[[loop_type]]=unique(the_model$trainingData$.outcome)
 number =paste(strsplit(strsplit(info, "_")[[1]][4],"-")[[1]][3:6],collapse="-")
 

    all_models_list_by_loop[[loop_type]][[number]]=the_model$result
    rm(the_model)
gc()
  #=strsplit(paste(strsplit(info, "_")[[1]][4]
}
all_result=lapply(all_models_list_by_loop,function(x){do.call(rbind,x)})
save_file("all_models_list_by_loop")

names(all_result[[1]])[1]=c("n.trees")
names(all_result[[1]])[4]="interaction.depth"
names(all_result[[1]])[5]="shrinkage"
names(all_result[[1]])[2]="n.minobsinnode"
all_result[[1]][,2]=rep(3,dim(all_result[[1]])[1])
#names_interested=names(all_result[[2]])[c(1:8,15,16,18,19)]
names_interested_o=names(all_result[[2]])[c(1:8)]
names_interested=c(names_interested_o, "Mean_Balanced_Accuracy", "logLossSD", "AUCSD", "AccuracySD",   "KappaSD","loop_type")
all_result=lapply(all_result,function(x){if(!"Mean_Balanced_Accuracy" %in% names(x)){names(x)[which(names(x)=="Balanced_Accuracy")]="Mean_Balanced_Accuracy" }; return(x)})

results=do.call(rbind,lapply(names(all_result),function(x){all_result[[x]]$loop_type=rep(x,dim(all_result[[x]])[1]); return(all_result[[x]][,names_interested])}))
ggplot(results)
results=results[!grepl("NA",rownames(results)),]
saveRDS(results,file="../results.rds")
gc()
rm(list())


results=readRDS("result.rds")
#overall_accuracy=readRDS("overall_accuracy.rds")
#blind_blast_accuracy=as.data.frame(do.call(rbind,overall_accuracy))
#blind_blast_accuracy$loop=rownames(blind_blast_accuracy)
#colnames(blind_blast_accuracy)=c("mean","sd","loop_type")
#blind_blast_accuracy=blind_blast_accuracy[!rownames(blind_blast_accuracy) %in% c("H1_14","H1_15","H1_10","L1_10","L2_12","L3_12"),]
# extract the 
all_results=results
results=all_results[all_results$shrinkage>0.001,]
results_smaller=all_results[all_results$shrinkage<=0.001,]

ggplot(results, aes(x=n.trees, y=Accuracy,group=interaction.depth, color=interaction.depth,ymin=Accuracy-AccuracySD/2, ymax=Accuracy+AccuracySD/2))+geom_errorbar( width=0.9)  + geom_line()+facet_wrap(~loop_type,scales="free" )+ggtitle("Gradient Boost Machine model complexity tuning") +theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(data = blind_blast_accuracy, aes(yintercept = mean))
file_name="/Users/xlong3/lab_work_data/machine_learning_cdr/proline_classifier/gbm_grid_search.png"

ggsave(file=file_name,  width = 300, height = 200,units = c( "mm"))
results_smaller_5=results_smaller[results_smaller$n.minobsinnode==5,]
results_smaller_2=results_smaller[results_smaller$n.minobsinnode==2,]

ggplot(results_smaller_2, aes(x=n.trees, y=Accuracy,group=interaction.depth, color=interaction.depth,ymin=Accuracy-AccuracySD/2, ymax=Accuracy+AccuracySD/2))+geom_errorbar( width=0.9)  + geom_line()+facet_wrap(~loop_type,scales="free" )+ggtitle("Gradient Boost Machine model complexity tuning") +theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(data = blind_blast_accuracy, aes(yintercept = mean))

ggplot(results_smaller_5, aes(x=n.trees, y=Accuracy,group=interaction.depth, color=interaction.depth,ymin=Accuracy-AccuracySD/2, ymax=Accuracy+AccuracySD/2))+geom_errorbar( width=0.9)  + geom_line()+facet_wrap(~loop_type,scales="free" )+ggtitle("Gradient Boost Machine model complexity tuning") +theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(data = blind_blast_accuracy, aes(yintercept = mean))



#get test result 
library("plot3D")
best_parameters_each_loop=list()
all_gird_search_results=list()
result_table=as.data.frame(matrix(nrow=length(names(all_models_list_by_loop)),ncol=length(c("branch_level","trees","shrinkage","min_node","accuracy","accuracy_sd"))))
colnames(result_table)=c("branch_level","trees","shrinkage","min_node","accuracy","accuracy_sd")
rownames(result_table)=names(all_models_list_by_loop)
for(each_loop in names(all_models_list_by_loop)){
  all_gird_search_results[[each_loop]]=do.call(rbind,all_models_list_by_loop[[each_loop]])
  
  # do a plot 
  x=unique( all_gird_search_results[[each_loop]]$n.trees)
  y=sort(unique( all_gird_search_results[[each_loop]]$interaction.depth))
  z=matrix(nrow=length(x),ncol=length(y))
  z=as.data.frame(z)

  rownames(z)=x; colnames(z)=y
  # only select the mino equal to 5
  sub_table=all_gird_search_results[[each_loop]][all_gird_search_results[[each_loop]]$n.minobsinnode==5,]
  for(ind in 1:dim(sub_table)[1]){
    z[as.character(sub_table[ind,"n.trees"]),as.character(sub_table[ind,"interaction.depth"])]=sub_table[ind,"Accuracy"]
  }
  

  max_indses=which(all_gird_search_results[[each_loop]]$Accuracy==max(all_gird_search_results[[each_loop]]$Accuracy))
  max_indses_index=which(all_gird_search_results[[each_loop]][max_indses,"Kappa"]==max(all_gird_search_results[[each_loop]][max_indses,"Kappa"]))
  final_max_ind=max_indses[max_indses_index]
  if(length(final_max_ind)!=1){
    which_i=which(all_gird_search_results[[each_loop]][final_max_ind,"n.trees"]==min(all_gird_search_results[[each_loop]][final_max_ind,"n.trees"]))
    final_max_ind=final_max_ind[which_i]
    if(length(final_max_ind)!=1){
      which_i=which(all_gird_search_results[[each_loop]][final_max_ind,"interaction.depth"]==min(all_gird_search_results[[each_loop]][final_max_ind,"interaction.depth"]))
      final_max_ind=final_max_ind[which_i]
      
    }
    
  }
  
  parameters= rownames(all_gird_search_results[[each_loop]])[final_max_ind]
  parameters= t(data.frame(strsplit(as.character(parameters), "-")))
  
  result_table[each_loop,c("branch_level","trees","shrinkage","min_node")]=parameters
  result_table[each_loop,c("accuracy","accuracy_sd")]=all_gird_search_results[[each_loop]][final_max_ind,c("Accuracy","AccuracySD")]
  best_parameters_each_loop[[each_loop]]=parameters
 # slice3D(all_gird_search_results[[each_loop]]$n.trees, all_gird_search_results[[each_loop]]$shrinkage, all_gird_search_results[[each_loop]]$interaction.depth,colvar=all_gird_search_results[[each_loop]]$Accuracy)
#  M <- mesh(all_gird_search_results[[each_loop]]$n.trees, all_gird_search_results[[each_loop]]$shrinkage, all_gird_search_results[[each_loop]]$interaction.depth)
}

rownames(result_table)=split_vector_and_replace(rownames(result_table),"_",1,2,"-")
result_table=result_table[order(order_factor_by_two_component(rownames(result_table),"-",1,2),rownames(result_table)),]
write.csv(result_table,file=paste(result_dir,"best_parameters.csv"),row.names = TRUE,col.names  = TRUE)
save_file("best_parameters_each_loop")
# print the err


