library("caret")
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
args=c(2, "L3_9", 4, 6, 1500, 0.01, 5)
split_index=args[1]
num_core=args[3]
loop_type=args[2]
interaction.depth=args[4]
n.trees=args[5]
shrinkage=args[6]
n.minobsinnode=args[7]
subsitution_matrix="PAM30"

registerDoMC(num_core)

#gbmGrid=splitted_grid[[split_index]]
gbmGrid=data.frame(interaction.depth=as.numeric(interaction.depth),n.trees=as.numeric(n.trees),shrinkage=as.numeric(shrinkage),n.minobsinnode=as.numeric(n.minobsinnode))

parameter_spe = paste(unlist(gbmGrid),collapse="-")
each_method="gbm_test"
cluster_dis="north"




  data_by_loop_type_list_unduplicated=readRDS(paste(c(overall_prefix,"machine_learning_cdr/proline_classifier/data_by_loop_type_list_unduplicated_no_filtering.rds"),collapse = ""))
  
  gbm_result_dir=paste(c("./proline_classifier/rmsd_cluster_hits_rmsd/"),collapse="")

to_save_file=paste(c(gbm_result_dir,"",loop_type,"_",paste(c(each_method,cluster_dis,parameter_spe),collapse="-"),"_trained_model_extra_test.rds"),collapse="")

if(file.exists(to_save_file)){
print("file_already_exists!")}else{
data=data_by_loop_type_list_unduplicated[[loop_type]][[1]]
data$cluster_type=sub("-","_",data$cluster_type)
data$cluster_type=as.factor(as.character(data$cluster_type))




#sequences$rmsd_cluster = as.character(sequences$rmsd_cluster)
sequences=data

features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
each_loop_length_data_feature_string=as.formula(paste(c("cluster_type ~ ",paste(features,collapse=" + ")),collapse=""))
all_cases =  sequences[,c(features,"id","cluster_type")]
the_levels=unique(unlist(data_by_loop_type_list_unduplicated[[loop_type]][[1]][,features]))
for(each_f in features){
  all_cases[,each_f]=factor(all_cases[,each_f],levels=the_levels)
}
all_cases=all_cases[complete.cases(all_cases), ]
all_cases$cluster_type=gsub("-","_",all_cases$cluster_type)
all_cases$cluster_type=gsub(",",".",all_cases$cluster_type)

all_cases$cluster_type=as.factor(as.character(all_cases$cluster_type))

# tune the parameter
trained_model = generic_train(each_loop_length_data_feature_string_rmsd,each_method,all_cases) 


pred_result = trained_model$pred # get the prediction vs observation
conf_table = table(all_cases$cluster_type, predict(trained_model,all_cases)) # get the confusion table and save it

print("finished the cluster prediction")
print(conf_table)

saveRDS(trained_model,file =to_save_file)



}
