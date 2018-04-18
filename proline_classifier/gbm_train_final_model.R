
args <- commandArgs(trailingOnly = TRUE)
setwd("..")
system("Rscript ./proline_classifier/0.load_function_and_data.R")
execute_training_rscript<-function(args){
  loop_type=args[1]
  num_core=args[2]
  interaction.depth=args[3]
  n.trees=args[4]
  shrinkage=args[5]
  n.minobsinnode=args[6]
  registerDoMC(num_core)
  
  gbmGrid=data.frame(interaction.depth=as.numeric(interaction.depth),n.trees=as.numeric(n.trees),shrinkage=as.numeric(shrinkage),n.minobsinnode=as.numeric(n.minobsinnode))
  
  parameter_spe = paste(unlist(gbmGrid),collapse="-")
  each_method="gbm_test"
  gbm_result_dir=paste(c("./proline_classifier/rmsd_cluster_hits_rmsd/"),collapse="")
  
  to_save_file=paste(c(gbm_result_dir,"",loop_type,"_",paste(c(each_method,cluster_dis,parameter_spe),collapse="-"),"_final_model.rds"),collapse="")
  
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
      all_cases$cluster_type=gsub("\\*","none",all_cases$cluster_type)
      
      all_cases$cluster_type=as.factor(as.character(all_cases$cluster_type))
      
      # tune the parameter

      trained_model = train_final_model(each_loop_length_data_feature_string,each_method,all_cases,gbmGrid) 
      saveRDS(trained_model,file =to_save_file)
      
      
      
    }
}
