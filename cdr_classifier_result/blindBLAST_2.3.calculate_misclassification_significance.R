# This script calculate the significance of any misclssification based on the cluster random assignment result. 
current_d=getwd()

if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

most_similar_score_result=list()
all_loop_sig_frame_list=list()
all_dists=lapply(all_dist_matrix,function(x){max(unlist(x),na.rm=TRUE)-min(unlist(x),na.rm=TRUE)})
all_neighborhood_dist=unlist(all_dists)/6



for(each_l in names(data_by_loop_type_list_unduplicated_modified)){
  #get sequences 
 #  g(sequences, this_simi, this_dis) %=% getting_similar_sequences_similarity_matrix_rmsd_matrix(each_l)
#if(!is.data.frame(sequences)){if(!sequences){next}}
    features = data_by_loop_type_list_unduplicated_modified[[each_l]][[4]]


  # extract the wrong cases of blindBLAST in this loop 
    this_loop_wrong_cases= conf_tables_all_loops_blindBLAST[[each_l]]
    this_loop_wrong_cases$Var1 = as.character(this_loop_wrong_cases$Var1)
    this_loop_wrong_cases$Var2 = as.character(this_loop_wrong_cases$Var2)
    
    this_loop_wrong_cases$misclassification = paste(this_loop_wrong_cases$Var1, this_loop_wrong_cases$Var2, sep = "")
    this_loop_wrong_cases$sd<-NULL
    this_loop_wrong_cases= this_loop_wrong_cases[this_loop_wrong_cases$Var1!=this_loop_wrong_cases$Var2,]
    sig_data_frame = as.data.frame(matrix(nrow = length(this_loop_wrong_cases$misclassification), ncol = 8))
    colnames(sig_data_frame) =c("error_type","Var1","Var2","significance","error_count","sd","mean_simu_error","effect_size")
    rownames(sig_data_frame) = this_loop_wrong_cases$misclassification
  for (misclassification in this_loop_wrong_cases$misclassification) {
    
    this_mis_info = this_loop_wrong_cases[this_loop_wrong_cases$misclassification ==
                                            misclassification, ]
    
    
    this_error = this_mis_info["Freq"]
    names(all_significance_simulation[[each_l]])= gsub("\\*","none",names(all_significance_simulation[[each_l]]))
    if(is.null(all_significance_simulation[[each_l]][[misclassification]])){next()}
    simu_error = all_significance_simulation[[each_l]][[misclassification]]
    aycdf <- ecdf(simu_error)
    sig_data_frame[misclassification, "error_type"] = misclassification
    sig_data_frame[misclassification, "Var1"] = this_mis_info["Var1"]
    sig_data_frame[misclassification, "Var2"] = this_mis_info["Var2"]
    sig_data_frame[misclassification, "significance"] = aycdf(this_error)
    sig_data_frame[misclassification, "error_count"] = this_error
    sig_data_frame[misclassification, "sd"] = sd(simu_error)
    sig_data_frame[misclassification, "mean_simu_error"] = mean(simu_error)
    sig_data_frame[misclassification, "effect_size"] = (this_error - mean(simu_error)) /
      sd(simu_error)
    
  }
    all_loop_sig_frame_list[[each_l]] = sig_data_frame
  sig_data_frame=""
  all_the_values = ""
  gc()


}
readRDS("all_loop_sig_frame_list")

