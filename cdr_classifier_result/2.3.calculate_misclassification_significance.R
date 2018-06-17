# This script calculate the significance of any misclssification based on the cluster random assignment result. 
current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

most_similar_score_result=list()
all_loop_sig_frame_list=list()
all_dists=lapply(all_dist_matrix,function(x){max(unlist(x),na.rm=TRUE)-min(unlist(x),na.rm=TRUE)})
all_neighborhood_dist=unlist(all_dists)/6



for(each_l in names(all_similarity_matrix)){
  #get sequences 
   g(sequences, this_simi, this_dis) %=% getting_similar_sequences_similarity_matrix_rmsd_matrix(each_l)
if(!is.data.frame(sequences)){if(!sequences){next}}
    features = data_by_loop_type_list_unduplicated[[each_l]][[4]]


  # extract the wrong cases of blindBLAST in this loop 
  this_loop_wrong_cases=count_number_for_misclassification(each_l)  
  if(is.null(this_loop_wrong_cases)){next}
  all_loop_sig_frame_list[[each_l]]=get_significance(this_loop_wrong_cases)
  
  sig_data_frame=""
  all_the_values = ""
  gc()


}
save_file("all_loop_sig_frame_list")
