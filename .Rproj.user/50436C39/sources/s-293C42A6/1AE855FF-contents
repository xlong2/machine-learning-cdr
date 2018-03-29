# get the summary by category 
simulated_prediction=list()
r <-3 # number of repeats
k <- 10 # number of folds

for(loop_type in names(data_by_loop_type_list_unduplicated)){
  
  # tryCatch({
  sequences = data_by_loop_type_list_unduplicated[[loop_type]][[1]]#load the sequences from a loop length
  seq_dim=dim(sequences)[1]
  features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  each_loop_length_data_feature_string = data_by_loop_type_list_unduplicated[[loop_type]][[2]]
  sequences$cluster_type=as.character(sequences$cluster_type) 
  all_cases =  sequences[,c(features,"cluster_type","identifier","id")]
  #make folds and separates folds in and folds out 
  training_cases=sequences
  training_cases=sequences[,c(features,"cluster_type","id","identifier")]
  all_results=list()
  print(loop_type)
  #for(sim_i in 1:num_it){
  #  print(sim_i)
    returned_results=make_3_10_cross_val(training_cases,r,k)
    divisions_out=returned_results[[1]][[1]]
    divisions_in=returned_results[[1]][[2]]
    this_iter_result=list()
    for(i in 1:30){
      divisions_out[[i]]
      divisions_in[[i]]
      result=sample(divisions_in[[i]],length(divisions_out[[i]]),replace =TRUE)
      result_clus=training_cases[result,"cluster_type"]
      this_iter_result[[i]]=data.frame(query_cluster=training_cases[divisions_out[[i]],"cluster_type"],result=result_clus)
    }
  #  all_results[[sim_i]]=do.call(rbind,this_iter_result)
#  }
    all_results_frame=do.call(rbind,this_iter_result)
  
 # all_results_frame=do.call(rbind,all_results)
  wrong_pred_all_results_frame=all_results_frame[all_results_frame[,1]!=all_results_frame[,2],]
  
  wrong_pred_all_results_frame_table=as.data.frame(table(wrong_pred_all_results_frame[,1],wrong_pred_all_results_frame[,2])/r)
  
  wrong_pred_all_results_frame_table=wrong_pred_all_results_frame_table[order(wrong_pred_all_results_frame_table$Freq),]
  wrong_pred_all_results_frame_table$ratio=wrong_pred_all_results_frame_table$Freq/seq_dim
  
  simulated_prediction[[loop_type]]=wrong_pred_all_results_frame_table
  gc()
}   

save_file("simulated_prediction")
  

# code the similarity comparison with x be the mean similarity between the cases in one cluster and its corresponding sequence pair of the most
