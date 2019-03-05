# find which CDR loops are wrongly classified by blindBLAST, and enforce template searching in the correct clusters for these CDRs
# for the enforced correct query-template CDR pairs, retrieve the corresponding RMSDs.


current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

get_accuracy_per_fold_enforcing_corrent_fold


enforcing_correct_rmsd_list=list()
for (loop_type in names(ten_foldcv_blindblastlist)[1]) {
  
  returned_list=getting_similar_sequences_similarity_matrix_rmsd_matrix(loop_type)
  if(!is.logical(returned_list)){
    g(sequences,similarity_matrix,rmsd_matrix)  %=%  returned_list
  }else{next}
  the_result_list=ten_foldcv_blindblastlist[[loop_type]]
  
  features=data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  this_loop_enforcing_correct_rmsds=lapply(1:length(the_result_list),get_accuracy_per_fold_enforcing_corrent_fold)
  enforcing_correct_rmsd_list[[loop_type]]=this_loop_enforcing_correct_rmsds
  if(is.null(enforcing_correct_rmsd_list[[loop_type]])){
    enforcing_correct_rmsd_list[[loop_type]]<-NULL
  }
}



save_file("enforcing_correct_rmsd_list")
