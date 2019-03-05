current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

each_method = "blindblast_just_get_alignment"
cluster_dis = "north"
subsitution_matrix_name = "wahtw"
subsitution_matrix = "PAM30"


# do blindBLAST template searching for each CDR loop 
overall_accuracy = list()
no_foldcv_blindblastlist = list()
the_loop_names=names(data_by_loop_type_list_unduplicated)[1:length(names(data_by_loop_type_list_unduplicated))]
left_loops=the_loop_names[!the_loop_names%in% names(no_foldcv_blindblastlist) ]
for (loop_type in the_loop_names) {
  returned_list=getting_similar_sequences_similarity_matrix_rmsd_matrix(loop_type)
  if(!is.logical(returned_list)){
    g(sequences,similarity_matrix,rmsd_matrix)  %=%  returned_list
  }else{print("jump to next");next;}
  sequences = data_by_loop_type_list_unduplicated[[loop_type]][[1]]
  each_loop_length_data_feature_string = data_by_loop_type_list_unduplicated[[loop_type]][[2]]
  features=data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  sequences$cluster_type = as.character(sequences$cluster_type)
  all_cases =  sequences[, c(features, "cluster_type", "identifier", "id")]
  #make folds and separates folds in and folds out
  r <- 3 # number of repeats
  k <- 10 # number of folds
  g(folds.list,folds.list.out)  %=%generate_folds_foldsout(sequences,r,k)
  
  all_result_list = mclapply(1:length(folds.list), get_accuracy_per_fold_overload_final, mc.cores =    4)
  
  no_foldcv_blindblastlist[[loop_type]] = do.call(rbind, all_result_list)
  
}

save_file("no_foldcv_blindblastlist")


# get the error count 
error_c_list = list()
for(each_loop in names(no_foldcv_blindblastlist)){
  fr = no_foldcv_blindblastlist[[each_loop]]
  query_c = sapply(strsplit(fr[,1],"\\."),"[[",1)
  pred_c = sapply(strsplit(fr[,2],"\\."),"[[",1)
  error_c_list[[each_loop]] =as.data.frame(table(query_c, pred_c))
}

