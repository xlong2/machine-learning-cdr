
all_loops_similarity_matrix_by_self=list()
all_loops=names(data_by_loop_type_list_unduplicated)[order(unlist(lapply(data_by_loop_type_list_unduplicated,function(x){dim(x[[1]])[1]})),decreasing=TRUE)]
all_similarity_matrix=mclapply(all_loops,get_pairwise_sequence_simi_self_cal,mc.cores=3)
for(loop_type in names(data_by_loop_type_list_unduplicated)){
  
  bb=get_pairwise_sequence_simi_self_cal(sequences,features)
  all_loops_similarity_matrix_by_self[[loop_type]]=bb
  #saveRDS(bb,file=paste(c("all_loops_similarity_matrix",loop_type,".rds"),collapse=""))
  gc()
}



all_loops_similarity_matrix=list()
for(loop_type in names(data_by_loop_type_list_unduplicated)){
  sequences = data_by_loop_type_list_unduplicated[[loop_type]][[1]]#load the sequences from a loop length
  features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  bb=get_pairwise_sequence_simi(sequences,features)
  saveRDS(bb,file=paste(c("all_loops_similarity_matrix",loop_type,".rds"),collapse=""))
  gc()
}


for(loop_type in names(data_by_loop_type_list_unduplicated)){
  
  all_loops_similarity_matrix[[loop_type]]=readRDS(file=paste(c("all_loops_similarity_matrix",loop_type,".rds"),collapse=""))
  
}
saveRDS(all_loops_similarity_matrix,file="../all_loops_similarity_matrix.rds")
