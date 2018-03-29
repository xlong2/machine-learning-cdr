all_loops=names(data_by_loop_type_list_unduplicated)
all_similarity_matrix=mclapply(all_loops,get_pairwise_sequence_simi_self_cal,mc.cores=3)
names(all_similarity_matrix)=all_loops
saveRDS(all_similarity_matrix,file="./proline_classifier/Data_processed/all_similarity_matrix.rds")
# calculate the misclassifications counts