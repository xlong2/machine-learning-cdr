# perform 10 fold 3 repeats random prediction simulation
current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
    source("0.load_function_and_data.R")
    
}

all_loops=names(data_by_loop_type_list_unduplicated)
all_similarity_matrix=mclapply(all_loops,get_pairwise_sequence_simi_self_cal,mc.cores=3)
names(all_similarity_matrix)=all_loops
saveRDS(all_similarity_matrix,file="./Data_processed/all_similarity_matrix.rds")
# calculate the misclassifications counts
