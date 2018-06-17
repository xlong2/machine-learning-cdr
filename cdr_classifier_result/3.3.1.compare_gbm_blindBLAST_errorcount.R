# This script compare the blindBLAST identification result to the GBM prediction result.

current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

#compare the blindBLAST to GBM specifically
diff_between_blast_gbm=list()
for(loop in names(ten_foldcv_blindblastlist)){
  if((!loop%in%names(conf_tables_all_loops_blindBLAST_diff)) |(!loop %in%names(conf_tables_all_loops_gbm_diff))){next}
  blind_r=conf_tables_all_loops_blindBLAST_diff[[loop]]
  colnames(blind_r)[4]="sds"
  gbm_r=conf_tables_all_loops_gbm_diff[[loop]]
  colnames(gbm_r)[4]="sds"
  
  if(!is.null(gbm_r)&!is.null(blind_r)){
    combined_f=merge(blind_r,gbm_r,all=TRUE,by=c("Var1","Var2"))
    combined_f[is.na(combined_f)] <- 0
    combined_f$gbm_improved=combined_f$Freq.x-combined_f$Freq.y
    combined_f=combined_f[order(combined_f$gbm_improved,decreasing = TRUE),]
    colnames(combined_f)=c("query_cluster","template_cluster","blindBLAST_error_count","blindBLAST_sd","gbm_error_count","gbm_sd","improvement_by_gbm")
    diff_between_blast_gbm[[loop]]=combined_f
  }
}
write_list_into_single_csv(diff_between_blast_gbm,paste(c(result_dir,"diff_between_blast_gbm_error_count.csv"),collapse = ""))
print(diff_between_blast_gbm)


# find clusters that have smaller than 50 cases
all_sparse_names=lapply(data_by_loop_type_list_unduplicated_for_blindBLAST,function(x){
  a=table(x[[1]]$cluster_type)
  names(a[a<50])
  })
results=lapply(names(all_sparse_names),function(x){
  print(x)
  fr=diff_between_blast_gbm[[x]]
  int_f=fr[fr$query_cluster%in%all_sparse_names[[x]],]
  print(int_f)
  if(!is.null(int_f)){
  aa=split(int_f,as.character(int_f$query_cluster))
  total_imp=lapply(aa,function(x){sum(x$improvement_by_gbm)})
  return(aa)}else{
    return(NULL)
    }
  
})

total_error_count_for_sparse_clusters=lapply(names(all_sparse_names),function(x){
  print(x)
  fr=diff_between_blast_gbm[[x]]
  int_f=fr[fr$query_cluster%in%all_sparse_names[[x]],]
  print(int_f)
  if(!is.null(int_f)){
    aa=split(int_f,as.character(int_f$query_cluster))
    total_imp=lapply(aa,function(x){sum(x$improvement_by_gbm)})
    return(total_imp)}else{
      return(NULL)
    }
  
})

names(total_error_count_for_sparse_clusters)=names(all_sparse_names)

total_error_count_for_sparse_clusters=total_error_count_for_sparse_clusters[!unlist(lapply(total_error_count_for_sparse_clusters,is.null))]
total_error_count_for_sparse_clusters=total_error_count_for_sparse_clusters[names(total_error_count_for_sparse_clusters)!="L1_16"]
error_count_for_sparse_clusters=unlist(total_error_count_for_sparse_clusters,recursive = FALSE)
barplot(table(as.numeric(error_count_for_sparse_clusters)))
length(which(error_count_for_sparse_clusters>0))/length(error_count_for_sparse_clusters)
length(which(error_count_for_sparse_clusters<=0))/length(error_count_for_sparse_clusters)

names(results)=names(all_sparse_names)

results=results[!unlist(lapply(results,is.null))]
results=results[names(results)!="L1_16"]
names(results)
results_unlist=unlist(results,recursive = FALSE)
write_list_into_single_csv_sp(results_unlist,"./Data_processed/sparse_cluster_error_counts.csv")
  



#compare blindBLAST result to GBM result 
options(digits=2)
save_file("diff_between_blast_gbm")

