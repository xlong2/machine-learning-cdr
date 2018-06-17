check_index<-function(split_index){
  old_files=list.files(pattern="*extra_test.rds",path="/home/xlong/machine_learning_cdr/proline_classifier/rmsd_cluster_hits_rmsd",full.names = FALSE)
  print(c(old_files,"line 176 "))
  old_files=grep(loop_type,old_files,value=TRUE)
  print(old_files)
  old_indexes=sapply(old_files,function(x){as.numeric(strsplit(strsplit(as.character(x),"-")[[1]][3],"_")[[1]][1])})
  print(c(old_indexes,"OLD_INDEXES "))
  split_index=as.numeric(split_index)
  if(split_index %in% old_indexes){ split_index=split_index+1; check_index(split_index)}else{
    print(to_save_file)
    saveRDS(trained_model,file =to_save_file)}  }


all_similarity_matrix_names=names(all_loops_similarity_matrix)
all_loops_similarity_matrix=lapply(names(all_loops_similarity_matrix),function(x){
  this_matrix=all_loops_similarity_matrix[[x]]
  this_matrix_clusters=split_vector_and_replace(colnames(this_matrix),"\\.",1,1,"")
  this_matrix_pdbs=split_vector_and_replace(colnames(this_matrix),"\\.",2,3,".")
  unique_clusters=unique(data_by_loop_type_list_unduplicated[[x]][[1]]$cluster_type)
  the_cluster=paste(c(split_vector_and_replace(x,"_",1,2,"-"),"none"),collapse="-")
  this_matrix_clusters[!this_matrix_clusters%in%unique_clusters]=rep(the_cluster,length(this_matrix_clusters[!this_matrix_clusters%in%unique_clusters]))
  colnames(this_matrix)=paste(this_matrix_clusters,this_matrix_pdbs,sep=".")
  rownames(this_matrix)=paste(this_matrix_clusters,this_matrix_pdbs,sep=".")
  return(this_matrix)
})
names(all_loops_similarity_matrix)=all_similarity_matrix_names
save_file("all_loops_similarity_matrix")



# extract legend
plot=ggplot(loop_distribution_total)+facet_wrap(~loop_type, scales = "free" )+ 
  geom_bar(aes(x=loop,y=case_number,fill=cluster_identifier,width=0.8),position =  position_dodge(width = 0.90),stat = "identity",size=0.2)+ 
  ggtitle("Canonical cdr loop cluster distribution")+theme(plot.title = element_text(hjust = 0.5))
