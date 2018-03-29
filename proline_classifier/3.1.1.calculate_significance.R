

# perform 10 fold 3 repeats random prediction simulation

all_significance_simulation=list()  # the list to record all the test count 
for(each_l in names(all_similarity_matrix)){
  sequences = data_by_loop_type_list_unduplicated[[each_l]][[1]]#load the sequences from a loop length
  features = data_by_loop_type_list_unduplicated[[each_l]][[4]]
  case_ids= sequences[,"identifier"]
  splitted_sequence=split(sequences,sequences$cluster_type)
  clusters=unique(sequences$cluster_type)
  prediction_it=list()
  iteration_n=1000
  for(it in 1:iteration_n){
    print(c(each_l,it/iteration_n))
    all_clu_re=list()
    for(cluster in clusters){
      ind_cl=which(sequences$cluster_type==cluster)
      this_ind_most_sim_ind_list=list()
      for(this_ind in ind_cl){
        relevant_indexes=(1:dim(sequences)[1])[1:dim(sequences)[1]!=ind_cl]
        
        this_ind_most_sim_ind=sample(relevant_indexes,1,replace=TRUE)
        this_ind_most_sim_ind_list=c(this_ind_most_sim_ind_list,this_ind_most_sim_ind)
      }
      final_table=data.frame(sequences[ind_cl,"cluster_type"],sequences[unlist(this_ind_most_sim_ind_list),"cluster_type"])
      all_clu_re[[cluster]]=final_table
      
    }
    prediction_it[[it]]=as.data.frame(table(do.call(rbind,all_clu_re)))
    all_clu_re=""
    gc()
  }# end of 1000 iteration for this loop
  
  prediction_all=do.call(rbind, prediction_it)
  prediction_all$error_type=paste(prediction_all[,1],prediction_all[,2],sep="")
  all_the_values=split(prediction_all[,"Freq"],prediction_all[,"error_type"])
  prediction_it=""
  
  all_significance_simulation[[each_l]]=all_the_values
  all_the_values=""
  gc()
}
save_file("all_significance_simulation")
