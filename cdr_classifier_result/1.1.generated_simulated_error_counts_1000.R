# perform  random cluster assignment 1000 times
current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

all_significance_simulation=list()  # the list to record all the test count 
for(each_l in names(all_similarity_matrix)[18:19]){
  tryCatch({
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
  },error=function(e){print(each_l)})
}

all_significance_simulation=lapply(all_significance_simulation,function(x){
  x_n=names(x); 
  c_x_n=gsub("\\*","none",x_n)
  names(x)=c_x_n;
  print(names(x))
  return(x)
  })


save_file("all_significance_simulation")

