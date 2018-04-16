average_accuracy_from_simulation=list()
sd_from_simulation=list()
for(loop in names(data_by_loop_type_list_unduplicated)){
  clusters=unique(data_by_loop_type_list_unduplicated[[loop]][[1]]$cluster_type)
  clusters=gsub("\\*","none",clusters)
  all_combn=expand.grid(clusters,clusters)
  all_combn=as.data.frame(t(all_combn))
  all_names=names(all_significance_simulation[[loop]])
  all_s=all_significance_simulation[[loop]]
  all_num=sum(unlist(lapply(all_significance_simulation[[loop]],function(x){x[1]})))
  similar=all_combn[,unlist(lapply(all_combn,function(x){x[1]==x[2]}))]
  correct_num=  lapply(similar,function(x){
    x=paste(unlist(x),collapse="")
    #print(possible_ids)
    if( any(x %in% all_names)){
      correct=all_significance_simulation[[loop]][[x]]
      return(correct)
     }
    
  })
  starter=correct_num[[1]]
  for(ind in 2:length(correct_num)){
    starter=starter+correct_num[[ind]]
  }
  
  accuracy_sd = sd(correct_num_sum/all_num)
  
  all_errors=lapply(all_combn,function(x){
x=paste(unlist(x),collapse="")
    #print(possible_ids)
    if( any(x %in% all_names)){
      the_mean=mean(all_significance_simulation[[loop]][[x]])
      return(the_mean)
    }else{
      return(NA)
    }
   
  })
  simulated_error_average=as.data.frame(cbind(t(all_combn),as.numeric(unlist(all_errors))))
  simulated_error_average$V3=as.numeric(as.character(simulated_error_average$V3))
  accuracy=sum(simulated_error_average[simulated_error_average$Var1==simulated_error_average$Var2,3])/sum(simulated_error_average[,3])
  average_accuracy_from_simulation[[loop]]=accuracy
  sd_from_simulation[[loop]]=accuracy_sd
}
random_acc_sd = as.data.frame(cbind(accuracy=average_accuracy_from_simulation,sd=sd_from_simulation))
random_acc_sd = random_acc_sd[!rownames(random_acc_sd)%in%c("H2_12","L1_17"),]
save_file("random_acc_sd")
