# getting result of blindBLAST in 3 repeat 10 fold cv
# Getting rmsds of blindBLAST template vs query CDR
# need to use generate_training_grids.R for passing parameters

current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")

}

each_method = "blindblast_just_get_alignment"
  cluster_dis = "north"
subsitution_matrix_name = "wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
subsitution_matrix = "PAM30"


# do blindBLAST template searching for each CDR loop 
overall_accuracy = list()
ten_foldcv_blindblastlist = list()
the_loop_names=names(data_by_loop_type_list_unduplicated)
left_loops=the_loop_names[!the_loop_names%in% names(ten_foldcv_blindblastlist) ]
for (loop_type in left_loops) {
  returned_list=getting_similar_sequences_similarity_matrix_rmsd_matrix(loop_type)
  if(!is.logical(returned_list)){
  g(sequences,similarity_matrix,rmsd_matrix)  %=%  returned_list
  }else{print("jump to next");next;}
  
  each_loop_length_data_feature_string = data_by_loop_type_list_unduplicated[[loop_type]][[2]]
  features=data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  sequences$cluster_type = as.character(sequences$cluster_type)
  all_cases =  sequences[, c(features, "cluster_type", "identifier", "id")]
  #make folds and separates folds in and folds out
  r <- 3 # number of repeats
  k <- 10 # number of folds
  g(folds.list,folds.list.out)  %=%generate_folds_foldsout(sequences,3,10)
  
  all_result_list = mclapply(1:length(folds.list), get_accuracy_per_fold, mc.cores =    4)  # run blindBLAST for each folds out 
  acc_result=calculate_accuracy_mean_std(all_result_list)
  overall_accuracy[[loop_type]] = c(acc_result[[1]], acc_result[[2]])
  ten_foldcv_blindblastlist[[loop_type]] = all_result_list

}# end of iterating all loop types



save_file("ten_foldcv_blindblastlist")




results_rbind_allloops=lapply(ten_foldcv_blindblastlist,function(each_r){
  result_l=chunk2(1:length(each_r),3)
#  for(each_r in ten_foldcv_blindblastlist){
  results=lapply(result_l,function(x){
  fr=do.call(rbind,each_r[x])
  fr$obs=split_vector_and_replace(fr[,1],"\\.",1,1,"")
  fr$pred=split_vector_and_replace(fr[,2],"\\.",1,1,"")
  #dim(fr[fr$obs==fr$pred,])[1]/dim(fr)[1]
  repeat_number=x[1]%/%10+1
  print(dim(fr)[1])
  print(each_r)
  fr$repeats=rep(repeat_number,dim(fr)[1])
  return(fr)
  })
  results_rbind=do.call(rbind,results)
  return(results_rbind)
  #}
  
  
})


# calculate accuracies for blindBLAST for each repeat
# calculate the std for accuracies for each repeat

all_r=do.call(rbind,results_rbind_allloops)
split_by_r=split(all_r,all_r$repeats)
accuracies_by_repeat_b=lapply(split_by_r,function(x){
  accuracy=dim(x[x$pred==x$obs,])[1]/dim(x)[1]
})
sds=sd(unlist(accuracies_by_repeat_b))

# generate blindBLAST summary 
#blind_blast_cv_result_summary = as.data.frame(do.call(rbind, overall_accuracy))
blind_blast_cv_result_summary=as.data.frame(cbind(blindBLAST_mean_accu,blindBLAST_accu_std ))


#rownames(blind_blast_cv_result_summary) = names(overall_accuracy)
colnames(blind_blast_cv_result_summary) = c("mean", "sd")
blind_blast_cv_result_summary$loop_type = rownames(blind_blast_cv_result_summary)
save_file("blind_blast_cv_result_summary")



# summarize result for blindBLAST
blind_blast_cv_result_summary$loop=split_vector_and_replace(blind_blast_cv_result_summary$loop_type,"_",1,1,"-")
blind_blast_cv_result_summary$length=split_vector_and_replace(blind_blast_cv_result_summary$loop_type,"_",2,2,"-")
blind_blast_cv_result_summary_reordered=reorder_factor(blind_blast_cv_result_summary,"loop","length")
#blindBLAST_accu_std

# parse result of random assignment
random_acc_sd$loop=split_vector_and_replace(rownames(random_acc_sd),"_",1,1,"")
random_acc_sd$length=split_vector_and_replace(rownames(random_acc_sd),"_",2,2,"")
random_acc_sd_reordered=reorder_factor(random_acc_sd,"loop","length")
blind_blast_cv_result_summary_reordered_selected=blind_blast_cv_result_summary_reordered[,c( "mean" ,    "sd" , "loop", "length")]
blind_blast_cv_result_summary_reordered_selected$method=rep("blindBLAST",dim(blind_blast_cv_result_summary_reordered_selected)[1])
colnames(random_acc_sd_reordered)[1:2]=c("mean","sd")
random_acc_sd_reordered$method=rep("random",dim(random_acc_sd_reordered)[1])

# combine blindBLAST result and that of random assignment 
combined_blindBLAST_random=as.data.frame(rbind(blind_blast_cv_result_summary_reordered_selected,random_acc_sd_reordered))
combined_blindBLAST_random$sd=as.numeric(combined_blindBLAST_random$sd)
combined_blindBLAST_random$mean=as.numeric(combined_blindBLAST_random$mean)
combined_blindBLAST_random$min=combined_blindBLAST_random$mean-combined_blindBLAST_random$sd/2
combined_blindBLAST_random$max=combined_blindBLAST_random$mean+combined_blindBLAST_random$sd/2
save_file("combined_blindBLAST_random")
# plot accuracy for both random assignment and blindBLAST
p1=ggplot(combined_blindBLAST_random,
       aes(
         x = length,
         y = mean,
         ymin = mean - sd / 2,
         ymax = mean + sd / 2,fill=method
       ))  +geom_bar(
         stat = "identity" ,position=position_dodge(width = 0.90))+
  geom_errorbar(position = position_dodge(0.9)) + facet_grid(~loop,scales="free",space="free") + theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position=c(.8, .8),axis.title.x=element_blank())
options(digits = 2)

save_file("p1")

