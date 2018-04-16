# getting result of blindBLAST in 3 repeat 10 fold cv
# getting rmsds of the wrong cases and rmsds when we enforce correct cluster template selection.



each_method = "blindblast_just_get_alignment"
  cluster_dis = "north"
subsitution_matrix_name = "wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
subsitution_matrix = "PAM30"


#do LOOCV
#loop through all sequences
#construct a ml model using all sequences except for the single query seq
# predict the class of the query seq
# construct a blast database using the predicted cluster but without the seq
# run blast search and at the same time calculate the prediction rmsd
# store all such prediction and rmsd into two list
# retrieve them after finish looping
# store the rds into a directory

# use gbm to do LOOCV by iterating all seqs
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
  
  all_result_list = mclapply(1:length(folds.list), get_accuracy_per_fold, mc.cores =
                               4)
  acc_result=calculate_accuracy_mean_std(all_result_list)
  overall_accuracy[[loop_type]] = c(acc_result[[1]], acc_result[[2]])
  ten_foldcv_blindblastlist[[loop_type]] = all_result_list

}# end of iterating all loop types

#for each loop find out which ones are not the real cluster and get the query and predicted pdbs








save_file("ten_foldcv_blindblastlist")

blind_blast_cv_result_summary = as.data.frame(do.call(rbind, overall_accuracy))
rownames(blind_blast_cv_result_summary) = names(overall_accuracy)
colnames(blind_blast_cv_result_summary) = c("mean", "sd")
blind_blast_cv_result_summary$loop_type = rownames(blind_blast_cv_result_summary)
save_file("blind_blast_cv_result_summary")




# get the blindblast 10 fold cv accuracies
blind_blast_cv_result_summary$loop=split_vector_and_replace(blind_blast_cv_result_summary$loop_type,"_",1,1,"-")
blind_blast_cv_result_summary$length=split_vector_and_replace(blind_blast_cv_result_summary$loop_type,"_",2,2,"-")
blind_blast_cv_result_summary_reordered=reorder_factor(blind_blast_cv_result_summary,"loop","length")


random_acc_sd$loop=split_vector_and_replace(rownames(random_acc_sd),"_",1,1,"")
random_acc_sd$length=split_vector_and_replace(rownames(random_acc_sd),"_",2,2,"")
random_acc_sd_reordered=reorder_factor(random_acc_sd,"loop","length")
blind_blast_cv_result_summary_reordered_selected=blind_blast_cv_result_summary_reordered[,c( "mean" ,    "sd" , "loop", "length")]
blind_blast_cv_result_summary_reordered_selected$method=rep("blindBLAST",dim(blind_blast_cv_result_summary_reordered_selected)[1])
colnames(random_acc_sd_reordered)[1:2]=c("mean","sd")
random_acc_sd_reordered$method=rep("random",dim(random_acc_sd_reordered)[1])
combined_blindBLAST_random=as.data.frame(rbind(blind_blast_cv_result_summary_reordered_selected,random_acc_sd_reordered))
combined_blindBLAST_random$sd=as.numeric(combined_blindBLAST_random$sd)
combined_blindBLAST_random$mean=as.numeric(combined_blindBLAST_random$mean)
combined_blindBLAST_random$min=combined_blindBLAST_random$mean-combined_blindBLAST_random$sd/2
combined_blindBLAST_random$max=combined_blindBLAST_random$mean+combined_blindBLAST_random$sd/2

p1=ggplot(combined_blindBLAST_random,
       aes(
         x = length,
         y = mean,
         ymin = mean - sd / 2,
         ymax = mean + sd / 2,colour=method
       )) + geom_errorbar(stat = "identity") + facet_grid(~loop,scales="free",space="free")+geom_point() + theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position=c(.8, .8),axis.title.x=element_blank())
options(digits = 2)


grids=list()
grids[[1]]=p1
grids[[2]]=the_plot
p=grid.arrange(arrangeGrob(grids[[1]],grids[[2]], ncol=1, nrow=2, heights=c(3,2) ,bottom=textGrob("wrong prediction type", gp=gpar(fontsize=14))  ))

save_figure_specific_size(p,"blindBLAST_corrected_with_accuracy_rmsdplot.pdf",7,7)


# selecting out the incorrect ones and enforcing the search within just the right cluster
# get what happens if I enforce the correct cluster choice.


