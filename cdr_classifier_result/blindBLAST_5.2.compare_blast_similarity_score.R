# compare the blast score 

each_method = "blindblast_just_get_alignment"
cluster_dis = "north"
subsitution_matrix_name = "wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
subsitution_matrix = "PAM30"




# use gbm to do LOOCV by iterating all seqs
the_loop_names=names(data_by_loop_type_list_unduplicated)

# selecting out small clusters; 
mean_bitscores_and_bitscore_std_comparison=list()
for (loop_type in the_loop_names) {
  returned_list=getting_similar_sequences_similarity_matrix_rmsd_matrix(loop_type)
  if(!is.logical(returned_list)){
    g(sequences,similarity_matrix,rmsd_matrix)  %=%  returned_list
  }else{print("jump to next");next;}
  
  each_loop_length_data_feature_string = data_by_loop_type_list_unduplicated[[loop_type]][[2]]
  features=data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  sequences$cluster_type = as.character(sequences$cluster_type)
  all_cases =  sequences[, c(features, "cluster_type", "identifier", "id")]
  #make folds and separates folds in and folds out
  r <- 1 # number of repeats
  k <- 10 # number of folds
  small_clusters=as.data.frame(table(sequences$cluster_type))
  sparse_clusters = as.character(small_clusters[small_clusters$Freq<=20,"Var1"])
  g(folds.list,folds.list.out)  %=%generate_folds_foldsout(sequences,r,k)
  
  sparse_cluster_ids=which(sequences$cluster_type%in% sparse_clusters)
  if(length(sparse_cluster_ids)==0){next}
  folds.list.out=lapply(folds.list.out,function(x){
    x[x%in%sparse_cluster_ids]
  })
  all_result_list = mclapply(1:length(folds.list), get_bit_scores, mc.cores =
                               1)
  # get the mean difference 
  all_result_list=do.call(rbind,all_result_list)
  
  all_result_list$mean=apply(all_result_list,1,function(x){print(x); mean(as.numeric(as.character(x[2:6])),na.rm=TRUE) })
  all_result_list[,1]=split_vector_and_replace(all_result_list$X1,"\\.",1,1,"")
  mean_by_cluster=split(all_result_list,all_result_list[,1])
  mean_by_cluster_list=lapply(mean_by_cluster,function(x){mean(x$mean,na.rm=TRUE)})
  #calculate mean accuracy
  
  # block out anything 
  # randomly choose 5 matches and get the std of the bit scores... 
  random_bit_scores = mclapply(1:length(folds.list), get_random_bit_scores, mc.cores =
                               4)
  random_bit_scores=do.call(rbind,random_bit_scores)
  random_bit_scores$std=apply(random_bit_scores,1,function(x){sd(x[2:6])})
  random_bit_scores[,1]=split_vector_and_replace(random_bit_scores$X1,"\\.",1,1,"")
  std_by_cluster=split(random_bit_scores,random_bit_scores[,1])
  mean_std_by_cluster=lapply(std_by_cluster,function(x){mean(x$std,na.rm=TRUE)})
  
  # split by the cluster name 
  compare_values=cbind(a=unlist(mean_by_cluster_list),b=unlist(mean_std_by_cluster))
  
  mean_bitscores_and_bitscore_std_comparison[[loop_type]]=compare_values
}#


combined_v=do.call(rbind, mean_bitscores_and_bitscore_std_comparison)
combined_v=as.data.frame(combined_v[complete.cases(combined_v), ])

point_d=as.data.frame(rbind(c(0,0),c(5,5)))
combined_v=combined_v[!grepl("none",rownames(combined_v)),]
figure_sparse=ggplot(combined_v,
       aes(
         x = a,
         y = b
       ))+geom_point() +
  geom_text_repel(aes(label=rownames(combined_v)),hjust=0, vjust=0)+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+xlab("")+
  theme(text = element_text(size=12), axis.text.x = element_text(angle=0, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 
figure_sparse=figure_sparse+ coord_fixed()+ylab("Bit Score Std")+xlab("Bit Score Mean")+geom_line(data=point_d,aes(V1, V2, group = 1))

ggsave("./Plots/bitscores_mean_std_sparse_cluster.pdf", plot =figure_sparse,
       width = 7, height = 7, units = c("in"))
