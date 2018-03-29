


rbind_diff_between_blast_gbm=do.call(rbind,diff_between_blast_gbm)

loops_with_cis=rbind_diff_between_blast_gbm[grepl("cis",rbind_diff_between_blast_gbm$query_cluster)|grepl("cis",rbind_diff_between_blast_gbm$template_cluster),]
loops_with_cis=loops_with_cis[loops_with_cis$improvement_by_gbm>=1|loops_with_cis$improvement_by_gbm<=-1,]

loops_with_cis$loop_length=lapply(strsplit(as.character(loops_with_cis$query_cluster),"-"),function(x){paste(x[1:2],collapse = "-")})
#loops_with_cis=loops_with_cis[order(loops_with_cis$improvement_by_gbm,decreasing = TRUE),]
query_clust_short=split_vector_and_replace(loops_with_cis$query_cluster,"-",3,4,"-")
template_clust_short=split_vector_and_replace(loops_with_cis$template_cluster,"-",3,4,"-")
loops_with_cis$misclassification=paste(query_clust_short,template_clust_short,sep="->")
melted_loops_with_cis=as.data.frame(melt(loops_with_cis[,c("blindBLAST_error_count", "gbm_error_count","loop_length", "misclassification")],id.vars=c("loop_length", "misclassification")))
melted_loops_with_cis$misclassification=gsub("\\.",",",melted_loops_with_cis$misclassification)
melted_loops_with_cis=as.data.frame(melted_loops_with_cis)
for(x in 1:4){melted_loops_with_cis[,x]=unlist(melted_loops_with_cis[,x]);}
melted_loops_with_cis$loop_length=factor(melted_loops_with_cis$loop_length,levels=order_factor_by_two_component(melted_loops_with_cis$loop_length,"-",1,2))
  
figure=plot_figure(melted_loops_with_cis,"misclassification","value","variable","loop_length",c("none"))



non_cis_rbind_diff_between_blast_gbm=rbind_diff_between_blast_gbm[!grepl("cis",rbind_diff_between_blast_gbm$query_cluster)&!grepl("cis",rbind_diff_between_blast_gbm$template_cluster),]
non_cis_rbind_diff_between_blast_gbm=non_cis_rbind_diff_between_blast_gbm[non_cis_rbind_diff_between_blast_gbm$improvement_by_gbm>3|non_cis_rbind_diff_between_blast_gbm$improvement_by_gbm<=-3,]
non_cis_query_clust_short=split_vector_and_replace(non_cis_rbind_diff_between_blast_gbm$query_cluster,"-",3,4,"-")
non_cis_rbind_diff_between_blast_gbm$loop_length=split_vector_and_replace(non_cis_rbind_diff_between_blast_gbm$template_cluster,"-",1,2,"-")
non_cis_template_clust_short=split_vector_and_replace(non_cis_rbind_diff_between_blast_gbm$template_cluster,"-",3,4,"-")
non_cis_rbind_diff_between_blast_gbm$misclassification=paste(non_cis_query_clust_short,non_cis_template_clust_short,sep="->")
non_cis_rbind_diff_between_blast_gbm_melt=as.data.frame(melt(non_cis_rbind_diff_between_blast_gbm[,c("blindBLAST_error_count", "gbm_error_count","loop_length", "misclassification")],id.vars=c("loop_length", "misclassification")))
non_cis_rbind_diff_between_blast_gbm_melt$loop_length=factor(non_cis_rbind_diff_between_blast_gbm_melt$loop_length,levels=order_factor_by_two_component(non_cis_rbind_diff_between_blast_gbm_melt$loop_length,"-",1,2))

figure_non_cis=plot_figure(non_cis_rbind_diff_between_blast_gbm_melt,"misclassification","value","variable","loop_length",c("none"))
figure=figure+scale_y_continuous(position = "right")
figure_non_cis=figure_non_cis+scale_y_continuous(position = "right")


grids=list()
grids[[1]]=figure
grids[[2]]=figure_non_cis
p=grid.arrange(arrangeGrob(grids[[1]],grids[[2]], ncol=1, nrow=2, heights=c(1,1) ,bottom=textGrob("wrong prediction type", gp=gpar(fontsize=14))  ))

plot_dir="./proline_classifier/Plots/"
save_figure_specific_size(p,"gbm_improved_from_blindBLAST.pdf",7,7)
  
