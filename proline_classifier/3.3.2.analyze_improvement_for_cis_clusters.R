# This script compare between blindBLAST and GBM in cis related misclassifications. 

current_d=getwd()
if(grepl("proline_classifier",current_d)){
  source("0.load_function_and_data.R")
  
}

rbind_diff_between_blast_gbm=do.call(rbind,diff_between_blast_gbm)

loops_with_cis=rbind_diff_between_blast_gbm[grepl("cis",rbind_diff_between_blast_gbm$query_cluster)|grepl("cis",rbind_diff_between_blast_gbm$template_cluster),]
loops_with_cis=loops_with_cis[loops_with_cis$improvement_by_gbm>=1|loops_with_cis$improvement_by_gbm<=-1,]

loops_with_cis$loop_length=lapply(strsplit(as.character(loops_with_cis$query_cluster),"-"),function(x){paste(x[1:2],collapse = "-")})
#loops_with_cis=loops_with_cis[order(loops_with_cis$improvement_by_gbm,decreasing = TRUE),]
query_clust_short=split_vector_and_replace(loops_with_cis$query_cluster,"-",3,4,"-")
template_clust_short=split_vector_and_replace(loops_with_cis$template_cluster,"-",3,4,"-")
loops_with_cis$misclassification=paste(query_clust_short,template_clust_short,sep="->")
melted_loops_with_cis=as.data.frame(melt(loops_with_cis[,c("blindBLAST_error_count", "gbm_error_count","loop_length", "misclassification","blindBLAST_sd","gbm_sd")],id.vars=c("loop_length", "misclassification","blindBLAST_sd","gbm_sd")))
melted_loops_with_cis$misclassification=gsub("\\.",",",melted_loops_with_cis$misclassification)
melted_loops_with_cis=as.data.frame(melted_loops_with_cis)
for(x in 1:4){melted_loops_with_cis[,x]=unlist(melted_loops_with_cis[,x]);}
melted_loops_with_cis[melted_loops_with_cis$variable=="blindBLAST_error_count","sd"]=melted_loops_with_cis[melted_loops_with_cis$variable=="blindBLAST_error_count","blindBLAST_sd"]
melted_loops_with_cis[melted_loops_with_cis$variable=="gbm_error_count","sd"]=melted_loops_with_cis[melted_loops_with_cis$variable=="gbm_error_count","gbm_sd"]
melted_loops_with_cis[,c("blindBLAST_sd",    "gbm_sd")]<-NULL
melted_loops_with_cis$variable=as.character(melted_loops_with_cis$variable)
melted_loops_with_cis$loop_length=factor(as.character(melted_loops_with_cis$loop_length),levels=order_factor_by_two_component(as.character(melted_loops_with_cis$loop_length),"-",1,2))
melted_loops_with_cis$low=melted_loops_with_cis$value-melted_loops_with_cis$sd/2
melted_loops_with_cis$high=melted_loops_with_cis$value+melted_loops_with_cis$sd/2

figure=plot_figure(melted_loops_with_cis,"misclassification","value","variable","loop_length",c("none"))



figure=ggplot(melted_loops_with_cis,aes(x=misclassification,y=value,ymin = low , ymax = high,fill=variable))+geom_bar(
  stat = "identity" ,position=position_dodge(width = 0.90))+
  geom_errorbar(position = position_dodge(0.9))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  xlab("")+facet_grid(~loop_length,scales="free",space="free")+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 



non_cis_rbind_diff_between_blast_gbm=rbind_diff_between_blast_gbm[!grepl("cis",rbind_diff_between_blast_gbm$query_cluster)&!grepl("cis",rbind_diff_between_blast_gbm$template_cluster),]
non_cis_rbind_diff_between_blast_gbm=non_cis_rbind_diff_between_blast_gbm[non_cis_rbind_diff_between_blast_gbm$improvement_by_gbm>3|non_cis_rbind_diff_between_blast_gbm$improvement_by_gbm<=-3,]
non_cis_query_clust_short=split_vector_and_replace(non_cis_rbind_diff_between_blast_gbm$query_cluster,"-",3,4,"-")
non_cis_rbind_diff_between_blast_gbm$loop_length=split_vector_and_replace(non_cis_rbind_diff_between_blast_gbm$template_cluster,"-",1,2,"-")
non_cis_template_clust_short=split_vector_and_replace(non_cis_rbind_diff_between_blast_gbm$template_cluster,"-",3,4,"-")
non_cis_rbind_diff_between_blast_gbm$misclassification=paste(non_cis_query_clust_short,non_cis_template_clust_short,sep="->")
melted_non_cis_rbind_diff_between_blast_gbm=as.data.frame(melt(non_cis_rbind_diff_between_blast_gbm[,c("blindBLAST_error_count", "gbm_error_count","loop_length", "misclassification","blindBLAST_sd","gbm_sd")],id.vars=c("loop_length", "misclassification","blindBLAST_sd","gbm_sd")))
melted_non_cis_rbind_diff_between_blast_gbm$misclassification=gsub("\\.",",",melted_non_cis_rbind_diff_between_blast_gbm$misclassification)
melted_non_cis_rbind_diff_between_blast_gbm=as.data.frame(melted_non_cis_rbind_diff_between_blast_gbm)
for(x in 1:4){melted_non_cis_rbind_diff_between_blast_gbm[,x]=unlist(melted_non_cis_rbind_diff_between_blast_gbm[,x]);}
melted_non_cis_rbind_diff_between_blast_gbm[melted_non_cis_rbind_diff_between_blast_gbm$variable=="blindBLAST_error_count","sd"]=melted_non_cis_rbind_diff_between_blast_gbm[melted_non_cis_rbind_diff_between_blast_gbm$variable=="blindBLAST_error_count","blindBLAST_sd"]
melted_non_cis_rbind_diff_between_blast_gbm[melted_non_cis_rbind_diff_between_blast_gbm$variable=="gbm_error_count","sd"]=melted_non_cis_rbind_diff_between_blast_gbm[melted_non_cis_rbind_diff_between_blast_gbm$variable=="gbm_error_count","gbm_sd"]
melted_non_cis_rbind_diff_between_blast_gbm[,c("blindBLAST_sd",    "gbm_sd")]<-NULL
melted_non_cis_rbind_diff_between_blast_gbm$variable=as.character(melted_non_cis_rbind_diff_between_blast_gbm$variable)
melted_non_cis_rbind_diff_between_blast_gbm$loop_length=factor(as.character(melted_non_cis_rbind_diff_between_blast_gbm$loop_length),levels=order_factor_by_two_component(as.character(melted_non_cis_rbind_diff_between_blast_gbm$loop_length),"-",1,2))
melted_non_cis_rbind_diff_between_blast_gbm$low=melted_non_cis_rbind_diff_between_blast_gbm$value-melted_non_cis_rbind_diff_between_blast_gbm$sd/2
melted_non_cis_rbind_diff_between_blast_gbm$high=melted_non_cis_rbind_diff_between_blast_gbm$value+melted_non_cis_rbind_diff_between_blast_gbm$sd/2



figure_non_cis=ggplot(melted_non_cis_rbind_diff_between_blast_gbm,aes(x=misclassification,y=value,ymin = low , ymax = high,fill=variable))+geom_bar(
  stat = "identity" ,position=position_dodge(width = 0.90))+
  geom_errorbar(position = position_dodge(0.9))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  xlab("")+facet_grid(~loop_length,scales="free",space="free")+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 

figure=figure+scale_y_continuous(position = "right")
figure_non_cis=figure_non_cis+scale_y_continuous(position = "right")










grids=list()
grids[[1]]=figure
grids[[2]]=figure_non_cis
p=grid.arrange(arrangeGrob(grids[[1]],grids[[2]], ncol=1, nrow=2, heights=c(1,1) ,bottom=textGrob("wrong prediction type", gp=gpar(fontsize=14))  ))

plot_dir="./Plots/"
print(c("figure plotted", "./Plots/gbm_improved_from_blindBLAST.pdf"))
save_figure_specific_size(p,"gbm_improved_from_blindBLAST.pdf",7,7)
  
