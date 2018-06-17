# This script plot cluster member sizes versus the blindBLAST accuracies 

current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

# get loop member size
 member_sizes= unlist(lapply(data_by_loop_type_list_unduplicated,function(x){dim(x[[1]])[1]}))
 cluster_number=unlist(lapply(data_by_loop_type_list_unduplicated,function(x){
   dim(x[[1]])[1]
   length(unique(x[[1]]$cluster_type))
   }))
 cluster_number=as.data.frame(cluster_number)
 cluster_number$loops=rownames(cluster_number)
 
 balance_ratio=unlist(lapply(data_by_loop_type_list_unduplicated,function(x){
   dim(x[[1]])[1]
   aa=x[[1]]$cluster_type
   a_frame= as.data.frame(table(aa))
   a_frame=a_frame[order(a_frame$Freq,decreasing = TRUE),]
   a_frame[1,"Freq"]/a_frame[2,"Freq"]
 }))
 balance_ratio=as.data.frame(balance_ratio)
 balance_ratio$loops=rownames(balance_ratio)
 
 member_sizes=as.data.frame(member_sizes)
 member_sizes$loops=rownames(member_sizes)
# get accuracies
 blindBLAST_acc=combined_blindBLAST_random[combined_blindBLAST_random$method=="blindBLAST",]
 blindBLAST_acc$loops=sapply(strsplit(rownames(blindBLAST_acc),"\\."),"[[",2)
 merged_acc=merge(member_sizes,blindBLAST_acc)
 acc_size_number= merge(merged_acc,cluster_number)
 acc_size_number_ratio= merge(acc_size_number,balance_ratio)
 acc_size_number_ratio$labels=paste( split_vector_and_replace(acc_size_number_ratio$loops,"_",1,2,"-"), acc_size_number_ratio$cluster_number,sep="   ")

 ggplot(acc_size_number_ratio,
                         aes(
                           x = member_sizes,
                           y = mean
                         ))+geom_point()+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                                                        position=position_dodge(.9)) +
     geom_label_repel(aes(label=labels),hjust=0, vjust=0)+
     theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+xlab("")+
     theme(text = element_text(size=12), axis.text.x = element_text(angle=0, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 

   
#   figure_r_list=list()
#   grids[[1]]=figure_error_count
#   grids[[2]]=fig
   
   

#ggsave("./Plots/loop_member_size_and_other_info_vs_accuracy.pdf", plot =figure_r,
#       width = 7, height = 7, units = c("in"))



# 
acc_size_number_ratio$lean_label=split_vector_and_replace(acc_size_number_ratio$loops,"_",1,2,"-")
melted_data=melt(acc_size_number_ratio,id.vars=c("loops","mean","loop","length","method","min","max","labels","sd","lean_label"))
split_info=split(melted_data,melted_data$variable)
figure_lists=list()
for(x in 1:length(split_info)){
  x_name= names(split_info)[x]
  figure_lists[[x]]=ggplot(split_info[[x]],
       aes(
         x = value,
         y = mean
       ))+geom_point()+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                                      position=position_dodge(.9)) + xlab(x_name)+
  #geom_text_repel(aes(label=lean_label),hjust=0, vjust=0)+
  #facet_wrap(~variable,scales="free",nrow=2,ncol=2)+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  theme(text = element_text(size=12), axis.text.x = element_text(angle=0, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 

}

melted_data_only_sd=melted_data[melted_data$variable=="member_sizes",]
figure_lists[[4]]=ggplot(melted_data_only_sd,
                         aes(
                           x = value,
                           y = sd
                         ))+geom_point() + xlab("member size")+
  #geom_text_repel(aes(label=lean_label),hjust=0, vjust=0)+
  #facet_wrap(~variable,scales="free",nrow=2,ncol=2)+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  theme(text = element_text(size=12), axis.text.x = element_text(angle=0, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 


p2=grid.arrange(arrangeGrob(figure_lists[[1]],figure_lists[[2]],
                            figure_lists[[3]],figure_lists[[4]], 
                            ncol=2, nrow=2, heights=c(2,2),widths=c(2,2)   ))

save_figure_specific_size(p2,"member_size_cluster_number_balance_ratio_vs_accuracy_splitted.pdf",7,7)
print(c("Plotted as ./Plots/member_size_cluster_number_balance_ratio_vs_accuracy_splitted.pdf"))
