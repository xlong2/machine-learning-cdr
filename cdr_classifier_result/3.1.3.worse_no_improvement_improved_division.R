# categorize the misclassification based on the significance score
current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

all_error_sig_info=do.call(rbind,lapply(names(all_loop_sig_frame_list),function(y){
  x=all_loop_sig_frame_list[[y]];x$loop=rep(y,dim(x)[1])
  # if the first column has the cluster one and the two clusters are not the ssame
  H1_related=x[x[,"Var1"]!=x[,"Var2"],]
  # most prevalent cluster 
  each_l=paste(unlist(strsplit(y,"_")[[1]]),collapse="-")
  if(each_l=="L3-9"){
    most_p="L3-9-cis7-1$"
  }else{
    most_p=paste(c(each_l,"-1$"),collapse="")
  }
  H1_cluster_one_recovery_fail=H1_related[grepl(most_p,H1_related$Var1),]
  H1_cluster_one_precision_fail=H1_related[grepl(most_p,H1_related$Var2),]
  H1_cluster_one_recovery_fail$error_type=rep("1_reco",dim(H1_cluster_one_recovery_fail)[1])
  H1_cluster_one_precision_fail$error_type=rep("1_prec",dim(H1_cluster_one_precision_fail)[1])
  other_types=H1_related[!(grepl(most_p,H1_related$Var1) |grepl(most_p,H1_related$Var2)),]
  other_types$error_type=rep("non_1",dim(other_types)[1])
  all_annotated=rbind(rbind(other_types,H1_cluster_one_precision_fail),H1_cluster_one_recovery_fail)
  return(all_annotated)
  # if the 
}))
all_error_sig_info=all_error_sig_info[all_error_sig_info$Var1!=all_error_sig_info$Var2,]
all_error_sig_info=all_error_sig_info[complete.cases(all_error_sig_info),]
all_error_sig_info$sig_id=rep("insignificant",dim(all_error_sig_info)[1])
all_error_sig_info[all_error_sig_info$significance<=0.025,"sig_id"]=rep("smaller",length(all_error_sig_info[all_error_sig_info$significance<=0.025,"sig_id"]))
all_error_sig_info[all_error_sig_info$significance>=0.975,"sig_id"]=rep("larger",length(all_error_sig_info[all_error_sig_info$significance>=0.975,"sig_id"]))





effect_size_large=all_error_sig_info[all_error_sig_info$effect_size>=2 &all_error_sig_info$sig_id=="larger",]
effect_size_large=effect_size_large[order(effect_size_large$error_count,decreasing=TRUE),]
effect_size_large_sig=effect_size_large[effect_size_large$error_count>3,]


no_improvement_worse=all_error_sig_info[all_error_sig_info$error_count>3 & all_error_sig_info$sig_id=="larger",]
no_improvement_worse=no_improvement_worse[,c("error_count","mean_simu_error","sd","Var1","Var2","significance")]
names(no_improvement_worse)=c("error_count","mean_simu_error_count","sd","query_cluster","template_cluster","significance")
names(no_improvement_worse)=c("error_count","mean_simu_error_count","sd","query_cluster","template_cluster","significance")
file.remove("./Data_processed/no_improvement_worse.csv")
no_improvement_worse=no_improvement_worse[,c("query_cluster","template_cluster","error_count", "mean_simu_error_count",   "sd","significance")]
write.csv(no_improvement_worse,file="./Data_processed/no_improvement_worse.csv")


not_much_improvement=all_error_sig_info[all_error_sig_info$error_count>3 &(all_error_sig_info$sig_id=="insignificant"),]
not_much_improvement= not_much_improvement[order(not_much_improvement$error_count),]
not_much_improvement_output=not_much_improvement[,c("error_count","mean_simu_error","sd","Var1","Var2","significance")]
names(not_much_improvement_output)=c("error_count","mean_simu_error_count","sd","query_cluster","template_cluster","significance")
not_much_improvement_output=not_much_improvement_output[,c("query_cluster","template_cluster","error_count", "mean_simu_error_count",   "sd","significance")]

write.csv(not_much_improvement_output,file="./Data_processed/not_much_improvement_output.csv")



improvement_but_still_many_errors=all_error_sig_info[all_error_sig_info$error_count>3 &(all_error_sig_info$sig_id=="smaller"),]
improvement_but_still_many_errors=improvement_but_still_many_errors[,c("error_count","mean_simu_error","sd","Var1","Var2","significance")]
names(improvement_but_still_many_errors)=c("error_count","mean_simu_error_count","sd","query_cluster","template_cluster","significance")
improvement_but_still_many_errors=improvement_but_still_many_errors[,c("query_cluster","template_cluster","error_count", "mean_simu_error_count",   "sd","significance")]

write.csv(improvement_but_still_many_errors,file="./Data_processed/improvement_but_still_many_errors.csv")

improvement_and_not_very_much_error=all_error_sig_info[all_error_sig_info$error_count<=3 &(all_error_sig_info$sig_id=="smaller"),]
improvement_and_not_very_much_error=improvement_and_not_very_much_error[,c("error_count","mean_simu_error","sd","Var1","Var2","significance")]
names(improvement_and_not_very_much_error)=c("error_count","mean_simu_error_count","sd","query_cluster","template_cluster","significance")
improvement_and_not_very_much_error=improvement_and_not_very_much_error[,c("query_cluster","template_cluster","error_count", "mean_simu_error_count",   "sd","significance")]
write.csv(improvement_and_not_very_much_error,file="./Data_processed/improvement_and_not_very_much_error.csv")

improvement_and_not_very_much_error$sig_stat=rep("better than random and little errors",dim(improvement_and_not_very_much_error)[1])
improvement_but_still_many_errors$sig_stat=rep("better than random but still many errors",dim(improvement_but_still_many_errors)[1])
not_much_improvement_output$sig_stat=rep("not different than random",dim(not_much_improvement_output)[1])
no_improvement_worse$sig_stat=rep("worse than random",dim(no_improvement_worse)[1])
all_frames_I_care=rbind(rbind(rbind(improvement_and_not_very_much_error,improvement_but_still_many_errors),not_much_improvement_output),no_improvement_worse)
all_frames_I_care$loop_type=unlist(lapply(strsplit(as.character(all_frames_I_care$query_cluster),"-"),function(x){paste(x[1:2],collapse="-")}))
all_frames_I_care$sig_stat=factor(all_frames_I_care$sig_stat,levels=c("worse than random","not different than random","better than random but still many errors","better than random and little errors"))
#generate colors
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector=col_vector[c(1:3,5:(n))]
col_vector=c(col_vector,"#bdbdbd")
pie(rep(1,n), col=col_vector)
col_vec=col_vector[1:length(unique(all_frames_I_care$sig_stat))]
saveRDS(all_frames_I_care,file="./Data_processed/all_frames_I_care.rds")

abc=ggplot(all_frames_I_care, aes(mean_simu_error_count, error_count,colour = factor(sig_stat)))+
  geom_point( position=position_dodge(width = 0.90),size = 2.5)+
  scale_colour_manual(breaks = all_frames_I_care$sig_stat, name="category",
                      values =col_vec)+
 xlab(" error count from random assignment")+ylab("blindBLAST error count")+   theme_classic()+theme_bw()+
  theme(strip.text.x = element_text(size = 12),legend.position = c(0.75, 0.9),
        axis.text.x = element_text(size = 11, colour = "black"),plot.title =element_text(hjust = 0.5,size=18,face="bold"),
        plot.margin=unit(c(2,2,2,2),"mm"),
        axis.text.y = element_text(size=11))
save_figure_specific_size(abc,"misclassifications_categorized_by_significance_1.pdf",7,7)


print(c("plotted ./Plots/misclassifications_categorized_by_significance_1.pdf"))


