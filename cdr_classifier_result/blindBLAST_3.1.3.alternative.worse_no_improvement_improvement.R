# categorize the misclassification based on the significance score
current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}


for(each_loop in names(all_loop_sig_list_chi_test_new)){
  each_l_f = all_loop_sig_list_chi_test_new[[each_loop]]
  new_v = conf_tables_all_loops_blindBLAST_diff[[each_loop]]
  for(x in 1:dim(each_l_f)[1]){
    t_1 = each_l_f[x,"Var1"]; t_2= each_l_f[x,"Var2"]
    if(length(new_v[as.character(new_v$Var1)==as.character(t_1) &as.character(new_v$Var2)==as.character(t_2),"Freq"])!=0){
      each_l_f[x,"error_count"]=new_v[as.character(new_v$Var1)==as.character(t_1) &as.character(new_v$Var2)==as.character(t_2),"Freq"]
    }
  }
  all_loop_sig_list_chi_test_new[[each_loop]] =each_l_f
  
}

all_error_sig_info=do.call(rbind,lapply(names(all_loop_sig_list_chi_test_new),function(y){
  x=all_loop_sig_list_chi_test_new[[y]];x$loop=rep(y,dim(x)[1])
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
#all_error_sig_info=all_error_sig_info[complete.cases(all_error_sig_info),]
all_error_sig_info$sig_id=rep("insignificant",dim(all_error_sig_info)[1])
which_smaller = all_error_sig_info$error_count <all_error_sig_info$mean_simu_error & all_error_sig_info$significance<=0.025
which_larger = all_error_sig_info$error_count >all_error_sig_info$mean_simu_error & all_error_sig_info$significance>=0.975
all_error_sig_info[which_smaller,"sig_id"] = rep("smaller",length(all_error_sig_info[which_smaller,"sig_id"]))
all_error_sig_info[which_larger,"sig_id"] = rep("larger",length(all_error_sig_info[which_larger,"sig_id"]))



no_improvement_worse = all_error_sig_info[ all_error_sig_info$sig_id=="larger",]
no_improvement_worse = no_improvement_worse[order(no_improvement_worse$error_count,decreasing=TRUE),]
no_improvement_worse = no_improvement_worse[,c("error_count","mean_simu_error","expected","sd","Var1","Var2","significance","chi_sq_seq","sig_id")]
names(no_improvement_worse)=c("error_count","mean_simu_error","expected","sd","query_cluster","template_cluster","significance","chi_sq_seq","sig_id")

file.remove("./Data_processed/no_improvement_worse.csv")
no_improvement_worse=no_improvement_worse[,c("query_cluster","template_cluster","mean_simu_error","error_count", "expected","significance","sig_id")]
no_improvement_worse <- apply(no_improvement_worse,2,as.character)

write.csv(no_improvement_worse,file="./Data_processed/no_improvement_worse.csv")


not_much_improvement=all_error_sig_info[ (all_error_sig_info$sig_id=="insignificant"),]
not_much_improvement= not_much_improvement[order(not_much_improvement$error_count,decreasing= TRUE),]
not_much_improvement_output=not_much_improvement[,c("error_count","mean_simu_error","expected","sd","Var1","Var2","significance","significance","sig_id")]

names(not_much_improvement_output)=c("error_count","mean_simu_error","expected","sd","query_cluster","template_cluster","significance","significance","sig_id")
not_much_improvement_output = as.data.frame(not_much_improvement_output[,c("query_cluster","mean_simu_error","template_cluster","error_count", "expected","significance","sig_id")])
not_much_improvement_output <- apply(not_much_improvement_output,2,as.character)

write.csv(not_much_improvement_output,file="./Data_processed/not_much_improvement_output.csv")



improvement_but_still_many_errors=all_error_sig_info[(all_error_sig_info$sig_id=="smaller"),]
improvement_but_still_many_errors= improvement_but_still_many_errors[order(improvement_but_still_many_errors$error_count,decreasing= TRUE),]
improvement_but_still_many_errors=improvement_but_still_many_errors[,c("error_count","mean_simu_error","expected","sd","Var1","Var2","significance","significance","sig_id")]

names(improvement_but_still_many_errors)=c("error_count","mean_simu_error","expected","sd","query_cluster","template_cluster","significance","significance","sig_id")
improvement_but_still_many_errors = as.data.frame(improvement_but_still_many_errors[,c("query_cluster","mean_simu_error","template_cluster","error_count", "expected","significance","sig_id")])
improvement_but_still_many_errors <- apply(improvement_but_still_many_errors,2,as.character)

write.csv(improvement_but_still_many_errors,file="./Data_processed/improvement_but_still_many_errors.csv")




all_frames_I_care=rbind(rbind(improvement_but_still_many_errors,not_much_improvement_output),no_improvement_worse)
write.csv(all_frames_I_care,"./Data_processed/work.csv")
all_frames_I_care= as.data.frame(all_frames_I_care)
all_frames_I_care$loop_type=unlist(lapply(strsplit(as.character(all_frames_I_care$query_cluster),"-"),function(x){paste(x[1:2],collapse="-")}))
#generate colors
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector=col_vector[c(1:3,5:(n))]
col_vector=c(col_vector,"#bdbdbd")
pie(rep(1,n), col=col_vector)
col_vec=col_vector[1:length(unique(all_frames_I_care$sig_id))]
saveRDS(all_frames_I_care,file="./Data_processed/all_frames_I_care.rds")




abc=ggplot(all_frames_I_care, aes(expected, error_count,colour = factor(sig_id)))+
  geom_point( position=position_dodge(width = 0.90),size = 2.5)+
  scale_colour_manual(breaks = all_frames_I_care$sig_id, name="category",
                      values =col_vec)+
  xlab(" error count from random assignment")+ylab("blindBLAST error count")+   theme_classic()+theme_bw()+
  theme(strip.text.x = element_text(size = 12),legend.position = c(0.75, 0.9),
        axis.text.x = element_text(size = 11, colour = "black"),plot.title =element_text(hjust = 0.5,size=18,face="bold"),
        plot.margin=unit(c(2,2,2,2),"mm"),
        axis.text.y = element_text(size=11))
save_figure_specific_size(abc,"misclassifications_categorized_by_significance_1.pdf",7,7)

all_frames_I_care$mean_simu_error = as.numeric(as.character(all_frames_I_care$mean_simu_error))
all_frames_I_care$error_count = as.numeric(as.character((all_frames_I_care$error_count)))
abc=ggplot(all_frames_I_care, aes(expected, error_count,colour = factor(sig_id)))+
  geom_point( position=position_dodge(width = 0.90),size = 2.5)+
  scale_colour_manual(breaks = all_frames_I_care$sig_id, name="category",
                      values =col_vec)+
  xlab(" error count from random assignment")+ylab("blindBLAST error count")+   theme_classic()+theme_bw()+
  theme(strip.text.x = element_text(size = 12),legend.position = c(0.75, 0.9),
        axis.text.x = element_text(size = 11, colour = "black"),plot.title =element_text(hjust = 0.5,size=18,face="bold"),
        plot.margin=unit(c(2,2,2,2),"mm"),
        axis.text.y = element_text(size=11))+geom_abline(intercept=0, slope=1)
save_figure_specific_size(abc,"suggested_changes/misclassifications_categorized_by_significance_1.pdf",7,7)




print(c("plotted ./Plots/misclassifications_categorized_by_significance_1.pdf"))


read.csv("./Data_processed/work.csv")

