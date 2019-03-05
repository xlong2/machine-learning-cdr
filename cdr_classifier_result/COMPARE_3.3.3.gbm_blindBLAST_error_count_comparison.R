# this script compare the misclassification error counts between GBM and blindBLAST CDR cluster prediction result
current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

conf_tables_all_loops_blindBLAST
conf_tables_all_loops_gbm = merged_testing_result
accuracy_list=list(); 
for(x_loop in names(conf_tables_all_loops_blindBLAST)){
  bl=conf_tables_all_loops_blindBLAST[[x_loop]]
  accuracy_list[["blindBLAST"]][[x_loop]]=calculate_accuracy(bl,c("Var1", "Var2"),"Freq")
  bgm=conf_tables_all_loops_gbm[[x_loop]]
  print(x_loop)
  print(bgm)
  accuracy_list[["gbm"]][[x_loop]]=calculate_accuracy(bgm,c("Var1", "Var2"),"Freq")
}
ori=accuracy_list
accuracy_list[["blindBLAST"]]=unlist(accuracy_list[["blindBLAST"]])
accuracy_list[["gbm"]]=unlist(accuracy_list[["gbm"]])
common=intersect(names(accuracy_list[["gbm"]]),names(accuracy_list[["blindBLAST"]]))
accuracy_list=lapply(accuracy_list,function(x){x[common]})
accuracy_gbm_blast=as.data.frame(accuracy_list)
accuracy_gbm_blast_remove_unknow=accuracy_gbm_blast[complete.cases(accuracy_gbm_blast), ]

accuracy_gbm_blast_remove_unknow$loop=split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",1,1,"-")
accuracy_gbm_blast_remove_unknow$length=as.numeric(split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",2,2,"-"))
accuracy_gbm_blast_remove_unknow=reorder_factor(accuracy_gbm_blast_remove_unknow,"loop","length")

accuracy_gbm_blast_remove_unknow_melt=melt(accuracy_gbm_blast_remove_unknow,id.vars = c("loop","length"))
accuracy_gbm_blast_remove_unknow_melt=as.data.frame(accuracy_gbm_blast_remove_unknow_melt)



gbm_folds_sd=gbm_per_loop_sd
sd_list=list(); sd_list[["blindBLAST"]]=unlist(blindBLAST_accu_std); sd_list[["gbm"]]=unlist(gbm_folds_sd)
common=intersect(names(sd_list[["gbm"]]),names(sd_list[["blindBLAST"]]))
sd_list=lapply(sd_list,function(x){x[common]})
sd_gbm_blast=as.data.frame(sd_list)


sd_gbm_blast$loop = rownames(sd_gbm_blast)
colnames(sd_gbm_blast)[1:2]=c("blindBLAST","gbm")




sd_gbm_blast_remove_unknow=sd_gbm_blast[complete.cases(sd_gbm_blast), ]
colnames(sd_gbm_blast_remove_unknow)[1:2]=c("blindBLAST","gbm")
rownames(sd_gbm_blast_remove_unknow)=sd_gbm_blast_remove_unknow$loop
sd_gbm_blast_remove_unknow$loop = sapply(strsplit(sd_gbm_blast_remove_unknow$loop,"_"),"[[",1)
sd_gbm_blast_remove_unknow$length= sapply(strsplit(rownames(sd_gbm_blast_remove_unknow),"_"),"[[",2)
sd_gbm_blast_remove_unknow=reorder_factor(sd_gbm_blast_remove_unknow,"loop","length")

sd_gbm_blast_remove_unknow_melt=melt(sd_gbm_blast_remove_unknow,id.vars = c("loop","length"))
sd_gbm_blast_remove_unknow_melt=as.data.frame(sd_gbm_blast_remove_unknow_melt)
colnames(sd_gbm_blast_remove_unknow_melt)[4]="sd"
colnames(accuracy_gbm_blast_remove_unknow_melt)[4]="accuracy"
accuracy_sd_gbm_blast_remove_unknow_melt=merge(accuracy_gbm_blast_remove_unknow_melt,sd_gbm_blast_remove_unknow_melt)
accuracy_sd_gbm_blast_remove_unknow_melt$low=accuracy_sd_gbm_blast_remove_unknow_melt$accuracy-accuracy_sd_gbm_blast_remove_unknow_melt$sd/2
accuracy_sd_gbm_blast_remove_unknow_melt$high=accuracy_sd_gbm_blast_remove_unknow_melt$accuracy+accuracy_sd_gbm_blast_remove_unknow_melt$sd/2
accuracy_sd_gbm_blast_remove_unknow_melt=as.data.frame(accuracy_sd_gbm_blast_remove_unknow_melt)
fig=ggplot(accuracy_sd_gbm_blast_remove_unknow_melt,aes(x=length,y=accuracy,ymin = low , ymax = pmin(high,1),fill=variable))+geom_bar(
  stat = "identity" ,position=position_dodge(width = 0.90))+
  geom_errorbar(position = position_dodge(0.9))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  xlab("")+facet_grid(~loop,scales="free",space="free")+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 






error_count_list=list(); 
gbm_errorcount_list= lapply(conf_tables_all_loops_gbm,function(x){
  sum(x[as.character(x$Var1)!=as.character(x$Var2),"Freq"])
})
blindBLAST_errorcount_lists
for(x_loop in names(conf_tables_all_loops_blindBLAST_diff)){
  print(x_loop)
  bl=conf_tables_all_loops_blindBLAST_diff[[x_loop]]
  values_bl=bl[as.character(bl$Var1)!=as.character(bl$Var2),"Freq"];
  if(!is.null(values_bl)){
    print(sum(values_bl))
    error_count_list[["blindBLAST"]][[x_loop]]=sum(values_bl)
  }else{
    error_count_list[["blindBLAST"]][[x_loop]]=NA
  }
  bgm=conf_tables_all_loops_gbm_diff[[x_loop]]
  values_gbm=bgm[as.character(bgm$Var1)!=as.character(bgm$Var2),"Freq"];
  if(!is.null(values_gbm)){
    print(sum(values_gbm))
    error_count_list[["gbm"]][[x_loop]]=sum(values_gbm)
  }else{
    error_count_list[["gbm"]][[x_loop]]=NA
  }
}
error_count_list_with_sd=list()
t1=as.data.frame(error_count_list[["blindBLAST"]])
t1$loops=rownames(t1)
t2=as.data.frame(t(as.data.frame(lapply(blindBLAST_errorcount_lists,sd)))); t2$loops=rownames(t2)
fr1=merge(t1,t2,by="loops")
colnames(fr1)=c("loop","error_count","sd")
rownames(fr1)=fr1$loop

t1=as.data.frame(error_count_list[["gbm"]])
t1$loops=rownames(t1)
t2=as.data.frame(t(as.data.frame(gbm_error_count_sd))); t2$loops=rownames(t2)
fr2=merge(t1,t2,by="loops")
colnames(fr2)=c("loop","error_count","sd")
rownames(fr2)=fr2$loop
common_n=intersect(rownames(fr1),rownames(fr2)); fr1=fr1[rownames(fr1)%in%common_n,]; fr2=fr2[rownames(fr2)%in%common_n,]
error_count_list_with_sd[[1]]=fr1
error_count_list_with_sd[[2]]=fr2
names(error_count_list_with_sd)=c("blindBLAST","gbm")
error_count_list_frame=do.call(rbind,lapply(names(error_count_list_with_sd),frame_for_plot))
colnames(error_count_list_frame)=c("loop","value","sd","length","variable")

figure_error_count=plot_figure(error_count_list_frame,"length","value","variable","loop",c(0.9,0.9))
error_count_list_frame$low=error_count_list_frame$value-error_count_list_frame$sd/2
error_count_list_frame$high=error_count_list_frame$value + error_count_list_frame$sd/2


figure_error_count=ggplot(error_count_list_frame,aes(x=length,y=value,ymin = low , ymax = high,fill=variable))+geom_bar(
  stat = "identity" ,position=position_dodge(width = 0.90))+
  geom_errorbar(position = position_dodge(0.9))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  xlab("")+facet_grid(~loop,scales="free",space="free")+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1),strip.text.y = element_text(angle =-90)) 


figure_error_count=figure_error_count+scale_y_continuous(position = "right")
fig=fig+scale_y_continuous(limits=c(0.4,1.05),oob = rescale_none,position = "right")
grids=list()
grids[[1]]=figure_error_count
grids[[2]]=fig
p=grid.arrange(arrangeGrob(grids[[1]],grids[[2]], ncol=1, nrow=2, heights=c(1,1) ,bottom=textGrob("wrong prediction type", gp=gpar(fontsize=14))  ))

save_figure_specific_size(p,"suggested_changes/gbm_blindBLAST_accuracy_and_error_count.pdf",7,7)

