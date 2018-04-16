library("scales")

conf_tables_all_loops_blindBLAST
conf_tables_all_loops_gbm
accuracy_list=list(); 
for(x_loop in names(conf_tables_all_loops_blindBLAST)){
  bl=conf_tables_all_loops_blindBLAST[[x_loop]]
  accuracy_list[["blindBLAST"]][[x_loop]]=calculate_accuracy(bl,c("Var1", "Var2"),"Freq")
  bgm=conf_tables_all_loops_gbm[[x_loop]]
  accuracy_list[["gbm"]][[x_loop]]=calculate_accuracy(bgm,c("Var1", "Var2"),"Freq")
}

accuracy_gbm_blast=as.data.frame(accuracy_list)
accuracy_gbm_blast_remove_unknow=accuracy_gbm_blast[complete.cases(accuracy_gbm_blast), ]
accuracy_gbm_blast_remove_unknow$loop=split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",1,1,"-")
accuracy_gbm_blast_remove_unknow$length=as.numeric(split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",2,2,"-"))
accuracy_gbm_blast_remove_unknow=reorder_factor(accuracy_gbm_blast_remove_unknow,"loop","length")

accuracy_gbm_blast_remove_unknow_melt=melt(accuracy_gbm_blast_remove_unknow,id.vars = c("loop","length"))
accuracy_gbm_blast_remove_unknow_melt=as.data.frame(accuracy_gbm_blast_remove_unknow_melt)




gbm_folds_sd=lapply(names(all_folds_sd_list),function(x){if(!x%in% names(gbm_folds_sd)){gbm_folds_sd[[x]]=NA};return(gbm_folds_sd[[x]])})
sd_list=list(); sd_list[["blindBLAST"]]=unlist(all_folds_sd_list); sd_list[["gbm"]]=unlist(gbm_folds_sd)
sd_gbm_blast=as.data.frame(sd_list)
sd_gbm_blast_remove_unknow=sd_gbm_blast[complete.cases(sd_gbm_blast), ]
sd_gbm_blast_remove_unknow$loop=split_vector_and_replace(rownames(sd_gbm_blast_remove_unknow),"_",1,1,"-")
sd_gbm_blast_remove_unknow$length=as.numeric(split_vector_and_replace(rownames(sd_gbm_blast_remove_unknow),"_",2,2,"-"))
sd_gbm_blast_remove_unknow=reorder_factor(sd_gbm_blast_remove_unknow,"loop","length")

sd_gbm_blast_remove_unknow_melt=melt(sd_gbm_blast_remove_unknow,id.vars = c("loop","length"))
sd_gbm_blast_remove_unknow_melt=as.data.frame(sd_gbm_blast_remove_unknow_melt)
colnames(sd_gbm_blast_remove_unknow_melt)[4]="sd"
accuracy_sd_gbm_blast_remove_unknow_melt=merge(accuracy_gbm_blast_remove_unknow_melt,sd_gbm_blast_remove_unknow_melt)
accuracy_sd_gbm_blast_remove_unknow_melt$low=accuracy_sd_gbm_blast_remove_unknow_melt$value-accuracy_sd_gbm_blast_remove_unknow_melt$sd/2
accuracy_sd_gbm_blast_remove_unknow_melt$high=accuracy_sd_gbm_blast_remove_unknow_melt$value+accuracy_sd_gbm_blast_remove_unknow_melt$sd/2
accuracy_sd_gbm_blast_remove_unknow_melt=as.data.frame(accuracy_sd_gbm_blast_remove_unknow_melt)
fig=ggplot(accuracy_sd_gbm_blast_remove_unknow_melt,aes(x=length,y=value,ymin = low , ymax = pmin(high,1),fill=variable))+geom_bar(
                                           stat = "identity" ,position=position_dodge(width = 0.90))+
  geom_errorbar(position = position_dodge(0.9))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  xlab("")+facet_grid(~loop,scales="free",space="free")+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 


error_count_list=list(); 

for(x_loop in names(conf_tables_all_loops_blindBLAST_diff)){
  bl=conf_tables_all_loops_blindBLAST_diff[[x_loop]]
  values_bl=bl[as.character(bl$Var1)!=as.character(bl$Var2),"Freq"];
  if(!is.null(values_bl)){
  error_count_list[["blindBLAST"]][[x_loop]]=sum(values_bl)
  }else{
    error_count_list[["blindBLAST"]][[x_loop]]=NA
  }
  bgm=conf_tables_all_loops_gbm_diff[[x_loop]]
  values_gbm=bgm[as.character(bgm$Var1)!=as.character(bgm$Var2),"Freq"];
  if(!is.null(values_gbm)){
  error_count_list[["gbm"]][[x_loop]]=sum(values_gbm)
  }else{
    error_count_list[["gbm"]][[x_loop]]=NA
  }
}
error_count_list_frame=frame_manipulating_for_ploting(error_count_list)
figure_error_count=plot_figure(error_count_list_frame,"length","value","variable","loop",c(0.9,0.9))
figure_error_count=figure_error_count+scale_y_continuous(position = "right")
fig=fig+scale_y_continuous(limits=c(0.5,1.05),oob = rescale_none,position = "right")
grids=list()
grids[[1]]=figure_error_count
grids[[2]]=fig
p=grid.arrange(arrangeGrob(grids[[1]],grids[[2]], ncol=1, nrow=2, heights=c(1,1) ,bottom=textGrob("wrong prediction type", gp=gpar(fontsize=14))  ))

save_figure_specific_size(p,"gbm_blindBLAST_accuracy_and_error_count.pdf",7,7)


