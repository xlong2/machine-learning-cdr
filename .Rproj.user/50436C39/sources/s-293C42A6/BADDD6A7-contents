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

fig=plot_figure(accuracy_gbm_blast_remove_unknow_melt,"length","value","variable","loop",c(0.9,0.9))

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
fig=fig+scale_y_continuous(position = "right")
grids=list()
grids[[1]]=figure_error_count
grids[[2]]=fig
p=grid.arrange(arrangeGrob(grids[[1]],grids[[2]], ncol=1, nrow=2, heights=c(1,1) ,bottom=textGrob("wrong prediction type", gp=gpar(fontsize=14))  ))

save_figure_specific_size(p,"gbm_blindBLAST_accuracy_and_error_count.pdf",7,7)


