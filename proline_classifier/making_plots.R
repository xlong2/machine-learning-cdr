# Copyright 2018 Xiyao Long <xlong2@jhu.edu>

#  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files, to deal 
#  in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#  the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

initial_message="A program for performing grid search for  Gradient Boosted Model training in x repeats n folds scheme and getting the best parameter;   The implementation features parallelization for not only different models but for different data category so that better parallelization  for datasets of different data sizes are possibe. The training process is also implemented with default sparse class upsampling method to balance data with imbalanced classes"

# by default ArgumentParser will add an help option

parser$add_argument("-o", "--output_dir", type="character",metavar="data_dir",default="Data_dir", help="Directory to place all intermediate data files and results ")
parser$add_argument("-p","--plot_dir",help="directory to place generated plots", default="Plot_dir")
args <- parser$parse_args()


result_dir=args$output_dir

setwd(c(result_dir))  # enter the specified directory 
data.sources = list.files(pattern="*.rds")
for(x in data.sources){
  print(x)
  tryCatch({
    x_name=strsplit(x,"\\.")[[1]][1]
    assign(x_name,readRDS(x));#print("successfully loaded")
  },error=function(e) {e})
}
setwd("../")
# comparison between the two methods
blindblast_by_loop=lapply(ten_foldcv_blindblastlist,function(x){if(n_repeats>1){split_c=chunk2(1:length(x),n_repeats)}else{split_c=list();split_c[[1]]=1:length(x)}; lapply(split_c,function(y){do.call(rbind,x[y])})})
blindblast_by_loop_bind=lapply(blindblast_by_loop,function(x){
  accuracies=lapply(x,function(y){
    all=y
    all[,1]=split_vector_and_replace(all[,1],"\\.",1,1,"")
    all[,2]=split_vector_and_replace(all[,2],"\\.",1,1,"")
    accu=dim(all[all[,1]==all[,2],])[1]/dim(all)[1]
  })
  accuracies=unlist(accuracies)
  
})

blindBLAST_mean_accu=lapply(blindblast_by_loop_bind,mean)
blindBLAST_accu_std = lapply(blindblast_by_loop_bind,sd)

# comparison between the two methods on misclassification specific error counts 

#blind_blast_cv_result_summary = as.data.frame(do.call(rbind, overall_accuracy))
blind_blast_cv_result_summary=as.data.frame(cbind(blindBLAST_mean_accu,blindBLAST_accu_std ))


#rownames(blind_blast_cv_result_summary) = names(overall_accuracy)
colnames(blind_blast_cv_result_summary) = c("mean", "sd")
blind_blast_cv_result_summary$loop_type = rownames(blind_blast_cv_result_summary)
save_file("blind_blast_cv_result_summary")



# summarize result for blindBLAST
blind_blast_cv_result_summary$loop=split_vector_and_replace(blind_blast_cv_result_summary$loop_type,"_",1,1,"-")
blind_blast_cv_result_summary$length=split_vector_and_replace(blind_blast_cv_result_summary$loop_type,"_",2,2,"-")
blind_blast_cv_result_summary_reordered=reorder_factor(blind_blast_cv_result_summary,"loop","length")
#blindBLAST_accu_std


blind_blast_cv_result_summary_reordered_selected=blind_blast_cv_result_summary_reordered[,c( "mean" ,    "sd" , "loop", "length")]
blind_blast_cv_result_summary_reordered_selected$method=rep("blindBLAST",dim(blind_blast_cv_result_summary_reordered_selected)[1])

















#colnames(random_acc_sd_reordered)[1:2]=c("mean","sd")
#random_acc_sd_reordered$method=rep("random",dim(random_acc_sd_reordered)[1])

# combine blindBLAST result and that of random assignment 
combined_blindBLAST_random=as.data.frame(rbind(blind_blast_cv_result_summary_reordered_selected,random_acc_sd_reordered))
combined_blindBLAST_random$sd=as.numeric(combined_blindBLAST_random$sd)
combined_blindBLAST_random$mean=as.numeric(combined_blindBLAST_random$mean)
combined_blindBLAST_random$min=combined_blindBLAST_random$mean-combined_blindBLAST_random$sd/2
combined_blindBLAST_random$max=combined_blindBLAST_random$mean+combined_blindBLAST_random$sd/2
save_file("combined_blindBLAST_random")





conf_tables_all_loops_blindBLAST
conf_tables_all_loops_gbm
accuracy_list=list(); 
for(x_loop in names(conf_tables_all_loops_blindBLAST)){
  bl=conf_tables_all_loops_blindBLAST[[x_loop]]
  accuracy_list[["blindBLAST"]][[x_loop]]=calculate_accuracy(bl,c("Var1", "Var2"),"Freq")
  bgm=conf_tables_all_loops_gbm[[x_loop]]
  accuracy_list[["gbm"]][[x_loop]]=calculate_accuracy(bgm,c("Var1", "Var2"),"Freq")
}
ori=accuracy_list
accuracy_list[["blindBLAST"]]=unlist(blindBLAST_mean_accu)
accuracy_list[["gbm"]]=unlist(gbm_accuracy_by_loops)
common=intersect(names(accuracy_list[["gbm"]]),names(accuracy_list[["blindBLAST"]]))
accuracy_list=lapply(accuracy_list,function(x){x[common]})
accuracy_gbm_blast=as.data.frame(accuracy_list)
accuracy_gbm_blast_remove_unknow=accuracy_gbm_blast[complete.cases(accuracy_gbm_blast), ]
accuracy_gbm_blast_remove_unknow$loop=split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",1,1,"-")
accuracy_gbm_blast_remove_unknow$length=as.numeric(split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",2,2,"-"))
accuracy_gbm_blast_remove_unknow=reorder_factor(accuracy_gbm_blast_remove_unknow,"loop","length")

accuracy_gbm_blast_remove_unknow_melt=melt(accuracy_gbm_blast_remove_unknow,id.vars = c("loop","length"))
accuracy_gbm_blast_remove_unknow_melt=as.data.frame(accuracy_gbm_blast_remove_unknow_melt)




gbm_folds_sd=lapply(names(all_folds_sd_list),function(x){if(!x%in% names(gbm_folds_sd)){gbm_folds_sd[[x]]=NA};return(gbm_folds_sd[[x]])})
sd_list=list(); sd_list[["blindBLAST"]]=unlist(blindBLAST_accu_std); sd_list[["gbm"]]=unlist(gbm_folds_sd)
sd_list[["gbm"]]=unlist(gbm_sd_by_loops)
common=intersect(names(sd_list[["gbm"]]),names(sd_list[["blindBLAST"]]))
sd_list=lapply(sd_list,function(x){x[common]})
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






