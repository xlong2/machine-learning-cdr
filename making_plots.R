# Copyright 2018 Xiyao Long <xlong2@jhu.edu>

#  MIT license
#  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files, to deal 
#  in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#  the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#load functions
file.sources = list.files(pattern="functions.R")
sapply(file.sources,source,.GlobalEnv)
file.sources = list.files(pattern="*function.R")
sapply(file.sources,source,.GlobalEnv)



# Install required R packages if they do not exist
list.of.packages <- c("Matrix", "grid","caret","RWebLogo","MLmetrics","parallel","pryr","protr",
                      "gbm","ggplot2","ggrepel","reshape2","gridExtra","doMC","RColorBrewer",
                      "scales","e1071","argparse")
local_package_dir="~/R_libs/"
specify_R_package_diretory=FALSE
install_and_load_packages(list.of.packages,local_package_dir,specify_R_package_diretory)


suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

initial_message="A program for performing grid search for  Gradient Boosted Model training in x repeats n folds scheme and getting the best parameter;   The implementation features parallelization for not only different models but for different data category so that better parallelization  for datasets of different data sizes are possibe. The training process is also implemented with default sparse class upsampling method to balance data with imbalanced classes"

# by default ArgumentParser will add an help option

parser$add_argument("-o", "--output_dir", type="character",metavar="data_dir",default="Data_dir", help="Directory to place all intermediate data files and results ")
parser$add_argument("-p","--plot_dir",help="directory to place generated plots", default="Plot_dir")

args <- parser$parse_args()
plot_dir = args$plot_dir
result_dir=args$output_dir

file.sources = list.files(pattern="functions.R")
sapply(file.sources,source,.GlobalEnv)
file.sources = list.files(pattern="*function.R")
sapply(file.sources,source,.GlobalEnv)


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





accuracy_sd_gbm_blast_remove_unknow_melt=process_results_for_accuracies(conf_tables_all_loops_blindBLAST, conf_tables_all_loops_gbm)
fig=ggplot(accuracy_sd_gbm_blast_remove_unknow_melt,aes(x=length,y=value,ymin = low , ymax = pmin(high,1),fill=variable))+geom_bar(
  stat = "identity" ,position=position_dodge(width = 0.90))+
  geom_errorbar(position = position_dodge(0.9))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  xlab("")+facet_grid(~loop,scales="free",space="free")+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 

save_figure_specific_size(fig,"gbm_BLAST_accuracy_and_error_count.pdf",7,7)



error_count_list_frame=process_error_counts(conf_tables_all_loops_blindBLAST_diff,conf_tables_all_loops_gbm_diff)



figure_error_count=ggplot(error_count_list_frame,aes(x=length,y=value,ymin = low , ymax = high,fill=variable))+geom_bar(
  stat = "identity" ,position=position_dodge(width = 0.90))+
  geom_errorbar(position = position_dodge(0.9))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  xlab("")+facet_grid(~loop,scales="free",space="free")+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 


figure_error_count=figure_error_count+scale_y_continuous(position = "right")
fig=fig+scale_y_continuous(limits=c(0.5,1.05),oob = rescale_none,position = "right")
grids=list()
grids[[1]]=figure_error_count
grids[[2]]=fig
p=grid.arrange(arrangeGrob(grids[[1]],grids[[2]], ncol=1, nrow=2, heights=c(1,1) ,bottom=textGrob("wrong prediction type", gp=gpar(fontsize=14))  ))
p
save_figure_specific_size(p,"gbm_BLAST_accuracy_and_error_count_comparison.pdf",7,7)


