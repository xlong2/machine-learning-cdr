# getting result of blindBLAST in 3 repeat 10 fold cv
# Getting rmsds of blindBLAST template vs query CDR
# need to use generate_training_grids.R for passing parameters
# this script generate machine learning model parameters and launch jobs for training models in parallel
# pass arguments for how many cores available for computation. 
# user should store the data file in the same directory
# load all the functions

#!/usr/bin/env Rscript
# Copyright 2012-2013 Trevor L Davis <trevor.l.davis@gmail.com>
# Copyright 2008 Allen Day
#
#  This file is free software: you may copy, redistribute and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation, either version 2 of the License, or (at your
#  option) any later version.
#
#  This file is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

initial_message="A program for performing grid search for  Gradient Boosted Model training in x repeats n folds scheme and getting the best parameter;   The implementation features parallelization for not only different models but for different data category so that better parallelization  for datasets of different data sizes are possibe. The training process is also implemented with default sparse class upsampling method to balance data with imbalanced classes"



# by default ArgumentParser will add an help option
parser$add_argument("-d", "--data ",type="character", default="data_table.tsv",  help="the data file")

parser$add_argument("-o", "--output_dir", type="character",metavar="data_dir",default="Data_dir", help="Directory to place all intermediate data files and results ")
parser$add_argument("-p","--plot_dir",help="directory to place generated plots", default="Plot_dir")
parser$add_argument( "-re", "--n_repeats",help=" the number of repeats for performing the blindBLAST prediction ", default=3, type="integer")
parser$add_argument( "-fold", "--n_folds",help=" the number of repeats for performing the blindBLAST prediction ", default=10, type ="integer" )
parser$add_argument("-c", "--n_cores",help=" number of cores available for computing",default=1,type="integer")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()



file.sources = list.files(pattern="functions.R")
sapply(file.sources,source,.GlobalEnv)
file.sources = list.files(pattern="*function.R")
sapply(file.sources,source,.GlobalEnv)
#setwd(direct)  # go back to previous directory 



# this script runs from the Rstudio environment after loading the R project 
list.of.packages <- c("Matrix", "grid","caret","RWebLogo","MLmetrics","parallel","pryr","protr",
                      "gbm","ggplot2","ggrepel","reshape2","gridExtra","doMC","RColorBrewer",
                      "scales","e1071","argparse")
local_package_dir="~/R_libs/"
specify_R_package_diretory=FALSE
install_and_load_packages(list.of.packages,local_package_dir,specify_R_package_diretory)





current_d = getwd()
if (!grepl("proline_classifier", current_d)) {
  setwd("./proline_classifier")
  
} 
cur_d=getwd()
prefix=strsplit(cur_d,"\\/")[[1]][length(strsplit(cur_d,"\\/")[[1]])]
data(AAPAM30)
AAPAM30_index=1:20
names(AAPAM30_index)=rownames(AAPAM30)




# parameters passed to Rscript for gbm model training
total_max_core=args$n_cores

n_repeats=args$n_repeats
n_folds=args$n_folds
plot_dir=args$plot_dir
output_dir=args$output_dir
command_mk=paste(c("mkdir", output_dir),collapse=" ")
system(command_mk)

data_file=args$data
g(the_method,each_method,cluster_dis,result_dir,subsitution_matrix,subsitution_matrix_name,plot_dir,file)%=% specify_variable_names(prefix, output_dir,  plot_dir ,data_file)

data_by_loop_type_list_unduplicated=read_in_data(data_file)
system("mkdir ./tmp")
save_file_sp("data_by_loop_type_list_unduplicated",result_dir,data_file)
data_save_f=paste(c(result_dir,data_file,".rds"),collapse="")



each_method = "blindblast_just_get_alignment"
cluster_dis = "north"
subsitution_matrix_name = "wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
subsitution_matrix = "PAM30"


# do blindBLAST template searching for each CDR loop 
overall_accuracy = list()
ten_foldcv_blindblastlist = list()
the_loop_names=names(data_by_loop_type_list_unduplicated)
left_loops=the_loop_names[!the_loop_names%in% names(ten_foldcv_blindblastlist) ]
for (loop_type in left_loops) {

  sequences=data_by_loop_type_list_unduplicated[[loop_type]][[1]]
  each_loop_length_data_feature_string = data_by_loop_type_list_unduplicated[[loop_type]][[2]]
  features=data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  sequences$cluster_type = as.character(sequences$cluster_type)
  all_cases =  sequences[, c(features, "cluster_type", "identifier", "id")]
  #make folds and separates folds in and folds out
  r <- n_repeats # number of repeats
  k <- n_folds # number of folds
  g(folds.list,folds.list.out)  %=%generate_folds_foldsout(sequences,r,k)
  
  all_result_list = mclapply(1:length(folds.list), get_accuracy_per_fold_overload_final, mc.cores =  total_max_core  )  # run blindBLAST for each folds out 
  acc_result=calculate_accuracy_mean_std(all_result_list)
  overall_accuracy[[loop_type]] = c(acc_result[[1]], acc_result[[2]])
  ten_foldcv_blindblastlist[[loop_type]] = all_result_list
  
}# end of iterating all loop types



save_file("ten_foldcv_blindblastlist")




results_rbind_allloops=lapply(ten_foldcv_blindblastlist,function(each_r){
  if(n_repeats>1){
  result_l=chunk2(1:length(each_r),n_repeats)
  }else{
    result_l=list()
    result_l[[1]]=1:length(each_r)
  }
  #  for(each_r in ten_foldcv_blindblastlist){
  results=lapply(result_l,function(x){
    fr=do.call(rbind,each_r[x])
    fr$obs=split_vector_and_replace(fr[,1],"\\.",1,1,"")
    fr$pred=split_vector_and_replace(fr[,2],"\\.",1,1,"")
    #dim(fr[fr$obs==fr$pred,])[1]/dim(fr)[1]
    repeat_number=x[1]%/%10+1
    print(dim(fr)[1])
    print(each_r)
    fr$repeats=rep(repeat_number,dim(fr)[1])
    return(fr)
  })
  results_rbind=do.call(rbind,results)
  return(results_rbind)
  #}
  
  
})

blindBLAST_summ=cal_blindBLAST_accuracy(ten_foldcv_blindblastlist,n_folds)
print(blindBLAST_summ)  
# calculate accuracies for blindBLAST for each repeat
# calculate the std for accuracies for each repeat
