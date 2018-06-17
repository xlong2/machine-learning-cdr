# this script runs the GBM model trainig with parameters specifying which loop and specific model hyperparamters
#args=c( "H1_13", 4, 6, 1500, 0.01, 5)


args <- commandArgs(trailingOnly = TRUE)
print(args)

# load the environment 
currect_d=getwd()
print(c("current directory is ",currect_d))
source("functions.R")
source("utility_function.R")
getwd()
file=paste(c("./tmp/","data_by_loop_type_list_unduplicated.rds"),collapse="")
data_by_loop_type_list_unduplicated=readRDS(file)

# load required packages 
list.of.packages <- c("Matrix","stats", "grid","caret","RWebLogo","MLmetrics","parallel","pryr","protr","gbm","ggplot2","ggrepel","reshape2","gridExtra","doMC","RColorBrewer","e1071")
local_package_dir="~/R_libs/"
specify_R_package_diretory=FALSE
install_and_load_packages(list.of.packages,local_package_dir,specify_R_package_diretory)

args_list=list()
args_list[[1]]=args
args_list[[2]]=data_by_loop_type_list_unduplicated

execute_training_rscript(args_list)
