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
parser$add_argument(  "-md", "--gbm_models_dir",help= "directory to place gbm trained models",default= "Gbm_models")
parser$add_argument(  "-comp", "--complexity", help="the complexity levels interested in testing separated by ':'",default="3:6")
parser$add_argument( "-tn", "--n_trees",help=" the number of trees interested in testing, separated by ':'", default="5:10")
parser$add_argument( "-re", "--n_repeats",help=" the number of repeats for performing the model training ", default=3, type="integer")
parser$add_argument( "-fold", "--n_folds",help=" the number of repeats for performing the model training ", default=10, type ="integer" )
parser$add_argument("-c", "--n_cores",help=" number of cores available for computing",default=1,type="integer")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
print_info<-function(x){
  
}


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
complexity = as.numeric(strsplit(args$complexity,":")[[1]])
gbm_models_dir=args$gbm_models_dir
command_mk=paste(c("mkdir", gbm_models_dir),collapse=" ")
system(command_mk)
trees = as.numeric(strsplit(args$n_trees,":")[[1]])

n_repeats=args$n_repeats
n_folds=args$n_folds
plot_dir=args$plot_dir
output_dir=args$output_dir
min_node_n = 5
data_file=args$data
g(the_method,each_method,cluster_dis,result_dir,subsitution_matrix,subsitution_matrix_name,plot_dir,file)%=% specify_variable_names(prefix, output_dir,  plot_dir ,data_file)

data_by_loop_type_list_unduplicated=read_in_data(data_file)
system("mkdir ./tmp")
save_file_sp("data_by_loop_type_list_unduplicated",result_dir,data_file)
data_save_f=paste(c(result_dir,data_file,".rds"),collapse="")

for (loop in names(data_by_loop_type_list_unduplicated)) {
  eta = generate_eta(loop)   # generate learning rate


  arguments=generate_arguments(loop, total_max_core, complexity, trees, eta, min_node_n)
  mclapply(arguments,gbm_train_script,mc.cores=1)
  

}



# read all models



all_models = list.files(pattern = "*extra_test.rds",
                        path = gbm_models_dir,
                        full.names = FALSE)
all_models=all_models[grepl("gbm",all_models)]
all_models_list_by_loop = list()
index_i = 1
all_t = list()
for (each_file in all_models) {
  tryCatch({
    each_file_r = paste(c(gbm_models_dir,"/", each_file), collapse =
                          "")
    the_model = readRDS(each_file_r)
    
    index_i = index_i + 1
    print(each_file)
    info = strsplit(each_file, "\\/")[[1]][length(strsplit(each_file, "\\/")[[1]])]
    loop_type = paste(strsplit(info, "_")[[1]][1:2], collapse = "_")
    all_t[[loop_type]] = unique(the_model$trainingData$.outcome)
    number = paste(strsplit(strsplit(info, "_")[[1]][4], "-")[[1]][3:6], collapse =
                     "-")
    model_result = the_model$result
    all_models_list_by_loop[[loop_type]][[number]] = model_result
    rm(the_model)
    gc()
  }, error = function(e) {
    print(e)
  })
}





all_result = lapply(all_models_list_by_loop, function(x) {
  z=lapply(x,function(y){
    y=as.data.frame(y)
    return(y)})
  w=do.call(rbind, z)

  return(w)
  
})



save_file("all_result")


# find the best parameter for each loop and length type
best_para_result_list=get_best_para(all_gird_search_results,all_result)
best_parameters_each_loop=best_para_result_list[[2]]
save_file("best_parameters_each_loop")

result_table=best_para_result_list[[1]]

print(result_table)


