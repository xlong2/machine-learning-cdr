# this script is a refactor of the Untitled.R in the R/ directory
# Purpose:
#    Read in a speficic loop type 
#    Do LOOCV for blind blast
#    Refer to the rmsd_table to get the rmsd from the template
#    Output the prediction result 
#    Do a blast search within the predicted cluster; Or a blast search outside of the clusters or globally if it is not identified as any cluster. 
#    Retrieve the rmsd between the predicted sequence and query sequence. Record it 
#    Also this script should enable blind blast search as currently employed in the rosetta
#    It should output the data format so that comparison between the prediction accuracy and rmsds are easy

# output:
#   File1:    Confusion table with rownames be the real cases and colnames be the prediction cases of all the predictions curated from LOOCV
#             rmsd_cluster_guided_blast_rmsd_",loop_type,"_",paste(c(each_method,cluster_dis),collapse="-"),"_conftable.rds"
#   File2:    Mean of the three best templated for each sequence and stored as a table  :
#               seq
#                H1_13-1.4312    1.003    1.0023  2.323  3.22
#                ...
#    
# Input:
#   
#   loop_type : "H2_9"
#   the_method  "blindblast"
#   cluster_dis :   "north"  / "rmsd"

list.of.packages <- c("Matrix", "grid","caret","MLmetrics","parallel","pryr","protr","gbm","ggplot2","reshape2","gridExtra","doMC","RColorBrewer","e1071")
local_package_dir="~/R_libs/"
install_and_load_packages<-function(package_list,local_package_dir){
new.packages <- list.of.packages[!(package_list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,lib=local_package_dir)
for(pack in package_list){
  library(pack,character.only = T)
}
}
install_and_load_packages(list.of.packages,local_package_dir)

#generate colors
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector=col_vector[c(1:3,5:(n))]
col_vector=c(col_vector,"#bdbdbd")
pie(rep(1,n), col=col_vector)


  the_method="blindblast_just_get_alignment"
  cluster_dis="north"





overall_prefix="./proline_classifier/Data_processed"

file=paste(c(overall_prefix,"/data_by_loop_type_list_unduplicated_for_blindBLAST.rds"),collapse = "")
data_by_loop_type_list_unduplicated=readRDS(file)
#data_by_loop_type_list_unduplicated=data_by_loop_type_list_unduplicated_for_blindBLAST

subsitution_matrix_name ="wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
subsitution_matrix="PAM30"
result_dir = "./proline_classifier/Data_processed/"
mkcommand=paste(c("mkdir ",result_dir), collapse=" ")
system(mkcommand)
methods= c(the_method) # specify the machine learning methods to be used 
each_method=the_method

# end of iterating all folds

the_method="somerandommethod"

data(AAPAM30)
AAPAM30_index=1:20
names(AAPAM30_index)=rownames(AAPAM30)
direct=getwd()
plot_dir="./proline_classifier/Plots/"

save_file1="./proline_classifier/Data_processed/overall_accuracy.rds"
save_file3 = paste(c(result_dir,"all_pred_tables_realcluster_list",".rds"),collapse="")
all_pred_tables_realcluster_list=readRDS(save_file3)
overall_accuracy=readRDS(save_file1)

setwd("./proline_classifier/")

file.sources = list.files(pattern="*functions.R")
sapply(file.sources,source,.GlobalEnv)
file.sources = list.files(pattern="*function.R")
sapply(file.sources,source,.GlobalEnv)
setwd(direct)
setwd("./proline_classifier/Data_processed/")
data.sources = list.files(pattern="*.rds")
for(x in data.sources){
  print(x)
  tryCatch({
  x_name=strsplit(x,"\\.")[[1]][1]
  assign(x_name,readRDS(x));print("successfully loaded")},error=function(e) {e})
}

lapply(data_by_loop_type_list_unduplicated,function(x){dim(x[[1]])})

# modify all all_similarity_matrix names to H1-13-

# modify the nofiltering to make the cluster_type column and also the identifier column to match to that of the 
setwd(direct)


