# Purpose:
#    Read in data table derived from Ipygclassify
#    Install all required R libraries that are required for other scripts as well. 


direct=getwd()
if(!grepl("cdr_classifier_result",direct)){
  setwd("./cdr_classifier_result/")  # enter directory ./cdr_classifier_result
  
}

# load all the functions
file.sources = list.files(pattern="functions.R")
sapply(file.sources,source,.GlobalEnv)
file.sources = list.files(pattern="*function.R")
sapply(file.sources,source,.GlobalEnv)
#setwd(direct)  # go back to previous directory 


# this script runs from the Rstudio environment after loading the R project 
list.of.packages <- c("Matrix", "grid","caret","RWebLogo","MLmetrics","parallel","pryr","protr",
                      "gbm","ggplot2","ggrepel","reshape2","gridExtra","doMC","RColorBrewer",
                      "scales","e1071")
local_package_dir="~/R_libs/"
specify_R_package_diretory=FALSE
install_and_load_packages(list.of.packages,local_package_dir,specify_R_package_diretory)


cur_d=getwd()
prefix=strsplit(cur_d,"\\/")[[1]][length(strsplit(cur_d,"\\/")[[1]])]
g(the_method,each_method,cluster_dis,result_dir,subsitution_matrix,subsitution_matrix_name,plot_dir,file)%=% specify_variable_names(prefix,"Data_processed","Plots","data_by_loop_type_list_unduplicated_for_blindBLAST")
gbm_models_dir="./rmsd_cluster_hits_rmsd/"

data_by_loop_type_list_unduplicated=readRDS(file)
data_by_loop_type_list_unduplicated=data_by_loop_type_list_unduplicated_for_blindBLAST
unique(data_by_loop_type_list_unduplicated_for_blindBLAST[[1]][[1]]$cluster_type)
data(AAPAM30)
AAPAM30_index=1:20
names(AAPAM30_index)=rownames(AAPAM30)
plot_dir="./Plots/"



# load all the data
setwd("./Data_processed/")  # enter the specified directory 
data.sources = list.files(pattern="*.rds")
for(x in data.sources){
 print(x)
  tryCatch({
  x_name=strsplit(x,"\\.")[[1]][1]
  assign(x_name,readRDS(x));#print("successfully loaded")
  },error=function(e) {e})
}
#data_by_loop_type_list_unduplicated=data_by_loop_type_list_unduplicated_no_filtering_original

print("line 46")
setwd("../")  # go back to cdr_classifier_result directory


