# this script generate machine learning model parameters and launch jobs for training models in parallel  


current_d=getwd()
if(!grepl("proline_classifier",current_d)){
  setwd("./proline_classifier")
  
}else{
}
source('./0.load_function_and_data.R')  

Rscript_dir = "./training_script/"
system(paste(c("mkdir ",Rscript_dir),collapse = " "))
executable_dir="./"
Rscript_name = "6.1.gbm_train_final_model.R"


cores=4;

a_new_script=TRUE

master_condor_file="/Users/longxiyao/condor_script_master.sh"   # not needed if not using condor
condor_script_dir="/Users/longxiyao/Google\ Drive/2018_spring/lab_work_data/machine_learning_cdr/proline_classifier/condor_script"  # also host shell scripts no needed either way
condor_submit=FALSE
   # get into the directory where the .sh file and also .con file resides 
parallel_install_max=1
count=0
para_script="run_final_models_master.sh"
for(loop in names(data_by_loop_type_list_unduplicated)){
  nl=dim(data_by_loop_type_list_unduplicated[[loop]][[1]])[1]
  eta=max (0.01, 0.1*min(1, nl/10000))
  paras=best_parameters_each_loop[[loop]]   # find the corresponding best parameters 
  args=c(loop,cores,paras[["interaction.depth"]],paras[["n.trees"]], eta,paras[["n.minobsinnode"]])
  
  
  exe_sh_file=generate_Rscript_command(args,Rscript_dir,Rscript_name)
  print(exe_sh_file)
  count= parallelize_jobs(parallel_install_max,count,a_new_script,exe_sh_file, para_script)
  a_new_script=FALSE
  if(!condor_submit){
#    system(paste(c("nohup  sh ",exe_sh_file," &"),collapse=""))   # run stuff in parallel
  }else{
    # write condor script
    setwd(condor_script_dir) 
    a= generate_condor_script(args,exe_sh_file,condor_script_dir,master_condor_file)
    
    system(a)
  }
}
print(c("running ", para_script))
system(paste(c("sh ",para_script),collapse=""))
