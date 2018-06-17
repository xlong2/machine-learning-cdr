setwd("..")
source('./cdr_classifier_result/0.load_function_and_data.R') 
cores=4;

executable_dir="./cdr_classifier_result/"  # where the Rscript that load data and run major computational task sits
master_condor_file="/Users/longxiyao/condor_script_master.sh"   # not needed if not using condor
condor_script_dir="./condor_script"  # also host shell scripts no needed either way
Rscript_dir = "./training_script/"
Rscript_name = "2.2.1.gbm_train_test_splitted_grid_H1.R"

condor_submit=FALSE
#setwd(condor_script_dir)    # get into the directory where the .sh file and also .con file resides 
batches=c(6,7)

for(batch in batches){
for(loop in c("H1_13","L3_9","H2_10")){
  nl=dim(data_by_loop_type_list_unduplicated[[loop]][[1]])[1]
  eta=max (0.01, 0.1*min(1, nl/10000))
  paras=best_parameters_each_loop[[loop]] #paras[["n.trees"]]
  args=c(loop,cores,paras[["interaction.depth"]],paras[["n.trees"]], eta,paras[["n.minobsinnode"]], batch)
  
  exe_sh_file=generate_Rscript_command_H1(args,Rscript_dir,Rscript_name)
  print(exe_sh_file)
  

  if(!condor_submit){
    system(paste(c("nohup  sh ",exe_sh_file," &"),collapse=""))   # run stuff in parallel
  }else{
    # write condor script
    a= generate_condor_script(args,exe_sh_file,condor_script_dir,master_condor_file)
    
    system(a)
  }
}
}

