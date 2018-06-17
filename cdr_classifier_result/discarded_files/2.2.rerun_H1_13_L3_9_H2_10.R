setwd("..")
source('./proline_classifier/0.load_function_and_data.R') 
core=4;

executable_dir="./proline_classifier/"  # where the Rscript that load data and run major computational task sits
master_condor_file="/Users/longxiyao/condor_script_master.sh"   # not needed if not using condor
condor_script_dir="/Users/longxiyao/Google\ Drive/2018_spring/lab_work_data/machine_learning_cdr/proline_classifier/condor_script"  # also host shell scripts no needed either way

condor_submit=FALSE
setwd(condor_script_dir)    # get into the directory where the .sh file and also .con file resides 

for(loop in c("H1_13", "H2_10", "L3_9")){
  nl=dim(data_by_loop_type_list_unduplicated[[loop]][[1]])[1]
  eta=max (0.01, 0.1*min(1, nl/10000))
  paras=best_parameters_each_loop[[loop]]
  args=c(loop,cores,paras[["interaction.depth"]],15, eta,5)
  
  exe_sh_file=generate_Rscript_command(args,Rscript_dir,Rscript_name)
  print(exe_sh_file)
  
  
  
  if(!condor_submit){
    system(paste(c("nohup  sh ",exe_sh_file," &"),collapse=""))   # run stuff in parallel
  }else{
    # write condor script
    a= generate_condor_script(args,exe_sh_file,condor_script_dir,master_condor_file)
    
    system(a)
  }
}
