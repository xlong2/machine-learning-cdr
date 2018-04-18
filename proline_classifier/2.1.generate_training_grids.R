setwd("..")
source('./proline_classifier/0.load_function_and_data.R') 
data_by_loop_type_list_unduplicated
cores=4;
complexity=c(3,6,9);
trees=c(100,300,500,1000,1500, 2000,2500,3000)
min_node_n=5
executable_dir="./proline_classifier/"
master_condor_file="/Users/longxiyao/condor_script_master.sh"
condor_script_dir="/Users/longxiyao/Google\ Drive/2018_spring/lab_work_data/machine_learning_cdr/proline_classifier/condor_script"
setwd(condor_script_dir)
for(loop in names(data_by_loop_type_list_unduplicated)){
  
  nl=dim(data_by_loop_type_list_unduplicated[[loop]][[1]])[1]
  eta=max (0.01, 0.1*min(1, nl/10000))
  args=c( loop, cores, complexity_n, the_tree_n, eta, min_node_n)
  


  
  executable=paste(args, collapse="_")
  
  exe_sh_file=paste(c(executable,"exe.sh"),collapse="_")  # customize shell script name
  
  Rscript_command_line=paste(c("Rscript 2.gbm_train_test_splitted_grid.R ",paste(args,collapse=" ")),collapse=" ")    # write shell script
  write(paste(c("cd ",Rscript_dir),collapse=" "), file = exe_sh_file )
  write(Rscript_command_line, file = exe_sh_file,append=TRUE )
  
  
  if(!condor_submit){
    system(paste(c("nohup  sh ",exe_sh_file," &"),collapse=""))   # run stuff in parallel
  }else{
    # write condor script
    con_file=paste(c(executable,".con"),collapse="_")  #customize condor script name
    system(paste(c("cp ", master_condor_file," ", con_file),collapse=""))   # copy the master file to script directory
    
    
    
    x <- readLines(con_file)
    y <- gsub( "\\$executable", exe_sh_file, x )
    z <- gsub( "\\$condor_script_dir", condor_script_dir, y )
    w <- gsub( "\\$cores", cores, z )
    cat(w, file=con_file, sep="\n")
    
    
    
    print(paste(c( condor_script_dir,"/",con_file),collapse=""))
    a=paste(c("condor_submit ", condor_script_dir,"/",con_file),collapse="")
    print(a)
    system(a)
  }
  
  
  
}


