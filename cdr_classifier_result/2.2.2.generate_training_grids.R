# this script generate machine learning model parameters and launch jobs for training models in parallel


current_d = getwd()
if (!grepl("cdr_classifier_result", current_d)) {
  setwd("./cdr_classifier_result")
  
} else{
  
}
source('./0.load_function_and_data.R')

Rscript_dir = "./training_script/"
system(paste(c("mkdir ",Rscript_dir),collapse = " "))
executable_dir = "./"


# parameters passed to Rscript for gbm model training
cores = 1
complexity = c(3, 6, 9)

trees = c(100, 300, 500, 1000, 1500, 2000, 2500, 3000)
trees = c(100)

min_node_n = 5

Rscript_name = "2.2.1.gbm_train_test_splitted_grid.R"
master_condor_file = "/Users/longxiyao/condor_script_master.sh"
condor_submit = FALSE
a_new_script = TRUE
parallel_script = "grid_training_parallel.sh"
condor_script_dir = "/Users/longxiyao/Google\ Drive/2018_spring/lab_work_data/machine_learning_cdr/cdr_classifier_result/condor_script"
count = 0
max_core = 1
for (loop in names(data_by_loop_type_list_unduplicated)) {
  eta = generate_eta(loop)   # generate learning rate
  argument_table = data.frame(expand.grid(loop, cores, complexity, trees, eta, min_node_n))
  
  for (index in 1:dim(argument_table)[1]) {   # iterate all parameter combination
    args = argument_table[index, ]
    
    exe_sh_file = generate_Rscript_command(args, Rscript_dir, Rscript_name)
    print(count)
    count = parallelize_jobs(max_core,
                             count,
                             a_new_script,
                             exe_sh_file,
                             parallel_script)
    a_new_script = FALSE
    
    
    
    if (!condor_submit) {
      
    } else{
      # write condor script
      setwd(condor_script_dir)
      
      a = generate_condor_script(args,
                                 exe_sh_file,
                                 condor_script_dir,
                                 master_condor_file)
      
      system(a)
    }
    
    
  }
  
}
print(parallel_script)
system(paste(c("  sh ", parallel_script), collapse = ""))   # run stuff in parallel



