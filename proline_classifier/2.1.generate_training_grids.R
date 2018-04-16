best_parameters_each_loop[["H1_13"]]
data_by_loop_type_list_unduplicated
core=4;
complexity=c(3,6,9);
trees=c(100,300,500,1000,1500, 2000,2500,3000)
min_node_n=5
for(loop in names(data_by_loop_type_list_unduplicated)){
  
  nl=dim(data_by_loop_type_list_unduplicated[[loop]][[1]])[1]
  eta=max (0.01, 0.1*min(1, nl/10000))
  args=c( loop, core, complexity_n, the_tree_n, eta, min_node_n)
  execute_training_rscript(args)
  
  
}


