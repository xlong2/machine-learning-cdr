library("e1071")
core=4;


for(loop in c("H1_13", "H2_10", "L3_9")){
nl=dim(data_by_loop_type_list_unduplicated[[loop]][[1]])[1]
eta=max (0.01, 0.1*min(1, nl/10000))
paras=best_parameters_each_loop[[loop]]
args=c(loop,core,paras[["interaction.depth"]],paras[["n.trees"]], eta,5)
execute_training_rscript(args)
}