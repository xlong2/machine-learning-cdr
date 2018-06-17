#This script get the best parameter 

current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

best_parameters_each_loop = list()
all_gird_search_results = list()
result_table = as.data.frame(matrix(nrow = length(names(
  all_models_list_by_loop
)), ncol = length(
  c(
    "branch_level",
    "trees",
    "shrinkage",
    "min_node",
    "accuracy",
    "accuracy_sd"
  )
)))
colnames(result_table) = c("branch_level",
                           "trees",
                           "shrinkage",
                           "min_node",
                           "accuracy",
                           "accuracy_sd")
rownames(result_table) = names(all_result)
for (each_loop in names(all_result)) {
  all_gird_search_results[[each_loop]] = all_result[[each_loop]]
  
  # do a plot
  x = unique(all_gird_search_results[[each_loop]]$n.trees)
  y = sort(unique(all_gird_search_results[[each_loop]]$interaction.depth))
  z = matrix(nrow = length(x), ncol = length(y))
  z = as.data.frame(z)
  
  rownames(z) = x
  colnames(z) = y
  # only select the mino equal to 5
  sub_table = all_gird_search_results[[each_loop]][all_gird_search_results[[each_loop]]$n.minobsinnode ==
                                                     5, ]
  for (ind in 1:dim(sub_table)[1]) {
    z[as.character(sub_table[ind, "n.trees"]), as.character(sub_table[ind, "interaction.depth"])] =
      sub_table[ind, "Accuracy"]
  }
  
  
  max_indses = which(all_gird_search_results[[each_loop]]$Accuracy == max(all_gird_search_results[[each_loop]]$Accuracy))
  max_indses_index = which(all_gird_search_results[[each_loop]][max_indses, "Kappa"] ==
                             max(all_gird_search_results[[each_loop]][max_indses, "Kappa"]))
  final_max_ind = max_indses[max_indses_index]
  if (length(final_max_ind) != 1) {
    which_i = which(all_gird_search_results[[each_loop]][final_max_ind, "n.trees"] ==
                      min(all_gird_search_results[[each_loop]][final_max_ind, "n.trees"]))
    final_max_ind = final_max_ind[which_i]
    if (length(final_max_ind) != 1) {
      which_i = which(all_gird_search_results[[each_loop]][final_max_ind, "interaction.depth"] ==
                        min(all_gird_search_results[[each_loop]][final_max_ind, "interaction.depth"]))
      final_max_ind = final_max_ind[which_i]
      
    }
    
  }
  
  parameters = rownames(all_gird_search_results[[each_loop]])[final_max_ind]
  parameters = t(data.frame(strsplit(as.character(parameters), "-")))
  
  result_table[each_loop, c("branch_level", "trees", "shrinkage", "min_node")] =
    parameters
  result_table[each_loop, c("accuracy", "accuracy_sd")] = all_gird_search_results[[each_loop]][final_max_ind, c("Accuracy", "AccuracySD")]
  best_parameters_each_loop[[each_loop]] = parameters
  # slice3D(all_gird_search_results[[each_loop]]$n.trees, all_gird_search_results[[each_loop]]$shrinkage, all_gird_search_results[[each_loop]]$interaction.depth,colvar=all_gird_search_results[[each_loop]]$Accuracy)
  #  M <- mesh(all_gird_search_results[[each_loop]]$n.trees, all_gird_search_results[[each_loop]]$shrinkage, all_gird_search_results[[each_loop]]$interaction.depth)
}

rownames(result_table) = split_vector_and_replace(rownames(result_table), "_", 1, 2, "-")
result_table = result_table[order(order_factor_by_two_component(rownames(result_table), "-", 1, 2),  rownames(result_table)), ]

print(result_table)
write.csv(
  result_table,
  file = paste(result_dir, "best_parameters.csv"),
  row.names = TRUE,
  col.names  = TRUE
)
best_parameters_each_loop=lapply(best_parameters_each_loop,function(x){
  print(x)
  x=as.data.frame(x)
  x=lapply(x,as.character)
  names(x)=c("interaction.depth","n.trees","eta","n.minobsinnode")
  x=as.data.frame(x)
  return(x)
})
best_parameters_each_loop[["L3_9"]][,1]=9
best_parameters_each_loop[["L3_9"]][,2]=1500

save_file("best_parameters_each_loop")
# print the err

