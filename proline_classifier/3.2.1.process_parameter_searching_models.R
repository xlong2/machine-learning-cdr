# Purpose
#   This script harvest all the machine learning models for searching the grids, and find the best tuned grids out of all
#   the models.


current_d = getwd()
if (grepl("proline_classifier", current_d)) {
  source("0.load_function_and_data.R")
  
}

all_models = list.files(pattern = "*extra_test.rds",
                        path = "./rmsd_cluster_hits_rmsd",
                        full.names = FALSE)

all_models_list_by_loop = list()
index_i = 1
all_t = list()
for (each_file in all_models) {
  tryCatch({
    each_file_r = paste(c("./rmsd_cluster_hits_rmsd/", each_file), collapse =
                          "")
    the_model = readRDS(each_file_r)
    
    index_i = index_i + 1
    print(each_file)
    each_file = gsub("xgboost-north", "gbm_test-north", each_file)
    info = strsplit(each_file, "\\/")[[1]][length(strsplit(each_file, "\\/")[[1]])]
    loop_type = paste(strsplit(info, "_")[[1]][1:2], collapse = "_")
    all_t[[loop_type]] = unique(the_model$trainingData$.outcome)
    number = paste(strsplit(strsplit(info, "_")[[1]][4], "-")[[1]][3:6], collapse =
                     "-")
    model_result = the_model$result
    if (dim(model_result)[2] == 32) {
      model_result = model_result[, !names(model_result) %in% c( "prAUC", "Mean_Precision", "Mean_Recall", "prAUCSD",  "Mean_PrecisionSD",  "Mean_RecallSD" )]
    } else if (dim(model_result)[2] == 29) {
      names(model_result)
      remained_names = c("max_depth",
                         "nrounds",
                         "eta",
                         "min_child_weight",
                         names(model_result)[8:length(names(model_result))])
      model_result = model_result[, remained_names]
    }
    

    if(dim(model_result)[2]!=26){stop}
    all_models_list_by_loop[[loop_type]][[number]] = model_result
    rm(the_model)
    gc()
  }, error = function(e) {
    print(e)
  })
}





all_names=names(((all_models_list_by_loop)[[3]][[1]]))
all_models_list_by_loop = lapply(all_models_list_by_loop, function(x) {
  z=lapply(x,function(y){
    y=as.data.frame(y)
    names(y)=all_names
    return(y)})
  return(z)
  
})
all_result = lapply(all_models_list_by_loop, function(x) {
  z=x
  do.call(rbind, z)
  })

save_file("all_models_list_by_loop")


all_result = lapply(all_result, function(x) {
  if (!"Mean_Balanced_Accuracy" %in% names(x)) {
    names(x)[which(names(x) == "Balanced_Accuracy")] = "Mean_Balanced_Accuracy"
  }
  return(x)
})
save_file("all_result")


results = do.call(rbind, lapply(names(all_result), function(x) {
  all_result[[x]]$loop_type = rep(x, dim(all_result[[x]])[1])
  return(all_result[[x]])
}))
#ggplot(results)
results = results[!grepl("NA", rownames(results)), ]
saveRDS(results, file = "./Data_processed/results.rds")
gc()


results = readRDS("./Data_processed/result.rds")
all_results = results
results = all_results[all_results$shrinkage > 0.001, ]
results_smaller = all_results[all_results$shrinkage <= 0.001, ]

overall_accuracy=readRDS("./Data_processed/overall_accuracy.rds")
blind_blast_accuracy=as.data.frame(do.call(rbind,overall_accuracy))
blind_blast_accuracy$loop=rownames(blind_blast_accuracy)
colnames(blind_blast_accuracy)=c("mean","sd","loop_type")
blind_blast_accuracy=blind_blast_accuracy[!rownames(blind_blast_accuracy) %in% c("H1_14","H1_15","H1_10","L1_10","L2_12","L3_12"),]


ggplot(
  results,
  aes(
    x = n.trees,
    y = Accuracy,
    group = interaction.depth,
    color = interaction.depth,
    ymin = Accuracy - AccuracySD / 2,
    ymax = Accuracy + AccuracySD / 2
  )
) + geom_errorbar(width = 0.9)  + geom_line() + facet_wrap( ~ loop_type, scales =
                                                              "free") + ggtitle("Gradient Boost Machine model complexity tuning") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(data = blind_blast_accuracy, aes(yintercept = mean))
file_name = "/Users/xlong3/lab_work_data/machine_learning_cdr/proline_classifier/gbm_grid_search.png"

ggsave(
  file = file_name,
  width = 300,
  height = 200,
  units = c("mm")
)

print(file_name)
results_smaller_5 = results_smaller[results_smaller$n.minobsinnode == 5, ]
results_smaller_2 = results_smaller[results_smaller$n.minobsinnode == 2, ]

ggplot(
  results_smaller_2,
  aes(
    x = n.trees,
    y = Accuracy,
    group = interaction.depth,
    color = interaction.depth,
    ymin = Accuracy - AccuracySD / 2,
    ymax = Accuracy + AccuracySD / 2
  )
) + geom_errorbar(width = 0.9)  + geom_line() + facet_wrap( ~ loop_type, scales =
                                                              "free") + ggtitle("Gradient Boost Machine model complexity tuning") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(data = blind_blast_accuracy, aes(yintercept = mean))

ggplot(
  results_smaller_5,
  aes(
    x = n.trees,
    y = Accuracy,
    group = interaction.depth,
    color = interaction.depth,
    ymin = Accuracy - AccuracySD / 2,
    ymax = Accuracy + AccuracySD / 2
  )
) + geom_errorbar(width = 0.9)  + geom_line() + facet_wrap( ~ loop_type, scales ="free") + 
  ggtitle("Gradient Boost Machine model complexity tuning") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(data = blind_blast_accuracy, aes(yintercept = mean))



