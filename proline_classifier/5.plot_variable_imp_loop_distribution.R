# This script plot the feature importance extracted from GBM models 

current_d=getwd()
if(grepl("proline_classifier",current_d)){
  source("0.load_function_and_data.R")
  
}



pattern = paste(c(
  "L3_9_gbm_test-north-",
  paste(unlist(best_parameters_each_loop[["L3_9"]][, c("interaction.depth", "n.trees", "shrinkage", "n.minobsinnode")]), collapse =
          "-"),
  "_trained_model_extra_test.rds"
),
collapse = "")
path = "./rmsd_cluster_hits_rmsd/"
model = readRDS(file = paste(c(path, pattern), collapse = ""))

gbmImp <- varImp(model, scale = TRUE)
#
var_import = as.data.frame(gbmImp[[1]])
var_import$var = rownames(var_import)
var_import = var_import[order(gbmImp[[1]][, 1], decreasing = TRUE), ]
var_import = var_import[1:20, ]
var_import = var_import[order(var_import[, 1], decreasing = FALSE), ]
var_import$var = factor(var_import$var, levels = var_import$var)
ggplot(var_import, aes(x = var, y = Overall, fill = "red")) + geom_bar(stat =
                                                                         "identity") + coord_flip() + ggtitle("Residue importance for classifying L3_9 using gbm") +
  theme(plot.title = element_text(hjust = 0.5))
gbmImp <- varImp(model, scale = TRUE)

var_import$var

