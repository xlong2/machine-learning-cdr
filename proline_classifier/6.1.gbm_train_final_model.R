# script for training the final model after finding the best hyperparameter. 
args <- commandArgs(trailingOnly = TRUE)
print(args)
getwd()
source("0.load_function_and_data.R")
execute_training_rscript_final_model(args)
