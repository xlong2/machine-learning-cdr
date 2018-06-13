library("caret")
library(parallel)
library(doMC)
args <- commandArgs(trailingOnly = TRUE)
split_index=args[1]
num_core=3
loop_type=args[2]
test=FALSE
test_jazz=FALSE
print(paste(c(loop_type,"splitted_grid.rds"),collapse="_"))
splitted_grid=readRDS(paste(c(loop_type,"splitted_grid.rds"),collapse="_"))
print(splitted_grid)
gbmGrid=splitted_grid[[split_index]]
each_method="gbm_test"
cluster_dis="north"

registerDoMC(num_core)

if(test){
  overall_prefix="/Volumes/lab/macbook/lab_work_data/"#"/Volumes/xlong/"
  
  load(paste(c(overall_prefix,"DTW_length_independent/data_by_loop_type_list_unduplicated.Rdata"),collapse = ""))
  result_dir=paste(c(overall_prefix,"/machine_learning_cdr/proline_classifier/"),collapse="")
  
}else{
  overall_prefix="/home/xlong/"
  
  load("/home/xlong/DTW_length_independent/data_by_loop_type_list_unduplicated.Rdata")
  result_dir=paste(c(overall_prefix,"machine_learning_cdr/proline_classifier/rmsd_cluster_hits_rmsd/"),collapse="")
}
data=data_by_loop_type_list_unduplicated[[loop_type]][[1]]
#data=data[data[,"V8"]=="P",]
#data$cluster_type=as.character(data$cluster_type)
#data[grepl("cis",data$cluster_type),"cluster_type"]=rep("cis",length(data[grepl("cis",data$cluster_type),"cluster_type"]))
#data[!grepl("cis",data$cluster_type),"cluster_type"]=rep("noncis",length(data[!grepl("cis",data$cluster_type),"cluster_type"]))
data$cluster_type=sub("-","_",data$cluster_type)
data$cluster_type=as.factor(as.character(data$cluster_type))





subsitution_matrix_name ="wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
subsitution_matrix="PAM30"


even_out_all_classes<-function(training_cases){
  
  
  splitted_class=split(training_cases,training_cases$cluster_type)
  proportions=unlist(lapply(split(training_cases,training_cases$cluster_type),function(x){dim(x)[1]}))
  max_class=names(proportions)[which(proportions==max(proportions))]
  remaining_class=names(proportions)[names(proportions)!=max_class]
  for(each_remaining in remaining_class){
    ratio=proportions[max_class]%/%proportions[each_remaining]
    if(ratio>1){
      ratio=ratio-1
    }
    to_be_added=do.call(rbind,replicate(ratio,splitted_class[[each_remaining]],simplify = FALSE))
    training_cases=rbind(training_cases,to_be_added)    
  }
  return(training_cases)
}



generic_train<-function(each_loop_length_data_feature_string_rmsd,each_method,training_cases){
  
  # to add classes to off set 
  training_cases=even_out_all_classes(training_cases)
  
  getModelInfo()$gbm$parameters
  # Max shrinkage for gbm
  nl = nrow(training_cases)
  max(0.01, 0.1*min(1, nl/10000))
  # Max Value for interaction.depth
  floor(sqrt(NCOL(training_cases)))
  single_grid <-  expand.grid(interaction.depth =c(6),#*c(6), #c(1, 3, 6, 9, 10),
                              n.trees =1000,#3000, #(5:60)*100, 
                              shrinkage =c(0.001), #c(0.001),#seq(.0005, .02,.0005),
                              n.minobsinnode = c(5)) # you can also put something   
  single_grid <-  expand.grid(interaction.depth =c(6),#*c(6), #c(1, 3, 6, 9, 10),
                              n.trees =3000,#3000, #(5:60)*100, 
                              shrinkage =c(0.01), #c(0.001),#seq(.0005, .02,.0005),
                              n.minobsinnode = c(2)) # you can also put something   
  
  single_grid_test <-  expand.grid(interaction.depth =c(2),#*c(6), #c(1, 3, 6, 9, 10),
                                   n.trees =100,#3000, #(5:60)*100,
                                   shrinkage =c(0.01), #c(0.001),#seq(.0005, .02,.0005),
                                   n.minobsinnode = c(2))

  fitControl <- trainControl(method = "repeatedcv",
                             repeats = 3,
                             preProcOptions = list(thresh = 0.95),
                             ## Estimate class probabilities
                             classProbs = TRUE,
                             ## Evaluate performance using
                             ## the following function
                             savePredictions="final",
                             summaryFunction = multiClassSummary)
  
  
  fitcontrol_test <- trainControl(method = "repeatedcv",
                                  repeats = 1,number=2,
                                  preProcOptions = list(thresh = 0.95),
                                  ## Estimate class probabilities
                                  classProbs = TRUE,
                                  ## Evaluate performance using
                                  ## the following function
                                  savePredictions="final",
                                  summaryFunction = multiClassSummary)
  # Method + Date + distribution
  set.seed(1)
  print(training_cases$weights)
  weights=training_cases$weights
  if(test_jazz){
    time= system.time(trained_model <- train(each_loop_length_data_feature_string_rmsd, data = training_cases,
                                             #distribution = "multinomial",
                                             method = "gbm", bag.fraction = 0.5,   # fold number 10
                                             #nTrain = round(nrow(training_cases) *.75),
                                             trControl = fitcontrol_test,
                                             verbose = TRUE,
                                             tuneGrid = single_grid_test,
                                             #savePredictions=final,
                                             ## Specify which metric to optimize
                                             metric = "kappa"))
    print(time)
  }else{
    time= system.time(trained_model <- train(each_loop_length_data_feature_string_rmsd, data = training_cases,
                                             #distribution = "adaboost",
                                             method = "gbm", bag.fraction = 0.5,   # fold number 10 
                                             #nTrain = round(nrow(training_cases) *.75),
                                             trControl = fitControl,
                                             verbose = TRUE,
                                             tuneGrid = gbmGrid,
                                             #savePredictions=final,
                                             ## Specify which metric to optimize
                                             metric = "kappa"))
  }
  print(time)
  
  return(trained_model)
}





#sequences$rmsd_cluster = as.character(sequences$rmsd_cluster)
sequences=data

each_loop_length_data_feature_string = data_by_loop_type_list_unduplicated[[loop_type]][[2]]
features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
features=features
each_loop_length_data_feature_string=as.formula(paste(c("cluster_type ~ ",paste(features,collapse=" + ")),collapse=""))
#if(cluster_dis=="north"){
each_loop_length_data_feature_string_rmsd = each_loop_length_data_feature_string
all_cases =  sequences[,c(features,"identifier","id","cluster_type")]
all_cases=all_cases[complete.cases(all_cases), ]
all_cases$cluster_type=gsub("-","_",all_cases$cluster_type)
all_cases$cluster_type=as.factor(as.character(all_cases$cluster_type))

# tune the parameter
trained_model = generic_train(each_loop_length_data_feature_string_rmsd,each_method,all_cases)  # do LOOCV using the selected machine learning method


pred_result = trained_model$pred # get the prediction vs observation
#pretty.gbm.tree(trained_model)
#pred_result_table= data.frame(obs = pred_result$obs, pred = pred_result$pred,identifier = all_cases$identifier,id = all_cases$id)
conf_table = table(all_cases$cluster_type, predict(trained_model,all_cases)) # get the confusion table and save it

print("finished the cluster prediction")
print(conf_table)
check_index<-function(split_index){
old_files=list.files(pattern="*extra_test.rds",path="/home/xlong/machine_learning_cdr/proline_classifier/rmsd_cluster_hits_rmsd",full.names = FALSE)
print(c(old_files,"line 176 "))
old_files=grep(loop_type,old_files,value=TRUE)
print(old_files)
old_indexes=sapply(old_files,function(x){as.numeric(strsplit(strsplit(as.character(x),"-")[[1]][3],"_")[[1]][1])})
print(c(old_indexes,"OLD_INDEXES "))
split_index=as.numeric(split_index)
if(split_index %in% old_indexes){ split_index=split_index+1; check_index(split_index)}else{
print("saved file")
print(paste(c(result_dir,"",loop_type,"_",paste(c(each_method,cluster_dis,split_index),collapse="-"),"_trained_model_extra_test.rds"),collapse=""))
saveRDS(trained_model,file =paste(c(result_dir,"",loop_type,"_",paste(c(each_method,cluster_dis,split_index),collapse="-"),"_trained_model_extra_test.rds"),collapse=""))}  }
check_index(split_index)
