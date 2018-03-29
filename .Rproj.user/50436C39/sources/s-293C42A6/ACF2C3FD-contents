library("caret")
print("line 2")
library(parallel)
library(doMC)
print("line 5")
args <- commandArgs(trailingOnly = TRUE)
#args=c(2, "L3_9", 4, 6, 1500, 0.01, 5)
split_index=args[1]
num_core=args[3]
loop_type=args[2]
interaction.depth=args[4]
n.trees=args[5]
shrinkage=args[6]
n.minobsinnode=args[7]
print("line 15")
test=FALSE
test_jazz= FALSE
#gbmGrid=splitted_grid[[split_index]]
gbmGrid=data.frame(interaction.depth= as.numeric(interaction.depth) , n.trees = as.numeric(n.trees), shrinkage= as.numeric(shrinkage), n.minobsinnode = as.numeric(n.minobsinnode))
print("this is gbmGrid")
print(gbmGrid)
print("line 22")
parameter_spe = paste(unlist(gbmGrid),collapse="-")
each_method="gbm_test"
cluster_dis="north"

registerDoMC(num_core)

if(test){
  setwd("/Users/xlong3/lab_work_data/machine_learning_cdr/proline_classifier")
  
 # overall_prefix="/Volumes/lab/macbook/lab_work_data/"#"/Volumes/xlong/"
  overall_prefix="/Users/xlong3/lab_work_data/machine_learning_cdr/proline_classifier"
  data_by_loop_type_list_unduplicated=readRDS(paste(c(overall_prefix,"/data_by_loop_type_list_unduplicated_longer_features.rds"),collapse = ""))
  lapply(data_by_loop_type_list_unduplicated,function(x){dim(x[[1]])})
  data_by_loop_type_list_unduplicated=readRDS(paste(c(overall_prefix,"/data_by_loop_type_list_unduplicated_no_filtering.rds"),collapse = ""))
  overall_prefix="/Users/xlong3/lab_work_data"
    result_dir=paste(c(overall_prefix,"/machine_learning_cdr/proline_classifier/"),collapse="")
  
}else{
  overall_prefix="/home/xlong/"
  
  load("/home/xlong/DTW_length_independent/data_by_loop_type_list_unduplicated.Rdata")
  data_by_loop_type_list_unduplicated=readRDS(paste(c(overall_prefix,"machine_learning_cdr/proline_classifier/data_by_loop_type_list_unduplicated_no_filtering.rds"),collapse = ""))
  
  result_dir=paste(c(overall_prefix,"machine_learning_cdr/proline_classifier/rmsd_cluster_hits_rmsd/"),collapse="")
}

to_save_file=paste(c(result_dir,"",loop_type,"_",paste(c(each_method,cluster_dis,parameter_spe),collapse="-"),"_trained_model_extra_test.rds"),collapse="")
print(to_save_file)
if(file.exists(to_save_file)){
print("file_already_exists!")}else{
data=data_by_loop_type_list_unduplicated[[loop_type]][[1]]
#data=data[data[,"V8"]=="P",]
#data$cluster_type=as.character(data$cluster_type)
#data[grepl("cis",data$cluster_type),"cluster_type"]=rep("cis",length(data[grepl("cis",data$cluster_type),"cluster_type"]))
#data[!grepl("cis",data$cluster_type),"cluster_type"]=rep("noncis",length(data[!grepl("cis",data$cluster_type),"cluster_type"]))
data$cluster_type=sub("-","_",data$cluster_type)
data$cluster_type=as.factor(as.character(data$cluster_type))





subsitution_matrix_name ="wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
subsitution_matrix="PAM30"


even_out_all_classes_non_cases_no_repeat<-function(training_cases){
  
  non_cases=training_cases[grepl("none",training_cases$cluster_type),]
  training_cases=training_cases[!grepl("none",training_cases$cluster_type),]
  training_cases$cluster_type=as.factor(as.character(training_cases$cluster_type))
  splitted_class=split(training_cases,training_cases$cluster_type)
  proportions=unlist(lapply(split(training_cases,training_cases$cluster_type),function(x){dim(x)[1]}))
  max_class=names(proportions)[which(proportions==max(proportions))]
  remaining_class=names(proportions)[names(proportions)!=max_class]
  for(each_remaining in remaining_class){
    ratio=proportions[max_class]%/%proportions[each_remaining]
    if(ratio>1){
      ratio=ratio-1
    }
    temp=splitted_class[[each_remaining]]
    to_be_added=temp[rep(seq_len(nrow(training_cases)), ratio), ] 
    #to_be_added=do.call(rbind,replicate(ratio,as.matrix(splitted_class[[each_remaining]]),simplify = FALSE))
    
    training_cases=rbind(training_cases,to_be_added)    
  }
  training_cases=rbind(training_cases,non_cases)
  training_cases$cluster_type=as.factor(as.character(training_cases$cluster_type))
  return(training_cases)
}



even_out_all_classes<-function(training_cases){
  
  training_cases$cluster_type=as.factor(as.character(training_cases$cluster_type))
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





make_3_10_cross_val<-function(training_cases,r,k){

  all_unique_ids=lapply(split(training_cases,training_cases$cluster_type),function(x){unique(x$id)})
  
  split_data=split(training_cases,training_cases$cluster_type)
  new_training_data=do.call(rbind,lapply(names(split_data),function(x){
    for(repeat_n in 1:r){
    all_unique_id_sample=lapply(all_unique_ids,function(y){if(length(y)>=k){sample(1:k,size=length(y),replace=TRUE)}else{sample(1:k,size=length(y),replace=FALSE)}})
    
    for(ind in 1:length(all_unique_ids[[x]])){
    all_unique_ids[[x]]
      repnum=dim(split_data[[x]][split_data[[x]]$id==all_unique_ids[[x]][ind],])[1]
    split_data[[x]][split_data[[x]]$id==all_unique_ids[[x]][ind],paste(c("fold.num",repeat_n),collapse="")]=rep(all_unique_id_sample[[x]][ind],repnum)
    }
    }
    return(split_data[[x]])
    }))

  
# Create folds and repeats here - you could create your own if you want #


folds.list.out <- list()
folds.list <- list()
list.counter <- 1
for (y in 1:r) {
  newcol <- paste('fold.num', y, sep='')
  for (z in 1:k) {
    out_rown= which(new_training_data[,newcol]==z)
folds_in =which(new_training_data[,newcol]!=z)
    a=new_training_data[folds_in,"id"]%in%  new_training_data[out_rown,"id"]
    print(a[a])

    sub=new_training_data[out_rown,"id"]
    out_rown=out_rown[which(!duplicated(sub))]
    a=new_training_data[folds_in,"id"]%in%  new_training_data[out_rown,"id"]
    print(a[a])
    if(length(a[a])!=0){
      print("The 1fold in and 9 fold out is not correct, check!")
      stop()
    }
    folds.list.out[[list.counter]] <- out_rown
    folds.list[[list.counter]] <- which(new_training_data[,newcol]!=z)
    list.counter <- list.counter + 1
  }
}
index_list=list()
index_list[[1]]=folds.list.out
index_list[[2]]=folds.list
returned_list=list()
returned_list[[1]]=index_list;
returned_list[[2]]=new_training_data
return(returned_list)

}


generic_train<-function(each_loop_length_data_feature_string_rmsd,each_method,training_cases){
   # to add classes to off set 
  training_cases=even_out_all_classes(training_cases)
  if(dim(training_cases)[1]*0.6*0.5 < as.numeric(as.character(gbmGrid[["n.minobsinnode"]]))){
    training_cases=rbind(training_cases,training_cases)
    training_cases=rbind(training_cases,training_cases)

  }
  
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


  

  
  # Method + Date + distribution
  set.seed(1)
  print(training_cases$weights)
  weights=training_cases$weights
  if(test_jazz){
    r <- 1 # number of repeats
    k <- 2 # number of folds
    returned_results=make_3_10_cross_val(training_cases,r,k)
    folds_spec=returned_results[[1]]
    training_cases=returned_results[[2]]
    folds.list.out=folds_spec[[1]]
    folds.list=folds_spec[[2]]

    fitcontrol_test <- trainControl(method = "repeatedcv",
                                    repeats = 1,number=2,
                                    preProcOptions = list(thresh = 0.95),
                                    ## Estimate class probabilities
                                    index=folds.list
                                    , indexOut=folds.list.out,
                                    classProbs = TRUE,
                                    ## Evaluate performance using
                                    ## the following function
                                    savePredictions="final",
                                    summaryFunction = multiClassSummary)
    
    fitcontrol_test_lean <- trainControl(method = "repeatedcv",
                                    repeats = 1,number=2,
                                    preProcOptions = list(thresh = 0.95),
                                    ## Estimate class probabilities

                                    classProbs = TRUE,
                                    ## Evaluate performance using
                                    ## the following function
                                    savePredictions="final",
                                    summaryFunction = multiClassSummary)
    
   trained_model <- train(each_loop_length_data_feature_string_rmsd, data = training_cases,
                                             #distribution = "multinomial",
                                             method = "gbm", bag.fraction = 0.5,   # fold number 10
                                             #nTrain = round(nrow(training_cases) *.75),
                                             trControl = fitcontrol_test,
                                             tuneGrid = single_grid_test,
                                              verbose=TRUE,
                                             ## Specify which metric to optimize
                                             metric = "kappa")
  
  }else{
    if(loop_type %in% c("H1_13", "H2_10", "L3_9" ) && n.trees>=1000 ){
       r=1}else{r=3}
    #r <- 3 # number of repeats
    k <- 10 # number of folds
    returned_results=make_3_10_cross_val(training_cases,r,k)
    folds_spec=returned_results[[1]]
    training_cases=returned_results[[2]]   # the training cases would have its mo
    folds.list.out=folds_spec[[1]]
    folds.list=folds_spec[[2]]

    fitControl <- trainControl(method = "repeatedcv",
                               repeats = r,
                               number=k,
                               preProcOptions = list(thresh = 0.95),
                               index=folds.list
                               , indexOut=folds.list.out,
                               ## Estimate class probabilities
                               classProbs = TRUE,
                               returnResamp="all",
                               ## Evaluate performance using
                               ## the following function
                               savePredictions="all",
                               summaryFunction = multiClassSummary)
    trained_model=""
   trained_model <- train(each_loop_length_data_feature_string_rmsd, data = training_cases,
                                             #distribution = "adaboost",
                                             method = "gbm", bag.fraction = 0.5,   # fold number 10 
                                             #nTrain = round(nrow(training_cases) *.75),
                                             trControl = fitControl,
                                             verbose = TRUE,
                                             tuneGrid = gbmGrid,

                                             ## Specify which metric to optimize
                                             metric = "kappa")
  }

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
all_cases =  sequences[,c(features,"id","cluster_type")]
H1_13features=data_by_loop_type_list_unduplicated[["H1_13"]][[4]]
the_levels=unique(unlist(data_by_loop_type_list_unduplicated[["H1_13"]][[1]][,H1_13features]))
for(each_f in features){
  all_cases[,each_f]=factor(all_cases[,each_f],levels=the_levels)
}
all_cases=all_cases[complete.cases(all_cases), ]
all_cases$cluster_type=gsub("-","_",all_cases$cluster_type)
all_cases$cluster_type=gsub(",",".",all_cases$cluster_type)

all_cases$cluster_type=as.factor(as.character(all_cases$cluster_type))

# tune the parameter
trained_model = generic_train(each_loop_length_data_feature_string_rmsd,each_method,all_cases)  # do LOOCV using the selected machine learning method


pred_result = trained_model$pred # get the prediction vs observation
#pretty.gbm.tree(trained_model)
#pred_result_table= data.frame(obs = pred_result$obs, pred = pred_result$pred,identifier = all_cases$identifier,id = all_cases$id)
conf_table = table(all_cases$cluster_type, predict(trained_model,all_cases)) # get the confusion table and save it

print("finished the cluster prediction")
print(conf_table)

saveRDS(trained_model,file =to_save_file)
print(to_save_file)

check_index<-function(split_index){
old_files=list.files(pattern="*extra_test.rds",path="/home/xlong/machine_learning_cdr/proline_classifier/rmsd_cluster_hits_rmsd",full.names = FALSE)
print(c(old_files,"line 176 "))
old_files=grep(loop_type,old_files,value=TRUE)
print(old_files)
old_indexes=sapply(old_files,function(x){as.numeric(strsplit(strsplit(as.character(x),"-")[[1]][3],"_")[[1]][1])})
print(c(old_indexes,"OLD_INDEXES "))
split_index=as.numeric(split_index)
if(split_index %in% old_indexes){ split_index=split_index+1; check_index(split_index)}else{
print(to_save_file)
saveRDS(trained_model,file =to_save_file)}  }



#check_index(split_index)

}
