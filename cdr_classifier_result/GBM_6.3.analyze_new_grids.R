path_files = list.files(path="./Data_processed/suggested_change/models_dir",full.names = TRUE)
path_files = path_files[grepl("trained_model_extra_test.rds", path_files)]

path_file_l = lapply(strsplit(path_files,"_"),length)
which_f = which(path_file_l==11)
for(which_f_name in which_f){
  file_n = path_file_l
new_name = paste(c(strsplit(path_files[which_f_name],"_")[[1]][1:7],"3",strsplit(path_files[which_f_name],"_")[[1]][8:11]),collapse="_")
system(paste(c("mv ", path_files[which_f_name], " ",new_name),collapse=" "))
}


path_files = list.files(path="./Data_processed/suggested_change/models_dir",full.names = TRUE)
path_files = path_files[grepl("trained_model_extra_test.rds", path_files)]

useful_info= lapply(strsplit(sapply(strsplit(path_files,"/"),"[[",length(strsplit(path_files[[1]],"/")[[1]])),"_"),function(x){x[c(1,2,4)]})
useful_info_parsed = lapply(useful_info,function(x){c(paste(c(x[1],x[2]),collapse="_"),  strsplit(x[3],"-")[[1]][3:7])    })
useful_info_parsed_data_frame = as.data.frame(cbind(path_files,as.data.frame(t(as.data.frame(useful_info_parsed,stringsAsFactors = FALSE)))),stringsAsFactors=FALSE)

splitte_result = split(useful_info_parsed_data_frame, apply(useful_info_parsed_data_frame,1, function(x){paste(x[c(2,7)],collapse="-")}))
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

validation_result_set = list()
for(each_set in names(splitte_result)){
  print(each_set)
  data_f = splitte_result[[each_set]]
  data_f$accuracy = rep(NA,dim(data_f)[1])
  data_f$sd = rep(NA,dim(data_f)[1])
  
  all_frames = list()
  for(ind in 1:dim(data_f)[1]){
    model_info = paste(unlist(data_f[ind,2:6]),collapse="-")
    print(as.character(unlist(data_f[ind,1])))
    model_result = readRDS(as.character(unlist(data_f[ind,1])))
    data_f[ind,"accuracy"] = model_result$results$Accuracy
    
    result = sd(unlist(lapply(chunk(unique(model_result$pred$Resample),3),function(x){
      rel = model_result$pred[model_result$pred$Resample%in% x, ]; 
      dim(rel[as.character(rel$obs)== as.character(rel$pred),])[1]/dim(rel)[1]
      
      })))
    data_f[ind,"sd"] = result
  }
  validation_result_set[[each_set]] = data_f
}



validation_result_set_split_repeat = lapply(validation_result_set,function(x){
  split(x,sapply(strsplit((as.character(x$path_files)),"_"),"[[",8))
})







max_f_new = lapply(validation_result_set_split_repeat,function(x){lapply(x,function(x){x[which(x[,"accuracy"]==max(x[,"accuracy"])),]})})
max_f = unlist(max_f_new,recursive = FALSE)
names(max_f) = gsub("\\.","-",names(max_f))
save_file("validation_result_set")
save_file("max_f")





testing_path_files = list.files(path="./Data_processed/suggested_change/testing_result",full.names = TRUE)
testing_path_files = testing_path_files[grepl("model_extra_test_result.tsv", testing_path_files)]

testing_path_files_l = lapply(strsplit(as.character(testing_path_files),"_"),length)
which_t_f = which(testing_path_files_l==12)
for(which_f_name in which_t_f){
  new_name = paste(c(strsplit(as.character(testing_path_files[which_f_name]),"_")[[1]][1:7],"3",strsplit(testing_path_files[which_f_name],"_")[[1]][8:12]),collapse="_")
  system(paste(c("mv ", testing_path_files[which_f_name], " ",new_name),collapse=" "))
}


testing_path_files = list.files(path="./Data_processed/suggested_change/testing_result",full.names = TRUE)
testing_path_files = testing_path_files[grepl("model_extra_test_result.tsv", testing_path_files)]





useful_info_testing= lapply(strsplit(sapply(strsplit(testing_path_files,"/"),"[[",length(strsplit(testing_path_files[[1]],"/")[[1]])),"_"),function(x){x[c(1,2,4,5)]})
useful_info_parsed_testing = lapply(useful_info_testing,function(x){c(paste(c(x[1],x[2]),collapse="_"),  strsplit(x[3],"-")[[1]][3:7],x[4] )    })
useful_info_parsed_testing_data_frame = as.data.frame(cbind(testing_path_files,as.data.frame(t(as.data.frame(useful_info_parsed_testing,stringsAsFactors = FALSE)))),stringsAsFactors=FALSE)


splitted_testing_result = split(useful_info_parsed_testing_data_frame, apply(useful_info_parsed_testing_data_frame,1, function(x){paste(x[c(2,7,8)],collapse="-")}))


testing_result_list = list()
for(each_set in names(max_f)){
  result = max_f[[each_set]]
  data_f_t = splitted_testing_result[[each_set]]
  testing_result = as.character( data_f_t[which(apply(data_f_t[,c("V2","V3","V4","V5")],1,function(x){ paste(x,collapse = "")})==paste(unlist(result[1,c(3:6)]),collapse = "")),1] )
  loop = strsplit(each_set,"-")[[1]][1]
  testing_ind = paste(strsplit(each_set,"-")[[1]][2:3],collapse="-")
  if(!loop%in% names(testing_result_list)){ testing_result_list[[loop]] = list()   }

  testing_result_list[[loop]][[testing_ind]] =  as.data.frame(table(as.character(read.table(testing_result)[,1]), as.character(read.table(testing_result)[,2])))
}
save_file("testing_result_list")

merge_c<-function(x,y){
  merge(x,y,by=c("Var1","Var2"),all=TRUE)
}

repeats_n =3
merged_testing_result = list()
gbm_per_loop_sd = list()
gbm_error_count_sd = list()
gbm_errorcount_list = list()
chunks_n_l = list()
all_records= list()
simply_merged = list()
for(loop in names(testing_result_list) ){
  data_l = testing_result_list[[loop]]
  chunks = split(names(data_l), sapply(strsplit(names(data_l),"-"),"[[",2))
  chunks_n_l[[loop]] = chunks
  splitted_d = lapply(chunks, function(x){data_l[x]})
  splitted_d_re = lapply(splitted_d, function(data_l){
    for(ind in 1:(length(data_l)-1)){
      if(ind==1){
        meg = data_l[[ind]]
      }
      print(meg)
      meg = merge_c(meg,data_l[[ind+1]])
    }
    meg$Freq = apply(meg[,3:dim(meg)[2]],1,function(x){sum(x,na.rm=TRUE)})
    meg = meg[,c("Var1","Var2","Freq")]
  })
  all_records[[loop]] = splitted_d_re
  loop_sd = sd(unlist( lapply(splitted_d_re,function(x){
      sum(x[as.character(x[,1])==as.character(x[,2]),3])/sum(x[,3])
    })))
  
  gbm_errorcount_list[[loop]] = unlist(lapply(splitted_d_re,function(x){
        sum(x[as.character(x[,1])!=as.character(x[,2]),3])
     }))
  gbm_per_loop_sd[[loop]] = loop_sd
  merged = merge(merge(splitted_d_re[[1]],splitted_d_re[[2]],by=c("Var1","Var2") ), splitted_d_re[[3]],by=c("Var1","Var2"))
  sd = apply(as.data.frame(merged),1,function(x){sd(unlist(x[3:5]))})
  average = apply(as.data.frame(merged),1,function(x){mean(as.numeric(unlist(x[3:5])),na.rm=TRUE)})
  gbm_error_count_sd[[loop]] = sd(unlist(lapply(merged[as.character(merged[,1])!=as.character(merged[,2]),3:5],function(x){sum(x)})))
  conf_t = merged[,c(1:2)]
  conf_t$Freq = average
  conf_t$sd = sd
  conf_t=conf_t[conf_t$Freq>0,]
  merged_testing_result[[loop]] = conf_t
  simply_merged[[loop]] = merged
  conf_t = conf_t[conf_t$Freq > 0 &
                    as.character(conf_t[, 1]) != as.character(conf_t[, 2]), ]
  
  
  


  
}
save_file("merged_testing_result")

accuracy_by_loop = lapply(merged_testing_result, function(x){
  sum(x[as.character(x$Var1)==as.character(x$Var2), "Freq"],na.rm = TRUE)/ sum(x[, "Freq"],na.rm = TRUE)
})

names = split(names(merged_testing_result),sapply(strsplit(names(merged_testing_result),"_"),"[[",1))
accuracy_by_l = lapply(names, function(j){x = do.call(rbind,merged_testing_result[j]);
sum(x[as.character(x$Var1)==as.character(x$Var2), "Freq"],na.rm = TRUE)/ sum(x[, "Freq"],na.rm = TRUE)

})

gbm_binded_by_loop = lapply(gbm_binded,function(x){
  split(x,sapply(strsplit(as.character(x[,1]),"_"),"[[",1))
})

save_file("accuracy_by_loop")
save_file("gbm_error_count_sd")
rbind_all_test_result = do.call(rbind, merged_testing_result)
save_file("rbind_all_test_result")
total_accu = sum(rbind_all_test_result[as.character(rbind_all_test_result$Var1)==as.character(rbind_all_test_result$Var2), "Freq"],na.rm = TRUE)/ sum(rbind_all_test_result[, "Freq"],na.rm = TRUE)



# calculating the accuracies for the method 
gbm_curated_by_loop = lapply(testing_result_list,function(check){
  split_check = split(check,sapply(strsplit(names(check),"-"),"[[",2))
  lapply(split_check,function(y){
    do.call(rbind,y)
  })
})
gbm_by_large_loop = split(gbm_curated_by_loop,sapply(strsplit(names(gbm_curated_by_loop),"_"),"[[",1))
gbm_by_large_loop_combine = lapply(gbm_by_large_loop,function(x){
  lapply(names(x[[1]]),function(y){
    do.call(rbind,lapply(x,function(g){g[[y]]}))
  })
})
gbm_by_large_loop_combine_result = lapply(gbm_by_large_loop_combine,function(x){lapply(x,function(y){
  
  sum(y[as.character(y[,1])==as.character(y[,2]),3])/sum(y[,3])  })
})

result = as.data.frame(t(as.data.frame(lapply(gbm_by_large_loop_combine_result,function(x){
  return(c(mean(unlist(x)),
  sd(unlist(x))))
}))))

comb_result_1 = lapply(1:3,function(x){do.call(rbind,lapply(gbm_by_large_loop_combine,function(y){y[[x]]}))})
mean(unlist(lapply(comb_result_1,function(x){sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])})))
sd(unlist(lapply(comb_result_1,function(x){sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])})))

comb_result_1 = lapply(comb_result_1,function(x){x[!grepl("none",x$Var1),]})



all_rbind = do.call(rbind,lapply(simply_merged,function(x){
  x=x[!grepl("none",x$Var1),]
  
}))
sd(unlist(lapply(3:5,function(x){y = all_rbind[,c(1:2,x)];  sum(y[as.character(y$Var1)==as.character(y$Var2),3])/sum(y[,3])  })))
all_rbind = do.call(rbind,simply_merged)


accuracy_l = data.frame(result)
write.table(accuracy_l, file = "accuracy_by_loop.csv")
