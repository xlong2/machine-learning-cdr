# make the chi square test 
expected_list = list()
all_loop_sig_list_chi_test = all_loop_sig_frame_list
all_loop_sig_list_chi_test_new= list()
for(each_loop in names(all_loop_sig_list_chi_test)){
  error_c = error_c_list[[each_loop]]/3
  
  data=data_by_loop_type_list_unduplicated[[each_loop]][[1]]
  head(data)
  each_cluster_re = table(data$cluster_type)
  
  result = as.data.frame(each_cluster_re)
  perc = (result$Freq/sum(result$Freq))
  expected_v = list()
  for(ind in 1:dim(result)[1]){
    each_count= result[ind, 2]
    
    expected_v[[ind ]]=each_count*perc
  }
  names(expected_v) = gsub("\\*","none",result$Var1)
  expected_v = as.data.frame(do.call(cbind,expected_v))
  result$Var1 = gsub("\\*","none",result$Var1)
  rownames(result)= result$Var1
  result = as.data.frame(cbind(result,expected_v))
  result$predicted=rownames(result)
  
  result$Var1<-NULL
  result$Freq<-NULL
  melted_result = as.data.frame(melt(result))
  expected_list[[each_loop]] = result

  
  frame = all_loop_sig_list_chi_test[[each_loop]]
  frame = frame[complete.cases(frame),]
  chi_v = rep(NA,dim(frame)[1])
  p_v = rep(NA,dim(frame)[1])
  exp_v = rep(NA,dim(frame)[1])
  actual_e_c = rep(NA, dim(frame)[1])
  for(each_ind in 1:dim(frame)[1]){
    query = frame[each_ind,"Var1"]
    pred = frame[each_ind, "Var2"]
    tryCatch({
    expected = result[pred,query]
    expected2 = sum(result[pred,!colnames(result)%in%c(pred, query,"predicted")])
    
    },error=function(e){})
    
    if(length(expected)==0){
      p_v[[each_ind]]=NA
      exp_v[[each_ind]] =  NA
      actual_e_c[[each_ind]] =  NA
      next()

    }
    #if(expected<10 |expected2<10){
    #  p_v[[each_ind]]=NA
    #  exp_v[[each_ind]] =  NA
    #  actual_e_c[[each_ind]] =  NA
    #  next()
    #}
    
    tryCatch({
    actual = error_c[error_c$query_c == query & error_c$pred_c==pred,"Freq"]
    actual2 = sum(error_c[error_c$query_c == query & !error_c$pred_c%in%c(query,pred),"Freq"])
    
    },error=function(e){})
    if(length(actual)==0){
      actual = 0

    }
    chi_sq = (actual - expected)^2/(expected) +(actual2 - expected2)^2/(expected2)
    p_v[[each_ind]]=pchisq(chi_sq, df=1, lower.tail=TRUE)
    exp_v[[each_ind]] =  expected
    actual_e_c[[each_ind]] =  actual
    expected2<=NULL
    expected<-NULL
  }
  frame$chi_sq_seq = p_v
  frame$expected = exp_v
  frame$error_count = actual_e_c
  frame = frame[!is.na(frame$chi_sq_seq),]
  if(dim(frame)[1]==0){
    all_loop_sig_list_chi_test_new[[each_loop]]<-NULL
  }else{
  all_loop_sig_list_chi_test_new[[each_loop]] = frame }
}
#calculate the dimension
save_file("all_loop_sig_list_chi_test_new")
