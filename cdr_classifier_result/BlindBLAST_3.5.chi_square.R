# make the chi square test
library("reshape2")


loop_data = readRDS("./Data_processed/loop_data.rds")
conf_tables_blindBLAST = readRDS("./Data_processed/conf_tables_blindBLAST.rds")
chi_square_test_list = list()

#For each loop and length type
for (each_l in names(loop_data)) {
  # wrongs cases counts 
  error_c = conf_tables_all_loops_blindBLAST[[each_l]][,1:3]
  colnames(error_c) = c("query_c",     "pred_c",      "Freq")
  error_c[,1:2] = lapply(error_c[,1:2],as.character)
  this_loop_wrong_cases = error_c[error_c$query_c != error_c$pred_c, ]
  
  # frame for recording date
  frame_colnames = c("query", "pred", "error_count","error_count_other", "expected", "expected_other","chi_sq","sig")
  sig_data_frame = as.data.frame(matrix(
    nrow = dim(this_loop_wrong_cases)[1],
    ncol = length(frame_colnames)
  ))
  colnames(sig_data_frame) = frame_colnames
  sig_data_frame[, 1:3] = this_loop_wrong_cases[, 1:3]
  
  
  #calculate expected error
  data = loop_data[[each_l]][[1]]
  head(data)
  result = as.data.frame(table(data$cluster_type))
  
  perc = (result$Freq / sum(result$Freq))
  expected_v=as.data.frame(matrix(nrow=nrow(result),ncol=nrow(result)))

  for (ind in 1:nrow(result)) {
    each_count = result[ind, 2]
    expected_v[ind,] = each_count * perc
  }
  colnames(expected_v)= as.character(gsub("\\*", "none", result$Var1))
  rownames(expected_v)= as.character(gsub("\\*", "none", result$Var1))
  expected_v = melt(as.matrix(expected_v))
  
  for (each_ind in 1:dim(sig_data_frame)[1]) {
    query = sig_data_frame[each_ind, "query"]
    pred = sig_data_frame[each_ind, "pred"]
    tryCatch({
      expected = expected_v[as.character(expected_v$Var1)==as.character(query)& as.character(expected_v$Var2)==as.character(pred),"value"]
      expected2 = sum(expected_v[expected_v$Var1==query & !expected_v$Var2%in%c(query), "value" ])
      sig_data_frame[each_ind, "expected"] = expected
      sig_data_frame[each_ind, "expected_other"] = expected2
      
    }, error = function(e) {
    })
    if(length(expected)==0 | length(expected2)==0){
      next()
    }
    
    #retrieve actual error
    tryCatch({
      actual = error_c[error_c$query_c == query &  error_c$pred_c == pred, "Freq"]
      actual2 = sum(error_c[error_c$query_c == query & !error_c$pred_c %in% c(query), "Freq"])
      sig_data_frame[each_ind, "error_count_other"] = actual2
    }, error = function(e) {
    })
    
    if (length(actual) == 0 ) {
      actual = 0  }
    if (length(actual2) == 0 ) {
      actual2 = 0  }
    
    # calculate chi square
    chi_sq = (actual - expected) ^ 2 / (expected) + (actual2 - expected2) ^
      2 / (expected2)
    # p value
    p_v=pchisq(chi_sq, df=1, lower.tail=TRUE)
    
    sig_data_frame[each_ind, "chi_sq"] = chi_sq
    sig_data_frame[each_ind, "sig"] = p_v

    #clean
    expected2 <= NULL
    expected <- NULL
    actual <- NULL
    actual2 <- NULL
  }
  
  
  # clean
  sig_data_frame = sig_data_frame[!is.na(sig_data_frame$sig), ]

  if(dim(sig_data_frame)[1]>=1){
    chi_square_test_list[[each_l]] = sig_data_frame
    
  }else{
    chi_square_test_list[[each_l]] <- NULL
    
  }
  

}








