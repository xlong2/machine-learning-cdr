

#compare the blindBLAST to GBM specifically
diff_between_blast_gbm=list()
for(loop in names(ten_foldcv_blindblastlist)){
  blind_r=conf_tables_all_loops_blindBLAST_diff[[loop]]
  gbm_r=conf_tables_all_loops_gbm_diff[[loop]]
  if(!is.null(gbm_r)&!is.null(blind_r)){
    combined_f=merge(blind_r,gbm_r,all=TRUE,by=c("Var1","Var2"))
    combined_f[is.na(combined_f)] <- 0
    combined_f$gbm_improved=combined_f$Freq.x-combined_f$Freq.y
    combined_f=combined_f[order(combined_f$gbm_improved,decreasing = TRUE),]
    colnames(combined_f)=c("query_cluster","template_cluster","blindBLAST_error_count","gbm_error_count","improvement_by_gbm")
    diff_between_blast_gbm[[loop]]=combined_f
  }
}
write_list_into_single_csv(diff_between_blast_gbm,paste(c(result_dir,"diff_between_blast_gbm_error_count.csv"),collapse = ""))

#compare blindBLAST result to GBM result 
save_file("diff_between_blast_gbm")




