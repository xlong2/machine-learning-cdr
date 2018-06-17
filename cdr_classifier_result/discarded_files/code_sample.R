result_dir="./"

data.sources = list.files(pattern="*.rds")
for(x in data.sources){
  print(x)
  tryCatch({
    x_name=strsplit(x,"\\.")[[1]][1]
    assign(x_name,readRDS(x));print("successfully loaded")},error=function(e){print(e)})
}


# write a list into csv file
write_list_into_single_csv<-function(table_list,output_file){
  count=1
  for(x in names(table_list)){
    if(count==1){write.table("",file=output_file)}
    write.table(x,file=output_file,append=TRUE)
    write.table(table_list[[x]],sep=",",file=output_file,append=TRUE,col.names = TRUE)
    write.table("\n",file=output_file,append=TRUE)
    count=count+1
  }
}

save_file<-function(var_name){
  saveRDS(get(var_name),file=paste(c(result_dir,var_name,".rds"),collapse=""))
}


#compare the blindBLAST to GBM 
diff_between_blast_gbm=list()
for(loop in names(conf_tables_all_loops_blindBLAST_diff)){
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
save_file("diff_between_blast_gbm")




