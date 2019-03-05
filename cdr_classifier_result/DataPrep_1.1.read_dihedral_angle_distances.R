# This script put the between CDR loop distance table together and save it as a data matrix. 
all_dist_matrix=list()
for(loop_ty in names(data_by_loop_type_list_unduplicated)){

  if(file.exists(paste(c(loop_ty,"_dist_matrix.rds"),collapse=""))){
    dist_matrix=readRDS(paste(c(loop_ty,"_dist_matrix.rds"),collapse=""))
    print(dist_matrix[1:5,1:5])
    all_dist_matrix[[loop_ty]]=dist_matrix
  }else{
    print(loop_ty)
    
  }
}
save_file("all_dist_matrix")
