# set the maximum number of cluster radius of be 

current_d=getwd()
if(grepl("proline_classifier",current_d)){
  source("0.load_function_and_data.R")
  
}
all_the_max_distances=list()

for(loop_n in names(data_by_loop_type_list_unduplicated)){
  data=data_by_loop_type_list_unduplicated[[loop_n]][[1]]
  #the center of clusters
  centers=data[data$center==1,]
  dis=all_dist_matrix[[loop_n]]
  if(is.null(dis)){ next}
  dis=forceSymmetric(dis)
  for(each_enter in 1:dim(centers)[1]){
    print(centers[each_enter,])
    id=centers[each_enter,"identifier"]
    cluster_d=as.character(centers[each_enter,"cluster_type"])
    cluster_data=data[data$fullcluster==cluster_d,]
    print(cluster_data)
    #get identifiers
    cluster_identifiers=cluster_data[,"identifier"]
    distances_list=c()
    for(x in cluster_identifiers){
      pdb1=substr(strsplit(id,"\\.")[[1]][2],1,4)
      pdb2=substr(strsplit(x,"\\.")[[1]][2],1,4)
      tryCatch({
        print(dis[pdb1,pdb2])
        distances_list=c(distances_list,dis[pdb1,pdb2])
      },error=function(e){})
      
    }
    this_max=max(distances_list,na.rm = TRUE)
    all_the_max_distances[[cluster_d]]=this_max
  }
  
}

save_file("all_the_max_distances")
  