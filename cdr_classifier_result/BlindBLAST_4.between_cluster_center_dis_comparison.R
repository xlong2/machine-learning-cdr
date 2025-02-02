# write the script for calculating the between cluster center dihedral distance
# getting the cluster center of each cluster and calculating the distances

current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

between_exemplar_distances=list()
for(each_l in names(data_by_loop_type_list_unduplicated)){
  #get sequences 
  print(each_l)
  this_dis=all_dist_matrix[[each_l]]
  if(is.null(this_dis)){ next }
  this_dis=forceSymmetric(this_dis)
  this_simi=all_similarity_matrix[[each_l]]
  sequences = data_by_loop_type_list_unduplicated[[each_l]][[1]]#load the sequences from a loop length
  # get the cluster exemplars
  exemplars_rows=sequences[sequences$center==1,]
  exemplars_rows = exemplars_rows[!duplicated(exemplars_rows$cluster_type),]
  if(dim(exemplars_rows)[1]<=1){next }
  exemplars_rows$cluster_type=as.character(exemplars_rows$cluster_type)
  all_combi=combn(1:dim(exemplars_rows)[1],2)
  exemplars=exemplars_rows$PDB
  exemplars_dist=this_dis[tolower(exemplars), tolower(exemplars)]
  all_combi=as.data.frame(all_combi)
  distance_info=lapply(all_combi,function(x){c(exemplars_rows[x[1],"cluster_type"], exemplars_rows[x[2],"cluster_type"],exemplars_dist[x[1],x[2]])})
  between_exemplar_distances[[each_l]]=as.data.frame(do.call(rbind,distance_info))
}


#get all misclassifications in loop types 
all_loop_types=split_vector_and_replace(names(data_by_loop_type_list_unduplicated),"_",1,2,"-")  
all_frames_I_care
good_ones=all_frames_I_care[which(grepl("little",all_frames_I_care$sig_stat)),]

bad_ones=all_frames_I_care[which(!grepl("little",all_frames_I_care$sig_stat)),]
bad_ones_chosen=bad_ones[grepl(paste(all_loop_types,collapse="|"),bad_ones$loop_type),]
good_ones_chosen=good_ones[grepl(paste(all_loop_types,collapse="|"),good_ones$loop_type),]
good_ones_chosen=good_ones_chosen[!(grepl("none",good_ones_chosen$query_cluster)|grepl("none",good_ones_chosen$template_cluster)),]
bad_ones_chosen=bad_ones_chosen[!(grepl("none",bad_ones_chosen$query_cluster)|grepl("none",bad_ones_chosen$template_cluster)),]


write.csv("",file="./Data_processed/this_loop_good.csv")

for(loop in all_loop_types){
  print(loop)
  this_loop_good=good_ones_chosen[good_ones_chosen$loop_type==loop,]
  if(dim(this_loop_good)==0) next
  this_loop_bad=bad_ones_chosen[bad_ones_chosen$loop_type==loop,]
  if(dim(this_loop_bad)==0) next
  
  loop_t=paste(strsplit(loop,"-")[[1]][1:2],collapse="_")
  this_loop_dis=between_exemplar_distances[[loop_t]]
  this_loop_dis_r=this_loop_dis
  this_loop_dis_r[,1]=this_loop_dis[,2]
  this_loop_dis_r[,2]=this_loop_dis[,1]
  this_loop_dis_c=rbind(this_loop_dis,this_loop_dis_r)
  colnames(this_loop_dis_c)=c("query_cluster","template_cluster","dist")
  merged_dis_info=merge(this_loop_good,this_loop_dis_c)
  merged_dis_info_s=merged_dis_info[,c("query_cluster","template_cluster","loop_type","dist")]
  write.table("good",sep=",",file="./Data_processed/this_loop_good.csv",append = TRUE)
  write.table(merged_dis_info_s,sep=",",file="./Data_processed/this_loop_good.csv",append = TRUE)


  merged_dis_info_bad=merge(this_loop_bad,this_loop_dis_c)
  merged_dis_info_bad_s=merged_dis_info_bad[,c("query_cluster","template_cluster","loop_type","dist")]
  write.table("bad",sep=",",file="./Data_processed/this_loop_good.csv",append = TRUE)
  write.table(merged_dis_info_bad_s,sep=",",file="./Data_processed/this_loop_good.csv",append = TRUE)
}

