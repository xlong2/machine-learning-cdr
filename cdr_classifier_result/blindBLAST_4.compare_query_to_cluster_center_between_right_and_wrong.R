#retrieve the number of structural neighbors for any query CDR inside its own cluster 

current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

all_dist_matrix=lapply(all_dist_matrix,function(x){forceSymmetric(x)})
wrong_and_right_cases_list=list()
number_filter=5   # set the minimum number of data points to be 5
for(loop in names(ten_foldcv_blindblastlist)){
  if(is.null(all_dist_matrix[[loop]])){next}
single_f=ten_foldcv_blindblastlist[[loop]]
single_f_r=do.call(rbind,single_f)
single_f_r$query_c=split_vector_and_replace(single_f_r$X1,"\\.",1,1,"")
single_f_r$template_c=split_vector_and_replace(single_f_r$X2,"\\.",1,1,"")
wrong_cases=single_f_r[single_f_r$query_c!=single_f_r$template_c,]
right_cases=single_f_r[single_f_r$query_c==single_f_r$template_c,]
right_cases=right_cases[!duplicated(right_cases),]
wrong_cases=wrong_cases[!duplicated(wrong_cases),]
if(dim(wrong_cases)[1]<number_filter |dim(right_cases)[1]<number_filter){print("small!"); next}


wrong_cases$dis=unlist(apply(wrong_cases,1,query_to_cluster_center_dis))
wrong_cases=add_label_to_vector(wrong_cases,"wrong","worng_or_right")
  
right_cases$dis=unlist(apply(right_cases,1,query_to_cluster_center_dis))
right_cases=add_label_to_vector(right_cases,"right","worng_or_right")
right_cases=right_cases[!is.na(right_cases$dis),]
wrong_cases=wrong_cases[!is.na(wrong_cases$dis),]
if(dim(wrong_cases)[1]<number_filter |dim(right_cases)[1]<number_filter){print("small again!"); next}


wrong_and_right_cases_list[[loop]]=rbind(wrong_cases,right_cases)
}

combined_data_r=do.call(rbind,wrong_and_right_cases_list)
combined_data_r$loop=split_vector_and_replace(combined_data_r$query_c,"-",1,2,"-")
combined_data_r$loop_d=split_vector_and_replace(combined_data_r$loop,"-",1,1,"-")
combined_data_r$length=as.numeric(split_vector_and_replace(combined_data_r$loop,"-",2,2,"-"))

combined_data_r=reorder_factor(combined_data_r,"loop_d","length")
#combined_data_r=combined_data_r[!(combined_data_r$worng_or_right=="right"& grepl("none",combined_data_r$query_c)),]

combined_data_r=combined_data_r[!grepl("none",combined_data_r$query_c),]
combined_data_r$worng_or_right=factor(combined_data_r$worng_or_right,levels=c("wrong","right"))
 the_plot=ggplot(combined_data_r, aes(x=dis,fill=worng_or_right)) +  geom_density(alpha = 0.43,color=NA)+xlim(0,12.5)+
   facet_wrap(~loop,scales="free")+theme_classic()+theme(legend.position = c(0.8,0.1))
 save_figure_specific_size(the_plot,"dihedral_angle_dis_cutoff_5.pdf",7,7)
   
 #the_plot_with_query=ggplot(combined_data_r, aes(x=dis,fill=worng_or_right)) +  geom_density(alpha = 0.5)+
#   facet_wrap(~query_c,scales="free")+theme_classic()+theme(legend.position = c(0.8,0.1))
 

 
 
 

 
