# This script calculate the number of structural neighbors of any query CDR in its own 
#clusters and make the density plot individually for the right cases and wrong cases when 
#using blindBLAST for the cluster prediction method. 

current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

wrong_right_cases_neighbor_number=list()
filter_number_cut=5
for(loop in names(ten_foldcv_blindblastlist)){
  if(is.null(all_dist_matrix[[loop]])){next}
  single_f=ten_foldcv_blindblastlist[[loop]]
  single_f_r=do.call(rbind,single_f)
  single_f_r$query_c=split_vector_and_replace(single_f_r$X1,"\\.",1,1,"")
  single_f_r$template_c=split_vector_and_replace(single_f_r$X2,"\\.",1,1,"")
  wrong_cases=single_f_r[single_f_r$query_c!=single_f_r$template_c,]
  right_cases=single_f_r[single_f_r$query_c==single_f_r$template_c,]
  right_cases=right_cases[!duplicated(right_cases),]
  right_cases=add_label_to_vector(right_cases,"right","worng_or_right")
  if(!dim(right_cases)[1]>filter_number_cut){next}
  right_cases$neighbor_number=unlist(apply(right_cases,1,query_structure_number_neighbor))
  
  
  
  
  wrong_cases=wrong_cases[!duplicated(wrong_cases),]
  wrong_cases=add_label_to_vector(wrong_cases,"wrong","worng_or_right")
  if(!dim(wrong_cases)[1]>filter_number_cut){next}
  wrong_cases$neighbor_number=unlist(apply(wrong_cases,1,query_structure_number_neighbor))
  
  tab=filter_by_case_number_by_comparison(right_cases,wrong_cases,filter_number_cut)
    
  tab=tab[!grepl("none",tab$query_c),]
  wrong_right_cases_neighbor_number[[loop]]=tab
}


combined_data_neighbor=do.call(rbind,wrong_right_cases_neighbor_number)
combined_data_neighbor$loop=split_vector_and_replace(combined_data_neighbor$query_c,"-",1,2,"-")
combined_data_neighbor$length=as.numeric(split_vector_and_replace(combined_data_neighbor$loop,"-",2,2,"-"))
combined_data_neighbor$loop_d=split_vector_and_replace(combined_data_neighbor$loop,"-",1,1,"-")

combined_data_neighbor=reorder_factor(combined_data_neighbor,"loop_d","length")
  
combined_data_neighbor$worng_or_right=factor(combined_data_neighbor$worng_or_right,levels=c("wrong", "right"))
the_plot=ggplot(combined_data_neighbor, aes(x=neighbor_number,fill=worng_or_right)) +  geom_density(alpha = 0.44,color =NA)+
  facet_wrap(~query_c,scales="free")+theme_classic()+theme(legend.position = c(0.8,0.1))

save_figure_specific_size(the_plot,"number_neighbor_dis_5_cut_off.pdf",7,7)
