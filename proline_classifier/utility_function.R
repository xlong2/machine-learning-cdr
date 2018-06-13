# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}


chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 


# write a list of dataframes into csv
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

#"L3-9-1"  --->  "L3-9"
split_vector_and_replace<-function(vector,separator,number_start,number_end,re_separator){
  mod_vector=unlist(lapply(strsplit( as.character(vector),separator),function(x){
  x=x[number_start:number_end]
  x=x[!is.na(x)]
  x=paste(x,collapse=re_separator)
  return(x)
  }))
  
  return(mod_vector)
}




#plot a figure with facet_grid  and specify the height and length you want to place the legend
#legend_pos is a vector
plot_figure<-function(data_frame,x_value,y_value,fill_value,facet_value,legend_pos){

  
figure_non_cis=ggplot(data_frame)+geom_bar(aes(x=get(x_value),y=get(y_value),fill=get(fill_value)),
                                           stat = "identity" ,position=position_dodge(width = 0.90))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = legend_pos)+
  xlab("")+facet_grid(~get(facet_value),scales="free",space="free")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 
return(figure_non_cis ) 
}



plot_geom_box_figure<-function(data_frame,x_value,y_value,fill_value,facet_value,legend_pos){
  
  figure_non_cis=ggplot(data_frame,
                        aes(
                          x = get(x_value),
                          y = get(y_value),fill=get(fill_value)
                        ))+geom_boxplot()+
    theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = legend_pos)+xlab("")+facet_grid(~get(facet_value),scales="free",space="free")+
    theme(text = element_text(size=12), axis.text.x = element_text(angle=0, hjust=1),axis.text.y = element_text(angle=90, hjust=1)) 
  return(figure_non_cis )
}


# the first component should be alphabet, the second component should be numeric
# exp L1-13 H1-15 H1-15 H1-14 L1-13  -> H1-14  H1-15 L1-13
order_factor_by_two_component<-function(the_vector,sep_c,first_n,second_n){
  
  splitted_v_1=sapply(strsplit(the_vector,sep_c),"[[",first_n)
  splitted_v_2=sapply(strsplit(the_vector,sep_c),"[[",second_n)
  the_f=data.frame(the_vector,splitted_v_1,splitted_v_2)
  split_by_1=split(the_f,splitted_v_1)
  sorted_v=lapply(names(split_by_1),function(x){
    y=split_by_1[[x]]
    values=as.character(unique(y[,3]))
    sorted_values=sort(as.numeric(unique(values)))
    sorted_values_with_id=paste(rep(x,length(sorted_values)),sorted_values,sep=sep_c)
    
    return(sorted_values_with_id) })
  factor_v=unlist(sorted_v)
  return(factor_v)
}

# assign the value to x
#asssign multiple factors to a data frame
reorder_factor<-function(data_fa,split_col_name,reorder_col_name){
  split_data_fa=split(data_fa,split_col_name)
  split_data_fa_trans=lapply(split_data_fa,function(x){
    values=unique(x[,reorder_col_name])
    values=values[order(as.numeric(values))]
    x[,reorder_col_name]=factor(x[,reorder_col_name],levels=values)
    return(x) })
  mod_split_data_fa=do.call(rbind,split_data_fa_trans)
  return(mod_split_data_fa)
}


# analyze variable importance
inch_width=7;inch_height=7;
save_figure_specific_size<-function(plots_list,file_name,inch_width,inch_height){
  ggsave(plots_list,file=paste(c(plot_dir,file_name),collapse=""),width=inch_width,height=inch_height,unit="in")
  
}

save_file<-function(var_name){
  saveRDS(get(var_name),file=paste(c(result_dir,var_name,".rds"),collapse=""))
}

calculate_accuracy<-function(data_f,def_col_names,freq_col_name){
  diff=sum(unlist(data_f[as.character(data_f[,def_col_names[1]])!=as.character(data_f[,def_col_names[2]]),freq_col_name]))
  all=sum(unlist(data_f[,freq_col_name]))
  acc_v=1-diff/all
  return(acc_v)
  
}


add_label_to_vector<-function(a_frame,label_value,label_name){
  a_frame[[label_name]]=rep(label_value,dim(a_frame)[1])
  return(a_frame)
}


install_and_load_packages<-function(package_list,local_package_dir,specify_R_package_diretory){
  
  
  new.packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){ 
    if(specify_R_package_diretory){install.packages(new.packages,lib=local_package_dir,dependencies=TRUE)}else{install.packages(new.packages,dependencies=TRUE)}
  }
  for(pack in package_list){
    print(pack)
    library(pack,character.only = T)
  }
  
}

Repeat1 <- function(d, n) {
  return(do.call("rbind", replicate(n, d, simplify = FALSE)))
}

Repeat2 <- function(d, n) {
  return(Reduce(rbind, list(d)[rep(1L, times=n)]))
}

Repeat3 <- function(d, n) {
  if ("data.table" %in% class(d)) return(d[rep(seq_len(nrow(d)), n)])
  return(d[rep(seq_len(nrow(d)), n), ])
}

