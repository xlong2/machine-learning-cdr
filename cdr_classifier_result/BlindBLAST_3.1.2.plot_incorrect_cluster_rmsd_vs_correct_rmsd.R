# figure of comparing RMSD of templates belonging to wrong clusters vs that of enforcing templates to be the right clusters

current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

all_corrected_ones=lapply(enforcing_correct_rmsd_list,function(x){this=do.call(rbind,x)
  
  })
all_corrected_ones=all_corrected_ones[!unlist(lapply(all_corrected_ones,is.null))]
all_corrected_ones_annot=lapply(names(all_corrected_ones),function(x){
  g=all_corrected_ones[[x]]
g=g[!grepl("none",g$X1),];  a=dim(g)[1]
  g$loop_type=rep(x,a); return(g) }); 
names(all_corrected_ones_annot)=names(all_corrected_ones)
write_list_into_single_csv(all_corrected_ones_annot,"./Data_processed/wrong_template_rmsd_vs_correct_rmsd.csv")


all_corrected_ones_annot_combined=do.call(rbind, all_corrected_ones_annot)
annot_combined <- melt(all_corrected_ones_annot_combined, id.vars=colnames(all_corrected_ones_annot_combined)[!colnames(all_corrected_ones_annot_combined) %in%c("rmsd","correct_cluster_rmsd")])
annot_combined$value=as.numeric(annot_combined$value)
annot_combined_sub=annot_combined[c(1:20,120:150,1000:1200),]
annot_combined_sub$value=as.numeric(annot_combined_sub$value)
annot_combined$loop=sapply(strsplit(annot_combined$loop_type,"_"),"[[",1)
annot_combined$loop_length=sapply(strsplit(annot_combined$loop_type,"_"),"[[",2)
the_ordered=lapply(split(annot_combined$loop_type,annot_combined$loop),function(x){
  a=unique(x);
  the_order=order(as.numeric(sapply(strsplit(a,"_"),"[[",2)));
  ordered_a=a[the_order]
  })
annot_combined$loop_type=factor(annot_combined$loop_type,levels=unlist(the_ordered))
annot_combined$loop

annot_combined=reorder_factor(annot_combined,"loop","loop_length")

the_plot=plot_geom_box_figure(annot_combined,"loop_length","value","variable","loop",c(0.9,0.9))
the_plot=the_plot+ylim(0, 5)



grids=list()
grids[[1]]=p1   # p1 plot object is generated from script 2.blind_blast.R
grids[[2]]=the_plot  
p=grid.arrange(arrangeGrob(grids[[1]],grids[[2]], ncol=1, nrow=2, heights=c(3,2) ,bottom=textGrob("wrong prediction type", gp=gpar(fontsize=14))  ))
print(c("plotted ./Plots/blindBLAST_corrected_with_accuracy_rmsdplot.pdf"))
save_figure_specific_size(p,"suggested_changes/blindBLAST_corrected_with_accuracy_rmsdplot.pdf",7,7)

for_computing=ten_foldcv_blindblastlist
for_computing[["H2_12"]]<-NULL
for_computing[["L1_17"]]<-NULL

for_computing_by_loop_all =lapply(for_computing, function(x){gg=chunk2(1:length(x),3); do.call(rbind,x[gg[[3]]])})
for_computing_by_loop_incorrect = lapply(for_computing_by_loop_all,function(x){x[sapply(strsplit(x[,1],"\\."),"[[",1)!=sapply(strsplit(x[,2],"\\."),"[[",1),]})
for_computing_by_loop =  lapply(for_computing_by_loop_all,function(x){x[sapply(strsplit(x[,1],"\\."),"[[",1)==sapply(strsplit(x[,2],"\\."),"[[",1),]})
accuracy_list
all_corrected_ones
final_result = list()
for(x in names(for_computing_by_loop)){
  
  accuracy_v = accuracy_list[["gbm"]][[x]]
  num = dim(data_by_loop_type_list_unduplicated_for_blindBLAST[[x]][[1]])[1]
  incorrect = for_computing_by_loop_incorrect[[x]]
  ori_correct = for_computing_by_loop[[x]]
  to_sample_num = ceiling(num*accuracy_v - dim(ori_correct)[1])
  if(to_sample_num<=0 | is.na(to_sample_num)){
    names(ori_correct)=c("V1","V2","V3","V4")
    final_result[[x]] = ori_correct
    next}
  corrected = all_corrected_ones[[x]]
  corrected[,2:3] = corrected[,5:6]
  corrected = corrected[,1:4]
  
  sampled_index = sample(1:dim(corrected)[1],to_sample_num)
  sampled_corrected = corrected[sampled_index,]
  names(ori_correct)=c("V1","V2","V3","V4")
  names(sampled_corrected) = c("V1", "V2", "V3", "V4")
  cbind_correct = rbind(ori_correct,sampled_corrected)
  incorrect_retained = incorrect[!incorrect$X1%in%sampled_corrected$V1,]
  names(incorrect_retained) = c("V1", "V2", "V3", "V4")
  final_cbind_result =rbind(cbind_correct,incorrect_retained)
  final_result[[x]] = final_cbind_result
  averaged_rmsd = mean(as.numeric(final_cbind_result[,3]),na.rm=TRUE)
  
}


loop_bind_rmsd = lapply(split(final_result,sapply(strsplit(names(final_result),"_"),"[[",1)),function(x){do.call(rbind,x)})
lapply(loop_bind_rmsd,function(x){mean(as.numeric(x[,3]),na.rm=TRUE)})
all_bind_re = do.call(rbind,final_result)
mean(as.numeric(all_bind_re[,3]),na.rm=TRUE)
1.030
sd(c(1.034002808,1.031711181,1.048322747))

loop_bind_rmsd = lapply(split(final_result,sapply(strsplit(names(final_result),"_"),"[[",1)),function(x){do.call(rbind,x)})
for_computing_by_loop_incorrect
num_greater_than1.5=lapply(for_computing_by_loop_incorrect,function(x){ dim(x[as.numeric(x[,3])>1.5,])[1]})

num_greater_than1.5=lapply(for_computing_by_loop_all,function(x){ dim(x[as.numeric(x[,3])>1.5,])[1]})
num_greater_than1.5_num = lapply(split(num_greater_than1.5,sapply(strsplit(names(num_greater_than1.5),"_"),"[[",1)),function(x){do.call(rbind,x)})
number = lapply(split(for_computing_by_loop_all,sapply(strsplit(names(for_computing_by_loop_all),"_"),"[[",1)),function(x){dim(do.call(rbind,x))[1]})

names = split(names(merged_testing_result),sapply(strsplit(names(merged_testing_result),"_"),"[[",1))

number = lapply(split(for_computing_by_loop_incorrect,sapply(strsplit(names(for_computing_by_loop_incorrect),"_"),"[[",1)),function(x){dim(do.call(rbind,x))[1]})
num =as.data.frame(lapply(data_by_loop_type_list_unduplicated_modified,function(x){dim(x[[1]])[1]}))
num2 = as.data.frame(num_greater_than1.5)
unlist(lapply(num_greater_than1.5_num,function(x){sum(x)}))/unlist(number)



