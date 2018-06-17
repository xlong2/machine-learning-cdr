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
save_figure_specific_size(p,"blindBLAST_corrected_with_accuracy_rmsdplot.pdf",7,7)


