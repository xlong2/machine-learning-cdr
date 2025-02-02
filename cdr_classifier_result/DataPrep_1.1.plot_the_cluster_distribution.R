# This script plot the CDR loop clusters  member size distributions
# A distribution is plotted for each  CDR loop type and whether the query CDRs that have  BLAST hits with correct cluster membership or not.
# The x axis is the number of CDRs that are within certain dihedral angle distance to the query CDR.
# The y axis is the density of query CDRs that have certain number of structural neighbors.
# The plots are generated to investigate whether there might be an effect between the number of structural neighbors a query CDR has and whether the query CDR cluster membership will match to that of the BLAST hit CDR cluster membership 


current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}



data_by_loop_type_list_unduplicated=readRDS("./Data_processed/data_by_loop_type_list_unduplicated_no_filtering.rds")

#make a plot of the the data
loop_distribution=lapply(data_by_loop_type_list_unduplicated,function(x){a=as.data.frame(table(x[[1]]$cluster_type)); 
k=a[order(a$Freq,decreasing=TRUE),]; gg=paste(strsplit(as.character(x[[1]][1,"loop_type"]),"_")[[1]][1:2],collapse="_"); k$loop=rep(gg,length(gg)); 
total=0
;return(k)})
loop_distribution_total=do.call(rbind,loop_distribution)



# plot the data distribution 
loop_distribution_total=melt(loop_distribution_total,id.var=c("Var1","loop"))
loop_distribution_total$Var1=as.character(loop_distribution_total$Var1)
names(loop_distribution_total)=c("cluster_identifier","loop","variable","case_number")
# plot the loop distribution
loop_distribution_total$loop_type=sapply(strsplit(loop_distribution_total$loop,"_"),"[[",1)
# parse  loop_type
loop_distribution_total=do.call(rbind,lapply(split(loop_distribution_total,loop_distribution_total$loop_type),function(x){  x$dim=rep(dim(x)[1],dim(x)[1]);return(x)}))
# parse cluster_identifier 
loop_distribution_total$cluster_identifier=as.factor(loop_distribution_total$cluster_identifier)
#parse  looplength
loop_distribution_total$looplength=as.numeric(sapply(strsplit(as.character(loop_distribution_total$loop),"_"),"[[",2))
# parse cluster_identifier 
loop_distribution_total$cluster_identifier= factor(loop_distribution_total$cluster_identifier,levels=c(levels(loop_distribution_total$cluster_identifier)[2:length(levels(loop_distribution_total$cluster_identifier))],"*"))


# separate data by loop_type
loop_distribution_total_list=list()
loop_distribution_total_list=split(loop_distribution_total,loop_distribution_total$loop_type)
figure_dim_list=list(list(66,45),list(65,45),list(120,45),list(46,45),list(78,45))
loop_list=as.list(names(loop_distribution_total_list))

# generate colors
n <- length(unique(loop_distribution_total$cluster_identifier))+1
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector=col_vector[c(1:3,5:(n))]
col_vector=c(col_vector,"#bdbdbd")
names(col_vector) <- c(as.character(unique(loop_distribution_total$cluster_identifier)),"other")
pie(rep(1,n), col=col_vector)






# plot cluster distribution for each loop individually 
for(index in 1:5){
  dim_f=figure_dim_list[[index]]
  the_data=loop_distribution_total_list[[index]]
  the_data=the_data[order(the_data$looplength),]
  the_data$position=the_data$case_number+7
  the_data$looplength=as.factor(the_data$looplength)
  the_data$fill_color=as.character(the_data$cluster_identifier);
  length1= !(the_data$fill_color==1 | grepl("cis",the_data$fill_color))
  the_data[length1,"fill_color"]=rep("other",length(which(length1==TRUE)))
  the_data$fill_color=factor(the_data$fill_color)

  #the_data$cluster_identifier=factor(the_data$cluster_identifier,levels=c(levels(the_data$cluster_identifier)[2:length(levels(the_data$cluster_identifier))],"*"))
  aa=ggplot(data=the_data) +
    geom_bar(stat = "identity", aes(x = cluster_identifier, y = case_number, fill=fill_color) ,position=position_dodge(width = 0.90))+
    geom_text(data=the_data,aes(x=cluster_identifier,y=position,label=cluster_identifier),size=3.5,vjust=0,label.size = 9) +
    facet_grid(~looplength,scales = "free_x", space="free")+ theme_classic()+
    theme(
      strip.text.x = element_text(size = 11, colour = "black"),
      plot.margin=unit(c(2,2,2,2),"mm"),axis.title.x=element_blank(),
      axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank(), axis.text.y = element_text(size=9)) +
    scale_fill_manual("Legend", values = col_vector)+ theme(legend.position = "none")#+ geom_text(aes(label=cluster_identifier), vjust=0) 
  system("mkdir ./Plots/")
  file_name=paste(c("./Plots/member_distribution_plot.",loop_list[[index]],".pdf"),collapse="")
  print(file_name)
  ggsave(aa,file=file_name, width=dim_f[[1]],height=dim_f[[2]],units = c( "mm"))
  
  }

girds=lapply(loop_distribution_total_list,function(x){
  ggplot(data=x) +
    geom_bar(stat = "identity", aes(x = cluster_identifier, y = case_number, fill=cluster_identifier) ,position=position_dodge(width = 0.90))+
    facet_grid(~loop,scales = "free_x", space="free")+ theme(plot.margin=unit(c(2,0,0,2),"mm"),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()) +scale_fill_manual("Legend", values = col_vector)+ theme(legend.position = "none")
})




# plot figure legend
loop_distribution_total$fill_color=loop_distribution_total$cluster_identifier
length1= !(loop_distribution_total$fill_color==1 | grepl("cis",loop_distribution_total$fill_color))
loop_distribution_total$fill_color=as.character(loop_distribution_total$cluster_identifier);
loop_distribution_total[length1,"fill_color"]=rep("other",length(which(length1==TRUE)))
the_data$fill_color=as.factor(the_data$fill_color)

plot=ggplot(data=loop_distribution_total) +
  geom_bar(stat = "identity", aes(x = loop, y = case_number, fill=fill_color,width=1 ),position=position_dodge(width = 0.90))+
  facet_grid(~loop_type,scales = "free", space="free")+scale_fill_manual("Legend", values = col_vector)+
  labs(cluster_identifier="loop&length type")+   guides(fill=guide_legend(title="Cluster identifier"))# add guide properties by aesthetic#Extract Legend 
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend <- g_legend(plot)

p=grid.arrange(arrangeGrob(legend, ncol=1, nrow=1),  widths=c(24),heights=c(101))
file_name="./Plots/member_distribution_figure_legend.pdf"
print(file_name)
ggsave(p,file=file_name)
