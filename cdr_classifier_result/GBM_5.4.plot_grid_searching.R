validation_result_set_for_plot = validation_result_set
validation_result_set_for_plot = lapply(validation_result_set_for_plot,function(x){
  colnames(x)[1:7] = c("file_name","loop_type","complexity","n.trees", "eta","min.node","testing.set"); return(x)})
# get a single result 
splitted_validation_result_set_for_plot = split(validation_result_set_for_plot,sapply(strsplit(names(validation_result_set_for_plot),"-"),"[[",2))
validationfor_plot_set_1 = splitted_validation_result_set_for_plot[[2]]
validationfor_plot_set_1_rbind = do.call(rbind, validationfor_plot_set_1)
rownames(validationfor_plot_set_1_rbind)=1:dim(validationfor_plot_set_1_rbind)[1]
validationfor_plot_set_1_rbind$min = validationfor_plot_set_1_rbind$accuracy - validationfor_plot_set_1_rbind$sd/2
#validationfor_plot_set_1_rbind$min = round(validationfor_plot_set_1_rbind$min, digits=3)
validationfor_plot_set_1_rbind$max = round(validationfor_plot_set_1_rbind$accuracy + validationfor_plot_set_1_rbind$sd/2,digits = 2)
#validationfor_plot_set_1_rbind$max = round(validationfor_plot_set_1_rbind$max, digits=3)

validationfor_plot_set_1_rbind$sd <-NULL
validationfor_plot_set_1_rbind[,4:5] = lapply(validationfor_plot_set_1_rbind[,4:5],function(x){as.numeric(as.character(x))})
validationfor_plot_set_1_rbind = as.data.frame(validationfor_plot_set_1_rbind)
validationfor_plot_set_1_rbind$loop_type = gsub("_","-",validationfor_plot_set_1_rbind$loop_type)
validationfor_plot_set_1_rbind$loop_type = factor(validationfor_plot_set_1_rbind$loop_type, levels = order_factor_by_two_component(as.character(validationfor_plot_set_1_rbind$loop_type),"-",1,2))
  
validationfor_plot_set_1_rbind$blindBLAST
conf_tables_blindBLAST = lapply(conf_tables_blindBLAST,function(x){
  colnames(x)=c("s","d","Count")
  return(x)
})
combined_blindBLAST = do.call(rbind,conf_tables_blindBLAST)
colnames(combined_blindBLAST)[1:3]=c("s","d","Count")
combined_blindBLAST$loop_length = paste(sapply(strsplit(as.character(combined_blindBLAST$s),"-"),"[[",1),sapply(strsplit(as.character(combined_blindBLAST$s),"-"),"[[",2),sep="-")
accuracy = lapply(conf_tables_blindBLAST,function(x){
  sum(x[as.character(x$s)==as.character(x$d),"Count"])/sum(x[,"Count"])
})
apply(validationfor_plot_set_1_rbind,1,function(x){
  
})
data = as.data.frame(t(as.data.frame(accuracy)))
data$loop_type = rownames(data)
data$loop_type = gsub("_","-",data$loop_type)
data = data[data$loop_type%in%validationfor_plot_set_1_rbind$loop_type,]
colnames(data)[1]="accu"
fig = ggplot(validationfor_plot_set_1_rbind,aes(x=n.trees,y=accuracy,ymin = min , ymax = max,color=complexity))  +geom_point(size=0.7,position=position_dodge(width=100),alpha=0.6)+
  geom_errorbar(position = position_dodge(100),alpha=0.6)+geom_hline(data = data, aes(yintercept = accu))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.9,0.9))+
  xlab("# Trees")+ylab("Validation Accuracy")+facet_wrap(~loop_type,scales="free")+
  theme(axis.title=element_text(size=14),  axis.text.x = element_text(angle=90, hjust=1),axis.text.y = element_text(angle=90, hjust=1,size=6),legend.position="bottom",legend.text = element_text(size = 13))
scaleFUN <- function(x) sprintf("%.2f", x)
fig = fig + scale_y_continuous(labels=scaleFUN)
fig




save_figure_specific_size(fig, "grid_searching_test.pdf", 7, 7)


