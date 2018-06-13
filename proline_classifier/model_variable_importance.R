
p=ggplot(values)+geom_bar(aes(x=loop,y=value,fill=methods),stat = "identity" ,position=position_dodge(width = 0.90))+
  ggtitle("GBM improved classification counts") +theme(plot.title = element_text(hjust = 0.5))+ylab("count")+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(p,file="error.pdf", width=4.5,height=3.5,units = c( "in"),limitsize = FALSE)



all_the_merged_comparision=list()
accu_summary=list()
bind_all_the_merged_comparision=do.call(rbind,all_the_merged_comparision)
#bind_all_the_merged_comparision=bind_all_the_merged_comparision[grepl("cis",bind_all_the_merged_comparision$wrong_prediction),]
#bind_all_the_merged_comparision_sub=bind_all_the_merged_comparision[bind_all_the_merged_comparision]
spllited_frames=split(bind_all_the_merged_comparision,paste(bind_all_the_merged_comparision$wrong_prediction,bind_all_the_merged_comparision$loop,sep=""))

for(name_a in names(spllited_frames)){
  x=spllited_frames[[name_a]]
  if(dim(as.data.frame(x))[1]==1){
    new_result=rbind(as.data.frame(x),as.data.frame(x))
    #print(new_result)
    new_result[2,"method"]=ifelse(new_result[1,"method"] == "blind-blast", "gbm", "blind-blast")
    new_result[2,"case_number"]=0
    new_result[2,"ratio"]=0
    print(new_result)
    spllited_frames[[name_a]]=new_result
  }else{
    print(x)
    
    
  }
}
spllited_frames_table=do.call(rbind,spllited_frames)
the_merged_comparision_melted=melt( (spllited_frames_table[,3:7]),id=c( "method"   ,     "wrong_prediction","loop"))


#the_merged_comparision_melted_L3_9=the_merged_comparision_melted[the_merged_comparision_melted$loop=="L3_9",]
#the_merged_comparision_melted=the_merged_comparision_melted[the_merged_comparision_melted$loop!="L3_9",]
ggplot(the_merged_comparision_melted, aes(x=wrong_prediction, y=value,shape=variable, 
                                          color=method))  + geom_point(aes(size = value)) +facet_wrap(loop~variable,scales="free" ,ncol=4)+ggtitle("gbm vs blindblast performance in L3_9 cluster prediction") +theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=11), axis.text.x = element_text(angle=90, hjust=1)) 


legend_plot=ggplot(data=the_merged_comparision_melted_cases) +
  geom_bar(stat = "identity", aes(x=wrong_prediction, y=value, fill=method) ,position=position_dodge(width = 0.90))+
  facet_grid(~loop,scales = "free_x", space="free")+ theme(plot.margin=unit(c(2,0,0,2),"mm"),axis.title.x=element_blank(),axis.title.y=element_text(angle=90),axis.text.x=element_blank(),axis.ticks.x = element_blank())+ylab('case number')
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend <- g_legend(legend_plot)

p1=ggplot(data=the_merged_comparision_melted_cases) +
  geom_bar(stat = "identity", aes(x=wrong_prediction, y=value, fill=method) ,position=position_dodge(width = 0.90))+
  facet_grid(~loop,scales = "free_x", space="free")+ theme(plot.margin=unit(c(2,0,0,2),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank())+ylab('case number')+ theme(legend.position = "none")


the_merged_comparision_melted_cases=the_merged_comparision_melted[the_merged_comparision_melted$variable=="case_number",]
the_merged_comparision_melted_cases$real=sapply(strsplit(the_merged_comparision_melted_cases$wrong_prediction,"_"),"[[",1)
the_merged_comparision_melted_cases$pred=sapply(strsplit(the_merged_comparision_melted_cases$wrong_prediction,"_"),"[[",2)
rec=which(the_merged_comparision_melted_cases$real=="1")
prec=which(the_merged_comparision_melted_cases$pred=="1")
others=(1:dim(the_merged_comparision_melted_cases)[1])[!(1:dim(the_merged_comparision_melted_cases)[1]) %in%c(rec,prec)]
the_merged_comparision_melted_cases[rec,"type"]=rep("cluster_1_reco",length(rec))
the_merged_comparision_melted_cases[prec,"type"]=rep("cluster_1_prec",length(prec))
the_merged_comparision_melted_cases[others,"type"]=rep("non_1",length(others))


split_t=split(the_merged_comparision_melted_cases,paste(the_merged_comparision_melted_cases$loop,the_merged_comparision_melted_cases$type,the_merged_comparision_melted_cases$method,sep=""))
summarized_data=do.call(rbind,lapply(split_t,function(x){y=x[1,]; 
y[["value"]]=sum(as.numeric(x$value));return(y)}))


p1=ggplot(data=summarized_data) +
  geom_bar(stat = "identity", aes(x=type, y=value, fill=method),position=position_dodge(width = 0.90))+
  facet_grid(~loop,scales = "free_x", space="free")+ 
  theme(plot.margin=unit(c(2,0,0,2),"mm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=14),
        axis.ticks.x = element_blank())+
  ylab('case number')+ theme(legend.position = "none")+ylab("case number")+
  theme(strip.text.x = element_text(size = 8))


p1=ggplot(data=the_merged_comparision_melted_cases) +
  geom_bar(stat = "identity", aes(x=type, y=value, fill=method),position=position_dodge(width = 0.90))+
  facet_grid(~loop,scales = "free_x", space="free")+ 
  theme(plot.margin=unit(c(2,0,0,2),"mm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=14),
        axis.ticks.x = element_blank())+
  ylab('case number')+ theme(legend.position = "none")+ylab("case number")+
  theme(strip.text.x = element_text(size = 8))




the_merged_comparision_melted_ratio=the_merged_comparision_melted[the_merged_comparision_melted$variable=="ratio",]
p2=ggplot(data=the_merged_comparision_melted_ratio) +
  geom_bar(stat = "identity", aes(x=wrong_prediction, y=value, fill=method) ,position=position_dodge(width = 0.90))+
  facet_grid(~loop,scales = "free_x", space="free")+ 
  theme(plot.margin=unit(c(2,0,0,2),"mm"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle=90),
        axis.text.x= element_text( angle=90,hjust=1) ,  
        axis.text.y = element_text(size=8))+
  ylab("ratio")+ ylim(0, 1)+
  theme(legend.position = "none",strip.text.x = element_text(size = 8))


gfile_name="/Users/xlong3/lab_work_data/machine_learning_cdr/proline_classifier/cis_cluster_rescue_summary.png"
tempfile_name="/Users/xlong3/lab_work_data/machine_learning_cdr/proline_classifier/cluster_rescue_summary.pdf"


ggsave(p,file=tempfile_name, width=100,height=30,units = c( "in"),limitsize = FALSE)


# plot the variable importance in helping finding the results
gbmImp <- varImp(model, scale = TRUE)
#
var_import=as.data.frame(gbmImp[[1]]); var_import$var=rownames(var_import)
var_import=var_import[order(gbmImp[[1]][,1],decreasing=TRUE),]
var_import=var_import[1:20,]
var_import=var_import[order(var_import[,1],decreasing=FALSE),]
var_import$var=factor(var_import$var,levels=var_import$var)
var_import$var=substring(var_import$var,2)
number=as.numeric(substr(var_import$var,1,nchar(var_import$var)-1))-11
res=substring(var_import$var,nchar(var_import$var))
var_import$var=paste(number,res,sep="")
var_import$var=factor(var_import$var,levels=var_import$var)

ggplot(var_import,aes(x=var,y=Overall, fill = "red"))+
  geom_bar(stat="identity")+coord_flip()+ 
  ggtitle("Residue importance for classifying L3_9 using gbm") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+xlab("relative important")+ylab("residue features")

file_name="/Users/xlong3/lab_work_data/machine_learning_cdr/proline_classifier/feature_importance_L3_9.png"
ggsave(file=file_name, width=200,height=200,units = c( "mm"))
gbmImp <- varImp(model, scale = TRUE)



var_import$var

