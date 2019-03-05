names_co = c("K", "T", "A", "R", "V", "S", "E", "Q", "P", "I", "L", "G", "N" ,"D" ,"Y" ,"F", "C", "M" ,"H", "W")
names(names_co) = as.character(0:19)
var_importance = list()
unique(sapply(strsplit(names(max_f),"-"),"[[",1))
for(loop in unique(sapply(strsplit(names(max_f),"-"),"[[",1))){
  model_files = max_f[grep(loop, names(max_f))]
  loop_l = as.numeric(strsplit(loop,"_")[[1]][2])
  counter=1 
  this_list=list()
  for(each_model in model_files){
    
    tryCatch({
    each_model = as.character(each_model[1,1])
    
    model = readRDS(file =each_model)
    object = model$finalModel
    
    scale.=TRUE; sort.=TRUE
    if (object$train.fraction < 1) {
      n.trees <- gbm.perf(object, method = "test", plot.it = FALSE)
    }
    
    get.rel.inf <- function(obj) {
      lapply(split(obj[[6]], obj[[1]]), sum)
    }
    temp <- unlist(lapply(object$trees[1:n.trees], get.rel.inf))
    rel.inf.compact <- unlist(lapply(split(temp, names(temp)),   sum))
    rel.inf.compact <- rel.inf.compact[names(rel.inf.compact) != 
                                         "-1"]
    rel.inf <- rep(0, length(object$var.names))
    i <- as.numeric(names(rel.inf.compact)) + 1
    rel.inf[i] <- rel.inf.compact
    names(rel.inf) <- object$var.names
    if (scale.) {
      rel.inf <- rel.inf/max(rel.inf)
    }
    if (sort.) {
      rel.inf <- rev(sort(rel.inf))
    }
    
    
    var_import = as.data.frame(rel.inf[1:200])
    var_import = var_import*100
    transfor_var=as.character(rownames(var_import))
    transfor_var = gsub("^V","", transfor_var)
    var_num = sapply(transfor_var,function(x){ the_c = substr(x, nchar(x), nchar(x)); print(x);if(!is.na(as.numeric(the_c))){print(the_c);the_c= names_co[[the_c]]}; pos= as.numeric(substr(x, 1, c(nchar(x)-1)))-11; print(pos);if(pos<=0){pos = pos-1; print(pos)};c( pos,the_c )})
    var_import$var = as.character(unlist(lapply(as.data.frame(var_num), function(x){ paste(unlist(x),collapse="")})))
    var_import = as.data.frame(var_import[sapply(as.data.frame(t(var_num))[,1],function(x){(as.numeric(as.character(x))+10)<=(loop_l+10)}),])
    
    var_import=var_import[!duplicated(as.character(var_import$var)),]
    names(var_import) = c("Overall", "var")
    var_import=var_import[sapply(gsub("-","",var_import$var),function(x){nchar(x)<=3}),]
    this_list[[counter]] = var_import
    if(counter==1){
      meg = var_import
    }else{
      meg = merge(meg,var_import,by="var",all=TRUE)
    }
    counter = counter +1
    },error=function(e){})
    
  }

  meg$average=apply(meg,1,function(x){mean(unlist(as.numeric(x[2:dim(meg)[2]])),na.rm=TRUE)})
  meg=meg[order(meg$average,decreasing = TRUE),]
  meg$sd = apply(meg, 1, function(x){sd(unlist(as.numeric(x[2:(dim(meg)[2]-1)])),na.rm=TRUE)})
  meg$var = factor(as.character(meg$var), levels = unique(as.character(meg$var)))
  var_importance[[loop]] = meg
  meg=""
  
}
var_importance_r = lapply(names(var_importance),function(x){fr = var_importance[[x]][1:20,]; fr$loop = rep(x,20); colnames(fr)[grepl("Overall",colnames(fr))]=paste("val",1:length(colnames(fr)[grepl("Overall",colnames(fr))]),sep=""); return(fr)})
names(var_importance_r) = names(var_importance)
melted_var_importance = lapply(var_importance_r, function(x){x=x[,grepl(paste(c("var" ,  "va", "loop"),collapse="|"),colnames(x))];melt(x,by.vars=c("var","loop"))})
melted_var_importance = lapply(melted_var_importance,function(x){x$var=factor(as.character(x$var), levels =rev( as.character(unique(x$var)))); return(x)})
loops = c("L3_9", "H2_10", "L3_10","L2_8")
var_list = list()
for(each_l in loops){
  loop_name = gsub("_","-",each_l)
p = ggplot(melted_var_importance[[each_l]], aes(x = var, y = value)) + geom_boxplot() + coord_flip()  +
xlab("")+ylab("") +theme_classic()+theme(plot.title = element_text(hjust = 0.5))+ggtitle(loop_name)
var_list[[each_l]] = p 
}
var_p = grid.arrange(arrangeGrob(var_list[[1]],var_list[[2]],var_list[[3]],var_list[[4]], ncol=2, nrow=2, heights=c(1,1) ,left = textGrob("Sequence Feature",rot=90, gp=gpar(fontsize=14)), bottom=textGrob("Relative Feature importance", gp=gpar(fontsize=14))  ))

save_figure_specific_size(var_p,"suggested_changes/variable_importance.pdf",7,7)

