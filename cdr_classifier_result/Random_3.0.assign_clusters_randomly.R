library(reshape2)
# comment library in if you want to write to disk

# load and process Xiyao's data
loop_data = readRDS("loop_data.rds")
conf_tables_blindBLAST = readRDS("conf_tables_blindBLAST.rds")

# these lists of lists are heavy, but we'll work with them
# access is ["loop"][[1]]
# and you always want 1, since there's no other useful list
ptm <- proc.time()
for (loop in names(loop_data)) {
  print(loop)
  # Xiyao bundles some clusters into "none"
  xiyao.clusters<-unique(c(as.character(conf_tables_blindBLAST[[loop]]$Var1), 
                           as.character(conf_tables_blindBLAST[[loop]]$Var2)))
  
  # temporarily extract into dataframe
  loop.df <- loop_data[[loop]][[1]]
  
  # reassign clusters to none if not in Xiyao's list
  loop.df$cluster[!loop.df$cluster_type %in% xiyao.clusters] <- "none"
  loop.df$cluster_type <- paste0(loop.df$CDR, "-", loop.df$length, "-", loop.df$cluster)
  
  # randomize 10000 times
  trials <- 10000
  random.assigment <- replicate(trials, table(data.frame(s=loop.df$cluster_type, d=sample(loop.df$cluster_type))))
  
  # first hack L1_12 since I don't think this none cluster should exist
  if (loop=="L1_12") {
    conf_tables_blindBLAST[[loop]] <- subset(conf_tables_blindBLAST[[loop]], Var1 != "L1-12-none" & Var2 != "L1-12-none")
  }
  # compare random assignments to actual
  # access all counts by cluster pair
  conf_tables_blindBLAST[[loop]]=   conf_tables_blindBLAST[[loop]][conf_tables_blindBLAST[[loop]]$Var1%in%colnames(random.assigment[,,1]),]
  conf_tables_blindBLAST[[loop]]=  conf_tables_blindBLAST[[loop]][conf_tables_blindBLAST[[loop]]$Var2%in%colnames(random.assigment[,,1]),]
  
  pv <- apply(conf_tables_blindBLAST[[loop]], 1, function(x) sum(random.assigment[x[1], x[2], ] >= as.double(x[3]))/trials )
  conf_tables_blindBLAST[[loop]]$pvalue <- pv
  # add the expect mean/sd for error count
  random.mean <- apply(random.assigment, c(1,2), mean)
  random.sd   <- apply(random.assigment, c(1,2), sd)
  conf_tables_blindBLAST[[loop]]$random.error <- apply(conf_tables_blindBLAST[[loop]], 1, function(x) random.mean[x[1], x[2]])
  conf_tables_blindBLAST[[loop]]$random.error.sd <- apply(conf_tables_blindBLAST[[loop]], 1, function(x) random.sd[x[1], x[2]])
  
  conf_tables_blindBLAST[[loop]]$diff <-apply(conf_tables_blindBLAST[[loop]], 1, function(x) as.numeric(x["Freq"])-as.numeric(x["random.error"]))
  # write the table to disk
  dir.create("tables", showWarnings = FALSE)
  # default style
  write.csv(file = paste0("tables/", loop, "_molten.csv"), x = conf_tables_blindBLAST[[loop]])
  
  # new style
  write.csv(file = paste0("tables/", loop, "_blindblast_counts.csv"), 
            x = dcast(conf_tables_blindBLAST[[loop]], formula = Var1 ~ Var2, value.var = "Freq"))
  write.csv(file = paste0("tables/", loop, "_blindblast_sd.csv"), 
            x = dcast(conf_tables_blindBLAST[[loop]], formula = Var1 ~ Var2, value.var = "sd"))
  
  # for the random results
  write.csv(file = paste0("tables/", loop, "_random_counts.csv"), 
            x = dcast(conf_tables_blindBLAST[[loop]], formula = Var1 ~ Var2, value.var = "random.error"))
  write.csv(file = paste0("tables/", loop, "_random_sd.csv"), 
            x = dcast(conf_tables_blindBLAST[[loop]], formula = Var1 ~ Var2, value.var = "random.error.sd"))
  
  write.csv(file = paste0("tables/", loop, "_random_sd.csv"), 
            x = dcast(conf_tables_blindBLAST[[loop]], formula = Var1 ~ Var2, value.var = "random.error.sd"))
  
  
  write.csv(file = paste0("tables/", loop, "_difference.csv"), 
            x = dcast(conf_tables_blindBLAST[[loop]], formula = Var1 ~ Var2, value.var = "diff"))
  
  
  
  # pvalues
  write.csv(file = paste0("tables/", loop, "_pvalues.csv"), 
            x = dcast(conf_tables_blindBLAST[[loop]], formula = Var1 ~ Var2, value.var = "pvalue"))
}
print(proc.time() - ptm)


conf_tables_blindBLAST_pro_list = list()
# group each loop result by the p value 
for(each in names(conf_tables_blindBLAST)){
  conf_tables_blindBLAST_pro_list[[each]] = conf_tables_blindBLAST[[each]]
  conf_tables_blindBLAST_pro_list[[each]]$sig_flag = (conf_tables_blindBLAST[[each]]$pvalue >0.975)
  between_exemplar_distances[[each]]$V1 = gsub("\\*", "none",between_exemplar_distances[[each]]$V1 )
  between_exemplar_distances[[each]]$V2 = gsub("\\*", "none",between_exemplar_distances[[each]]$V2 )
  if(length(between_exemplar_distances[[each]][[1]])==0){
    conf_tables_blindBLAST_pro_list[[each]]$dis= 
      rep(NA,dim(conf_tables_blindBLAST[[each]])[1]);next()}

  dis_f = dcast(as.data.frame(between_exemplar_distances[[each]]),formula =V1~V2,value.var="V3") 
  rownames(dis_f) = dis_f$V1; dis_f$V1<-NULL
  conf_tables_blindBLAST_pro_list[[each]]$dis = apply(conf_tables_blindBLAST_pro_list[[each]],1,function(x){print(x[["Var1"]]);
    dis_f[x[["Var1"]],x[["Var2"]]]
  })
  exp = unlist(lapply(conf_tables_blindBLAST_pro_list[[each]][["dis"]],function(x){is.null(x)}))
  conf_tables_blindBLAST_pro_list[[each]][["dis"]][exp] = rep(NA,length( conf_tables_blindBLAST_pro_list[[each]][["dis"]][exp] ))
  conf_tables_blindBLAST_pro_list[[each]][["dis"]] = c(conf_tables_blindBLAST_pro_list[[each]][["dis"]], recursive = TRUE)
  
}

merged_fr = do.call(rbind, conf_tables_blindBLAST_pro_list)
merged_fr= merged_fr[!grepl("none", merged_fr$Var1),]
merged_fr= merged_fr[!grepl("none", merged_fr$Var2),]

merged_fr = merged_fr[!is.na(merged_fr$dis),]
merged_fr = merged_fr[merged_fr$random.error>3,]
merged_fr$pro = merged_fr$Freq/merged_fr$random.error
merged_fr$pro = as.numeric(merged_fr$pro)
merged_fr$dis = as.numeric(merged_fr$dis)
merged_fr$loop = sapply(strsplit(rownames(merged_fr),"\\."),"[[",1)
#merged_fr = merged_fr[merged_fr$loop %in% c("H1_13", "H2_10", "H2_9","L3_9"),]
p = ggplot(merged_fr,aes(pro, dis ,color=as.factor(sig_flag)))+geom_point()+theme_classic()+theme(legend.position = c(0.8, 0.8))
ggsave("plot.pdf",p, width=7, height=7, unit="in")





result = list()
for(each in c("H1_13","H2_10","L3_9")){
  if(is.null(each)){next}
  print(each)
  if(is.null(x[["sig_flag"]])){next}
  x = conf_tables_blindBLAST[[each]]
  x=as.data.frame(x)
  
  x[["dis"]][exp] = rep(NA,length( x[["dis"]][exp] ))
  x$dis = c(x[["dis"]], recursive = TRUE)
 a =  x[x$sig_flag==TRUE & !is.na(x$dis) ,]
  b = x[x$sig_flag==FALSE & !is.na(x$dis) ,]
  sub_a = a[a$random.error>3,]
  sub_a = sub_a[!grepl("none",sub_a$Var1),]
  sub_a = sub_a[!grepl("none",sub_a$Var2),]
  
  sub_b = b[b$random.error>3,]
  sub_b = sub_b[!grepl("none",sub_b$Var2),]
  sub_b = sub_b[!grepl("none",sub_b$Var1),]
  
  result[[each]] = list()
  result[[each]][["good"]]= sub_a
  result[[each]][["bad"]]= sub_b
  
}
for(each_l in names(result)){
  write.csv(file = paste0("", each_l, "_between_center_dis_good.csv"), 
           result[[each_l]][["good"]])
  write.csv(file = paste0("", each_l, "_between_center_dis_bad.csv"), 
            result[[each_l]][["bad"]])
}
conf_tables_all_loops_blindBLAST

lapply(conf_tables_blindBLAST,function(x){
  x=x[!grepl("none",x$Var1),]
  sum(x[x$Var1!=x$Var2,"Freq"])/sum(x[,"Freq"])
})

bind_f = do.call(rbind, conf_tables_blindBLAST)
bind_f$loop = sapply(strsplit(row.names(bind_f),"_"),"[[",1)
ddply(bind_f,.(loop),.fun =function(x){
  x=x[!grepl("none",x$Var1),]
  sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])
})

ddply(bind_f,.(loop),.fun =function(x){
  sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])
})
bind_f$loop_length = unlist(lapply(strsplit(as.character(bind_f$Var1),"-"),function(x){paste(x[1:2],collapse="_")}))

ddply(bind_f,.(loop),.fun =function(x){
  sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])
})


ddply(bind_f,.(loop_length),.fun =function(x){
  sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])
})
rbind_all_test_result

rbind_all_test_result$loop = sapply(strsplit(rownames(rbind_all_test_result),"_"),"[[",1)
ddply(rbind_all_test_result,.(loop),.fun =function(x){
  sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])
})

rbind_all_test_result$all = rep("all",dim(rbind_all_test_result)[1])

ddply(rbind_all_test_result,.(loop),.fun =function(x){
  x=x[!grepl("none",x$Var1),]
  
  sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])
})
bind_f$all =  rep("all",dim(bind_f)[1])
ddply(bind_f,.(all),.fun =function(x){
  x=x[!grepl("none",x$Var1),]
  
  sum(x[as.character(x$Var1)==as.character(x$Var2),"Freq"])/sum(x[,"Freq"])
})




merged_fr = do.call(rbind, conf_tables_blindBLAST_pro_list)
merged_fr = merged_fr[as.character(merged_fr$Var1)!=as.character(merged_fr$Var2),]
merged_fr$sig_stat = "random"
merged_fr[merged_fr$pvalue<0.025,"sig_stat"] = rep("worse_than_random",length(merged_fr[merged_fr$pvalue<0.025,"sig_stat"]) )
merged_fr[merged_fr$pvalue>0.975,"sig_stat"] = rep("better_than_random",length(merged_fr[merged_fr$pvalue>0.975,"sig_stat"]) )

abc=ggplot(merged_fr, aes(random.error, Freq,colour = factor(sig_stat)))+
  geom_point( position=position_dodge(width = 0.90),size = 2.5)+
  scale_colour_manual(breaks = merged_fr$sig_stat, name="category",
                      values =col_vec)+
  xlab(" error count from random assignment")+ylab("blindBLAST error count")+   theme_classic()+theme_bw()+
  theme(strip.text.x = element_text(size = 12),legend.position = c(0.75, 0.9),
        axis.text.x = element_text(size = 11, colour = "black"),plot.title =element_text(hjust = 0.5,size=18,face="bold"),
        plot.margin=unit(c(2,2,2,2),"mm"),
        axis.text.y = element_text(size=11))+geom_abline(intercept=0, slope=1)
save_figure_specific_size(abc,"suggested_changes/misclassifications_categorized_by_significance_1.pdf",7,7)

