
aab=ggplot(all_frames_I_care, aes(mean_simu_error_count, error_count,colour = factor(loop_type)))+
  geom_point( position=position_dodge(width = 0.90))+
  scale_colour_manual(breaks = all_frames_I_care$loop_type, name="loop type",
                      values =col_vec)+
  facet_wrap(~sig_stat)+xlab(" error count from random assignment")+ylab("blindBLAST error count")+
  ggtitle("Misclassifications categorized by the significance test")+   
  theme_classic()+theme_bw()+
  theme(strip.text.x = element_text(size = 11),
        axis.text.x = element_text(size = 11, colour = "black"),plot.title =element_text(hjust = 0.5,size=18,face="bold"),
        plot.margin=unit(c(2,2,2,2),"mm"),
        axis.text.y = element_text(size=11)) 
ggsave(aab,units = c("in"),width=9,height=8,file="./proline_classifier/Plots/misclassifications_categorized_by_significance.pdf")
ggsave(aab,units = c("in"),dpi=300,width=9,height=8,file="./proline_classifier/Plots/misclassifications_categorized_by_significance.png")




#make a plot of the the data
loop_distribution = lapply(data_by_loop_type_list_unduplicated, function(x) {
    a = as.data.frame(table(x[[1]]$cluster_type))
    k = a[order(a$Freq, decreasing = TRUE), ]
    gg = paste(strsplit(as.character(k[1, 1]), "-")[[1]][1:2], collapse = "_")
    k$loop = rep(gg, length(gg))
    
    total = 0
    
    
    return(k)
})
loop_distribution_total = do.call(rbind, loop_distribution)



loop_distribution_total = melt(loop_distribution_total, id.var = c("Var1", "loop"))
loop_distribution_total$Var1 = as.character(loop_distribution_total$Var1)
ll = sapply(strsplit(as.character(loop_distribution_total$Var1), "-"), function(x) {
    paste(x[3:length(x)], collapse = "-")
})
loop_distribution_total$Var1 = ll

loop_distribution_total$Var1 =
factor(
loop_distribution_total$Var1,
levels = c(
"1" ,
"2" ,
"3",
"4",
"5" ,
"6" ,
"7",
"8",
"9",
"10",
"11",
"cis6-1",
"cis7-1",
"cis7-2",
"cis7-3" ,
"cis7,8-1" ,
"cis8-1" ,
"cis9-1"  ,
"none"
)
)
names(loop_distribution_total) = c("cluster_identifier", "loop", "variable", "case_number")


# plot the loop distribution

ggplot(loop_distribution_total,
aes(x = loop, y = case_number, fill = cluster_identifier)) + geom_bar(stat = "identity") +
ggtitle("Canonical cdr loop cluster distribution") + theme(plot.title = element_text(hjust = 0.5))



conf = as.data.frame.matrix(table(model$pred$obs, model$pred$pred))
conf_ratio = conf
apply(conf, 1, function(x) {
    x / sum(unlist(x))
})
for (x in 1:dim(conf_ratio)[1]) {
    conf_ratio[x, ] = conf_ratio[x, ] / sum(unlist(conf_ratio[x, ]))
}
conf_number = melt(as.matrix(conf))
conf_number = conf_number[conf_number$Var1 != conf_number$Var2, ]

melted_conf = melt(as.matrix(conf_ratio))
selected_conf = melted_conf[melted_conf[, 1] != melted_conf[, 2], ]
selected_conf_merged = merge(selected_conf, conf_number, by = c("Var1", "Var2"))
selected_conf_merged = selected_conf_merged[order(selected_conf_merged[, 3], decreasing =
TRUE), ]
selected_conf_merged = selected_conf_merged[selected_conf_merged$value.x !=
0, ]
selected_conf_merged$Var1 = gsub("_", "-", selected_conf_merged$Var1)
selected_conf_merged$Var2 = gsub("_", "-", selected_conf_merged$Var2)

conf_tables_all_loops_blindBLAST[["L3_9"]]
selected_conf_merged_blind_blast$method = "blind-blast"
selected_conf_merged$method = "gbm"
selected_conf_merged_blind_blast = selected_conf_merged_blind_blast[selected_conf_merged_blind_blast$value.y >
4, ]
selected_conf_merged = selected_conf_merged[apply(selected_conf_merged[, c("Var1", "Var2")], 1, function(x) {
    paste(x, collapse = "")
}) %in% apply(selected_conf_merged_blind_blast[, c("Var1", "Var2")], 1, function(x) {
    paste(x, collapse = "")
}), ]

the_merged_comparision = rbind(selected_conf_merged_blind_blast, selected_conf_merged)
the_merged_comparision$wrong_prediction = paste(
sapply(the_merged_comparision$Var1, function(x) {
    substr(as.character(x), 6, nchar(as.character(x)))
}),
sapply(the_merged_comparision$Var2, function(x) {
    substr(as.character(x), 6, nchar(as.character(x)))
}),
sep = "_"
)

names(the_merged_comparision)[3:4] = c("ratio", "case_number")
the_merged_comparision_melted = melt((the_merged_comparision[, 3:6]), id =
c("method"   ,     "wrong_prediction"))
ggplot(
the_merged_comparision_melted,
aes(
x = wrong_prediction,
y = value,
shape = variable,
color = method
)
)  + geom_point() + facet_wrap( ~ variable, scales = "free") + ggtitle("gbm vs blindblast performance in L3_9 cluster prediction") +
theme(plot.title = element_text(hjust = 0.5))
