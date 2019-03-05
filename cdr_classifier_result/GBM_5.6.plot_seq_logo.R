# This script plot seq logo for CDR loops belonging to some specific clusters 
current_d=getwd()
if(grepl("cdr_classifier_result",current_d)){
  source("0.load_function_and_data.R")
  
}

install.packages("ggseqlogo")
library("ggseqlogo")
loop_l=10
H2_10=data_by_loop_type_list_unduplicated[["H2_10"]][[1]]
H2_10_split=split(H2_10,H2_10$cluster_type)
H2_10_split_1=H2_10_split[["H2-10-1"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_1=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/H2_10_1.pdf", plot = H2_10_1,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)

H2_10_split_1=H2_10_split[["H2-10-2"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_2=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/H2_10_2.pdf", plot = H2_10_2,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)


H2_10_split_1=H2_10_split[["H2-10-6"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_6=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/H2_10_6.pdf", plot =H2_10_6,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)

















loop_l=10
H2_10=data_by_loop_type_list_unduplicated[["L3_10"]][[1]]
H2_10_split=split(H2_10,H2_10$cluster_type)
H2_10_split_1=H2_10_split[["L3-10-none"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_1=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3_10_none.pdf", plot = H2_10_1,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)

H2_10_split_1=H2_10_split[["L3-10-cis8-1"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_2=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3_10_cis8-1.pdf", plot = H2_10_2,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)


H2_10_split_1=H2_10_split[["L3-10-1"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_6=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3_10_1.pdf", plot =H2_10_6,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)

H2_10_split_1=H2_10_split[["L3-10-cis7,8-1"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_6=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3_10_cis7,8-1.pdf", plot =H2_10_6,
       width = 3.5, height = 2.5, units = c("in"), dpi = 1000)


loop_l=9

H2_10=data_by_loop_type_list_unduplicated[["L3_9"]][[1]]
H2_10_split=split(H2_10,H2_10$cluster_type)
H2_10_split_1=H2_10_split[["L3-9-1"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_1=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3-9-1.pdf", plot = H2_10_1,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)

H2_10_split_1=H2_10_split[["L3-9-2"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_1=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3-9-2.pdf", plot = H2_10_1,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)


H2_10_split_1=H2_10_split[["L3-9-cis7-1"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_2=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3-9-cis7-1.pdf", plot = H2_10_2,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)


H2_10_split_1=H2_10_split[["L3-9-cis7-2"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_6=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3-9-cis7-2.pdf", plot =H2_10_6,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)


H2_10_split_1=H2_10_split[["L3-9-cis7-3"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_6=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3-9-cis7-3.pdf", plot =H2_10_6,
       width = 3.5, height = 2.5, units = c("in"),dpi = 1000)

H2_10_split_1=H2_10_split[["L3-10-cis7,8-1"]]
H2_10_seqs_f=H2_10_split_1[,paste("V",(2+10):(2+10+loop_l-1),sep="")]
H2_10_seqs=sapply(as.data.frame(t(H2_10_seqs_f)),function(x){paste(x,collapse="")})
H2_10_6=ggseqlogo(H2_10_seqs,method="prob")
ggsave("./Plots/L3_10_cis7,8-1.pdf", plot =H2_10_6,
       width = 3.5, height = 2.5, units = c("in"), dpi = 1000)

