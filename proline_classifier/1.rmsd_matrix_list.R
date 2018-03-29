
L1_H1=read.table("/Users/xlong3/lab_work_data/DTW_project/L1_H1.distmat")
L1_H1_names=read.table("/Users/xlong3/lab_work_data/DTW_project/L1_H1.file.list",header=FALSE)
rownames(L1_H1)=sapply(strsplit(as.character(L1_H1_names[,1]),"\\/"),"[[",7)
colnames(L1_H1)=sapply(strsplit(as.character(L1_H1_names[,1]),"\\/"),"[[",7)

L2_H2=read.table("/Users/xlong3/lab_work_data/DTW_project/L2_H2.distmat")
L2_H2_names=read.table("/Users/xlong3/lab_work_data/DTW_project/L2_H2.file.list",header=FALSE)
rownames(L2_H2)=sapply(strsplit(as.character(L2_H2_names[,1]),"\\/"),"[[",7)
colnames(L2_H2)=sapply(strsplit(as.character(L2_H2_names[,1]),"\\/"),"[[",7)
L3_L3=read.table("/Users/xlong3/lab_work_data/DTW_project/L3_L3.distmat")
L3_L3_names=read.table("/Users/xlong3/lab_work_data/DTW_project/L3_L3.file.list",header=FALSE)
rownames(L3_L3)=sapply(strsplit(as.character(L3_L3_names[,1]),"\\/"),"[[",7)
colnames(L3_L3)=sapply(strsplit(as.character(L3_L3_names[,1]),"\\/"),"[[",7)