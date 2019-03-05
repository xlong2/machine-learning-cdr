# load all the functions
file.sources = list.files(pattern="functions.R")
sapply(file.sources,source,.GlobalEnv)
file.sources = list.files(pattern="*function.R")
sapply(file.sources,source,.GlobalEnv)
#setwd(direct)  # go back to previous directory


# this script runs from the Rstudio environment after loading the R project
list.of.packages <- c("Matrix", "grid","caret","RWebLogo","MLmetrics","parallel","pryr","protr",
"gbm","ggplot2","ggrepel","reshape2","gridExtra","doMC","RColorBrewer",
"scales","e1071")
local_package_dir="~/R_libs/"
specify_R_package_diretory=FALSE
install_and_load_packages(list.of.packages,local_package_dir)



result=readRDS("Data_processed/all_data_splitted.rds")
files=list.files("./renumbered_pdbs_Fv/",pattern = "*.pdb")
aho_anno=list()
aho_anno[["H1"]]= 24:42
aho_anno[["H2"]]=57:69
aho_anno[["L1"]]=24:42
aho_anno[["L2"]]=57:72
aho_anno[["L3"]]=107:138
amino_cods = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","M")
names(amino_cods) =c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
             "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE"
             
             )


results_list = list()
list_counter=1
for(each_f in files[1:length(files)]){
  each_f=each_f[1]
  print(each_f)
 # readLines(each_f)
  pdb_o="x"
  tryCatch({
  pdb_o = read.pdb(each_f)},error=function(e){})
  if(pdb_o=="x"){next()}
  file_table=pdb_o$atom
  #file_table = read.table(each_f,fill=TRUE)
  change_counter=0
  loop_sign="X"
  number=0
  
  for(each in 1:dim(file_table)[1]){
    #print(each)
    this_line_loop = file_table[each,"chain"]
    if(this_line_loop !=loop_sign){
      print(this_line_loop)
      loop_sign=this_line_loop
      change_counter = change_counter+1
      print(change_counter)
    }
    if(change_counter==3 |each==dim(file_table)[1]){
      number = each-1
      break()
    }
  }
  file_table = file_table[1:number, ]
  pdb_name = substr(each_f,1,4)
  for(each_name in names(aho_anno)){
    loop = strsplit(each_name,"_")[[1]][1]
    sub_file_table =file_table[which(as.character(file_table$chain)==substr(loop,1,1) ),]
    if(dim(sub_file_table)[1]==0){next()}
    sub_file_table= sub_file_table[!duplicated(sub_file_table$resno),]
    which_loop_corres= which(sub_file_table$resno%in% aho_anno[[each_name]] )
    if(length(which_loop_corres)==0){next}
    loop_length = length(which_loop_corres)
    resi_lines = sub_file_table[which_loop_corres,"resid" ]
    if(any(!resi_lines %in% names(amino_cods))  ){next()}
    seq1=paste(sapply(resi_lines,function(x){ amino_cods[[x]]  }),collapse="")
    
    
    before_ind = which_loop_corres[1]-10
    after_ind = which_loop_corres[length(which_loop_corres)]+10
    if(before_ind<1 | after_ind>dim(sub_file_table)[1]){next}
    this_loop_lines = sub_file_table[before_ind:after_ind,"resid" ]
    if(any(!this_loop_lines %in% names(amino_cods))){next()}
    seq=paste(sapply(this_loop_lines,function(x){ amino_cods[[x]]  }),collapse = "")
    
    results_list[[list_counter]]= c(pdb_name, each_name, loop,loop_length, seq1,seq)
    list_counter = list_counter+1
  }
}

sequence_frame = as.data.frame(do.call(rbind,results_list ))
sequence_frame_v = sequence_frame[!duplicated(sequence_frame[,1:4]),]
sequence_frame_v$ids = apply(sequence_frame_v[,c("V1","V2","V4")],1,function(x){paste(x,collapse="_")})

data_by_loop_type_list_unduplicated_modified = data_by_loop_type_list_unduplicated
for(each_loop in names(data_by_loop_type_list_unduplicated)){
  loop_data = data_by_loop_type_list_unduplicated[[each_loop]]
  frame = loop_data[[1]]
  to_be_delete=c()
  for(ind in 1:dim(frame)[1]){
    ids = unlist(frame[ind,c( "PDB","CDR","length")])
    ids[1]=tolower(ids[1])
    ids=paste(ids,collapse="_")
    which_id =which(sequence_frame_v$ids==ids) 
    if(length(which_id)==0){to_be_delete = c(to_be_delete,ind);next()}
    cdrseq=strsplit(as.character(sequence_frame_v[which_id,"V6"]),"")[[1]]
    feature_n = paste("V",2:(length(cdrseq)+1),sep="")
    frame[ind,feature_n] = cdrseq
  }
  if(length(to_be_delete)>0){
  frame=frame[-to_be_delete,]}
  loop_data[[1]] = frame
  loop_data[[4]]=feature_n
  loop_data[[2]] = as.formula(paste(c("cluster_type ~ ",paste(feature_n,collapse=" + ")),collapse=""))
  #frame   sequence_frame
  data_by_loop_type_list_unduplicated_modified[[each_loop]]=loop_data
}
save_file("data_by_loop_type_list_unduplicated_modified")
