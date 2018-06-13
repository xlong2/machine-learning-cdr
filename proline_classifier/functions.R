# make template database for blindBLAST 
make_reference_database <- function(member_seqs,each_fold) {
  #make a blast database
  seq_vec = apply(member_seqs[, features] , 1 , paste , collapse = "")
  names = member_seqs$identifier   # if the clustering scheme is by the torsion angles 
  cluster_types = member_seqs$cluster_type
  file = paste(c("./blast/",loop_type,each_fold,"some_random_file_wont_check_again.fasta"),collapse = "")
  file.remove(file)
  for (index in 1:length(names)) {
    write(paste(">", names[index], each_fold,collapse = ""), file, append = TRUE)
    write(seq_vec[[index]], file, append = TRUE)
  }
  db_name = paste(c("some_db",loop_type,the_method,cluster_dis,each_fold),collapse="_")
  arg = paste(c(' -in ', file, ' -out ', db_name, ' -dbtype prot -hash_index '),
              collapse = " ")
  command = paste("makeblastdb ", arg)
  system(command)
  return(list(db_name))
}

make_reference_database_with_f <- function(member_seqs,each_fold,features) {
  #make a blast database
  seq_vec = apply(member_seqs[, features] , 1 , paste , collapse = "")
  names = member_seqs$identifier   # if the clustering scheme is by the torsion angles 
  cluster_types = member_seqs$cluster_type
  file = paste(c("./blast/",loop_type,each_fold,"some_random_file_wont_check_again.fasta"),collapse = "")
  file.remove(file)
  for (index in 1:length(names)) {
    write(paste(">", names[index], each_fold,collapse = ""), file, append = TRUE)
    write(seq_vec[[index]], file, append = TRUE)
  }
  db_name = paste(c("some_db",loop_type,the_method,cluster_dis,each_fold),collapse="_")
  arg = paste(c(' -in ', file, ' -out ', db_name, ' -dbtype prot -hash_index '),
              collapse = " ")
  command = paste("makeblastdb ", arg)
  system(command)
  return(list(db_name))
}


find_the_best_similarity<-function(identifier,member_seqs_pdbs){
  print(identifier)
  if(length(member_seqs_pdbs)!=1){
    identifier_pdb=substr(split_vector_and_replace(identifier,"\\.",2,2,""),1,4)
    pdbs=substr(split_vector_and_replace(member_seqs_pdbs,"\\.",2,2,""),1,4)
      simi_col_4=substr(as.character(split_vector_and_replace(colnames(similarity_matrix),"\\.",2,2,"")),1,4)
    index=which(simi_col_4 %in%pdbs)  
    which_id=which(simi_col_4==identifier_pdb)
  best_similarity_template_id=colnames(similarity_matrix[,index])[
    which(similarity_matrix[which_id,index]==max(unlist(similarity_matrix[which_id,index]),na.rm=TRUE)[1])
    ]
 }else{best_similarity_template_id=member_seqs_pdbs};
  return(best_similarity_template_id[1])
}

#run blast to search template for a single fold in blindBLAST
runblast_and_retrive_rmsd <-function( seq,each_fold,member_seqs_pdbs,db_name){
   #the blastdbase is the directory of the database
  print(the_method)
  print(loop_type)
  print(features)
  print(cluster_dis)
  print(db_name)
  
  out_file = paste(c("./blast/blast_his_file_need_to_parse_",loop_type,the_method,cluster_dis,each_fold,".txt"),collapse = "")  # the file to record hits result
  file.remove(out_file);
  testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
  file = paste(c("./blast/testing_file_wont_check_again_",loop_type,the_method,cluster_dis,each_fold,".txt"),collapse=""); file.remove(file)
  write(paste(">", seq$identifier, collapse = ""), file, append =TRUE)
  write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
  find_flag=FALSE
  choice_index=1
  choices=c("PAM30","BLOSUM62")
  while(!find_flag){
    print("inside while")
   # remove(hit)
    if(choice_index!=3){
  sub_choice=choices[[choice_index]]
  command = paste(
        c(
          'blastp' ,
          " -db ",
          db_name,
          " -query ",
          file ,
          " -out" ,
          out_file ,
          "-outfmt 6 -word_size 2  -max_target_seqs 1024 -num_threads 4 ",
          '  -evalue 2000 -matrix ',
          "PAM30"
        ),
        collapse = " "
      )
      system(command)
      tryCatch({
      hit = read.table(out_file, stringsAsFactors = FALSE)
      find_flag=TRUE
      print(c("find flag now is ",find_flag))
      returned_id=hit[1,"V2"]
      },error=function(e){},finally = {
         print(find_flag)
      }  )#this table contains all the result of sequences and and the sequence corresponded hits
      choice_index=choice_index+1
    }else{
      identifier=seq$identifier
      returned_id=find_the_best_similarity(identifier,member_seqs_pdbs)
      if(is.character(returned_id)){
      find_flag=TRUE}
      choice_index=choice_index+1
    }
  }

  return(list(returned_id,choice_index-1) )
}

# make data division for k folds and r repeats croo validation
make_3_10_cross_val<-function(training_cases,r,k){
  
  all_unique_ids=lapply(split(training_cases,training_cases$cluster_type),function(x){unique(x$id)})
  
  split_data=split(training_cases,training_cases$cluster_type)
  new_training_data=do.call(rbind,lapply(names(split_data),function(x){
    for(repeat_n in 1:r){
      all_unique_id_sample=lapply(all_unique_ids,function(y){if(length(y)>=k){sample(1:k,size=length(y),replace=TRUE)}else{sample(1:k,size=length(y),replace=FALSE)}})
      
      for(ind in 1:length(all_unique_ids[[x]])){
        all_unique_ids[[x]]
        repnum=dim(split_data[[x]][split_data[[x]]$id==all_unique_ids[[x]][ind],])[1]
        split_data[[x]][split_data[[x]]$id==all_unique_ids[[x]][ind],paste(c("fold.num",repeat_n),collapse="")]=rep(all_unique_id_sample[[x]][ind],repnum)
      }
    }
    return(split_data[[x]])
  }))
  
  
  # Create folds and repeats here - you could create your own if you want #
  
  
  folds.list.out <- list()
  folds.list <- list()
  list.counter <- 1
  for (y in 1:r) {
    newcol <- paste('fold.num', y, sep='')
    for (z in 1:k) {
      out_rown= which(new_training_data[,newcol]==z)  # find the folds in 
      folds_in =which(new_training_data[,newcol]!=z)
      a=new_training_data[folds_in,"id"]%in%  new_training_data[out_rown,"id"]
      print(a[a])
      
      sub=new_training_data[out_rown,"id"]  # ids corresponding to the fold
      out_rown=out_rown[which(!duplicated(sub))]
      a=new_training_data[folds_in,"id"]%in%  new_training_data[out_rown,"id"]
      print(a[a])
      if(length(a[a])!=0){
        print("The 1fold in and 9 fold out is not correct, check!")
        stop()
      }
      folds.list.out[[list.counter]] <- out_rown
      folds.list[[list.counter]] <- which(new_training_data[,newcol]!=z)
      list.counter <- list.counter + 1
    }
  }
  index_list=list()
  index_list[[1]]=folds.list.out
  index_list[[2]]=folds.list
  returned_list=list()
  returned_list[[1]]=index_list;
  returned_list[[2]]=new_training_data
  return(returned_list)
}


#Retrieve the rmsd between the predicted sequence and query sequence by their ids. Record it 
get_rmsd_by_ids<-function(id1,id2){
  pdb1=substr(strsplit(as.character(id1),"\\.")[[1]][2],1,4)
  pdb2=substr(strsplit(as.character(id2),"\\.")[[1]][2],1,4)
    rmsd=rmsd_matrix[pdb1,pdb2]
  return(rmsd)
}

#
get_accuracy_per_fold<-function(each_fold){
  # get reference database
  print(c("the each_fold is ",each_fold))
  #print(sequences)
  #print(folds.list.out)
  print(each_fold)
  member_seqs = sequences[folds.list[[each_fold]],]
  returned_db = make_reference_database(member_seqs,each_fold) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
  member_seqs_pdbs=(member_seqs$identifier)
    fold_out_cases = sequences[folds.list.out[[each_fold]],]
  rmsd_list=data.frame(matrix(nrow=dim(fold_out_cases)[1],ncol=4))
  system("mkdir ./blast/")
  for(each_ind in 1:dim(fold_out_cases)[1]){
    print(each_ind)
    #tryCatch({
    seq=fold_out_cases[each_ind,]
    case_id = seq["identifier"]
    print(seq)
    found_template= runblast_and_retrive_rmsd( seq,each_fold,member_seqs_pdbs,returned_db[[1]])
    query_template_rmsd= get_rmsd_by_ids(as.character(case_id),found_template[[1]])
    system("mkdir ./garbage ")      
    system(paste(c("mv ", returned_db[[1]], "* " , "./garbage/"),collapse=""))   

    rmsd_list[each_ind,1:4]=c(as.character(case_id),as.character(found_template[[1]]),query_template_rmsd,found_template[[2]])
  }# end of iterating through the fold out cases for a single fold 
  return(rmsd_list)
} 





runblast_and_retrive_bitscore <-function( seq,each_fold,member_seqs_pdbs,db_name){
  #the blastdbase is the directory of the database
  print(the_method)
  print(loop_type)
  print(features)
  print(cluster_dis)
  print(db_name)
  
  out_file = paste(c("./blast/blast_his_file_need_to_parse_",loop_type,the_method,cluster_dis,each_fold,".txt"),collapse = "")  # the file to record hits result
  file.remove(out_file);
  testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
  file = paste(c("./blast/testing_file_wont_check_again_",loop_type,the_method,cluster_dis,each_fold,".txt"),collapse=""); file.remove(file)
  write(paste(">", seq$identifier, collapse = ""), file, append =TRUE)
  write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
  find_flag=FALSE
  choice_index=1
  choices=c("PAM30","BLOSUM62")
  while(!find_flag){
    print("inside while")
    # remove(hit)
    if(choice_index!=3){
      sub_choice=choices[[choice_index]]
      command = paste(
        c(
          'blastp' ,
          " -db ",
          db_name,
          " -query ",
          file ,
          " -out" ,
          out_file ,
          "-outfmt 6 -word_size 2  -max_target_seqs 1024 -num_threads 4 ",
          '  -evalue 2000 -matrix ',
          "PAM30"
        ),
        collapse = " "
      )
      system(command)
      tryCatch({
        hit = read.table(out_file, stringsAsFactors = FALSE)
        find_flag=TRUE
        print(c("find flag now is ",find_flag))
        # retrive 5 correct,  
        query_c=split_vector_and_replace(hit$V1,"\\.",1,1,"")
        template_c=split_vector_and_replace(hit$V2,"\\.",1,1,"")
        correct_ns=which(query_c==template_c)
        incorrect_ns=which(query_c!=template_c)
        if(length(correct_ns)>=5){correct_ns=correct_ns[1:5]}else{correct_ns=correct_ns[1:length(correct_ns)]}
        if(length(incorrect_ns)>=5){incorrect_ns=incorrect_ns[1:5]}else{incorrect_ns=incorrect_ns[1:length(incorrect_ns)]}
        correct_bit_scores=hit[correct_ns,"V12"]
        incorrect_bit_scores=hit[incorrect_ns,"V12"]
         # find the shorter one
        which_shorter=which.min(c(length(correct_bit_scores),length(incorrect_bit_scores)))
        min_len=min(c(length(correct_bit_scores),length(incorrect_bit_scores)))
        max_len=max(c(length(correct_bit_scores),length(incorrect_bit_scores)))
        
        diff_r=c()
        for(i in 1:max_len){
          #when j is not where I'm 
          if(i<=min_len){
            j=i}else{
              j=min_len
              }
          if(which_shorter==1){
            diff_r=c(diff_r, correct_bit_scores[j]-incorrect_bit_scores[i])
          }else{
            diff_r=c(diff_r, correct_bit_scores[i]-incorrect_bit_scores[j])
          }
          
          
      }
        return(diff_r ) 
      },error=function(e){},finally = {
        print(find_flag)
      }  )#this table contains all the result of sequences and and the sequence corresponded hits
      choice_index=choice_index+1
    }else{

      return(NULL)
    }
  }
  
  
}





runblast_for_unmatched_random_and_retrive_bitscore <-function( seq,each_fold,member_seqs_pdbs,db_name){
  #the blastdbase is the directory of the database
  print(the_method)
  print(loop_type)
  print(features)
  print(cluster_dis)
  print(db_name)
  
  out_file = paste(c("./blast/blast_his_file_need_to_parse_",loop_type,the_method,cluster_dis,each_fold,".txt"),collapse = "")  # the file to record hits result
  file.remove(out_file);
  testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
  file = paste(c("./blast/testing_file_wont_check_again_",loop_type,the_method,cluster_dis,each_fold,".txt"),collapse=""); file.remove(file)
  write(paste(">", seq$identifier, collapse = ""), file, append =TRUE)
  write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
  find_flag=FALSE
  choice_index=1
  choices=c("PAM30","BLOSUM62")
  while(!find_flag){
    print("inside while")
    # remove(hit)
    if(choice_index!=3){
      sub_choice=choices[[choice_index]]
      command = paste(
        c(
          'blastp' ,
          " -db ",
          db_name,
          " -query ",
          file ,
          " -out" ,
          out_file ,
          "-outfmt 6 -word_size 2  -max_target_seqs 1024 -num_threads 4 ",
          '  -evalue 2000 -matrix ',
          "PAM30"
        ),
        collapse = " "
      )
      system(command)
      tryCatch({
        hit = read.table(out_file, stringsAsFactors = FALSE)
        find_flag=TRUE
        print(c("find flag now is ",find_flag))
        # retrive 5 correct,  
        query_c=split_vector_and_replace(hit$V1,"\\.",1,1,"")
        template_c=split_vector_and_replace(hit$V2,"\\.",1,1,"")
        correct_ns=which(query_c==template_c)
        incorrect_ns=which(query_c!=template_c)
        if(length(incorrect_ns)>=5){t_n=5}else{t_n=length(incorrect_ns)}
        sampled_incorrect_ns=sample(incorrect_ns,t_n)
        incorrect_bit_scores_sampled=hit[sampled_incorrect_ns,"V12"]
        # find the shorter one
        std_v=c(incorrect_bit_scores_sampled)

        return(std_v )
      },error=function(e){},finally = {
        print(find_flag)
      }  )#this table contains all the result of sequences and and the sequence corresponded hits
      choice_index=choice_index+1
    }else{
      
      return(NULL)
    }
  }
  
  
}




get_bit_scores<-function(each_fold){
  # get reference database
  print(c("the each_fold is ",each_fold))
  print(each_fold)
  member_seqs = sequences[folds.list[[each_fold]],]
  returned_db = make_reference_database(member_seqs,each_fold) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
  member_seqs_pdbs=(member_seqs$identifier)
  fold_out_cases = sequences[folds.list.out[[each_fold]],]
  rmsd_list=data.frame(matrix(nrow=dim(fold_out_cases)[1],ncol=6))
  if(dim(rmsd_list)[1]<1){return(NULL)}
  for(each_ind in 1:dim(fold_out_cases)[1]){
    print(each_ind)
    #tryCatch({
    seq=fold_out_cases[each_ind,]
    case_id = seq["identifier"]
    print(seq)
    bit_scores_1_to_5= runblast_and_retrive_bitscore( seq,each_fold,member_seqs_pdbs,returned_db[[1]])
    rmsd_list[each_ind,1:(length(bit_scores_1_to_5)+1)]=c(as.character(case_id),unlist(bit_scores_1_to_5))
  }# end of iterating through the fold out cases for a single fold 
  return(rmsd_list)
} 






get_random_bit_scores<-function(each_fold){
  # get reference database
  print(c("the each_fold is ",each_fold))
  print(each_fold)
  member_seqs = sequences[folds.list[[each_fold]],]
  returned_db = make_reference_database(member_seqs,each_fold) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
  member_seqs_pdbs=(member_seqs$identifier)
  fold_out_cases = sequences[folds.list.out[[each_fold]],]
  rmsd_list=data.frame(matrix(nrow=dim(fold_out_cases)[1],ncol=6))
  if(dim(rmsd_list)[1]<1){return(NULL)}
  
  for(each_ind in 1:dim(fold_out_cases)[1]){
    print(each_ind)
    #tryCatch({
    seq=fold_out_cases[each_ind,]
    case_id = seq["identifier"]
    print(seq)
    random_bit_scores_1_5= runblast_for_unmatched_random_and_retrive_bitscore( seq,each_fold,member_seqs_pdbs,returned_db[[1]])
    rmsd_list[each_ind,1:(length(random_bit_scores_1_5)+1)]=c(as.character(case_id),random_bit_scores_1_5)
  }# end of iterating through the fold out cases for a single fold 
  return(rmsd_list)
}


# an overriden version if the member_seqs and fold_out_cases are given by 
get_accuracy_per_fold_overload<-function(each_fold,member_seqs,fold_out_cases,features){
  # get reference database
  
  print(c("the each_fold is ",each_fold))
  print(sequences)
  print(member_seqs)
  #print(folds.list.out)
  print(each_fold)
  returned_db = make_reference_database_with_f(member_seqs,each_fold,features) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
  member_seqs_pdbs=(member_seqs$identifier)
  rmsd_list=data.frame(matrix(nrow=dim(fold_out_cases)[1],ncol=4))
  
  for(each_ind in 1:dim(fold_out_cases)[1]){
    #tryCatch({
    seq=fold_out_cases[each_ind,]
    case_id = seq["identifier"]
    
    found_template= runblast_and_retrive_rmsd( seq,each_fold,member_seqs_pdbs,returned_db[[1]])
    query_template_rmsd= get_rmsd_by_ids(as.character(case_id),found_template[[1]])

    rmsd_list[each_ind,1:4]=c(as.character(case_id),as.character(found_template[[1]]),query_template_rmsd,found_template[[2]])
  }# end of iterating through the fold out cases for a single fold 
  
  system(paste(c("mv ", returned_db[[1]], "* " , "./garbage/"),collapse=""))   
  
  
  return(rmsd_list)
} 


get_accuracy_per_fold_enforcing_corrent_fold<-function(each_fold){
 
 # get reference database
  print(c("the each_fold is ",each_fold))
  print(sequences)
  #print(folds.list.out)
  print(each_fold)
  
  theresult_this_fold=the_result_list[[each_fold]]

  in_fold_member_seqs=sequences[!sequences$identifier%in%theresult_this_fold$X1,]
  theresult_this_fold$query_cluster=sapply(strsplit(theresult_this_fold$X1,"\\."),"[[",1)
  theresult_this_fold$template_cluster=sapply(strsplit(theresult_this_fold$X2,"\\."),"[[",1)
  unmatched=which(theresult_this_fold$query_cluster!=theresult_this_fold$template_cluster)
  if(length(unmatched)==0){return(NULL)}else{
  unmatched_cases=theresult_this_fold[unmatched,]
  matched_cases=theresult_this_fold[!unmatched,]
  fold_out_cases = sequences[sequences$identifier%in%unmatched_cases$X1,]
  #split the wrong fold out cases by query clusters
  splited_fold_out_cases=split(fold_out_cases,as.character(fold_out_cases$cluster_type))
  accuracy_per_fold_lists=list()
  for(each_cluster in names(splited_fold_out_cases)){
    this_cluster_fold_out_cases=splited_fold_out_cases[[each_cluster]]

  member_seqs= in_fold_member_seqs[in_fold_member_seqs$cluster_type %in% each_cluster,]
  if(dim(member_seqs)[1]!=0){
  accuracy_per_fold_lists[[each_cluster]]=get_accuracy_per_fold_overload(each_fold,member_seqs,this_cluster_fold_out_cases,features)
  }
  }
  rmsd_list=do.call(rbind,accuracy_per_fold_lists)
  colnames(rmsd_list)[3:4]=c("correct_cluster_rmsd","enforcing_method_choice")
  colnames(unmatched_cases)[3:4]=c("rmsd","method_choice")
  merged_result=merge(unmatched_cases[,1:4],rmsd_list,by=("X1"),all.x = TRUE)
  return(merged_result)}
} 




runblast_and_retrive_similarity <-function(returned_db, seq, features){
  blastdbase_name = returned_db[[1]]  #the blastdbase is the directory of the database
  out_file = paste(c("./blast/blast_his_file_need_to_parse_",loop_type,the_method,cluster_dis,"all",".txt"),collapse = "")  # the file to record hits result
  testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
  file = paste(c("./blast/testing_file_wont_check_again_",loop_type,the_method,cluster_dis,"all",".txt"),collapse=""); file.remove(file)
  write(paste(">", seq$identifier, collapse = ""), file, append =TRUE)
  write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
  if (grepl("PAM", subsitution_matrix)|grepl("BLO", subsitution_matrix)) {
    command = paste(
      c(
        'blastp' ,
        " -db ",
        blastdbase_name,
        " -query ",
        file ,
        " -out" ,
        out_file ,
        "-outfmt 6 -word_size 2  -max_target_seqs 10000 -num_threads 4 ",
        '  -evalue 200000 -matrix ',
        subsitution_matrix
      ),
      collapse = " "
    )
  }
  system(command)
  hit = read.table(out_file, stringsAsFactors = FALSE)  #this table contains all the result of sequences and and the sequence corresponded hits
  if(dim(hit)[1]==0){
    tryCatch({
      command = paste(
        c(
          'blastp' ,
          " -db ",
          blastdbase_name,
          " -query ",
          file ,
          " -out" ,
          out_file ,
          "-outfmt 6 -word_size 2  -max_target_seqs 1024 -num_threads 4 ",
          '  -evalue 2000 -matrix ',
          "BLOSUM62"
        ),
        collapse = " "
      )
      system(command)
      hit = read.table(out_file, stringsAsFactors = FALSE)  #this table contains all the result of sequences and and the sequence corresponded hits
      
    }, error=function(e){})
  }else{
    print("do nothing")
  }
  return(hit)
}


make_database <- function(member_seqs) {
  #make a blast database
  seq_vec = apply(member_seqs[, features] , 1 , paste , collapse = "")
  names = member_seqs$identifier   # if the clustering scheme is by the torsion angles 
  cluster_types = member_seqs$cluster_type
  file = paste(c("./blast/",loop_type,"all","some_random_file_wont_check_again.fasta"),collapse = "")
  file.remove(file)
  for (index in 1:length(names)) {
    write(paste(">", names[index], "all",collapse = ""), file, append = TRUE)
    write(seq_vec[[index]], file, append = TRUE)
  }
  db_name = paste(c("some_db",loop_type,the_method,cluster_dis,"all"),collapse="_")
  arg = paste(c(' -in ', file, ' -out ', db_name, ' -dbtype prot -hash_index '),
              collapse = " ")
  command = paste("makeblastdb ", arg)
  system(command)
  return(list(db_name))
}



get_pairwise_sequence_simi_self_cal<-function(loop_type){
  sequences = data_by_loop_type_list_unduplicated[[loop_type]][[1]]#load the sequences from a loop length
  features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  case_ids= sequences[,"identifier"]
  similarity_matrix=matrix(nrow=dim(sequences)[1],ncol=dim(sequences)[1])
  colnames(similarity_matrix)=case_ids
  rownames(similarity_matrix)=case_ids
  
  
  for(each_ind in 1:(length(case_ids)-1)){
    print(each_ind)
    a=c()
    for(j in (each_ind+1):length(case_ids)){
        a_seq = sequences[each_ind,features]
        b_seq =sequences[j,features]
        seq=rbind(a_seq,b_seq)
        similarity=sum(unlist(lapply(seq,function(x){AAPAM30[AAPAM30_index[[x[1]]],AAPAM30_index[[x[2]]]]})))
        
        similarity_matrix[each_ind,j]=similarity
    }
    print(similarity_matrix[each_ind,])
    
  }
  #forceSymmetric(similarity_matrix)
  return(similarity_matrix)
}


get_pairwise_sequence_simi<-function(sequences,features){
  
  case_ids= sequences[,"identifier"]
  similarity_matrix=matrix(nrow=dim(sequences)[1],ncol=dim(sequences)[1])
  colnames(similarity_matrix)=case_ids
  rownames(similarity_matrix)=case_ids
  
  for(each_ind in 1:length(case_ids)){
    tryCatch({
      member_seqs = sequences[-each_ind,]
      returned_db = make_database(member_seqs, features) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
      seq=sequences[each_ind,]
      hit_table = runblast_and_retrive_similarity(returned_db, seq,features)
      for(x in 1:dim(hit_table)[1]){
        similarity_matrix[each_ind,colnames(similarity_matrix)==hit_table[x,2]]=hit_table[x,12]
      } # finish recording all the available similarity values captured
      returned_db=""
    },error=function(e){})
    print(each_ind/length(case_ids))
  }
  return(similarity_matrix)
}


sig_info_function<-function(y){
    x=all_loop_sig_frame_list[[y]];x$loop=rep(y,dim(x)[1])
    # if the first column has the cluster one and the two clusters are not the ssame
    H1_related=x[x[,"Var1"]!=x[,"Var2"],]
    # most prevalent cluster 
    each_l=paste(unlist(strsplit(y,"_")[[1]]),collapse="-")
    if(each_l=="L3-9"){
      most_p="L3-9-cis7-1$"
    }else{
      most_p=paste(c(each_l,"-1$"),collapse="")
    }
    H1_cluster_one_recovery_fail=H1_related[grepl(most_p,H1_related$Var1),]
    H1_cluster_one_precision_fail=H1_related[grepl(most_p,H1_related$Var2),]
    H1_cluster_one_recovery_fail$error_type=rep("1_reco",dim(H1_cluster_one_recovery_fail)[1])
    H1_cluster_one_precision_fail$error_type=rep("1_prec",dim(H1_cluster_one_precision_fail)[1])
    other_types=H1_related[!(grepl(most_p,H1_related$Var1) |grepl(most_p,H1_related$Var2)),]
    other_types$error_type=rep("non_1",dim(other_types)[1])
    all_annotated=rbind(rbind(other_types,H1_cluster_one_precision_fail),H1_cluster_one_recovery_fail)
    return(all_annotated)
    # if the 
  
}

another_sig_info_function<-function(y){
  x=all_loop_sig_frame_list[[y]];x$loop=rep(y,dim(x)[1])
  # if the first column has the cluster one and the two clusters are not the ssame
  H1_related=x[x[,"Var1"]==x[,"Var2"],]
  # most prevalent cluster 
  each_l=paste(unlist(strsplit(y,"_")[[1]]),collapse="-")
  if(each_l=="L3-9"){
    most_p="L3-9-cis7-1$"
  }else{
    most_p=paste(c(each_l,"-1$"),collapse="")
  }
  H1_cluster_one_recovery_fail=H1_related[grepl(most_p,H1_related$Var1),]
  H1_cluster_one_recovery_fail$error_type=rep("cluster_1",dim(H1_cluster_one_recovery_fail)[1])
  other_types=H1_related[!(grepl(most_p,H1_related$Var1)) ,]
  other_types$error_type=rep("non_1",dim(other_types)[1])
  all_annotated=rbind(other_types,H1_cluster_one_recovery_fail)
  return(all_annotated)
  # if the 
}


calculate_accuracy_mean_std<-function(all_result_list){
  all_result=do.call(rbind,all_result_list)
  all_result$query_cluster=sapply(strsplit(all_result[,1],"\\."),"[[",1)
  all_result$template_cluster=sapply(strsplit(all_result[,2],"\\."),"[[",1)
  accuracy=dim(all_result[all_result$template_cluster==all_result$query_cluster,])[1]/dim(all_result)[1]
  
  #std
    all_result_list_split_by_repeat=chunk2(all_result_list,3)


   accuracies_all_repeats=lapply(all_result_list_split_by_repeat,function(x){
    x=do.call(rbind,x)
    x$query_cluster=sapply(strsplit(x[,1],"\\."),"[[",1)
    x$template_cluster=sapply(strsplit(x[,2],"\\."),"[[",1)
    accuracy=dim(x[x$template_cluster==x$query_cluster,])[1]/dim(x)[1]

    print(accuracy)
  })
  std=sd(unlist(accuracies_all_repeats))
  accuracy=mean(unlist(accuracies_all_repeats))
  return(list(accuracy,std))
}


# calculate the accuracies by repeats
calculate_accuracy_mean_std_by_repeats<-function(all_result_list,repeat_n){
  all_result_list_rbind=lapply(all_result_list,function(model_re){
    model_re$repeats= ceiling(as.numeric(gsub("Resample","",model_re[,"Resample"]))/10)
    model_re$Resample=as.numeric(gsub("Resample","",model_re[,"Resample"]))
    model_re=model_re[,c("pred","obs","Resample","repeats")]
    
    last_n=model_re$Resample[length(model_re$Resample)]
    model_re_by_resample=split(model_re,model_re$Resample)
    model_re_by_repeats=model_re
    loopa=paste(strsplit(model_re[1,1],"_")[[1]][1:2],collapse="_")
    lista=list()
    if(loopa%in% c("H1_13","H2_10","L3_9")){

      
    }else{
      resample_chunks=chunk2(1:length(model_re_by_resample),3)
      rbind_r=lapply(resample_chunks,function(x){
        print(x)
        do.call(rbind,model_re_by_resample[x])
      })
      
    }
    if(length(rbind_r)!=3){
      stop()
    }
    return(rbind_r)
  }
    
  )
  lapply(1:3,function(x){
    frames=lapply(all_result_list_rbind,function(y){
      y[[x]]
    })
    b_f=do.call(rbind,frames)
    dim(b_f[b_f[,1]==b_f[,2],])[1]/dim(b_f)[1]
  })
  
  all_result=do.call(rbind,all_result_list_rbind)
  all_result$query_cluster=sapply(strsplit(as.character(all_result[,1]),"\\."),"[[",1)
  all_result$template_cluster=sapply(strsplit(as.character(all_result[,2]),"\\."),"[[",1)
  split_by_repeats=split(all_result,all_result$repeats)
  accuracies_by_repeat=lapply(split_by_repeats,function(x){
    accuracy=dim(x[x$template_cluster==x$query_cluster,])[1]/dim(x)[1]
  })
  sds=sd(unlist(accuracies_by_repeat))

  result_list=list()
  result_list[[1]]=accuracies_by_repeat
  result_list[[2]]=sds
  return(result_list)
}





generate_folds_foldsout<-function(sequences,r,k){
  training_cases = sequences
  training_cases = sequences[, c(features, "cluster_type", "id", "identifier")]
  returned_results = make_3_10_cross_val(training_cases, r, k)
  folds_spec = returned_results[[1]]
  sequences = returned_results[[2]]   # the training cases would have its mo
  folds.list.out = folds_spec[[1]]
  folds.list = folds_spec[[2]]
  counted_folds = c(1:length(folds.list))[unlist(lapply(folds.list.out, length)) !=
                                            0]
  folds.list=folds.list[counted_folds];folds.list.out=folds.list.out[counted_folds]
  if(any(sequences[folds.list[[1]], "id"] %in% sequences[folds.list.out[[1]], "id"])){stop()}
  return(list(folds.list,folds.list.out))
}

getting_similar_sequences_similarity_matrix_rmsd_matrix<-function(loop_type){
  pass=FALSE
  tryCatch({
    sequences = data_by_loop_type_list_unduplicated[[loop_type]][[1]]#load the sequences from a loop length
    features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
    loop_type_type=strsplit(loop_type,"_")[[1]][1]
    similarity_matrix=forceSymmetric(all_similarity_matrix[[loop_type]])
    rmsd_matrix=forceSymmetric(as.matrix(all_rmsd_list[[loop_type_type]]))
    
    
    seq_pdbs=substr(sapply(strsplit(sequences$identifier,"\\."),"[[",2),1,4)
    similarity_matrix_pdbs=substr(sapply(strsplit(colnames(similarity_matrix),"\\."),"[[",2),1,4)
    rmsd_matrix_pdbs=colnames(rmsd_matrix)
    
    common_pdbs=Reduce(intersect, list(seq_pdbs,similarity_matrix_pdbs,rmsd_matrix_pdbs))
    sequences=sequences[seq_pdbs%in%common_pdbs,]
    which_i=which(similarity_matrix_pdbs%in%common_pdbs)
    similarity_matrix=similarity_matrix[which_i, which_i]
    rmsd_matrix=rmsd_matrix[rmsd_matrix_pdbs %in% common_pdbs, rmsd_matrix_pdbs %in% common_pdbs]
    pass=TRUE
  },error=function(e){})
  if(!pass){
    return(FALSE)
  }else(
    return(list(sequences,similarity_matrix,rmsd_matrix))
  )
}




# get template by similarity score

get_template_by_similarity_score<-function(cluster){
  ind_cl=which(sequences$cluster_type==cluster)
  this_ind_most_sim_ind_list=list()
  for(this_ind in ind_cl){
    # print(this_ind)
    relevant_indexes=(1:dim(sequences)[1])[(1:dim(sequences)[1])!=this_ind]
    the_max=max(unlist(this_simi[this_ind,]),na.rm=TRUE)
    this_ind_most_sim_ind=which(unlist(this_simi[this_ind,])==the_max)[1]
    # if(length(unlist(this_ind_most_sim_ind))>1){
    #    stop()
    #  }
    print(as.data.frame(c(this_ind, this_ind_most_sim_ind)))
    this_ind_most_sim_ind_list=c(this_ind_most_sim_ind_list,this_ind_most_sim_ind)
  }
  final_table=data.frame(sequences[ind_cl,"cluster_type"],sequences[unlist(this_ind_most_sim_ind_list),"cluster_type"])
  match_result=data.frame(query=ind_cl,template=unlist(this_ind_most_sim_ind_list))
  tem=cbind(match_result,final_table)
  tem$out_of_cluster_similarity=rep(NA,dim(tem)[1])
  tem$within_cluster_similarity=rep(NA,dim(tem)[1])
  all_disin_query_cluster=sequences[ind_cl,"dis"]
  percentile <- ecdf(all_disin_query_cluster)
  
  for(x in 1:dim(tem)[1]){
    query_id=tem[x,"query"]
    template_id=tem[x,2]
    query_cluster=tem[x,3]
    template_cluster=tem[x,4]
    temcluster_ids=which(sequences$cluster_type==template_cluster)
    template_dis_s=sequences[temcluster_ids,"dis"]
    temp_percentile=ecdf(template_dis_s)
    non_querycluster_cluster=which(sequences$cluster_type!=query_cluster)
    print(c("query_id ", query_id))
    within_querycluster_ind=ind_cl[ind_cl!=query_id]
    within_most_sim=max(unlist(this_simi[query_id,within_querycluster_ind]),na.rm=TRUE)
    query_to_querycluster_dis=sequences[query_id,"dis"]
    #query_dis_percentile = percentile(query_to_querycluster_dis)
    template_to_tempcluster_dis=sequences[template_id,"dis"]
    query_dis_percentile=percentile(query_to_querycluster_dis)
    temp_dis_percentile=temp_percentile(template_to_tempcluster_dis)
    out_most_sim=max(unlist(this_simi[query_id,non_querycluster_cluster]),na.rm=TRUE)
    #number_nei=which()
    tem[x,"out_of_cluster_similarity"]=out_most_sim
    tem[x,"within_cluster_similarity"]=within_most_sim
    tem[x,"query_to_cluster_distance"]= query_to_querycluster_dis
    tem[x,"query_to_clustercen"]=query_dis_percentile
    tem[x,"template_to_clustercen"]=temp_dis_percentile
    tem[x,"number_neighborhood_structure"]=length(which(unlist(this_dis[tem[x,"query"],])<all_neighborhood_dist[[each_l]]))
    
  }

  return(list(tem,final_table))
}# end of iterating all clusters in this loop type 



# get the number of wrong cases for each misclassification in blindBLAST, averaged by the number of repeat
count_number_for_misclassification<-function(each_l){
  this_loop_wrong_cases=enforcing_correct_rmsd_list[[each_l]]
  wrong_cases=do.call(rbind,this_loop_wrong_cases)
  if(is.null(wrong_cases)){return(NULL)}
  wrong_cases$template_c= sapply(strsplit(wrong_cases$X1,"\\."),"[[", 1)
  wrong_cases$query_c=sapply(strsplit(wrong_cases$X2.x,"\\."),"[[", 1)
  conf_t=as.data.frame(table(wrong_cases$template_c,wrong_cases$query_c))
  conf_t=conf_t[conf_t$Freq!=0,]
  conf_t$Freq=conf_t$Freq/3
  conf_t$misclassification=paste(conf_t$Var1,conf_t$Var2,sep="")
  conf_t$Var1=as.character(conf_t$Var1)
  conf_t$Var2=as.character(conf_t$Var2)
  return(conf_t)
}



#use the all_significance_simulation from the envrionment variables, and calculate the significance of every misclassification observed in blindBLAST
get_significance<-function(this_loop_wrong_cases){
sig_data_frame=data.frame(matrix(nrow=length(this_loop_wrong_cases$misclassification),ncol=8))
colnames(sig_data_frame)=c("error_type","Var1","Var2","significance","error_count","sd","mean_simu_error","effect_size")
rownames(sig_data_frame)=this_loop_wrong_cases$misclassification
for (misclassification in this_loop_wrong_cases$misclassification ){  # start iterating over the error types in this loop to assess its significance
  print(misclassification %in% names( all_significance_simulation[[each_l]]))
  if(!misclassification %in% names( all_significance_simulation[[each_l]])){next}
  this_mis_info=this_loop_wrong_cases[this_loop_wrong_cases$misclassification==misclassification,]
  this_error = this_mis_info["Freq"]
  simu_error=all_significance_simulation[[each_l]][[misclassification]]
  aycdf <- ecdf(simu_error)
  sig_data_frame[misclassification,"error_type"]=misclassification
  sig_data_frame[misclassification,"Var1"]=this_mis_info["Var1"]
  sig_data_frame[misclassification,"Var2"]=this_mis_info["Var2"]
  sig_data_frame[misclassification,"significance"]= aycdf(this_error)
  sig_data_frame[misclassification,"error_count"] =this_error
  sig_data_frame[misclassification,"sd"]=sd(simu_error)
  sig_data_frame[misclassification,"mean_simu_error"]=mean(simu_error)
  sig_data_frame[misclassification,"effect_size"]=(this_error-mean(simu_error))/sd(simu_error)
  
}# end of iterating over all errors in this this loop type
#sig_data_frame=sig_data_frame[sig_data_frame$mean_simu_error>1 | (sig_data_frame$significance>=0.975|sig_data_frame$significance<=0.025),]
return(sig_data_frame)

}


frame_manipulating_for_ploting<-function(accuracy_list){
  
  accuracy_gbm_blast=as.data.frame(accuracy_list)
  accuracy_gbm_blast_remove_unknow=accuracy_gbm_blast[complete.cases(accuracy_gbm_blast), ]
  accuracy_gbm_blast_remove_unknow$loop=split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",1,1,"-")
  accuracy_gbm_blast_remove_unknow$length=as.numeric(split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",2,2,"-"))
  accuracy_gbm_blast_remove_unknow=reorder_factor(accuracy_gbm_blast_remove_unknow,"loop","length")
  
  accuracy_gbm_blast_remove_unknow_melt=melt(accuracy_gbm_blast_remove_unknow,id.vars = c("loop","length"))
  accuracy_gbm_blast_remove_unknow_melt=as.data.frame(accuracy_gbm_blast_remove_unknow_melt)
  return(accuracy_gbm_blast_remove_unknow_melt)
}

frame_for_plot<-function(method_n){
  err_frame=error_count_list_with_sd[[method_n]]
  
  accuracy_gbm_blast_remove_unknow=err_frame[complete.cases(err_frame), ]
  accuracy_gbm_blast_remove_unknow$loop=split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",1,1,"-")
  accuracy_gbm_blast_remove_unknow$length=as.numeric(split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow),"_",2,2,"-"))
  accuracy_gbm_blast_remove_unknow=reorder_factor(accuracy_gbm_blast_remove_unknow,"loop","length")
  accuracy_gbm_blast_remove_unknow$methods=rep(method_n,dim(accuracy_gbm_blast_remove_unknow)[1])
  return(accuracy_gbm_blast_remove_unknow)
}


query_to_cluster_center_dis<-function(x){
  query_c=as.character(x["query_c"])
  print(c("query_c",query_c))
  loop_t=paste(strsplit(as.character(query_c),"-")[[1]][1:2],collapse="_")
  data=data_by_loop_type_list_unduplicated[[loop_t]][[1]]
  center=which(data$center==1&as.character(data$cluster_type)==query_c)
  if(length(center)==1){
  cluster_center_id=substr(strsplit(data[data$center==1&as.character(data$cluster_type)==query_c,"identifier"],"\\.")[[1]][2],1,4)
  
  query_id=substr(strsplit(as.character(x["X1"]),"\\.")[[1]][2],1,4)
  
  dis=all_dist_matrix[[loop_t]]
  this_dis=dis[rownames(dis)==cluster_center_id,colnames(dis)==query_id]}else{
    this_dis=NA
  }
  return(this_dis[1])
  # get the dihedral angle distance
  
}


#get number of neighboring structures
query_structure_number_neighbor<-function(x){
  query_c=as.character(unlist(x[["query_c"]]))
  this_pdb=substr(tolower(strsplit(x[["X1"]],"\\.")[[1]][2]),1,4)
  print(c("query_c",query_c))
  loop_t=paste(strsplit(as.character(unlist(query_c)),"-")[[1]][1:2],collapse="_")
  print(loop_t)
  data_for_loop=data_by_loop_type_list_unduplicated[[loop_t]][[1]]
  dis=forceSymmetric(all_dist_matrix[[loop_t]])
  #split the data by the clusters
  data_per_cluster=split(data_for_loop,as.character(data_for_loop$cluster_type))

  data_this_cluster= data_per_cluster[[query_c]]

    # get he radius
    all_pdbs_this_cluster=tolower(data_this_cluster$PDB)
    number=dim(data_this_cluster)[1]
    radius=1.5 #all_the_max_distances[[query_c]]/15
    if(length(radius)==0){
      return(NA)
    }else if(is.na(radius)){
      return(NA)
    }
    if(!this_pdb %in% colnames(dis) ){ return(NA)}
    pdbs=names(dis[this_pdb,])
    pdbs=pdbs[pdbs!=this_pdb]
    if(length(which(all_pdbs_this_cluster%in%pdbs))<=1){return(NA)}
    all_pdbs_this_cluster=all_pdbs_this_cluster[all_pdbs_this_cluster%in%pdbs]
    all_dis=dis[this_pdb,all_pdbs_this_cluster]
    all_dis=all_dis[!is.na(all_dis)]
    if(length(all_dis)<10){
      return(NA)
    }
    the_number=length(which(all_dis<radius))
    # count the number 
    return(the_number)

  }


# filtering the clusters with certain filter of cases numbers
filter_by_case_number_by_comparison<-function(right_cases,wrong_cases,filter_case_n){
  count_the_cluster=FALSE
  splitted_right=split(right_cases,right_cases$query_c)
  splitted_wrong=split(wrong_cases,wrong_cases$query_c)
  splitted_right_filtered=lapply(names(splitted_right),function(x){
    y=splitted_right[[x]];
print(y)
    real_y_n=dim(y[!is.na(y$neighbor_number),])[1]
    if(real_y_n<filter_case_n) return(NULL) else return(y)
  })
  names(splitted_right_filtered)=names(splitted_right)
  
  splitted_wrong_filtered=lapply(names(splitted_wrong),function(x){
    y=splitted_wrong[[x]];
    print(y)
    real_y_n=dim(y[!is.na(y$neighbor_number),])[1]
    if(real_y_n<filter_case_n) return(NULL) else return(y)
  })
  names(splitted_wrong_filtered)=names(splitted_wrong)
  
  
  splitted_right_filtered=splitted_right_filtered[!unlist(lapply(splitted_right_filtered,is.null))]
  splitted_wrong_filtered=splitted_wrong_filtered[!unlist(lapply(splitted_wrong_filtered,is.null))]
  
  
  splitted_right_filtered_row_binded=lapply(names(splitted_right_filtered),function(x){
    if(!is.null(splitted_wrong_filtered[[x]])){
    binded_info=rbind(splitted_right_filtered[[x]],splitted_wrong_filtered[[x]])
    return(binded_info)
    }else{
      return(NULL)
    }
  })
  right_filtered_row_binded=do.call(rbind,splitted_right_filtered_row_binded)
  return(right_filtered_row_binded)
}


even_out_all_classes<-function(training_cases){
  
  
  splitted_class=split(training_cases,training_cases$cluster_type)
  proportions=unlist(lapply(split(training_cases,training_cases$cluster_type),function(x){dim(x)[1]}))
  max_class=names(proportions)[which(proportions==max(proportions))]
  remaining_class=names(proportions)[names(proportions)!=max_class]
  for(each_remaining in remaining_class){
    ratio=proportions[max_class]%/%proportions[each_remaining]
    if(ratio>1){
      ratio=ratio-1
    }
    to_be_added=do.call(rbind,replicate(ratio,splitted_class[[each_remaining]],simplify = FALSE))
    training_cases=rbind(training_cases,to_be_added)    
  }
  return(training_cases)
}


# getting the training result
generic_train<-function(each_loop_length_data_feature_string_rmsd,each_method,training_cases,gbmGrid,loop_type){
  # to add classes to off set 
  training_cases=even_out_all_classes(training_cases)
  if(dim(training_cases)[1]*0.6*0.5 < as.numeric(as.character(gbmGrid[["n.minobsinnode"]]))){
    training_cases=rbind(training_cases,training_cases)
    training_cases=rbind(training_cases,training_cases)
  }
  n.trees=gbmGrid[["n.trees"]];
  

  set.seed(sample(1:100000,1))
  r <-3 # number of repeats
  k <- 10 # number of folds
  if(loop_type %in% c("H1_13", "H2_10", "L3_9" ) && n.trees>=1000 ){
    r=1}
  returned_results=make_3_10_cross_val(training_cases,r,k)  # make division
  folds_spec=returned_results[[1]]
  training_cases=returned_results[[2]]   # the training cases would have its mo
  folds.list.out=folds_spec[[1]]
  folds.list=folds_spec[[2]]
  fitControl <- trainControl(method = "repeatedcv",
                             repeats = r,
                             number=k,
                             #preProcOptions = list(thresh = 0.95),
                             index=folds.list
                             , indexOut=folds.list.out,
                             ## Estimate class probabilities
                             classProbs = TRUE,
                             returnResamp="all",
                             ## Evaluate performance using
                             ## the following function
                             savePredictions="final",
                             summaryFunction = multiClassSummary)
  trained_model=""
  trained_model <- train(each_loop_length_data_feature_string_rmsd, data = training_cases,
                         #distribution = "adaboost",
                         method = "gbm", bag.fraction = 0.5,   # fold number 10 
                         #nTrain = round(nrow(training_cases) *.75),
                         trControl = fitControl,
                         tuneGrid = gbmGrid,
                         verbose = TRUE,
                         
                         ## Specify which metric to optimize
                         metric = "kappa")
  return(trained_model)
  
}

convert_none<-function(model_re,blindBLAST_unique_clusters){
  obs=gsub("_","-",as.character(model_re$obs))
  pred=gsub("_","-",as.character(model_re$pred))
  
  obs_which=which(!obs%in%blindBLAST_unique_clusters)
  print(obs_which)
  pred_which=which(!pred%in%blindBLAST_unique_clusters)
  print(pred_which)
  pred[pred_which]=rep(paste(c(strsplit(loop,"_")[[1]],"none"),collapse="-"),length(pred_which))
  obs[obs_which]=rep(paste(c(strsplit(loop,"_")[[1]],"none"),collapse="-"),length(obs_which))
  returned_l=list()
  returned_l[[1]]=obs
  returned_l[[2]]=pred
  return(returned_l)
  
  }


train_final_model<-function(each_loop_length_data_feature_string_rmsd,each_method,training_cases,gbmGrid){
  # to add classes to off set 
  training_cases=even_out_all_classes(training_cases)
  if(dim(training_cases)[1]*0.6*0.5 < as.numeric(as.character(gbmGrid[["n.minobsinnode"]]))){
    training_cases=rbind(training_cases,training_cases)
    training_cases=rbind(training_cases,training_cases)
  }
  

  
  r <-1 # number of repeats
  k <- 1 # number of folds

  returned_results=make_3_10_cross_val(training_cases,r,k)  # make division
  folds_spec=returned_results[[1]]
  training_cases=returned_results[[2]]   # the training cases would have its mo

  fitControl <- trainControl(method = "none",
                             ## Estimate class probabilities
                             classProbs = TRUE,
                             returnResamp="all",
                             ## Evaluate performance using
                             ## the following function
                             savePredictions="final",
                             summaryFunction = multiClassSummary)
  trained_model=""
  trained_model <- train(each_loop_length_data_feature_string_rmsd, data = training_cases,
                         #distribution = "adaboost",
                         method = "gbm", bag.fraction = 0.5,   # fold number 10 
                         #nTrain = round(nrow(training_cases) *.75),
                         trControl = fitControl,
                         tuneGrid = gbmGrid,
                         verbose = TRUE,
                         
                         ## Specify which metric to optimize
                         metric = "kappa")
  return(trained_model)
  
}

# get accuracy from blindBLAST result with consideration of the n-cv-m-repeats scheme
get_conf_from_blindBLAST<-function(model_re_com,unique_cluster,repeats){
  model_re_com[,1]=sapply(strsplit(model_re_com[,1],"\\."),"[[",1)
  model_re_com[,2]=sapply(strsplit(model_re_com[,2],"\\."),"[[",1)
  which_1=which(!model_re_com[,1]%in%unique_cluster)
  model_re_com[which_1,1]=rep(paste(c(strsplit(loop,"_")[[1]],"none"),collapse="-"),length(which_1))
  which_2=which(!model_re_com[,2]%in%unique_cluster)
  model_re_com[which_2,2]=rep(paste(c(strsplit(loop,"_")[[1]],"none"),collapse="-"),length(which_2))
  
  
  conf_t=as.data.frame(table(model_re_com[,1],model_re_com[,2]))
  conf_t$Freq=conf_t$Freq/repeats
  accuracy_result=sum(conf_t[as.character(conf_t$Var1)==as.character(conf_t$Var2),"Freq"])/sum(conf_t[,"Freq"])
  
  return(accuracy_result)
}



generate_Rscript_command<-function(args,Rscript_dir,Rscript_n){
  
  args=unlist(lapply(args,as.character))
  executable=paste(args, collapse="_")

  exe_sh_file=paste(c(Rscript_dir,executable,"exe.sh"),collapse="_")  # customize shell script name
  Rscript_command_line=paste(c("Rscript ", Rscript_n,paste(args,collapse=" ")),collapse=" ")    # write shell script
 # write("cd ..", file = exe_sh_file,append=TRUE )
  write(Rscript_command_line, file = exe_sh_file,append=TRUE )
  return(exe_sh_file)
}

generate_Rscript_command_H1<-function(args,Rscript_dir,Rscript_n){
  
  args=unlist(lapply(args,as.character))
  executable=paste(args, collapse="_")
  exe_sh_file=paste(c(Rscript_dir,executable,"exe.sh"),collapse="_")  # customize shell script name
  file.remove(exe_sh_file)
  
  Rscript_command_line=paste(c("Rscript ", Rscript_n,paste(args,collapse=" ")),collapse=" ")    # write shell script
  # write("cd ..", file = exe_sh_file,append=TRUE )
  write(Rscript_command_line, file = exe_sh_file,append=TRUE )
  return(exe_sh_file)
}

generate_condor_script<-function(args,exe_sh_file,condor_script_dir,master_condor_file){
  executable=paste(args, collapse="_")
  con_file=paste(c(executable,".con"),collapse="_")  #customize condor script name
  system(paste(c("cp ", master_condor_file," ", con_file),collapse=""))   # copy the master file to script directory
  
  x <- readLines(con_file)
  y <- gsub( "\\$executable", exe_sh_file, x )
  z <- gsub( "\\$condor_script_dir", condor_script_dir, y )
  w <- gsub( "\\$cores", cores, z )
  cat(w, file=con_file, sep="\n")
  
  print(paste(c( condor_script_dir,"/",con_file),collapse=""))
  a=paste(c("condor_submit ", condor_script_dir,"/",con_file),collapse="")
  print(a)
  return(a)
}




preprocess_data<-function(loop_type,data_by_loop_type_list_unduplicated){
  print("1093")
  data=data_by_loop_type_list_unduplicated[[loop_type]][[1]]
  data$cluster_type=sub("-","_",data$cluster_type)
  data$cluster_type=as.factor(as.character(data$cluster_type))
  
  
  
  print("line 1100")
  #sequences$rmsd_cluster = as.character(sequences$rmsd_cluster)
  sequences=data
  
  print("line 1104")
  features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  print("line 1106")
  print(features)
  formu_str=paste(c("cluster_type ~ ",paste(features,collapse=" + ")),collapse="")
  print(formu_str)
  
  each_loop_length_data_feature_string=as.formula(formu_str)
  print("line 1108")
  all_cases =  sequences[,c(features,"id","cluster_type")]
  the_levels=unique(unlist(data_by_loop_type_list_unduplicated[[loop_type]][[1]][,features]))
  for(each_f in features){
    all_cases[,each_f]=factor(all_cases[,each_f],levels=the_levels)
  }
  print("line 1112")
  all_cases=all_cases[complete.cases(all_cases), ]
  all_cases$cluster_type=gsub("-","_",all_cases$cluster_type)
  all_cases$cluster_type=gsub(",",".",all_cases$cluster_type)
  all_cases$cluster_type=gsub("\\*","none",all_cases$cluster_type)
  
  all_cases$cluster_type=as.factor(as.character(all_cases$cluster_type))
  return(list(all_cases,each_loop_length_data_feature_string))
}

generate_eta<-function(loop){
  nl=dim(data_by_loop_type_list_unduplicated[[loop]][[1]])[1]
  eta=max (0.01, 0.1*min(1, nl/10000))
  return(eta)
}



execute_training_rscript<-function(args_list){   # function for 
  args=args_list[[1]]
  print("line 1141")
  print(args)
  
#  data_by_loop_type_list_unduplicated=args_list[[2]]
  getwd()
print('line 1129')
  # specify parameter as the args 
  print(args)
  print(args[[2]])
  print(args[2])
  loop_type=as.character(args[[1]])
  num_core=as.numeric(args[[2]])
  interaction.depth=as.numeric(args[3])
  n.trees=as.numeric(args[4])
  shrinkage=as.numeric(args[5])
  n.minobsinnode=as.numeric(args[6])
  registerDoMC(num_core)
  print("line 1141")

  # naming the model file for latter easy parsing 
  gbmGrid=
    data.frame(interaction.depth=as.numeric(interaction.depth),
               n.trees=as.numeric(n.trees),shrinkage=as.numeric(shrinkage),
               n.minobsinnode=as.numeric(n.minobsinnode))
  
  print(gbmGrid)
  parameter_spe = paste(unlist(gbmGrid),collapse="-")
  each_method="gbm_test"
  cluster_dis="north"
  gbm_result_dir=paste(c("./rmsd_cluster_hits_rmsd/"),collapse="")   # build a subdirectory in the current directory to store trained models
  print("line 1150")
  
  to_save_file=paste(c(gbm_result_dir,"",loop_type,"_",paste(c(each_method,cluster_dis,parameter_spe),collapse="-"),"_trained_model_extra_test.rds"),collapse="")
  
  # read pyigclassify file 
  print(getwd())

  print("line 1158")
  
  # perform gbm training 
  if(file.exists(to_save_file)){
    print("file_already_exists!")}else{
      print(loop_type)
      print(loop_type)
      preprocess_result=  preprocess_data(loop_type,data_by_loop_type_list_unduplicated)
      all_cases=preprocess_result[[1]]
      each_loop_length_data_feature_string=preprocess_result[[2]]
      print("line 1165")
      
      trained_model = generic_train(each_loop_length_data_feature_string,each_method,all_cases,gbmGrid,loop_type) 
      print("line 1168")
      
      saveRDS(trained_model,file =to_save_file)
      print(to_save_file)

      
      
    }
}










execute_training_rscript_comp<-function(args_list){   # function for 
  args=args_list[[1]]
  print("line 1141")
  print(args)
  
  #  data_by_loop_type_list_unduplicated=args_list[[2]]
  getwd()
  print('line 1129')
  # specify parameter as the args 
  print(args)
  print(args[[2]])
  print(args[2])
  loop_type=as.character(args[[1]])
  num_core=as.numeric(args[[2]])
  interaction.depth=as.numeric(args[3])
  n.trees=as.numeric(args[4])
  shrinkage=as.numeric(args[5])
  n.minobsinnode=as.numeric(args[6])
  batch=args_list[[3]]
  registerDoMC(num_core)
  print("line 1141")
  
  # naming the model file for latter easy parsing 
  gbmGrid=
    data.frame(interaction.depth=as.numeric(interaction.depth),
               n.trees=as.numeric(n.trees),shrinkage=as.numeric(shrinkage),
               n.minobsinnode=as.numeric(n.minobsinnode))
  
  print(gbmGrid)
  parameter_spe = paste(unlist(gbmGrid),collapse="-")
  each_method="gbm_test"
  cluster_dis="north"
  gbm_result_dir=paste(c("./rmsd_cluster_hits_rmsd/"),collapse="")   # build a subdirectory in the current directory to store trained models
  print("line 1150")
  
  to_save_file=paste(c(gbm_result_dir,"",loop_type,"_",paste(c(each_method,cluster_dis,parameter_spe,batch),collapse="-"),"_trained_model_extra_test.rds"),collapse="")
  
  # read pyigclassify file 
  print(getwd())
  
  print("line 1158")
  
  # perform gbm training 
  if(file.exists(to_save_file)){
    print("file_already_exists!")}else{
      print(loop_type)
      print(loop_type)
      preprocess_result=  preprocess_data(loop_type,data_by_loop_type_list_unduplicated)
      all_cases=preprocess_result[[1]]
      each_loop_length_data_feature_string=preprocess_result[[2]]
      print("line 1165")
      
      trained_model = generic_train(each_loop_length_data_feature_string,each_method,all_cases,gbmGrid,loop_type) 
      print("line 1168")
      
      saveRDS(trained_model,file =to_save_file)
      print(to_save_file)
      
      
      
    }
}


# load some naming variables into the script environment 
specify_variable_names<-function(cur_dire){
  #the_method="blindblast_just_get_alignment"
  cluster_dis="north"
  subsitution_matrix_name ="wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
  subsitution_matrix="PAM30"
  current_d=getwd()
  if(!grepl(cur_dire, current_d)){result_dir = paste(c("./",cur_dire,"/Data_processed/"),collapse="");
  overall_prefix=paste(c("./",cur_dire,"/Data_processed/"),collapse=""); plot_dir=paste(c("./",cur_dire,"/Plots/"),collapse="")
  }else{result_dir = "./Data_processed/";  plot_dir="./Plots/"; overall_prefix="./Data_processed";
 }
  mkcommand=paste(c("mkdir ",result_dir), collapse=" ")
  system(mkcommand)
  #methods= c(the_method) # specify the machine learning methods to be used 
  the_method="somerandommethod"
  each_method=the_method
  
  # end of iterating all folds
  file=paste(c(overall_prefix,"/data_by_loop_type_list_unduplicated_for_blindBLAST.rds"),collapse = "")
  data(AAPAM30)
  #AAPAM30_index=1:20
  #names(AAPAM30_index)=rownames(AAPAM30)
  system(paste(c("mkdir ",plot_dir),collapse=""))
  return(list(the_method,each_method,cluster_dis,result_dir,subsitution_matrix,subsitution_matrix_name,plot_dir,file))
}


execute_training_rscript_final_model<-function(args){
  print("at line 1230")
  loop_type=as.character(args[1])
  num_core=args[2]
  interaction.depth=args[3]
  n.trees=args[4]
  shrinkage=args[5]
  n.minobsinnode=args[6]
  registerDoMC(num_core)
  
  gbmGrid=data.frame(interaction.depth=as.numeric(interaction.depth),n.trees=as.numeric(n.trees),shrinkage=as.numeric(shrinkage),n.minobsinnode=as.numeric(n.minobsinnode))
  
  parameter_spe = paste(unlist(gbmGrid),collapse="-")
  each_method="gbm_test"
  gbm_result_dir=paste(c("./final_models/"),collapse="")
  system(paste(c("mkdir ", gbm_result_dir),collapse=""))
  
  to_save_file=paste(c(gbm_result_dir,"",loop_type,"_",paste(c(each_method,cluster_dis,parameter_spe),collapse="-"),"_final_model.rds"),collapse="")
  
  if(file.exists(to_save_file)){
    print("file_already_exists!")}else{
      print("at line 1249")
      preprocess_result=  preprocess_data(loop_type,data_by_loop_type_list_unduplicated)
      all_cases=preprocess_result[[1]]
      each_loop_length_data_feature_string=preprocess_result[[2]]
      # tune the parameter
      print("line 1254")
      trained_model = train_final_model(each_loop_length_data_feature_string,each_method,all_cases,gbmGrid) 
      saveRDS(trained_model,file =to_save_file)
      
      
      
    }
}

parallelize_jobs<-function(n,count,a_new_script,exe_sh_file,parallel_script){
  if(a_new_script){cat("\n",file=parallel_script); a_new_script=FALSE}
  cat(paste(c(" sh ",exe_sh_file," & \n"),collapse=""),file=parallel_script,append=TRUE )   # run stuff in parallel
  
  if(count==n  ){
    count=0
    cat("wait \n",file=parallel_script,append=TRUE)
    
  }
  count = count +1
  return(count)
}


# calculate dihedral angle distances
calc_dist<-function(x){
  loop_data=x
  dist_matrix=matrix(nrow=dim(loop_data)[1],ncol=dim(loop_data)[1])
  
  for(i in 1:(dim(dist_matrix)[1]-1)){
    for(j in (i+1):dim(dist_matrix)[1]){
      print(c(i,j))
      
      # used in directional statistics 
      aa=strsplit(as.character(loop_data[i,"dihedrals"]),":")[[1]]
      bb=strsplit(as.character(loop_data[j,"dihedrals"]),":")[[1]]
      n=as.numeric(as.character(loop_data[i,"length"]))
      chunk <- function(x,n) split(x, ceiling(seq_along(x)/3))
      aa_phi=as.numeric(sapply(chunk(aa,n),"[[",1))
      aa_psi=as.numeric(sapply(chunk(aa,n),"[[",2))
      bb_phi=as.numeric(sapply(chunk(bb,n),"[[",1))
      bb_psi=as.numeric(sapply(chunk(bb,n),"[[",2))
      d=sum(2*(1-cos(pi/180*(aa_phi-bb_phi))))+sum(2*(1-cos(pi/180*(aa_psi-bb_psi))))
      dist_matrix[i,j]=d
      average=acos(1-d/2)
      average_dis_matrix[i,j]=average
    }
  }
  original_dist=dist_matrix
  dist_matrix<- as.SparseSimilarityMatrix(forceSymmetric(dist_matrix))
  return(dist_matrix)
}


write_list_into_single_csv_sp<-function(table_list,output_file){
  count=1
  for(x in names(table_list)){
    if(count==1){write.table("",file=output_file)}
    write.table(x,file=output_file,append=TRUE)
    write.table(table_list[[x]],sep=",",file=output_file,col.names=FALSE,row.names  = FALSE,append=TRUE)
    count=count+1
  }
}

