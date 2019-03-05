# Copyright 2018 Xiyao Long <xlong2@jhu.edu>

#  MIT license
#  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files, to deal 
#  in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#  the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



# make template database for blindBLAST
make_reference_database <- function(member_seqs, each_fold) {
  #make a blast database
  seq_vec = apply(member_seqs[, features] , 1 , paste , collapse = "")
  names = member_seqs$identifier   # if the clustering scheme is by the torsion angles
  cluster_types = member_seqs$cluster_type
  file = paste(
    c(
      "./blast/",
      loop_type,
      each_fold,
      "some_random_file_wont_check_again.fasta"
    ),
    collapse = ""
  )
  file.remove(file)
  for (index in 1:length(names)) {
    write(paste(">", names[index], each_fold, collapse = ""), file, append = TRUE)
    write(seq_vec[[index]], file, append = TRUE)
  }
  db_name = paste(c("some_db", loop_type, the_method, cluster_dis, each_fold),
                  collapse = "_")
  arg = paste(c(' -in ', file, ' -out ', db_name, ' -dbtype prot -hash_index '),
              collapse = " ")
  command = paste("makeblastdb ", arg)
  system(command)
  return(list(db_name))
}

make_reference_database_with_f <-
  function(member_seqs, each_fold, features) {
    #make a blast database
    seq_vec = apply(member_seqs[, features] , 1 , paste , collapse = "")
    names = member_seqs$identifier   # if the clustering scheme is by the torsion angles
    cluster_types = member_seqs$cluster_type
    file = paste(
      c(
        "./blast/",
        loop_type,
        each_fold,
        "some_random_file_wont_check_again.fasta"
      ),
      collapse = ""
    )
    file.remove(file)
    for (index in 1:length(names)) {
      write(paste(">", names[index], each_fold, collapse = ""), file, append = TRUE)
      write(seq_vec[[index]], file, append = TRUE)
    }
    db_name = paste(c("some_db", loop_type, the_method, cluster_dis, each_fold),
                    collapse = "_")
    arg = paste(c(' -in ', file, ' -out ', db_name, ' -dbtype prot -hash_index '),
                collapse = " ")
    command = paste("makeblastdb ", arg)
    system(command)
    return(list(db_name))
  }


find_the_best_similarity <- function(identifier, member_seqs_pdbs) {
  print(identifier)
  if (length(member_seqs_pdbs) != 1) {
    identifier_pdb = substr(split_vector_and_replace(identifier, "\\.", 2, 2, ""), 1, 4)
    pdbs = substr(split_vector_and_replace(member_seqs_pdbs, "\\.", 2, 2, ""),
                  1,
                  4)
    simi_col_4 = substr(as.character(split_vector_and_replace(
      colnames(similarity_matrix), "\\.", 2, 2, ""
    )), 1, 4)
    index = which(simi_col_4 %in% pdbs)
    which_id = which(simi_col_4 == identifier_pdb)
    best_similarity_template_id = colnames(similarity_matrix[, index])[which(similarity_matrix[which_id, index] ==
                                                                               max(unlist(similarity_matrix[which_id, index]), na.rm = TRUE)[1])]
  } else{
    best_similarity_template_id = member_seqs_pdbs
  }
  
  return(best_similarity_template_id[1])
}

#run blast to search template for a single fold in blindBLAST
runblast_and_retrive_rmsd <-
  function(seq, each_fold, member_seqs_pdbs, db_name) {
    #the blastdbase is the directory of the database
    print(the_method)
    print(loop_type)
    print(features)
    print(cluster_dis)
    print(db_name)
    
    out_file = paste(
      c(
        "./blast/blast_his_file_need_to_parse_",
        loop_type,
        the_method,
        cluster_dis,
        each_fold,
        ".txt"
      ),
      collapse = ""
    )  # the file to record hits result
    file.remove(out_file)
    
    testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
    file = paste(
      c(
        "./blast/testing_file_wont_check_again_",
        loop_type,
        the_method,
        cluster_dis,
        each_fold,
        ".txt"
      ),
      collapse = ""
    )
    file.remove(file)
    write(paste(">", seq$identifier, collapse = ""), file, append = TRUE)
    write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
    find_flag = FALSE
    choice_index = 1
    choices = c("PAM30", "BLOSUM62")
    while (!find_flag) {
      print("inside while")
      # remove(hit)
      if (choice_index != 3) {
        sub_choice = choices[[choice_index]]
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
          find_flag = TRUE
          print(c("find flag now is ", find_flag))
          returned_id = hit[1, "V2"]
        }, error = function(e) {
        }, finally = {
          print(find_flag)
        })#this table contains all the result of sequences and and the sequence corresponded hits
        choice_index = choice_index + 1
      } else{
        identifier = seq$identifier
        returned_id = find_the_best_similarity(identifier, member_seqs_pdbs)
        if (is.character(returned_id)) {
          find_flag = TRUE
        }
        choice_index = choice_index + 1
      }
    }
    
    return(list(returned_id, choice_index - 1))
  }


runblast_and_retrive_rmsd_final <-
  function(seq, each_fold, member_seqs_pdbs, db_name) {
    #the blastdbase is the directory of the database
    print(the_method)
    print(loop_type)
    print(features)
    print(cluster_dis)
    print(db_name)
    
    out_file = paste(
      c(
        "./blast/blast_his_file_need_to_parse_",
        loop_type,
        the_method,
        cluster_dis,
        each_fold,
        ".txt"
      ),
      collapse = ""
    )  # the file to record hits result
    file.remove(out_file)
    
    testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
    file = paste(
      c(
        "./blast/testing_file_wont_check_again_",
        loop_type,
        the_method,
        cluster_dis,
        each_fold,
        ".txt"
      ),
      collapse = ""
    )
    file.remove(file)
    write(paste(">", seq$identifier, collapse = ""), file, append = TRUE)
    write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
    find_flag = FALSE
    choice_index = 1
    choices = c("PAM30", "BLOSUM62")
    while (!find_flag) {
      print("inside while")
      # remove(hit)
      if (choice_index != 3) {
        sub_choice = choices[[choice_index]]
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
          find_flag = TRUE
          print(c("find flag now is ", find_flag))
          returned_id = hit[1, "V2"]
        }, error = function(e) {
        }, finally = {
          print(find_flag)
        })#this table contains all the result of sequences and and the sequence corresponded hits
        choice_index = choice_index + 1
      } else{
        identifier = seq$identifier
        returned_id = NA
        find_flag = TRUE
        
        choice_index = choice_index + 1
      }
    }
    
    return(list(returned_id, choice_index - 1))
  }




# make data division for k folds and r repeats croo validation
make_3_10_cross_val <- function(training_cases, r, k) {
  all_unique_ids = lapply(split(training_cases, training_cases$cluster_type), function(x) {
    unique(x$id)
  })
  
  split_data = split(training_cases, training_cases$cluster_type)
  new_training_data = do.call(rbind, lapply(names(split_data), function(x) {
    for (repeat_n in 1:r) {
      all_unique_id_sample = lapply(all_unique_ids, function(y) {
        if (length(y) >= k) {
          sample(1:k, size = length(y), replace = TRUE)
        } else{
          sample(1:k, size = length(y), replace = FALSE)
        }
      })
      
      for (ind in 1:length(all_unique_ids[[x]])) {
        all_unique_ids[[x]]
        repnum = dim(split_data[[x]][split_data[[x]]$id == all_unique_ids[[x]][ind], ])[1]
        split_data[[x]][split_data[[x]]$id == all_unique_ids[[x]][ind], paste(c("fold.num", repeat_n), collapse =
                                                                                "")] = rep(all_unique_id_sample[[x]][ind], repnum)
      }
    }
    return(split_data[[x]])
  }))
  
  
  # Create folds and repeats here - you could create your own if you want #
  
  
  folds.list.out <- list()
  folds.list <- list()
  list.counter <- 1
  for (y in 1:r) {
    newcol <- paste('fold.num', y, sep = '')
    for (z in 1:k) {
      out_rown = which(new_training_data[, newcol] == z)  # find the folds in
      folds_in = which(new_training_data[, newcol] != z)
      a = new_training_data[folds_in, "id"] %in%  new_training_data[out_rown, "id"]
      print(a[a])
      
      sub = new_training_data[out_rown, "id"]  # ids corresponding to the fold
      out_rown = out_rown[which(!duplicated(sub))]
      a = new_training_data[folds_in, "id"] %in%  new_training_data[out_rown, "id"]
      print(a[a])
      if (length(a[a]) != 0) {
        print("The 1fold in and 9 fold out is not correct, check!")
        stop()
      }
      folds.list.out[[list.counter]] <- out_rown
      folds.list[[list.counter]] <-
        which(new_training_data[, newcol] != z)
      list.counter <- list.counter + 1
    }
  }
  index_list = list()
  index_list[[1]] = folds.list.out
  index_list[[2]] = folds.list
  returned_list = list()
  returned_list[[1]] = index_list
  
  returned_list[[2]] = new_training_data
  return(returned_list)
}


#Retrieve the rmsd between the predicted sequence and query sequence by their ids. Record it
get_rmsd_by_ids <- function(id1, id2) {
  pdb1 = substr(strsplit(as.character(id1), "\\.")[[1]][2], 1, 4)
  pdb2 = substr(strsplit(as.character(id2), "\\.")[[1]][2], 1, 4)
  rmsd = rmsd_matrix[pdb1, pdb2]
  return(rmsd)
}

#
get_accuracy_per_fold <- function(each_fold) {
  # get reference database
  print(c("the each_fold is ", each_fold))
  #print(sequences)
  #print(folds.list.out)
  print(each_fold)
  member_seqs = sequences[folds.list[[each_fold]], ]
  returned_db = make_reference_database(member_seqs, each_fold) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
  member_seqs_pdbs = (member_seqs$identifier)
  fold_out_cases = sequences[folds.list.out[[each_fold]], ]
  rmsd_list = data.frame(matrix(nrow = dim(fold_out_cases)[1], ncol = 4))
  system("mkdir ./blast/")
  for (each_ind in 1:dim(fold_out_cases)[1]) {
    print(each_ind)
    #tryCatch({
    seq = fold_out_cases[each_ind, ]
    case_id = seq["identifier"]
    print(seq)
    found_template = runblast_and_retrive_rmsd(seq, each_fold, member_seqs_pdbs, returned_db[[1]])
    query_template_rmsd = get_rmsd_by_ids(as.character(case_id), found_template[[1]])
    system("mkdir ./garbage ")
    system(paste(c("mv ", returned_db[[1]], "* " , "./garbage/"), collapse =
                   ""))
    
    rmsd_list[each_ind, 1:4] = c(
      as.character(case_id),
      as.character(found_template[[1]]),
      query_template_rmsd,
      found_template[[2]]
    )
  }# end of iterating through the fold out cases for a single fold
  return(rmsd_list)
}






get_accuracy_per_fold_overload_final <- function(each_fold) {
  # get reference database
  
  print(c("the each_fold is ", each_fold))
  #print(folds.list.out)
  print(each_fold)
  member_seqs = sequences[folds.list[[each_fold]], ]
  fold_out_cases = sequences[folds.list.out[[each_fold]], ]
  returned_db = make_reference_database_with_f(member_seqs, each_fold, features) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
  member_seqs_pdbs = (member_seqs$identifier)
  rmsd_list = data.frame(matrix(nrow = dim(fold_out_cases)[1], ncol = 3))
  
  for (each_ind in 1:dim(fold_out_cases)[1]) {
    #tryCatch({
    seq = fold_out_cases[each_ind, ]
    case_id = seq["identifier"]
    
    found_template = runblast_and_retrive_rmsd_final(seq, each_fold, member_seqs_pdbs, returned_db[[1]])
    
    rmsd_list[each_ind, 1:3] = c(as.character(case_id),
                                 as.character(found_template[[1]]),
                                 found_template[[2]])
  }# end of iterating through the fold out cases for a single fold
  
  system(paste(c("mv ", returned_db[[1]], "* " , "./garbage/"), collapse =
                 ""))
  
  
  return(rmsd_list)
}









runblast_and_retrive_bitscore <-
  function(seq, each_fold, member_seqs_pdbs, db_name) {
    #the blastdbase is the directory of the database
    print(the_method)
    print(loop_type)
    print(features)
    print(cluster_dis)
    print(db_name)
    
    out_file = paste(
      c(
        "./blast/blast_his_file_need_to_parse_",
        loop_type,
        the_method,
        cluster_dis,
        each_fold,
        ".txt"
      ),
      collapse = ""
    )  # the file to record hits result
    file.remove(out_file)
    
    testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
    file = paste(
      c(
        "./blast/testing_file_wont_check_again_",
        loop_type,
        the_method,
        cluster_dis,
        each_fold,
        ".txt"
      ),
      collapse = ""
    )
    file.remove(file)
    write(paste(">", seq$identifier, collapse = ""), file, append = TRUE)
    write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
    find_flag = FALSE
    choice_index = 1
    choices = c("PAM30", "BLOSUM62")
    while (!find_flag) {
      print("inside while")
      # remove(hit)
      if (choice_index != 3) {
        sub_choice = choices[[choice_index]]
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
          find_flag = TRUE
          print(c("find flag now is ", find_flag))
          # retrive 5 correct,
          query_c = split_vector_and_replace(hit$V1, "\\.", 1, 1, "")
          template_c = split_vector_and_replace(hit$V2, "\\.", 1, 1, "")
          correct_ns = which(query_c == template_c)
          incorrect_ns = which(query_c != template_c)
          if (length(correct_ns) >= 5) {
            correct_ns = correct_ns[1:5]
          } else{
            correct_ns = correct_ns[1:length(correct_ns)]
          }
          if (length(incorrect_ns) >= 5) {
            incorrect_ns = incorrect_ns[1:5]
          } else{
            incorrect_ns = incorrect_ns[1:length(incorrect_ns)]
          }
          correct_bit_scores = hit[correct_ns, "V12"]
          incorrect_bit_scores = hit[incorrect_ns, "V12"]
          # find the shorter one
          which_shorter = which.min(c(
            length(correct_bit_scores),
            length(incorrect_bit_scores)
          ))
          min_len = min(c(
            length(correct_bit_scores),
            length(incorrect_bit_scores)
          ))
          max_len = max(c(
            length(correct_bit_scores),
            length(incorrect_bit_scores)
          ))
          
          diff_r = c()
          for (i in 1:max_len) {
            #when j is not where I'm
            if (i <= min_len) {
              j = i
            } else{
              j = min_len
            }
            if (which_shorter == 1) {
              diff_r = c(diff_r,
                         correct_bit_scores[j] - incorrect_bit_scores[i])
            } else{
              diff_r = c(diff_r,
                         correct_bit_scores[i] - incorrect_bit_scores[j])
            }
            
            
          }
          return(diff_r)
        }, error = function(e) {
        }, finally = {
          print(find_flag)
        })#this table contains all the result of sequences and and the sequence corresponded hits
        choice_index = choice_index + 1
      } else{
        return(NULL)
      }
    }
    
    
  }





runblast_for_unmatched_random_and_retrive_bitscore <-
  function(seq, each_fold, member_seqs_pdbs, db_name) {
    #the blastdbase is the directory of the database
    print(the_method)
    print(loop_type)
    print(features)
    print(cluster_dis)
    print(db_name)
    
    out_file = paste(
      c(
        "./blast/blast_his_file_need_to_parse_",
        loop_type,
        the_method,
        cluster_dis,
        each_fold,
        ".txt"
      ),
      collapse = ""
    )  # the file to record hits result
    file.remove(out_file)
    
    testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
    file = paste(
      c(
        "./blast/testing_file_wont_check_again_",
        loop_type,
        the_method,
        cluster_dis,
        each_fold,
        ".txt"
      ),
      collapse = ""
    )
    file.remove(file)
    write(paste(">", seq$identifier, collapse = ""), file, append = TRUE)
    write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
    find_flag = FALSE
    choice_index = 1
    choices = c("PAM30", "BLOSUM62")
    while (!find_flag) {
      print("inside while")
      # remove(hit)
      if (choice_index != 3) {
        sub_choice = choices[[choice_index]]
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
          find_flag = TRUE
          print(c("find flag now is ", find_flag))
          # retrive 5 correct,
          query_c = split_vector_and_replace(hit$V1, "\\.", 1, 1, "")
          template_c = split_vector_and_replace(hit$V2, "\\.", 1, 1, "")
          correct_ns = which(query_c == template_c)
          incorrect_ns = which(query_c != template_c)
          if (length(incorrect_ns) >= 5) {
            t_n = 5
          } else{
            t_n = length(incorrect_ns)
          }
          sampled_incorrect_ns = sample(incorrect_ns, t_n)
          incorrect_bit_scores_sampled = hit[sampled_incorrect_ns, "V12"]
          # find the shorter one
          std_v = c(incorrect_bit_scores_sampled)
          
          return(std_v)
        }, error = function(e) {
        }, finally = {
          print(find_flag)
        })#this table contains all the result of sequences and and the sequence corresponded hits
        choice_index = choice_index + 1
      } else{
        return(NULL)
      }
    }
    
    
  }




get_bit_scores <- function(each_fold) {
  # get reference database
  print(c("the each_fold is ", each_fold))
  print(each_fold)
  member_seqs = sequences[folds.list[[each_fold]], ]
  returned_db = make_reference_database(member_seqs, each_fold) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
  member_seqs_pdbs = (member_seqs$identifier)
  fold_out_cases = sequences[folds.list.out[[each_fold]], ]
  rmsd_list = data.frame(matrix(nrow = dim(fold_out_cases)[1], ncol = 6))
  if (dim(rmsd_list)[1] < 1) {
    return(NULL)
  }
  for (each_ind in 1:dim(fold_out_cases)[1]) {
    print(each_ind)
    #tryCatch({
    seq = fold_out_cases[each_ind, ]
    case_id = seq["identifier"]
    print(seq)
    bit_scores_1_to_5 = runblast_and_retrive_bitscore(seq, each_fold, member_seqs_pdbs, returned_db[[1]])
    rmsd_list[each_ind, 1:(length(bit_scores_1_to_5) + 1)] = c(as.character(case_id), unlist(bit_scores_1_to_5))
  }# end of iterating through the fold out cases for a single fold
  return(rmsd_list)
}






get_random_bit_scores <- function(each_fold) {
  # get reference database
  print(c("the each_fold is ", each_fold))
  print(each_fold)
  member_seqs = sequences[folds.list[[each_fold]], ]
  returned_db = make_reference_database(member_seqs, each_fold) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
  member_seqs_pdbs = (member_seqs$identifier)
  fold_out_cases = sequences[folds.list.out[[each_fold]], ]
  rmsd_list = data.frame(matrix(nrow = dim(fold_out_cases)[1], ncol = 6))
  if (dim(rmsd_list)[1] < 1) {
    return(NULL)
  }
  
  for (each_ind in 1:dim(fold_out_cases)[1]) {
    print(each_ind)
    #tryCatch({
    seq = fold_out_cases[each_ind, ]
    case_id = seq["identifier"]
    print(seq)
    random_bit_scores_1_5 = runblast_for_unmatched_random_and_retrive_bitscore(seq, each_fold, member_seqs_pdbs, returned_db[[1]])
    rmsd_list[each_ind, 1:(length(random_bit_scores_1_5) + 1)] = c(as.character(case_id), random_bit_scores_1_5)
  }# end of iterating through the fold out cases for a single fold
  return(rmsd_list)
}


# an overriden version if the member_seqs and fold_out_cases are given by
get_accuracy_per_fold_overload <-
  function(each_fold,
           member_seqs,
           fold_out_cases,
           features) {
    # get reference database
    
    print(c("the each_fold is ", each_fold))
    print(sequences)
    print(member_seqs)
    #print(folds.list.out)
    print(each_fold)
    returned_db = make_reference_database_with_f(member_seqs, each_fold, features) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
    member_seqs_pdbs = (member_seqs$identifier)
    rmsd_list = data.frame(matrix(nrow = dim(fold_out_cases)[1], ncol = 4))
    
    for (each_ind in 1:dim(fold_out_cases)[1]) {
      #tryCatch({
      seq = fold_out_cases[each_ind, ]
      case_id = seq["identifier"]
      
      found_template = runblast_and_retrive_rmsd(seq, each_fold, member_seqs_pdbs, returned_db[[1]])
      query_template_rmsd = get_rmsd_by_ids(as.character(case_id), found_template[[1]])
      
      rmsd_list[each_ind, 1:4] = c(
        as.character(case_id),
        as.character(found_template[[1]]),
        query_template_rmsd,
        found_template[[2]]
      )
    }# end of iterating through the fold out cases for a single fold
    
    system(paste(c("mv ", returned_db[[1]], "* " , "./garbage/"), collapse =
                   ""))
    
    
    return(rmsd_list)
  }


get_accuracy_per_fold_enforcing_corrent_fold <- function(each_fold) {
  # get reference database
  print(c("the each_fold is ", each_fold))
  print(sequences)
  #print(folds.list.out)
  print(each_fold)
  
  theresult_this_fold = the_result_list[[each_fold]]
  
  in_fold_member_seqs = sequences[!sequences$identifier %in% theresult_this_fold$X1, ]
  theresult_this_fold$query_cluster = sapply(strsplit(theresult_this_fold$X1, "\\."), "[[", 1)
  theresult_this_fold$template_cluster = sapply(strsplit(theresult_this_fold$X2, "\\."), "[[", 1)
  unmatched = which(theresult_this_fold$query_cluster != theresult_this_fold$template_cluster)
  if (length(unmatched) == 0) {
    return(NULL)
  } else{
    unmatched_cases = theresult_this_fold[unmatched, ]
    matched_cases = theresult_this_fold[!unmatched, ]
    fold_out_cases = sequences[sequences$identifier %in% unmatched_cases$X1, ]
    #split the wrong fold out cases by query clusters
    splited_fold_out_cases = split(fold_out_cases,
                                   as.character(fold_out_cases$cluster_type))
    accuracy_per_fold_lists = list()
    for (each_cluster in names(splited_fold_out_cases)) {
      this_cluster_fold_out_cases = splited_fold_out_cases[[each_cluster]]
      
      member_seqs = in_fold_member_seqs[in_fold_member_seqs$cluster_type %in% each_cluster, ]
      if (dim(member_seqs)[1] != 0) {
        accuracy_per_fold_lists[[each_cluster]] = get_accuracy_per_fold_overload(each_fold,
                                                                                 member_seqs,
                                                                                 this_cluster_fold_out_cases,
                                                                                 features)
      }
    }
    rmsd_list = do.call(rbind, accuracy_per_fold_lists)
    colnames(rmsd_list)[3:4] = c("correct_cluster_rmsd", "enforcing_method_choice")
    colnames(unmatched_cases)[3:4] = c("rmsd", "method_choice")
    merged_result = merge(unmatched_cases[, 1:4],
                          rmsd_list,
                          by = ("X1"),
                          all.x = TRUE)
    return(merged_result)
  }
}




runblast_and_retrive_similarity <-
  function(returned_db, seq, features) {
    blastdbase_name = returned_db[[1]]  #the blastdbase is the directory of the database
    out_file = paste(
      c(
        "./blast/blast_his_file_need_to_parse_",
        loop_type,
        the_method,
        cluster_dis,
        "all",
        ".txt"
      ),
      collapse = ""
    )  # the file to record hits result
    testing_seq_vec = apply(seq[, features] , 1 , paste , collapse = "")  # turn the query dataframe into strings
    file = paste(
      c(
        "./blast/testing_file_wont_check_again_",
        loop_type,
        the_method,
        cluster_dis,
        "all",
        ".txt"
      ),
      collapse = ""
    )
    file.remove(file)
    write(paste(">", seq$identifier, collapse = ""), file, append = TRUE)
    write(testing_seq_vec[[1]], file, append = TRUE)      #write the sequence of the predicted
    if (grepl("PAM", subsitution_matrix) |
        grepl("BLO", subsitution_matrix)) {
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
    if (dim(hit)[1] == 0) {
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
        
      }, error = function(e) {
      })
    } else{
      print("do nothing")
    }
    return(hit)
  }


make_database <- function(member_seqs) {
  #make a blast database
  seq_vec = apply(member_seqs[, features] , 1 , paste , collapse = "")
  names = member_seqs$identifier   # if the clustering scheme is by the torsion angles
  cluster_types = member_seqs$cluster_type
  file = paste(
    c(
      "./blast/",
      loop_type,
      "all",
      "some_random_file_wont_check_again.fasta"
    ),
    collapse = ""
  )
  file.remove(file)
  for (index in 1:length(names)) {
    write(paste(">", names[index], "all", collapse = ""), file, append = TRUE)
    write(seq_vec[[index]], file, append = TRUE)
  }
  db_name = paste(c("some_db", loop_type, the_method, cluster_dis, "all"),
                  collapse = "_")
  arg = paste(c(' -in ', file, ' -out ', db_name, ' -dbtype prot -hash_index '),
              collapse = " ")
  command = paste("makeblastdb ", arg)
  system(command)
  return(list(db_name))
}



get_pairwise_sequence_simi_self_cal <- function(loop_type) {
  sequences = data_by_loop_type_list_unduplicated[[loop_type]][[1]]#load the sequences from a loop length
  features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
  case_ids = sequences[, "identifier"]
  similarity_matrix = matrix(nrow = dim(sequences)[1], ncol = dim(sequences)[1])
  colnames(similarity_matrix) = case_ids
  rownames(similarity_matrix) = case_ids
  
  
  for (each_ind in 1:(length(case_ids) - 1)) {
    print(each_ind)
    a = c()
    for (j in (each_ind + 1):length(case_ids)) {
      a_seq = sequences[each_ind, features]
      b_seq = sequences[j, features]
      seq = rbind(a_seq, b_seq)
      similarity = sum(unlist(lapply(seq, function(x) {
        AAPAM30[AAPAM30_index[[x[1]]], AAPAM30_index[[x[2]]]]
      })))
      
      similarity_matrix[each_ind, j] = similarity
    }
    print(similarity_matrix[each_ind, ])
    
  }
  #forceSymmetric(similarity_matrix)
  return(similarity_matrix)
}


get_pairwise_sequence_simi <- function(sequences, features) {
  case_ids = sequences[, "identifier"]
  similarity_matrix = matrix(nrow = dim(sequences)[1], ncol = dim(sequences)[1])
  colnames(similarity_matrix) = case_ids
  rownames(similarity_matrix) = case_ids
  
  for (each_ind in 1:length(case_ids)) {
    tryCatch({
      member_seqs = sequences[-each_ind, ]
      returned_db = make_database(member_seqs, features) #construct a blast database with the predicted cluster and find the best three hits , extract their rmsd
      seq = sequences[each_ind, ]
      hit_table = runblast_and_retrive_similarity(returned_db, seq, features)
      for (x in 1:dim(hit_table)[1]) {
        similarity_matrix[each_ind, colnames(similarity_matrix) == hit_table[x, 2]] =
          hit_table[x, 12]
      } # finish recording all the available similarity values captured
      returned_db = ""
    }, error = function(e) {
    })
    print(each_ind / length(case_ids))
  }
  return(similarity_matrix)
}


sig_info_function <- function(y) {
  x = all_loop_sig_frame_list[[y]]
  x$loop = rep(y, dim(x)[1])
  # if the first column has the cluster one and the two clusters are not the ssame
  H1_related = x[x[, "Var1"] != x[, "Var2"], ]
  # most prevalent cluster
  each_l = paste(unlist(strsplit(y, "_")[[1]]), collapse = "-")
  if (each_l == "L3-9") {
    most_p = "L3-9-cis7-1$"
  } else{
    most_p = paste(c(each_l, "-1$"), collapse = "")
  }
  H1_cluster_one_recovery_fail = H1_related[grepl(most_p, H1_related$Var1), ]
  H1_cluster_one_precision_fail = H1_related[grepl(most_p, H1_related$Var2), ]
  H1_cluster_one_recovery_fail$error_type = rep("1_reco", dim(H1_cluster_one_recovery_fail)[1])
  H1_cluster_one_precision_fail$error_type = rep("1_prec", dim(H1_cluster_one_precision_fail)[1])
  other_types = H1_related[!(grepl(most_p, H1_related$Var1) |
                               grepl(most_p, H1_related$Var2)), ]
  other_types$error_type = rep("non_1", dim(other_types)[1])
  all_annotated = rbind(
    rbind(other_types, H1_cluster_one_precision_fail),
    H1_cluster_one_recovery_fail
  )
  return(all_annotated)
  # if the
  
}

another_sig_info_function <- function(y) {
  x = all_loop_sig_frame_list[[y]]
  x$loop = rep(y, dim(x)[1])
  # if the first column has the cluster one and the two clusters are not the ssame
  H1_related = x[x[, "Var1"] == x[, "Var2"], ]
  # most prevalent cluster
  each_l = paste(unlist(strsplit(y, "_")[[1]]), collapse = "-")
  if (each_l == "L3-9") {
    most_p = "L3-9-cis7-1$"
  } else{
    most_p = paste(c(each_l, "-1$"), collapse = "")
  }
  H1_cluster_one_recovery_fail = H1_related[grepl(most_p, H1_related$Var1), ]
  H1_cluster_one_recovery_fail$error_type = rep("cluster_1", dim(H1_cluster_one_recovery_fail)[1])
  other_types = H1_related[!(grepl(most_p, H1_related$Var1)) , ]
  other_types$error_type = rep("non_1", dim(other_types)[1])
  all_annotated = rbind(other_types, H1_cluster_one_recovery_fail)
  return(all_annotated)
  # if the
}


calculate_accuracy_mean_std <- function(all_result_list) {
  all_result = do.call(rbind, all_result_list)
  all_result$query_cluster = sapply(strsplit(all_result[, 1], "\\."), "[[", 1)
  all_result$template_cluster = sapply(strsplit(all_result[, 2], "\\."), "[[", 1)
  accuracy = dim(all_result[all_result$template_cluster == all_result$query_cluster, ])[1] /
    dim(all_result)[1]
  
  #std
  all_result_list_split_by_repeat = chunk2(all_result_list, 3)
  
  
  accuracies_all_repeats = lapply(all_result_list_split_by_repeat, function(x) {
    x = do.call(rbind, x)
    x$query_cluster = sapply(strsplit(x[, 1], "\\."), "[[", 1)
    x$template_cluster = sapply(strsplit(x[, 2], "\\."), "[[", 1)
    accuracy = dim(x[x$template_cluster == x$query_cluster, ])[1] / dim(x)[1]
    
    print(accuracy)
  })
  std = sd(unlist(accuracies_all_repeats))
  accuracy = mean(unlist(accuracies_all_repeats))
  return(list(accuracy, std))
}


# calculate the accuracies by repeats
calculate_accuracy_mean_std_by_repeats <-
  function(all_result_list, repeat_n) {
    all_result_list_rbind = lapply(all_result_list, function(model_re) {
      model_re$repeats = ceiling(as.numeric(gsub("Resample", "", model_re[, "Resample"])) /
                                   10)
      model_re$Resample = as.numeric(gsub("Resample", "", model_re[, "Resample"]))
      model_re = model_re[, c("pred", "obs", "Resample", "repeats")]
      
      last_n = model_re$Resample[length(model_re$Resample)]
      model_re_by_resample = split(model_re, model_re$Resample)
      model_re_by_repeats = model_re
      loopa = paste(strsplit(model_re[1, 1], "_")[[1]][1:2], collapse = "_")
      lista = list()
      if (loopa %in% c("H1_13", "H2_10", "L3_9")) {
        
      } else{
        resample_chunks = chunk2(1:length(model_re_by_resample), 3)
        rbind_r = lapply(resample_chunks, function(x) {
          print(x)
          do.call(rbind, model_re_by_resample[x])
        })
        
      }
      if (length(rbind_r) != 3) {
        stop()
      }
      return(rbind_r)
    })
    lapply(1:3, function(x) {
      frames = lapply(all_result_list_rbind, function(y) {
        y[[x]]
      })
      b_f = do.call(rbind, frames)
      dim(b_f[b_f[, 1] == b_f[, 2], ])[1] / dim(b_f)[1]
    })
    
    all_result = do.call(rbind, all_result_list_rbind)
    all_result$query_cluster = sapply(strsplit(as.character(all_result[, 1]), "\\."), "[[", 1)
    all_result$template_cluster = sapply(strsplit(as.character(all_result[, 2]), "\\."), "[[", 1)
    split_by_repeats = split(all_result, all_result$repeats)
    accuracies_by_repeat = lapply(split_by_repeats, function(x) {
      accuracy = dim(x[x$template_cluster == x$query_cluster, ])[1] / dim(x)[1]
    })
    sds = sd(unlist(accuracies_by_repeat))
    
    result_list = list()
    result_list[[1]] = accuracies_by_repeat
    result_list[[2]] = sds
    return(result_list)
  }





generate_folds_foldsout <- function(sequences, r, k) {
  training_cases = sequences
  training_cases = sequences[, c(features, "cluster_type", "id", "identifier")]
  returned_results = make_3_10_cross_val(training_cases, r, k)
  folds_spec = returned_results[[1]]
  sequences = returned_results[[2]]   # the training cases would have its mo
  folds.list.out = folds_spec[[1]]
  folds.list = folds_spec[[2]]
  counted_folds = c(1:length(folds.list))[unlist(lapply(folds.list.out, length)) !=
                                            0]
  folds.list = folds.list[counted_folds]
  folds.list.out = folds.list.out[counted_folds]
  if (any(sequences[folds.list[[1]], "id"] %in% sequences[folds.list.out[[1]], "id"])) {
    stop()
  }
  return(list(folds.list, folds.list.out))
}

getting_similar_sequences_similarity_matrix_rmsd_matrix <-
  function(loop_type) {
    pass = FALSE
    tryCatch({
      sequences = data_by_loop_type_list_unduplicated[[loop_type]][[1]]#load the sequences from a loop length
      features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
      loop_type_type = strsplit(loop_type, "_")[[1]][1]
      similarity_matrix = forceSymmetric(all_similarity_matrix[[loop_type]])
      rmsd_matrix = forceSymmetric(as.matrix(all_rmsd_list[[loop_type_type]]))
      
      
      seq_pdbs = substr(sapply(strsplit(sequences$identifier, "\\."), "[[", 2), 1, 4)
      similarity_matrix_pdbs = substr(sapply(strsplit(colnames(
        similarity_matrix
      ), "\\."), "[[", 2), 1, 4)
      rmsd_matrix_pdbs = colnames(rmsd_matrix)
      
      common_pdbs = Reduce(intersect,
                           list(seq_pdbs, similarity_matrix_pdbs, rmsd_matrix_pdbs))
      sequences = sequences[seq_pdbs %in% common_pdbs, ]
      which_i = which(similarity_matrix_pdbs %in% common_pdbs)
      similarity_matrix = similarity_matrix[which_i, which_i]
      rmsd_matrix = rmsd_matrix[rmsd_matrix_pdbs %in% common_pdbs, rmsd_matrix_pdbs %in% common_pdbs]
      pass = TRUE
    }, error = function(e) {
    })
    if (!pass) {
      return(FALSE)
    } else
      (
        return(list(sequences, similarity_matrix, rmsd_matrix))
      )
  }




# get template by similarity score

get_template_by_similarity_score <- function(cluster) {
  ind_cl = which(sequences$cluster_type == cluster)
  this_ind_most_sim_ind_list = list()
  for (this_ind in ind_cl) {
    # print(this_ind)
    relevant_indexes = (1:dim(sequences)[1])[(1:dim(sequences)[1]) != this_ind]
    the_max = max(unlist(this_simi[this_ind, ]), na.rm = TRUE)
    this_ind_most_sim_ind = which(unlist(this_simi[this_ind, ]) == the_max)[1]
    # if(length(unlist(this_ind_most_sim_ind))>1){
    #    stop()
    #  }
    print(as.data.frame(c(this_ind, this_ind_most_sim_ind)))
    this_ind_most_sim_ind_list = c(this_ind_most_sim_ind_list, this_ind_most_sim_ind)
  }
  final_table = data.frame(sequences[ind_cl, "cluster_type"], sequences[unlist(this_ind_most_sim_ind_list), "cluster_type"])
  match_result = data.frame(query = ind_cl,
                            template = unlist(this_ind_most_sim_ind_list))
  tem = cbind(match_result, final_table)
  tem$out_of_cluster_similarity = rep(NA, dim(tem)[1])
  tem$within_cluster_similarity = rep(NA, dim(tem)[1])
  all_disin_query_cluster = sequences[ind_cl, "dis"]
  percentile <- ecdf(all_disin_query_cluster)
  
  for (x in 1:dim(tem)[1]) {
    query_id = tem[x, "query"]
    template_id = tem[x, 2]
    query_cluster = tem[x, 3]
    template_cluster = tem[x, 4]
    temcluster_ids = which(sequences$cluster_type == template_cluster)
    template_dis_s = sequences[temcluster_ids, "dis"]
    temp_percentile = ecdf(template_dis_s)
    non_querycluster_cluster = which(sequences$cluster_type != query_cluster)
    print(c("query_id ", query_id))
    within_querycluster_ind = ind_cl[ind_cl != query_id]
    within_most_sim = max(unlist(this_simi[query_id, within_querycluster_ind]), na.rm =
                            TRUE)
    query_to_querycluster_dis = sequences[query_id, "dis"]
    #query_dis_percentile = percentile(query_to_querycluster_dis)
    template_to_tempcluster_dis = sequences[template_id, "dis"]
    query_dis_percentile = percentile(query_to_querycluster_dis)
    temp_dis_percentile = temp_percentile(template_to_tempcluster_dis)
    out_most_sim = max(unlist(this_simi[query_id, non_querycluster_cluster]), na.rm =
                         TRUE)
    #number_nei=which()
    tem[x, "out_of_cluster_similarity"] = out_most_sim
    tem[x, "within_cluster_similarity"] = within_most_sim
    tem[x, "query_to_cluster_distance"] = query_to_querycluster_dis
    tem[x, "query_to_clustercen"] = query_dis_percentile
    tem[x, "template_to_clustercen"] = temp_dis_percentile
    tem[x, "number_neighborhood_structure"] = length(which(unlist(this_dis[tem[x, "query"], ]) <
                                                             all_neighborhood_dist[[each_l]]))
    
  }
  
  return(list(tem, final_table))
}# end of iterating all clusters in this loop type



# get the number of wrong cases for each misclassification in blindBLAST, averaged by the number of repeat
count_number_for_misclassification <- function(each_l) {
  this_loop_wrong_cases = enforcing_correct_rmsd_list[[each_l]]
  wrong_cases = do.call(rbind, this_loop_wrong_cases)
  if (is.null(wrong_cases)) {
    return(NULL)
  }
  wrong_cases$template_c = sapply(strsplit(wrong_cases$X1, "\\."), "[[", 1)
  wrong_cases$query_c = sapply(strsplit(wrong_cases$X2.x, "\\."), "[[", 1)
  conf_t = as.data.frame(table(wrong_cases$template_c, wrong_cases$query_c))
  conf_t = conf_t[conf_t$Freq != 0, ]
  conf_t$Freq = conf_t$Freq / 3
  conf_t$misclassification = paste(conf_t$Var1, conf_t$Var2, sep = "")
  conf_t$Var1 = as.character(conf_t$Var1)
  conf_t$Var2 = as.character(conf_t$Var2)
  return(conf_t)
}



#use the all_significance_simulation from the envrionment variables, and calculate the significance of every misclassification observed in blindBLAST
get_significance <- function(this_loop_wrong_cases) {
  sig_data_frame = data.frame(matrix(
    nrow = length(this_loop_wrong_cases$misclassification),
    ncol = 8
  ))
  colnames(sig_data_frame) = c(
    "error_type",
    "Var1",
    "Var2",
    "significance",
    "error_count",
    "sd",
    "mean_simu_error",
    "effect_size"
  )
  rownames(sig_data_frame) = this_loop_wrong_cases$misclassification
  for (misclassification in this_loop_wrong_cases$misclassification) {
    # start iterating over the error types in this loop to assess its significance
    print(misclassification %in% names(all_significance_simulation[[each_l]]))
    if (!misclassification %in% names(all_significance_simulation[[each_l]])) {
      next
    }
    this_mis_info = this_loop_wrong_cases[this_loop_wrong_cases$misclassification ==
                                            misclassification, ]
    this_error = this_mis_info["Freq"]
    simu_error = all_significance_simulation[[each_l]][[misclassification]]
    aycdf <- ecdf(simu_error)
    sig_data_frame[misclassification, "error_type"] = misclassification
    sig_data_frame[misclassification, "Var1"] = this_mis_info["Var1"]
    sig_data_frame[misclassification, "Var2"] = this_mis_info["Var2"]
    sig_data_frame[misclassification, "significance"] = aycdf(this_error)
    sig_data_frame[misclassification, "error_count"] = this_error
    sig_data_frame[misclassification, "sd"] = sd(simu_error)
    sig_data_frame[misclassification, "mean_simu_error"] = mean(simu_error)
    sig_data_frame[misclassification, "effect_size"] = (this_error - mean(simu_error)) /
      sd(simu_error)
    
  }# end of iterating over all errors in this this loop type
  #sig_data_frame=sig_data_frame[sig_data_frame$mean_simu_error>1 | (sig_data_frame$significance>=0.975|sig_data_frame$significance<=0.025),]
  return(sig_data_frame)
  
}


frame_manipulating_for_ploting <- function(accuracy_list) {
  accuracy_gbm_blast = as.data.frame(accuracy_list)
  accuracy_gbm_blast_remove_unknow = accuracy_gbm_blast[complete.cases(accuracy_gbm_blast),]
  accuracy_gbm_blast_remove_unknow$loop = split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow), "_", 1, 1, "-")
  accuracy_gbm_blast_remove_unknow$length = as.numeric(split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow), "_", 2, 2, "-"))
  accuracy_gbm_blast_remove_unknow = reorder_factor(accuracy_gbm_blast_remove_unknow, "loop", "length")
  
  accuracy_gbm_blast_remove_unknow_melt = melt(accuracy_gbm_blast_remove_unknow,
                                               id.vars = c("loop", "length"))
  accuracy_gbm_blast_remove_unknow_melt = as.data.frame(accuracy_gbm_blast_remove_unknow_melt)
  return(accuracy_gbm_blast_remove_unknow_melt)
}

frame_for_plot <- function(method_n) {
  err_frame = error_count_list_with_sd[[method_n]]
  
  accuracy_gbm_blast_remove_unknow = err_frame[complete.cases(err_frame),]
  accuracy_gbm_blast_remove_unknow$loop = split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow), "_", 1, 1, "-")
  accuracy_gbm_blast_remove_unknow$length = as.numeric(split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow), "_", 2, 2, "-"))
  accuracy_gbm_blast_remove_unknow = reorder_factor(accuracy_gbm_blast_remove_unknow, "loop", "length")
  accuracy_gbm_blast_remove_unknow$methods = rep(method_n, dim(accuracy_gbm_blast_remove_unknow)[1])
  return(accuracy_gbm_blast_remove_unknow)
}


query_to_cluster_center_dis <- function(x) {
  query_c = as.character(x["query_c"])
  print(c("query_c", query_c))
  loop_t = paste(strsplit(as.character(query_c), "-")[[1]][1:2], collapse =
                   "_")
  data = data_by_loop_type_list_unduplicated[[loop_t]][[1]]
  center = which(data$center == 1 &
                   as.character(data$cluster_type) == query_c)
  if (length(center) == 1) {
    cluster_center_id = substr(strsplit(data[data$center == 1 &
                                               as.character(data$cluster_type) == query_c, "identifier"], "\\.")[[1]][2], 1, 4)
    
    query_id = substr(strsplit(as.character(x["X1"]), "\\.")[[1]][2], 1, 4)
    
    dis = all_dist_matrix[[loop_t]]
    this_dis = dis[rownames(dis) == cluster_center_id, colnames(dis) == query_id]
  } else{
    this_dis = NA
  }
  return(this_dis[1])
  # get the dihedral angle distance
  
}


#get number of neighboring structures
query_structure_number_neighbor <- function(x) {
  query_c = as.character(unlist(x[["query_c"]]))
  this_pdb = substr(tolower(strsplit(x[["X1"]], "\\.")[[1]][2]), 1, 4)
  print(c("query_c", query_c))
  loop_t = paste(strsplit(as.character(unlist(query_c)), "-")[[1]][1:2], collapse =
                   "_")
  print(loop_t)
  data_for_loop = data_by_loop_type_list_unduplicated[[loop_t]][[1]]
  dis = forceSymmetric(all_dist_matrix[[loop_t]])
  #split the data by the clusters
  data_per_cluster = split(data_for_loop, as.character(data_for_loop$cluster_type))
  
  data_this_cluster = data_per_cluster[[query_c]]
  
  # get he radius
  all_pdbs_this_cluster = tolower(data_this_cluster$PDB)
  number = dim(data_this_cluster)[1]
  radius = 1.5 #all_the_max_distances[[query_c]]/15
  if (length(radius) == 0) {
    return(NA)
  } else if (is.na(radius)) {
    return(NA)
  }
  if (!this_pdb %in% colnames(dis)) {
    return(NA)
  }
  pdbs = names(dis[this_pdb, ])
  pdbs = pdbs[pdbs != this_pdb]
  if (length(which(all_pdbs_this_cluster %in% pdbs)) <= 1) {
    return(NA)
  }
  all_pdbs_this_cluster = all_pdbs_this_cluster[all_pdbs_this_cluster %in%
                                                  pdbs]
  all_dis = dis[this_pdb, all_pdbs_this_cluster]
  all_dis = all_dis[!is.na(all_dis)]
  if (length(all_dis) < 10) {
    return(NA)
  }
  the_number = length(which(all_dis < radius))
  # count the number
  return(the_number)
  
}


# filtering the clusters with certain filter of cases numbers
filter_by_case_number_by_comparison <-
  function(right_cases, wrong_cases, filter_case_n) {
    count_the_cluster = FALSE
    splitted_right = split(right_cases, right_cases$query_c)
    splitted_wrong = split(wrong_cases, wrong_cases$query_c)
    splitted_right_filtered = lapply(names(splitted_right), function(x) {
      y = splitted_right[[x]]
      
      print(y)
      real_y_n = dim(y[!is.na(y$neighbor_number), ])[1]
      if (real_y_n < filter_case_n)
        return(NULL)
      else
        return(y)
    })
    names(splitted_right_filtered) = names(splitted_right)
    
    splitted_wrong_filtered = lapply(names(splitted_wrong), function(x) {
      y = splitted_wrong[[x]]
      
      print(y)
      real_y_n = dim(y[!is.na(y$neighbor_number), ])[1]
      if (real_y_n < filter_case_n)
        return(NULL)
      else
        return(y)
    })
    names(splitted_wrong_filtered) = names(splitted_wrong)
    
    
    splitted_right_filtered = splitted_right_filtered[!unlist(lapply(splitted_right_filtered, is.null))]
    splitted_wrong_filtered = splitted_wrong_filtered[!unlist(lapply(splitted_wrong_filtered, is.null))]
    
    
    splitted_right_filtered_row_binded = lapply(names(splitted_right_filtered), function(x) {
      if (!is.null(splitted_wrong_filtered[[x]])) {
        binded_info = rbind(splitted_right_filtered[[x]], splitted_wrong_filtered[[x]])
        return(binded_info)
      } else{
        return(NULL)
      }
    })
    right_filtered_row_binded = do.call(rbind, splitted_right_filtered_row_binded)
    return(right_filtered_row_binded)
  }


even_out_all_classes <- function(training_cases) {
  splitted_class = split(training_cases, training_cases$cluster_type)
  proportions = unlist(lapply(split(training_cases, training_cases$cluster_type), function(x) {
    dim(x)[1]
  }))
  max_class = names(proportions)[which(proportions == max(proportions))]
  remaining_class = names(proportions)[names(proportions) != max_class]
  for (each_remaining in remaining_class) {
    ratio = proportions[max_class] %/% proportions[each_remaining]
    if (ratio > 1) {
      ratio = ratio - 1
    }
    to_be_added = do.call(rbind, replicate(ratio, splitted_class[[each_remaining]], simplify = FALSE))
    training_cases = rbind(training_cases, to_be_added)
  }
  return(training_cases)
}


# getting the training result
generic_train <-
  function(each_loop_length_data_feature_string_rmsd,
           each_method,
           training_cases,
           gbmGrid,
           loop_type) {
    # to add classes to off set
    training_cases = even_out_all_classes(training_cases)
    if (dim(training_cases)[1] * 0.6 * 0.5 < as.numeric(as.character(gbmGrid[["n.minobsinnode"]]))) {
      training_cases = rbind(training_cases, training_cases)
      training_cases = rbind(training_cases, training_cases)
    }
    n.trees = gbmGrid[["n.trees"]]
    
    
    
    set.seed(sample(1:100000, 1))
    r <- n_repeats # number of repeats
    k <- n_folds # number of folds
    
    returned_results = make_3_10_cross_val(training_cases, r, k)  # make division
    folds_spec = returned_results[[1]]
    training_cases = returned_results[[2]]   # the training cases would have its mo
    folds.list.out = folds_spec[[1]]
    folds.list = folds_spec[[2]]
    fitControl <- trainControl(
      method = "repeatedcv",
      repeats = r,
      number = k,
      #preProcOptions = list(thresh = 0.95),
      index = folds.list
      ,
      indexOut = folds.list.out,
      ## Estimate class probabilities
      classProbs = TRUE,
      returnResamp = "all",
      ## Evaluate performance using
      ## the following function
      savePredictions = "final",
      summaryFunction = multiClassSummary
    )
    trained_model = ""
    trained_model <-
      train(
        each_loop_length_data_feature_string_rmsd,
        data = training_cases,
        #distribution = "adaboost",
        method = "gbm",
        bag.fraction = 0.5,
        # fold number 10
        #nTrain = round(nrow(training_cases) *.75),
        trControl = fitControl,
        tuneGrid = gbmGrid,
        verbose = TRUE,
        
        ## Specify which metric to optimize
        metric = "kappa"
      )
    return(trained_model)
    
  }

convert_none <- function(model_re,
                         blindBLAST_unique_clusters,
                         loop) {
  obs = gsub("_", "-", as.character(model_re$obs))
  pred = gsub("_", "-", as.character(model_re$pred))
  
  obs_which = which(!obs %in% blindBLAST_unique_clusters)
  print(obs_which)
  pred_which = which(!pred %in% blindBLAST_unique_clusters)
  print(pred_which)
  pred[pred_which] = rep(paste(c(strsplit(loop, "_")[[1]], "none"), collapse =
                                 "-"), length(pred_which))
  obs[obs_which] = rep(paste(c(strsplit(loop, "_")[[1]], "none"), collapse =
                               "-"), length(obs_which))
  returned_l = list()
  returned_l[[1]] = obs
  returned_l[[2]] = pred
  return(returned_l)
  
}


train_final_model <-
  function(each_loop_length_data_feature_string_rmsd,
           each_method,
           training_cases,
           gbmGrid) {
    # to add classes to off set
    training_cases = even_out_all_classes(training_cases)
    if (dim(training_cases)[1] * 0.6 * 0.5 < as.numeric(as.character(gbmGrid[["n.minobsinnode"]]))) {
      training_cases = rbind(training_cases, training_cases)
      training_cases = rbind(training_cases, training_cases)
    }
    
    
    
    r <- 1 # number of repeats
    k <- 1 # number of folds
    
    returned_results = make_3_10_cross_val(training_cases, r, k)  # make division
    folds_spec = returned_results[[1]]
    training_cases = returned_results[[2]]   # the training cases would have its mo
    
    fitControl <- trainControl(
      method = "none",
      ## Estimate class probabilities
      classProbs = TRUE,
      returnResamp = "all",
      ## Evaluate performance using
      ## the following function
      savePredictions = "final",
      summaryFunction = multiClassSummary
    )
    trained_model = ""
    trained_model <-
      train(
        each_loop_length_data_feature_string_rmsd,
        data = training_cases,
        #distribution = "adaboost",
        method = "gbm",
        bag.fraction = 0.5,
        # fold number 10
        #nTrain = round(nrow(training_cases) *.75),
        trControl = fitControl,
        tuneGrid = gbmGrid,
        verbose = TRUE,
        
        ## Specify which metric to optimize
        metric = "kappa"
      )
    return(trained_model)
    
  }

# get accuracy from blindBLAST result with consideration of the n-cv-m-repeats scheme
get_conf_from_blindBLAST <-
  function(model_re_com,
           unique_cluster,
           repeats,
           loop) {
    print(model_re_com)
    print("line 1248")
    model_re_com[, 1] = sapply(strsplit(model_re_com[, 1], "\\."), "[[", 1)
    model_re_com[, 2] = sapply(strsplit(model_re_com[, 2], "\\."), "[[", 1)
    print("line 1251")
    which_1 = which(!model_re_com[, 1] %in% unique_cluster)
    print("line 1253")
    model_re_com[which_1, 1] = rep(paste(c(strsplit(loop, "_")[[1]], "none"), collapse =
                                           "-"), length(which_1))
    print("line 1255")
    which_2 = which(!model_re_com[, 2] %in% unique_cluster)
    print("line 1257")
    model_re_com[which_2, 2] = rep(paste(c(strsplit(loop, "_")[[1]], "none"), collapse =
                                           "-"), length(which_2))
    print("line 1259")
    
    conf_t = as.data.frame(table(model_re_com[, 1], model_re_com[, 2]))
    conf_t$Freq = conf_t$Freq / repeats
    accuracy_result = sum(conf_t[as.character(conf_t$Var1) == as.character(conf_t$Var2), "Freq"]) /
      sum(conf_t[, "Freq"])
    print("line 1264")
    
    return(accuracy_result)
  }



generate_Rscript_command <- function(args, Rscript_dir, Rscript_n) {
  args = unlist(lapply(args, as.character))
  executable = paste(args, collapse = "_")
  
  exe_sh_file = paste(c(Rscript_dir, executable, "exe.sh"), collapse = "_")  # customize shell script name
  Rscript_command_line = paste(c("Rscript ", Rscript_n, paste(args, collapse =
                                                                " ")), collapse = " ")    # write shell script
  # write("cd ..", file = exe_sh_file,append=TRUE )
  write(Rscript_command_line, file = exe_sh_file, append = TRUE)
  return(exe_sh_file)
}



generate_Rscript_command_with_function <- function(args, Rscript_dir) {
  args = unlist(lapply(args, as.character))
  executable = paste(args, collapse = "_")
  
  exe_sh_file = paste(c(Rscript_dir, executable, "exe.sh"), collapse = "_")  # customize shell script name
  Rscript_command_line = paste(c("Rscript ", Rscript_n, paste(args, collapse =
                                                                " ")), collapse = " ")    # write shell script
  # write("cd ..", file = exe_sh_file,append=TRUE )
  write(Rscript_command_line, file = exe_sh_file, append = TRUE)
  return(exe_sh_file)
}




generate_Rscript_command_H1 <- function(args, Rscript_dir, Rscript_n) {
  args = unlist(lapply(args, as.character))
  executable = paste(args, collapse = "_")
  exe_sh_file = paste(c(Rscript_dir, executable, "exe.sh"), collapse = "_")  # customize shell script name
  file.remove(exe_sh_file)
  
  Rscript_command_line = paste(c("Rscript ", Rscript_n, paste(args, collapse =
                                                                " ")), collapse = " ")    # write shell script
  # write("cd ..", file = exe_sh_file,append=TRUE )
  write(Rscript_command_line, file = exe_sh_file, append = TRUE)
  return(exe_sh_file)
}

generate_condor_script <-
  function(args,
           exe_sh_file,
           condor_script_dir,
           master_condor_file) {
    executable = paste(args, collapse = "_")
    con_file = paste(c(executable, ".con"), collapse = "_")  #customize condor script name
    system(paste(c("cp ", master_condor_file, " ", con_file), collapse = ""))   # copy the master file to script directory
    
    x <- readLines(con_file)
    y <- gsub("\\$executable", exe_sh_file, x)
    z <- gsub("\\$condor_script_dir", condor_script_dir, y)
    w <- gsub("\\$cores", cores, z)
    cat(w, file = con_file, sep = "\n")
    
    print(paste(c(condor_script_dir, "/", con_file), collapse = ""))
    a = paste(c("condor_submit ", condor_script_dir, "/", con_file),
              collapse = "")
    print(a)
    return(a)
  }




preprocess_data <-
  function(loop_type,
           data_by_loop_type_list_unduplicated) {
    print("1093")
    data = data_by_loop_type_list_unduplicated[[loop_type]][[1]]
    data$cluster_type = sub("-", "_", data$cluster_type)
    data$cluster_type = as.factor(as.character(data$cluster_type))
    
    
    
    print("line 1100")
    #sequences$rmsd_cluster = as.character(sequences$rmsd_cluster)
    sequences = data
    
    print("line 1104")
    features = data_by_loop_type_list_unduplicated[[loop_type]][[4]]
    print("line 1106")
    print(features)
    formu_str = paste(c("cluster_type ~ ", paste(features, collapse = " + ")), collapse =
                        "")
    print(formu_str)
    
    each_loop_length_data_feature_string = as.formula(formu_str)
    print("line 1108")
    all_cases =  sequences[, c(features, "id", "cluster_type")]
    the_levels = unique(unlist(data_by_loop_type_list_unduplicated[[loop_type]][[1]][, features]))
    for (each_f in features) {
      all_cases[, each_f] = factor(all_cases[, each_f], levels = the_levels)
    }
    print("line 1112")
    all_cases = all_cases[complete.cases(all_cases),]
    all_cases$cluster_type = gsub("-", "_", all_cases$cluster_type)
    all_cases$cluster_type = gsub(",", ".", all_cases$cluster_type)
    all_cases$cluster_type = gsub("\\*", "none", all_cases$cluster_type)
    
    all_cases$cluster_type = as.factor(as.character(all_cases$cluster_type))
    return(list(all_cases, each_loop_length_data_feature_string))
  }

generate_eta <- function(loop) {
  nl = dim(data_by_loop_type_list_unduplicated[[loop]][[1]])[1]
  eta = max (0.01, 0.1 * min(1, nl / 10000))
  return(eta)
}



execute_training_rscript <- function(args_list) {
  # function for
  args = args_list[[1]]
  print("line 1141")
  print(args)
  
  #  data_by_loop_type_list_unduplicated=args_list[[2]]
  getwd()
  
  loop_type = as.character(args[[1]])
  num_core = as.numeric(args[[2]])
  interaction.depth = as.numeric(args[3])
  n.trees = as.numeric(args[4])
  shrinkage = as.numeric(args[5])
  n.minobsinnode = as.numeric(args[6])
  registerDoMC(num_core)
  print("line 1141")
  
  # naming the model file for latter easy parsing
  gbmGrid =
    data.frame(
      interaction.depth = as.numeric(interaction.depth),
      n.trees = as.numeric(n.trees),
      shrinkage = as.numeric(shrinkage),
      n.minobsinnode = as.numeric(n.minobsinnode)
    )
  
  print(gbmGrid)
  parameter_spe = paste(unlist(gbmGrid), collapse = "-")
  each_method = "gbm_test"
  cluster_dis = "north"
  gbm_result_dir = paste(c(gbm_models_dir), collapse = "")   # build a subdirectory in the current directory to store trained models
  print("line 1150")
  
  to_save_file = paste(
    c(
      gbm_result_dir,
      "/",
      loop_type,
      "_",
      paste(c(
        each_method, cluster_dis, parameter_spe
      ), collapse = "-"),
      "_trained_model_extra_test.rds"
    ),
    collapse = ""
  )
  print(to_save_file)
  # read pyigclassify file
  print(getwd())
  
  print("line 1158")
  
  # perform gbm training
  if (file.exists(to_save_file)) {
    print("file_already_exists!")
  } else{
    print(loop_type)
    print(loop_type)
    preprocess_result =  preprocess_data(loop_type, data_by_loop_type_list_unduplicated)
    all_cases = preprocess_result[[1]]
    each_loop_length_data_feature_string = preprocess_result[[2]]
    print("line 1165")
    
    trained_model = generic_train(
      each_loop_length_data_feature_string,
      each_method,
      all_cases,
      gbmGrid,
      loop_type
    )
    print("line 1168")
    print(to_save_file)
    saveRDS(trained_model, file = to_save_file)
    print(c("model file saved as ", to_save_file))
    
    
    
  }
}












execute_training_rscript_comp <-
  function(args_list) {
    # function for
    args = args_list[[1]]
    print("line 1141")
    print(args)
    
    #  data_by_loop_type_list_unduplicated=args_list[[2]]
    getwd()
    print('line 1129')
    # specify parameter as the args
    print(args)
    print(args[[2]])
    print(args[2])
    loop_type = as.character(args[[1]])
    num_core = as.numeric(args[[2]])
    interaction.depth = as.numeric(args[3])
    n.trees = as.numeric(args[4])
    shrinkage = as.numeric(args[5])
    n.minobsinnode = as.numeric(args[6])
    batch = args_list[[3]]
    registerDoMC(num_core)
    print("line 1141")
    
    # naming the model file for latter easy parsing
    gbmGrid =
      data.frame(
        interaction.depth = as.numeric(interaction.depth),
        n.trees = as.numeric(n.trees),
        shrinkage = as.numeric(shrinkage),
        n.minobsinnode = as.numeric(n.minobsinnode)
      )
    
    print(gbmGrid)
    parameter_spe = paste(unlist(gbmGrid), collapse = "-")
    each_method = "gbm_test"
    cluster_dis = "north"
    gbm_result_dir = paste(c("./rmsd_cluster_hits_rmsd/"), collapse = "")   # build a subdirectory in the current directory to store trained models
    print("line 1150")
    
    to_save_file = paste(
      c(
        gbm_result_dir,
        "",
        loop_type,
        "_",
        paste(c(
          each_method, cluster_dis, parameter_spe, batch
        ), collapse = "-"),
        "_trained_model_extra_test.rds"
      ),
      collapse = ""
    )
    
    # read pyigclassify file
    print(getwd())
    
    print("line 1158")
    
    # perform gbm training
    if (file.exists(to_save_file)) {
      print("file_already_exists!")
    } else{
      print(loop_type)
      print(loop_type)
      preprocess_result =  preprocess_data(loop_type, data_by_loop_type_list_unduplicated)
      all_cases = preprocess_result[[1]]
      each_loop_length_data_feature_string = preprocess_result[[2]]
      print("line 1165")
      
      trained_model = generic_train(
        each_loop_length_data_feature_string,
        each_method,
        all_cases,
        gbmGrid,
        loop_type
      )
      print("line 1168")
      
      saveRDS(trained_model, file = to_save_file)
      print(to_save_file)
      
      
      
    }
  }


# load some naming variables into the script environment
specify_variable_names <-
  function(cur_dire,
           data_dir,
           plot_dir,
           data_table_file) {
    #the_method="blindblast_just_get_alignment"
    cluster_dis = "north"
    subsitution_matrix_name = "wahtw"  #"/Volumes/lab/macbook/lab_work_data/vall_rmsd/loop_sub_matrix.csv"
    subsitution_matrix = "PAM30"
    current_d = getwd()
    if (!grepl(cur_dire, current_d)) {
      result_dir = paste(c("./", cur_dire, "/", data_dir, "/"), collapse = "")
      
      overall_prefix = paste(c("./", cur_dire, "/", data_dir, "/"), collapse =
                               "")
      plot_dir = paste(c("./", cur_dire, "/Plots/"), collapse = "")
    } else{
      result_dir = paste(c("./", data_dir, "/"), collapse = "")
      plot_dir = paste(c("./", plot_dir, "/"), collapse = "")
      overall_prefix = paste(c("./", data_dir), collapse = "")
      
    }
    mkcommand = paste(c("mkdir ", result_dir), collapse = " ")
    system(mkcommand)
    #methods= c(the_method) # specify the machine learning methods to be used
    the_method = "somerandommethod"
    each_method = the_method
    
    # end of iterating all folds
    file = paste(c(overall_prefix, "/", data_table_file, ".rds"), collapse = "")
    data(AAPAM30)
    #AAPAM30_index=1:20
    #names(AAPAM30_index)=rownames(AAPAM30)
    system(paste(c("mkdir ", plot_dir), collapse = ""))
    return(
      list(
        the_method,
        each_method,
        cluster_dis,
        result_dir,
        subsitution_matrix,
        subsitution_matrix_name,
        plot_dir,
        file
      )
    )
  }


execute_training_rscript_final_model <- function(args) {
  print("at line 1230")
  loop_type = as.character(args[1])
  num_core = args[2]
  interaction.depth = args[3]
  n.trees = args[4]
  shrinkage = args[5]
  n.minobsinnode = args[6]
  registerDoMC(num_core)
  
  gbmGrid = data.frame(
    interaction.depth = as.numeric(interaction.depth),
    n.trees = as.numeric(n.trees),
    shrinkage = as.numeric(shrinkage),
    n.minobsinnode = as.numeric(n.minobsinnode)
  )
  
  parameter_spe = paste(unlist(gbmGrid), collapse = "-")
  each_method = "gbm_test"
  gbm_result_dir = paste(c("./final_models/"), collapse = "")
  system(paste(c("mkdir ", gbm_result_dir), collapse = ""))
  
  to_save_file = paste(c(
    gbm_result_dir,
    "",
    loop_type,
    "_",
    paste(c(
      each_method, cluster_dis, parameter_spe
    ), collapse = "-"),
    "_final_model.rds"
  ),
  collapse = "")
  
  if (file.exists(to_save_file)) {
    print("file_already_exists!")
  } else{
    print("at line 1249")
    preprocess_result =  preprocess_data(loop_type, data_by_loop_type_list_unduplicated)
    all_cases = preprocess_result[[1]]
    each_loop_length_data_feature_string = preprocess_result[[2]]
    # tune the parameter
    print("line 1254")
    trained_model = train_final_model(each_loop_length_data_feature_string,
                                      each_method,
                                      all_cases,
                                      gbmGrid)
    saveRDS(trained_model, file = to_save_file)
    
    
    
  }
}

parallelize_jobs <-
  function(n,
           count,
           a_new_script,
           exe_sh_file,
           parallel_script) {
    if (a_new_script) {
      cat("\n", file = parallel_script)
      a_new_script = FALSE
    }
    cat(paste(c(" sh ", exe_sh_file, " & \n"), collapse = ""), file = parallel_script, append =
          TRUE)   # run stuff in parallel
    
    if (count == n) {
      count = 0
      cat("wait \n", file = parallel_script, append = TRUE)
      
    }
    count = count + 1
    return(count)
  }


# calculate dihedral angle distances
calc_dist <- function(x) {
  loop_data = x
  dist_matrix = matrix(nrow = dim(loop_data)[1], ncol = dim(loop_data)[1])
  
  for (i in 1:(dim(dist_matrix)[1] - 1)) {
    for (j in (i + 1):dim(dist_matrix)[1]) {
      print(c(i, j))
      
      # used in directional statistics
      aa = strsplit(as.character(loop_data[i, "dihedrals"]), ":")[[1]]
      bb = strsplit(as.character(loop_data[j, "dihedrals"]), ":")[[1]]
      n = as.numeric(as.character(loop_data[i, "length"]))
      chunk <- function(x, n)
        split(x, ceiling(seq_along(x) / 3))
      aa_phi = as.numeric(sapply(chunk(aa, n), "[[", 1))
      aa_psi = as.numeric(sapply(chunk(aa, n), "[[", 2))
      bb_phi = as.numeric(sapply(chunk(bb, n), "[[", 1))
      bb_psi = as.numeric(sapply(chunk(bb, n), "[[", 2))
      d = sum(2 * (1 - cos(pi / 180 * (
        aa_phi - bb_phi
      )))) + sum(2 * (1 - cos(pi / 180 * (
        aa_psi - bb_psi
      ))))
      dist_matrix[i, j] = d
      average = acos(1 - d / 2)
      average_dis_matrix[i, j] = average
    }
  }
  original_dist = dist_matrix
  dist_matrix <-
    as.SparseSimilarityMatrix(forceSymmetric(dist_matrix))
  return(dist_matrix)
}


write_list_into_single_csv_sp <- function(table_list, output_file) {
  count = 1
  for (x in names(table_list)) {
    if (count == 1) {
      write.table("", file = output_file)
    }
    write.table(x, file = output_file, append = TRUE)
    write.table(
      table_list[[x]],
      sep = ",",
      file = output_file,
      col.names = FALSE,
      row.names  = FALSE,
      append = TRUE
    )
    count = count + 1
  }
}


read_in_data <- function(file_name) {
  d_table = read.table(file_name, header = TRUE)
  if (!all(c("CDR", "length", "cluster_type", "seq") %in% names(d_table))) {
    print("Not enough information")
    stop()
  } else if (!"identifier" %in% names(d_table)) {
    d_table$identifier = paste(d_table$cluster_type, rownames(d_table), sep =
                                 ".")
  }
  # check for feature length
  d_table$type_length = paste(d_table$CDR, d_table$length, sep = "_")
  d_table_s = split(d_table, d_table$type_length)
  d_table_s_listed = lapply(d_table_s, function(single_d_table) {
    feature_splitted = as.data.frame(do.call(rbind, apply(single_d_table, 1, function(x) {
      sp_c = as.list(strsplit(x[["seq"]], "")[[1]])
      
    })))
    feature_s = cbind(single_d_table, feature_splitted)
    return(feature_s)
  })
  data_list = lapply(names(d_table_s_listed), function(x) {
    this_data_list = list()
    this_data_list[[1]] = d_table_s_listed[[x]]
    this_data_list[[4]] = grep("V", names(d_table_s_listed[[x]]), value = TRUE)
    this_data_list[[2]] = as.formula(paste("cluster_type ~ ", paste(this_data_list[[4]], collapse =
                                                                      "+")))
    return(this_data_list)
  })
  names(data_list) = names(d_table_s_listed)
  return(data_list)
}


gbm_train_script <- function(args) {
  print("line 4")
  
  # load the environment
  currect_d = getwd()
  print(c("current directory is ", currect_d))
  #  source("functions.R")
  #  source("utility_function.R")
  cw = getwd()
  print(paste(c("current working directory is ", cw), collapse = "  "))
  data_by_loop_type_list_unduplicated = readRDS(data_save_f)
  
  # load required packages
  
  args_list = list()
  args_list[[1]] = args
  args_list[[2]] = data_by_loop_type_list_unduplicated
  
  execute_training_rscript(args_list)
}




get_best_para <- function(all_result) {
  print("line 1679")
  all_gird_search_results = list()
  
  
  col_names =  c("branch_level",
                 "trees",
                 "shrinkage",
                 "min_node",
                 "accuracy",
                 "accuracy_sd")
  print("line 1690")
  result_table = data.frame(matrix(nrow = length(names(all_result)), ncol = length(col_names)))
  colnames(result_table) = col_names
  rownames(result_table) = names(all_result)
  best_parameters_each_loop = list()
  for (each_loop in names(all_result)) {
    all_gird_search_results[[each_loop]] = all_result[[each_loop]]
    
    # do a plot
    x = unique(all_gird_search_results[[each_loop]]$n.trees)
    y = sort(unique(all_gird_search_results[[each_loop]]$interaction.depth))
    z = matrix(nrow = length(x), ncol = length(y))
    z = as.data.frame(z)
    
    rownames(z) = x
    colnames(z) = y
    # only select the mino equal to 5
    sub_table = all_gird_search_results[[each_loop]]
    for (ind in 1:dim(sub_table)[1]) {
      z[as.character(sub_table[ind, "n.trees"]), as.character(sub_table[ind, "interaction.depth"])] =
        sub_table[ind, "Accuracy"]
    }
    
    
    max_indses = which(all_gird_search_results[[each_loop]]$Accuracy == max(all_gird_search_results[[each_loop]]$Accuracy))
    max_indses_index = which(all_gird_search_results[[each_loop]][max_indses, "Kappa"] ==
                               max(all_gird_search_results[[each_loop]][max_indses, "Kappa"]))
    final_max_ind = max_indses[max_indses_index]
    if (length(final_max_ind) != 1) {
      which_i = which(all_gird_search_results[[each_loop]][final_max_ind, "n.trees"] ==
                        min(all_gird_search_results[[each_loop]][final_max_ind, "n.trees"]))
      final_max_ind = final_max_ind[which_i]
      if (length(final_max_ind) != 1) {
        which_i = which(all_gird_search_results[[each_loop]][final_max_ind, "interaction.depth"] ==
                          min(all_gird_search_results[[each_loop]][final_max_ind, "interaction.depth"]))
        final_max_ind = final_max_ind[which_i]
        
      }
      
    }
    
    parameters = rownames(all_gird_search_results[[each_loop]])[final_max_ind]
    parameters = t(data.frame(strsplit(as.character(parameters), "-")))
    
    result_table[each_loop, c("branch_level", "trees", "shrinkage", "min_node")] =
      parameters
    result_table[each_loop, c("accuracy", "accuracy_sd")] = all_gird_search_results[[each_loop]][final_max_ind, c("Accuracy", "AccuracySD")]
    best_parameters_each_loop[[each_loop]] = parameters
    # slice3D(all_gird_search_results[[each_loop]]$n.trees, all_gird_search_results[[each_loop]]$shrinkage, all_gird_search_results[[each_loop]]$interaction.depth,colvar=all_gird_search_results[[each_loop]]$Accuracy)
    #  M <- mesh(all_gird_search_results[[each_loop]]$n.trees, all_gird_search_results[[each_loop]]$shrinkage, all_gird_search_results[[each_loop]]$interaction.depth)
  }
  
  rownames(result_table) = split_vector_and_replace(rownames(result_table), "_", 1, 2, "-")
  result_table = result_table[order(order_factor_by_two_component(rownames(result_table), "-", 1, 2),
                                    rownames(result_table)),]
  write.csv(
    result_table,
    file = paste(result_dir, "best_parameters.csv"),
    row.names = TRUE,
    col.names  = TRUE
  )
  best_parameters_each_loop = lapply(best_parameters_each_loop, function(x) {
    print(x)
    x = as.data.frame(x)
    x = lapply(x, as.character)
    names(x) = c("interaction.depth", "n.trees", "eta", "n.minobsinnode")
    x = as.data.frame(x)
    return(x)
  })
  print(result_table)
  return_list = list()
  return_list[[1]] = result_table
  return_list[[2]] = best_parameters_each_loop
  return(return_list)
}




generate_arguments <-
  function(loop,
           total_max_core,
           complexity,
           trees,
           eta,
           min_node_n) {
    cores = 1
    argument_table = data.frame(expand.grid(loop, cores, complexity, trees, eta, min_node_n))
    
    cores = ceiling(total_max_core / dim(argument_table)[1])
    argument_table = data.frame(expand.grid(loop, cores, complexity, trees, eta, min_node_n))
    
    argument_table[, 2:6] = lapply(argument_table[, 2:6], as.numeric)
    arguments = as.data.frame(do.call(cbind, apply(argument_table, 1, as.list)))
    return(arguments)
  }


cal_blindBLAST_accuracy <-
  function(ten_foldcv_blindblastlist, n_folds) {
    #ten_foldcv_blindblastlist
    conf_tables_all_loops_blindBLAST = list()
    conf_tables_all_loops_blindBLAST_diff = list()
    new_accuracy_list = list()
    blindBLAST_errorcount_lists = list()
    
    all_folds_sd_list = list()
    for (loop in names(ten_foldcv_blindblastlist)) {
      print("line 1781")
      print(loop)
      model_re = ten_foldcv_blindblastlist[[loop]]
      model_re_by_repeats = lapply(chunk2(1:length(model_re), 3), function(x) {
        do.call(rbind, model_re[x])
      })
      print("line 1786")
      n_repeats = ceiling(length(model_re) %/% n_folds)
      print(length(model_re))
      print(n_repeats)
      print("line 1790")
      unique_cluster = gsub("\\*",
                            "none",
                            unique(data_by_loop_type_list_unduplicated[[loop]][[1]]$cluster_type))
      re_by_repeats = lapply(model_re_by_repeats, function(x) {
        # get_conf_from_blindBLAST(x, unique_cluster, 1)
        print("line 1791")
        print(x)
        conf_t_r = as.data.frame(table(sapply(strsplit(x[, 1], "\\."), "[[", 1),
                                       sapply(strsplit(x[, 2], "\\."), "[[", 1)))
      })
      error_by_repeats = lapply(re_by_repeats, function(x) {
        x = as.data.frame(x)
        error_c = x[as.character(x[, 1]) != as.character(x[, 2]),]
        sum(error_c$Freq)
      })
      
      blindBLAST_errorcount_lists[[loop]] = unlist(error_by_repeats)
      
      frame = re_by_repeats[[1]]
      for (ind in 2:(length(re_by_repeats))) {
        print(re_by_repeats[[ind]])
        frame = as.data.frame(merge(frame, re_by_repeats[[ind]], by = c("Var1", "Var2")))
      }
      
      this_loop_acc_by_repeats = lapply(3:dim(frame)[2], function(x) {
        sum(frame[as.character(frame[, 1]) == as.character(frame[, 2]), x]) / sum(frame[, x])
      })
      accuracy_av = mean(unlist(this_loop_acc_by_repeats))
      sd_by_repeats = sd(unlist(this_loop_acc_by_repeats))
      
      each_classification_sds = apply(frame, 1, function(x) {
        sd(unlist(x[3:5]))
      })
      sd_f = cbind(frame[, 1:2], each_classification_sds)
      #individual accuracy
      all_folds_acc = unlist(lapply(model_re, function(ind_re) {
        get_conf_from_blindBLAST(ind_re, unique_cluster, 1, loop)
      }))
      all_folds_acc = all_folds_acc[!is.na(all_folds_acc)]
      all_folds_sd_list[[loop]] = sd(all_folds_acc)
      # total accu
      
      model_re_com = do.call(rbind, model_re)
      
      accuracy_result = get_conf_from_blindBLAST(model_re_com, unique_cluster, n_repeats, loop)
      new_accuracy_list[[loop]] = accuracy_result
      print(loop)
      
      conf_t = as.data.frame(table(sapply(
        strsplit(model_re_com[, 1], "\\."), "[[", 1
      ), sapply(
        strsplit(model_re_com[, 2], "\\."), "[[", 1
      )) / n_repeats)
      conf_t = merge(conf_t, sd_f, by = c("Var1", "Var2"))
      conf_tables_all_loops_blindBLAST[[loop]] = conf_t
      
      conf_t = conf_t[conf_t$Freq > 0 &
                        as.character(conf_t[, 1]) != as.character(conf_t[, 2]), ]
      conf_tables_all_loops_blindBLAST_diff[[loop]] = conf_t
    }
    conf_tables_all_loops_blindBLAST_diff
    
    blindblast_by_loop = lapply(ten_foldcv_blindblastlist, function(x) {
      if (n_repeats > 1) {
        split_c = chunk2(1:length(x), n_repeats)
      } else{
        split_c = list()
        split_c[[1]] = 1:length(x)
      }
      lapply(split_c, function(y) {
        do.call(rbind, x[y])
      })
    })
    blindblast_by_loop_bind = lapply(blindblast_by_loop, function(x) {
      accuracies = lapply(x, function(y) {
        all = y
        all[, 1] = split_vector_and_replace(all[, 1], "\\.", 1, 1, "")
        all[, 2] = split_vector_and_replace(all[, 2], "\\.", 1, 1, "")
        accu = dim(all[all[, 1] == all[, 2],])[1] / dim(all)[1]
      })
      accuracies = unlist(accuracies)
      
    })
    
    blindBLAST_mean_accu = lapply(blindblast_by_loop_bind, mean)
    saveRDS(blindBLAST_mean_accu, file = paste(c(
      result_dir, "/", "blindBLAST_mean_accu.rds"
    ), collapse = ""))
    blindBLAST_accu_std = lapply(blindblast_by_loop_bind, sd)
    saveRDS(blindBLAST_accu_std, file = paste(c(
      result_dir, "/", "blindBLAST_accu_std.rds"
    ), collapse = ""))
    saveRDS(all_folds_sd_list, file = paste(c(
      result_dir, "/", "all_folds_sd_list.rds"
    ), collapse = ""))
    summary_info = as.data.frame(cbind(unlist(blindBLAST_mean_accu), unlist(blindBLAST_accu_std)))
    print(conf_tables_all_loops_blindBLAST)
    print(conf_tables_all_loops_blindBLAST_diff)
    saveRDS(conf_tables_all_loops_blindBLAST, file = paste(
      c(result_dir, "/", "conf_tables_all_loops_blindBLAST.rds"),
      collapse = ""
    ))
    saveRDS(conf_tables_all_loops_blindBLAST_diff,
            file = paste(
              c(
                result_dir,
                "/",
                "conf_tables_all_loops_blindBLAST_diff.rds"
              ),
              collapse = ""
            ))
    saveRDS(blindBLAST_errorcount_lists, file = paste(
      c(result_dir, "/", "blindBLAST_errorcount_lists.rds"),
      collapse = ""
    ))
    return(summary_info)
  }




cal_GBM_accuracy <-
  function(all_models_list_by_loop,
           best_parameters_each_loop ,
           all_models,
           data_by_loop_type_list_unduplicated) {
    conf_tables_all_loops_gbm_diff = list()
    conf_tables_all_loops_gbm = list()
    all_gbm_pred_result = list()
    gbm_folds_sd = list()
    gbm_errorcount_list = list()
    gbm_sd_by_loops = list()
    gbm_accuracy_by_loops = list()
    gbm_errorcount_sd_list = list()
    
    for (loop in names(all_models_list_by_loop)) {
      paras = best_parameters_each_loop[[loop]]
      paras = lapply(paras, as.character)
      best_para = paste(unlist(paras), collapse = "-")
      file_n = grep(loop, grep(best_para, all_models, value = TRUE), value =
                      TRUE)[1]
      model_file = paste(c("./", gbm_models_dir, "/", file_n, sep = ""), collapse =
                           "")
      model = readRDS(model_file)
      model_results_list = lapply(model_file, function(x) {
        model = readRDS(x)
        model_result = model$pred
        re = model_result[, c("pred", "obs", "Resample")]
        #print(re)
        names(re) = c("pred", "obs", "Resample")
        return(re)
      })
      
      if (length(model_results_list) == 1) {
        model_re = do.call(rbind, model_results_list)
        model_re_resample = split(model_re, model_re$Resample)
        resample_chunks = chunk2(1:length(model_re_resample), n_repeats)
        results_by_repeats = lapply(resample_chunks, function(x) {
          result = do.call(rbind, model_re_resample[x])
          
        })
        
      } else{
        # process irregular cases, overwrite the repeats numbering
        results_by_repeats = model_results_list
        
      }
      
      
      for (l in 1:length(results_by_repeats)) {
        the_f = results_by_repeats[[l]]
        the_f$repeats  = rep(l, dim(the_f)[1])
        results_by_repeats[[l]] = the_f
      }
      model_re = do.call(rbind, results_by_repeats)
      
      # split the results into 3 repeats
      
      gbm_folds_sd[[loop]] = model$results$AccuracySD
      
      
      
      blindBLAST_unique_clusters = gsub("\\ * ",
                                        "none",
                                        unique(data_by_loop_type_list_unduplicated[[loop]][[1]]$cluster_type))
      
      returned_l = convert_none(model_re, blindBLAST_unique_clusters, loop)
      obs = returned_l[[1]]
      pred = returned_l[[2]]
      obs = factor(obs, levels = blindBLAST_unique_clusters)
      pred = factor(pred, levels = blindBLAST_unique_clusters)
      
      model_re[, 1] = obs
      model_re[, 2] = pred
      model_re_resplit = split(model_re, model_re$repeats)
      all_gbm_pred_result[[loop]] = model_re
      
      accuracies = lapply(model_re_resplit, function(x) {
        dim(x[x[, 1] == x[, 2],])[1] / dim(x)[1]
      })
      
      errorcount_sd = lapply(model_re_resplit, function(x) {
        returned_l = convert_none(x, blindBLAST_unique_clusters, loop)
        obs = returned_l[[1]]
        pred = returned_l[[2]]
        obs = factor(obs, levels = blindBLAST_unique_clusters)
        pred = factor(pred, levels = blindBLAST_unique_clusters)
        conf_t = as.data.frame(table(obs, pred))
        conf_t = conf_t[conf_t$Freq > 0 &
                          conf_t[, 1] != conf_t[, 2],]
        conf_t[, 1] = gsub("_", "-", conf_t[, 1])
        conf_t[, 2] = gsub("_", "-", conf_t[, 2])
        sum(conf_t$Freq)
      })
      
      errorcount_by_mis_sd = lapply(model_re_resplit, function(x) {
        print(loop)
        returned_l = convert_none(x, blindBLAST_unique_clusters, loop)
        obs = returned_l[[1]]
        pred = returned_l[[2]]
        obs = factor(obs, levels = blindBLAST_unique_clusters)
        pred = factor(pred, levels = blindBLAST_unique_clusters)
        conf_t = as.data.frame(table(obs, pred))
        conf_t = conf_t[conf_t$Freq > 0 &
                          conf_t[, 1] != conf_t[, 2],]
        conf_t[, 1] = gsub("_", "-", conf_t[, 1])
        conf_t[, 2] = gsub("_", "-", conf_t[, 2])
        return(conf_t)
      })
      errorcount_by_mis_fr = errorcount_by_mis_sd[[1]]
      for (x in 2:length(errorcount_by_mis_sd)) {
        errorcount_by_mis_fr = merge(errorcount_by_mis_fr,
                                     errorcount_by_mis_sd[[x]],
                                     by = c("obs", "pred"))
      }
      sds = apply(errorcount_by_mis_fr, 1, function(x) {
        sd(as.numeric(unlist(x[3:(3 + length(errorcount_by_mis_sd) - 1)])))
      })
      errorcount_by_mis_fr = cbind(errorcount_by_mis_fr[, 1:2], sds)
      
      gbm_errorcount_list[[loop]] = unlist(errorcount_sd)
      gbm_errorcount_sd_list[[loop]] = sd(unlist(errorcount_sd))
      print(accuracies)
      gbm_sd_by_loops[[loop]] = sd(unlist(accuracies))
      gbm_accuracy_by_loops[[loop]] = mean(unlist(accuracies))
      conf_t = as.data.frame(table(obs, pred))
      colnames(conf_t)[1:2] = c("Var1", "Var2")
      conf_t$Freq = conf_t$Freq / 3   # average by the repeat number
      
      conf_tables_all_loops_gbm[[loop]] = conf_t
      conf_t = conf_t[conf_t$Freq > 0 &
                        conf_t[, 1] != conf_t[, 2],]
      colnames(errorcount_by_mis_fr) = c("Var1", "Var2" , "sds")
      conf_t = merge(conf_t, errorcount_by_mis_fr)
      conf_t[, 1] = gsub("_", "-", conf_t[, 1])
      conf_t[, 2] = gsub("_", "-", conf_t[, 2])
      conf_tables_all_loops_gbm_diff[[loop]] = conf_t
    }
    print(conf_tables_all_loops_gbm)
    saveRDS(conf_tables_all_loops_gbm, file = paste(
      c(result_dir, "/", "conf_tables_all_loops_gbm.rds"),
      collapse = ""
    ))
    saveRDS(conf_tables_all_loops_gbm_diff, file = paste(
      c(result_dir, "/", "conf_tables_all_loops_gbm_diff.rds"),
      collapse = ""
    ))
    saveRDS(gbm_sd_by_loops, file = paste(c(result_dir, "/", "gbm_sd_by_loops.rds"), collapse =
                                            ""))
    
    gbm_sd_by_loops
    gbm_accuracy_by_loops
    saveRDS(gbm_accuracy_by_loops, file = paste(c(
      result_dir, "/", "gbm_accuracy_by_loops.rds"
    ), collapse = ""))
    saveRDS(gbm_folds_sd, file = paste(c(result_dir, "/", "gbm_folds_sd.rds"), collapse =
                                         ""))
    saveRDS(gbm_errorcount_list, file = paste(c(
      result_dir, "/", "gbm_errorcount_list.rds"
    ), collapse = ""))
    result_f = data.frame(mean = unlist(gbm_accuracy_by_loops),
                          sd = unlist(gbm_sd_by_loops))
    return(result_f)
  }



process_results_for_accuracies <-
  function(conf_tables_all_loops_blindBLAST,
           conf_tables_all_loops_gbm) {
    conf_tables_all_loops_blindBLAST
    conf_tables_all_loops_gbm
    accuracy_list = list()
    
    for (x_loop in names(conf_tables_all_loops_blindBLAST)) {
      bl = conf_tables_all_loops_blindBLAST[[x_loop]]
      accuracy_list[["blindBLAST"]][[x_loop]] = calculate_accuracy(bl, c("Var1", "Var2"), "Freq")
      bgm = conf_tables_all_loops_gbm[[x_loop]]
      accuracy_list[["gbm"]][[x_loop]] = calculate_accuracy(bgm, c("Var1", "Var2"), "Freq")
    }
    ori = accuracy_list
    accuracy_list[["blindBLAST"]] = unlist(blindBLAST_mean_accu)
    accuracy_list[["gbm"]] = unlist(gbm_accuracy_by_loops)
    common = intersect(names(accuracy_list[["gbm"]]), names(accuracy_list[["blindBLAST"]]))
    accuracy_list = lapply(accuracy_list, function(x) {
      x[common]
    })
    accuracy_gbm_blast = as.data.frame(accuracy_list)
    accuracy_gbm_blast_remove_unknow = accuracy_gbm_blast[complete.cases(accuracy_gbm_blast),]
    accuracy_gbm_blast_remove_unknow$loop = split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow), "_", 1, 1, "-")
    accuracy_gbm_blast_remove_unknow$length = as.numeric(split_vector_and_replace(rownames(accuracy_gbm_blast_remove_unknow), "_", 2, 2, "-"))
    accuracy_gbm_blast_remove_unknow = reorder_factor(accuracy_gbm_blast_remove_unknow, "loop", "length")
    
    accuracy_gbm_blast_remove_unknow_melt = melt(accuracy_gbm_blast_remove_unknow,
                                                 id.vars = c("loop", "length"))
    accuracy_gbm_blast_remove_unknow_melt = as.data.frame(accuracy_gbm_blast_remove_unknow_melt)
    
    
    
    
    gbm_folds_sd = lapply(names(all_folds_sd_list), function(x) {
      if (!x %in% names(gbm_folds_sd)) {
        gbm_folds_sd[[x]] = NA
      }
      return(gbm_folds_sd[[x]])
    })
    sd_list = list()
    sd_list[["blindBLAST"]] = unlist(blindBLAST_accu_std)
    sd_list[["gbm"]] = unlist(gbm_folds_sd)
    sd_list[["gbm"]] = unlist(gbm_sd_by_loops)
    common = intersect(names(sd_list[["gbm"]]), names(sd_list[["blindBLAST"]]))
    sd_list = lapply(sd_list, function(x) {
      x[common]
    })
    sd_gbm_blast = as.data.frame(sd_list)
    sd_gbm_blast_remove_unknow = sd_gbm_blast[complete.cases(sd_gbm_blast),]
    sd_gbm_blast_remove_unknow$loop = split_vector_and_replace(rownames(sd_gbm_blast_remove_unknow), "_", 1, 1, "-")
    sd_gbm_blast_remove_unknow$length = as.numeric(split_vector_and_replace(rownames(sd_gbm_blast_remove_unknow), "_", 2, 2, "-"))
    sd_gbm_blast_remove_unknow = reorder_factor(sd_gbm_blast_remove_unknow, "loop", "length")
    
    sd_gbm_blast_remove_unknow_melt = melt(sd_gbm_blast_remove_unknow, id.vars = c("loop", "length"))
    sd_gbm_blast_remove_unknow_melt = as.data.frame(sd_gbm_blast_remove_unknow_melt)
    colnames(sd_gbm_blast_remove_unknow_melt)[4] = "sd"
    accuracy_sd_gbm_blast_remove_unknow_melt = merge(accuracy_gbm_blast_remove_unknow_melt,
                                                     sd_gbm_blast_remove_unknow_melt)
    accuracy_sd_gbm_blast_remove_unknow_melt$low = accuracy_sd_gbm_blast_remove_unknow_melt$value -
      accuracy_sd_gbm_blast_remove_unknow_melt$sd / 2
    accuracy_sd_gbm_blast_remove_unknow_melt$high = accuracy_sd_gbm_blast_remove_unknow_melt$value +
      accuracy_sd_gbm_blast_remove_unknow_melt$sd / 2
    accuracy_sd_gbm_blast_remove_unknow_melt = as.data.frame(accuracy_sd_gbm_blast_remove_unknow_melt)
    
    return(accuracy_sd_gbm_blast_remove_unknow_melt)
    
    
  }





process_error_counts <-
  function(conf_tables_all_loops_blindBLAST_diff,
           conf_tables_all_loops_gbm_diff) {
    error_count_list = list()
    
    gbm_errorcount_list
    blindBLAST_errorcount_lists
    for (x_loop in names(conf_tables_all_loops_blindBLAST_diff)) {
      print(x_loop)
      bl = conf_tables_all_loops_blindBLAST_diff[[x_loop]]
      values_bl = bl[as.character(bl$Var1) != as.character(bl$Var2), "Freq"]
      
      if (!is.null(values_bl)) {
        print(sum(values_bl))
        error_count_list[["blindBLAST"]][[x_loop]] = sum(values_bl)
      } else{
        error_count_list[["blindBLAST"]][[x_loop]] = NA
      }
      bgm = conf_tables_all_loops_gbm_diff[[x_loop]]
      values_gbm = bgm[as.character(bgm$Var1) != as.character(bgm$Var2), "Freq"]
      
      if (!is.null(values_gbm)) {
        print(sum(values_gbm))
        error_count_list[["gbm"]][[x_loop]] = sum(values_gbm)
      } else{
        error_count_list[["gbm"]][[x_loop]] = NA
      }
    }
    error_count_list_with_sd = list()
    t1 = as.data.frame(error_count_list[["blindBLAST"]])
    t1$loops = rownames(t1)
    t2 = as.data.frame(t(as.data.frame(
      lapply(blindBLAST_errorcount_lists, sd)
    )))
    t2$loops = rownames(t2)
    fr1 = merge(t1, t2, by = "loops")
    colnames(fr1) = c("loop", "error_count", "sd")
    rownames(fr1) = fr1$loop
    
    t1 = as.data.frame(error_count_list[["gbm"]])
    t1$loops = rownames(t1)
    t2 = as.data.frame(t(as.data.frame(lapply(
      gbm_errorcount_list, sd
    ))))
    t2$loops = rownames(t2)
    fr2 = merge(t1, t2, by = "loops")
    colnames(fr2) = c("loop", "error_count", "sd")
    rownames(fr2) = fr2$loop
    common_n = intersect(rownames(fr1), rownames(fr2))
    fr1 = fr1[rownames(fr1) %in% common_n, ]
    fr2 = fr2[rownames(fr2) %in% common_n, ]
    error_count_list_with_sd[[1]] = fr1
    error_count_list_with_sd[[2]] = fr2
    names(error_count_list_with_sd) = c("blindBLAST", "gbm")
    
    error_count_list_frame = do.call(rbind, lapply(names(error_count_list_with_sd), frame_for_plot))
    colnames(error_count_list_frame) = c("loop", "value", "sd", "length", "variable")
    
    figure_error_count = plot_figure(error_count_list_frame,
                                     "length",
                                     "value",
                                     "variable",
                                     "loop",
                                     c(0.9, 0.9))
    error_count_list_frame$low = error_count_list_frame$value - error_count_list_frame$sd /
      2
    error_count_list_frame$high = error_count_list_frame$value + error_count_list_frame$sd /
      2
    return(error_count_list_frame)
  }



read_gbm_models <- function(all_models, gbm_models_dir) {
  all_models = all_models[grepl("gbm", all_models)]
  all_models_list_by_loop = list()
  index_i = 1
  all_t = list()
  for (each_file in all_models) {
    tryCatch({
      each_file_r = paste(c(gbm_models_dir, "/", each_file), collapse =
                            "")
      the_model = readRDS(each_file_r)
      
      index_i = index_i + 1
      print(each_file)
      info = strsplit(each_file, "\\/")[[1]][length(strsplit(each_file, "\\/")[[1]])]
      loop_type = paste(strsplit(info, "_")[[1]][1:2], collapse = "_")
      all_t[[loop_type]] = unique(the_model$trainingData$.outcome)
      number = paste(strsplit(strsplit(info, "_")[[1]][4], "-")[[1]][3:6], collapse =
                       "-")
      model_result = the_model$result
      all_models_list_by_loop[[loop_type]][[number]] = model_result
      rm(the_model)
      gc()
    }, error = function(e) {
      print(e)
    })
  }
  return(all_models_list_by_loop)
  
}


