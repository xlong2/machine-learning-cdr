# read in all rmsd files for rmsd retrieval

current_d = getwd()
if (grepl("cdr_classifier_result", current_d)) {
  source("0.load_function_and_data.R")
  
}

all_rmsd_list = list()
type_list = c("L1_H1", "L2_H2", "L3_L3")
for (type in type_list) {
  L1_H1 = read.table(paste(c("/../DTW_project/", type, ".distmat"), collapse =
                             ""))
  name_file = paste(c("/../DTW_project/", type, ".file.list"), collapse =
                      "")
  L1_H1_names = sapply(strsplit(as.character(unlist(
    read.table(name_file)
  )), "\\/"), "[[", 7)
  marks = sapply(strsplit(L1_H1_names, "_"), "[[", 3)
  pdbs = substr(sapply(strsplit(L1_H1_names, "\\.|-"), "[[", 1), 1, 4)
  
  split_list = unique(strsplit(type, "_")[[1]])
  for (x_type in split_list) {
    L1_index = which(grepl(x_type, marks))
    L1_names = pdbs[L1_index]
    
    L1_H1_names[H1_index][grepl("3rkd", L1_H1_names[H1_index])]
    
    
    L1_matrix = L1_H1[L1_index, L1_index]
    
    undup_L1_names_ind = !duplicated(L1_names)
    L1_matrix = L1_matrix[undup_L1_names_ind, undup_L1_names_ind]
    colnames(L1_matrix) = L1_names[undup_L1_names_ind]
    rownames(L1_matrix) = L1_names[undup_L1_names_ind]
    all_rmsd_list[[x_type]] = L1_matrix
  }
}

saveRDS(all_rmsd_list, file = "./Data_processed/all_rmsd_list.rds")
