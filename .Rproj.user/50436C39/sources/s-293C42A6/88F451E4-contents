common_loops=intersect(names(ten_foldcv_blindblastlist),names(conf_tables_all_loops_gbm))
ten_foldcv_blindblastlist_common_l=ten_foldcv_blindblastlist[common_loops]
ten_foldcv_blindblastlist_common_l_rbind=do.call(rbind,lapply(ten_foldcv_blindblastlist_common_l,function(x){do.call(rbind,x)}))
ten_foldcv_blindblastlist_common_l_rbind$query_t=split_vector_and_replace(ten_foldcv_blindblastlist_common_l_rbind$X1,"\\.",1,1,"")
ten_foldcv_blindblastlist_common_l_rbind$tem_t=split_vector_and_replace(ten_foldcv_blindblastlist_common_l_rbind$X2,"\\.",1,1,"")
a=dim(ten_foldcv_blindblastlist_common_l_rbind[ten_foldcv_blindblastlist_common_l_rbind$query_t!=ten_foldcv_blindblastlist_common_l_rbind$tem_t,])[1]
b=dim(ten_foldcv_blindblastlist_common_l_rbind)[1]
# this comes from the 3 repeats 


c=1-a/b

conf_tables_all_loops_gbm_common_l=conf_tables_all_loops_gbm[common_loops]
conf_tables_all_loops_gbm_common_l_rbind=do.call(rbind,conf_tables_all_loops_gbm)
d=sum(conf_tables_all_loops_gbm_common_l_rbind[conf_tables_all_loops_gbm_common_l_rbind$Var1!=conf_tables_all_loops_gbm_common_l_rbind$Var2,3])
e=sum(conf_tables_all_loops_gbm_common_l_rbind[,3])

f=1-d/e
