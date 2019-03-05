all_data=do.call(rbind,lapply(data_by_loop_type_list_unduplicated,function(x){x[[1]]}))



spe_t=table(all_data$GSpecies)
spe_t[["Hu"]]/sum(unlist(spe_t))
spe_t[["Mo"]]/sum(unlist(spe_t))


