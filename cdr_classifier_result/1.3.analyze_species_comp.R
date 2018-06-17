all_data=do.call(rbind,lapply(data_by_loop_type_list_unduplicated,function(x){x[[1]]}))



spe_t=table(all_data$GSpecies)
spe_t[["Hu"]]/sum(unlist(spe_t))
spe_t[["Mo"]]/sum(unlist(spe_t))


0.3+0.4+0.3+0.4+2.8+0.3+0.3+0.3+0.3+0.2+0.2
