library(GenomicRanges)
library(traseR)
library(tidyverse)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
tbl_in_fn<-function(tmp_file){
  tmp_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  return(tmp_tbl)
}
cl_reduce_coord_fn<-function(hmec_dagger_01_tbl,tmp_res,res_num){
  
  tmp_bin_tbl<-hmec_dagger_01_tbl %>% filter(res==tmp_res) %>%unnest(cols = "bins")
  tmp_bin_tbl<-tmp_bin_tbl%>%mutate(end = as.numeric(bins) + res_num[res]-1)
  tmp_bin_tbl<-tmp_bin_tbl%>%distinct(chr,bins,end)
  inter_cl_Grange<-   GRanges(seqnames=tmp_bin_tbl$chr,
                              ranges = IRanges(start=as.numeric(tmp_bin_tbl$bins),
                                               end=tmp_bin_tbl$end,
                                               names=paste("cl_inter",1:nrow(tmp_bin_tbl),sep='_')
                              ))
  inter_cl_Grange<-IRanges::reduce(inter_cl_Grange)
  return(inter_cl_Grange)
  
}

#-------------------------------------------------------------------------------------------------------
# table with cluster of interest
union_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/GM12878_union_trans_res_dagger_tbl.Rda"
union_cl_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_GM12878_pval_tbl.Rda"

cl_union_tbl<-tbl_in_fn(union_cl_file)
hmec_dagger_01_tbl<-tbl_in_fn(union_hub_file) #%>% filter(res == "5kb" | res == "10kb")
hmec_dagger_01_tbl<-hmec_dagger_01_tbl %>% left_join(.,cl_union_tbl %>% dplyr::rename(node=cl)%>% dplyr::select(chr,node,res,bins))
#-------------------------------------------------------------------------------------------------------
# Spectral clustering results

hub_GRange<-do.call("c",lapply(unique(hmec_dagger_01_tbl$res),function(f){
  tmp<-cl_reduce_coord_fn(hmec_dagger_01_tbl,f,res_num)
  mcols(tmp)<-tibble(res=f)
  return(tmp)
}))

data(taSNP)
res_obj<-traseR(taSNP,IRanges::reduce(hub_GRange),rankby="pvalue",test.method="nonparametric")
