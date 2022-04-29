library(GenomicRanges)
library(traseR)
library(tidyverse)
library(valr)
library(furrr)
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

hub_sets_GRange_build_fn<-function(union_cl_file,union_hub_file){
  cl_union_tbl<-tbl_in_fn(union_cl_file)
  hmec_dagger_01_tbl<-tbl_in_fn(union_hub_file) #%>% filter(res == "5kb" | res == "10kb")
  hmec_dagger_01_tbl<-hmec_dagger_01_tbl %>% left_join(.,cl_union_tbl %>% dplyr::rename(node=cl)%>% dplyr::select(chr,node,res,bins))
  rm(cl_union_tbl)
  
  hub_GRange<-do.call("c",lapply(unique(hmec_dagger_01_tbl$res),function(f){
    tmp<-cl_reduce_coord_fn(hmec_dagger_01_tbl,f,res_num)
    mcols(tmp)<-tibble(res=f)
    return(tmp)
  }))
  return(hub_GRange)
}
#-------------------------------------------------------------------------------------------------------
# table with cluster of interest
union_hub_files<-c(H1="~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/H1_union_top_trans_res_dagger_tbl.Rda",
                   HMEC="~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/HMEC_union_top_trans_res_dagger_tbl.Rda",
                   GM12878="~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda")

union_cl_files<-c(H1="~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_H1_pval_tbl.Rda",
                  HMEC="~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda",
                  GM12878="~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_GM12878_pval_tbl.Rda")

#-------------------------------------------------------------------------------------------------------
# Spectral clustering results

hub_GRanges_l<-lapply(names(union_hub_files),function(i){
  message(i)
  tmp_union_cl_file<-union_cl_files[i]
  tmp_union_hub_file<-union_hub_files[i]
  return(IRanges::reduce(hub_sets_GRange_build_fn(tmp_union_cl_file,tmp_union_hub_file)))
})
names(hub_GRanges_l)<-names(union_hub_files)

hg19_coord <- read_delim("~/Documents/multires_bhicect/data/hg19.genome", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")

tmp_cl_tbl<-hub_GRanges_l[[1]] %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)

data(taSNP)

traseR_res_l<-lapply(1:10,function(x){
  rn_pol<-bed_shuffle(tmp_cl_tbl,genome = hg19_coord,max_tries=1e8,within=T)
  rn_GRange<-GRanges(seqnames=rn_pol$chrom,
                     ranges = IRanges(start=rn_pol$start,
                                      end=rn_pol$end))
  #length(unique(subjectHits(findOverlaps(rn_GRange,taSNP))))
  tmp<-traseR(taSNP,rn_GRange,rankby="pvalue",test.method="binomial",alternative = "greater")
  return(tmp$tb.all$p.value)
})

tmp_cl_tbl<-hub_GRanges_l[[1]] %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)
plan(multisession,workers=4)
traseR_SNP_count<-future_map_int(1:5e2,function(x){
  rn_pol<-bed_shuffle(tmp_cl_tbl,genome = hg19_coord,max_tries=1e8,within=T)
  rn_GRange<-GRanges(seqnames=rn_pol$chrom,
                     ranges = IRanges(start=rn_pol$start,
                                      end=rn_pol$end))
  return(length(unique(subjectHits(findOverlaps(rn_GRange,taSNP)))))
})
plan(sequential)

obs_count<-length(unique(subjectHits(findOverlaps(hub_GRanges_l[[1]],taSNP))))
pnorm(obs_count,mean = mean(traseR_SNP_count),sd = sd(traseR_SNP_count),lower.tail = F)
