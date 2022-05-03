library(GenomicRanges)
library(vroom)
library(furrr)
library(valr)
library(rtracklayer)
library(tidyverse)
library(formattable)
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

GO_set_enrich_fn<-function(cl_set_gene,cage_active_genes_vec,GOBP_set){
  fn_env<-environment()
  
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    print('node ready')
  })
  clusterExport(cl,c('cl_set_gene','cage_active_genes_vec'),envir = fn_env)
  go_pval<-parLapply(cl,GOBP_set,function(tmp_set){
    hitInSample<-sum(cl_set_gene %in% tmp_set)
    sampleSize<-length(cl_set_gene)
    hitInPop<-sum(cage_active_genes_vec %in% tmp_set)
    failInPop<-length(cage_active_genes_vec) - hitInPop
    p_val<-phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    OR_GO<-(hitInSample/sampleSize)/(hitInPop/length(cage_active_genes_vec))
    return(tibble(p.val=p_val,OR=OR_GO,in.gene=hitInSample))
  })
  stopCluster(cl)
  rm(cl)
  path_tbl<-do.call(bind_rows,go_pval)%>%mutate(Gene.Set=names(go_pval),FDR=p.adjust(p.val,method='fdr'))%>%dplyr::select(Gene.Set,FDR,OR,in.gene)
  return(path_tbl)
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
#-------------------------------------------------------------------------------------------------------
#Load EBI-GWAS data
ebi_gwas<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/EMBL_EBI/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv")
ebi_gwas_coord_tbl<-ebi_gwas %>% 
  filter(!(is.na(CHR_POS))) %>% 
  mutate(CHR_ID=paste0("chr",CHR_ID)) %>% 
  distinct(CHR_ID,CHR_POS,SNPS)

ebi_Grange<-   GRanges(seqnames=ebi_gwas_coord_tbl$CHR_ID,
                       ranges = IRanges(start=as.numeric(ebi_gwas_coord_tbl$CHR_POS),
                                        end=as.numeric(ebi_gwas_coord_tbl$CHR_POS),
                       ))
mcols(ebi_Grange)<-tibble(SNP=ebi_gwas_coord_tbl$SNPS)
ch = import.chain("~/Documents/multires_bhicect/data/epi_data/hg38ToHg19.over.chain")
ebi_Grange_hg19<-unlist(liftOver(ebi_Grange, ch))


pheno_set<-unique(unlist(map(ebi_gwas$MAPPED_TRAIT_URI,function(x){
  unlist(strsplit(x,", "))
})))

plan(multisession,workers=5)
pheno_snp_set_l<-future_map(pheno_set,function(x){
  ebi_gwas %>% 
    filter(!(is.na(CHR_POS))) %>% 
    filter( MAPPED_TRAIT_URI %in% x) %>% 
    dplyr::select(SNPS) %>% distinct %>% unlist
})
plan(sequential)
names(pheno_snp_set_l)<-pheno_set

tot_snp_set<-unique(unlist(pheno_snp_set_l))

snp_enrich_l<-lapply(seq_along(hub_GRanges_l),function(x){
  message(names(hub_GRanges_l)[x])
  hub_snp<-mcols(ebi_Grange_hg19)$SNP[unique(subjectHits(findOverlaps(hub_GRanges_l[[x]],ebi_Grange_hg19)))]
  return(GO_set_enrich_fn(hub_snp,tot_snp_set,pheno_snp_set_l))
  
})

formattable(snp_enrich_l[[3]] %>% 
  filter(FDR<=0.01) %>% arrange(FDR) %>% 
  left_join(ebi_gwas %>% distinct(MAPPED_TRAIT_URI,MAPPED_TRAIT),by=c("Gene.Set"="MAPPED_TRAIT_URI")) %>% 
    dplyr::select(MAPPED_TRAIT,FDR,OR,in.gene) %>% 
    arrange(FDR) %>% 
    mutate(FDR = round(FDR,digits = 20)) %>% 
    slice_head(n=15))

#-----------------------------------------------------------
# Enrichment by parent pheno categories
EFO_map<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/EMBL_EBI/gwas-efo-trait-mappings.tsv")
pheno_categ<-unique(EFO_map$`Parent term`)
efo_set_l<-map(pheno_categ,function(x){
return(EFO_map %>% 
  filter(`Parent term` == x) %>% 
  dplyr::select(`EFO URI`) %>% ungroup %>% distinct %>% unlist)
})

plan(multisession,workers=5)
pheno_efo_categ_snp_set_l<-future_map(efo_set_l,function(x){
  ebi_gwas %>% 
    filter(!(is.na(CHR_POS))) %>% 
    filter( MAPPED_TRAIT_URI %in% x) %>% 
    dplyr::select(SNPS) %>% distinct %>% unlist
})
plan(sequential)

names(pheno_efo_categ_snp_set_l)<-pheno_categ

tot_snp_set<-unique(unlist(pheno_efo_categ_snp_set_l))

snp_categ_enrich_l<-lapply(seq_along(hub_GRanges_l),function(x){
  message(names(hub_GRanges_l)[x])
  hub_snp<-mcols(ebi_Grange_hg19)$SNP[unique(subjectHits(findOverlaps(hub_GRanges_l[[x]],ebi_Grange_hg19)))]
  return(GO_set_enrich_fn(hub_snp,tot_snp_set,pheno_efo_categ_snp_set_l))
  
})

formattable(snp_categ_enrich_l[[1]] %>% 
              arrange(FDR) %>% 
              mutate(FDR = round(FDR,digits = 20)) %>% 
              filter(FDR<=0.01) %>% arrange(FDR))


#-----------------------------------------------------------
hub_snp_l<-lapply(hub_GRanges_l,function(x){
  mcols(ebi_Grange_hg19)$SNP[unique(subjectHits(findOverlaps(x,ebi_Grange_hg19)))]
})
names(hub_snp_l)<-names(hub_GRanges_l)
library(UpSetR)
UpSetR::upset(fromList(hub_snp_l),order.by = "freq")
univ_hub_snp<-fromList(hub_snp_l) %>%
  mutate(SNP=unique(unlist(hub_snp_l))) %>% 
  filter(H1 ==1 & HMEC == 1 & GM12878 == 1) %>% 
  dplyr::select(SNP) %>% 
  unlist

h1_hub_snp<-fromList(hub_snp_l) %>%
  mutate(SNP=unique(unlist(hub_snp_l))) %>% 
  filter(H1 ==1 & HMEC == 0 & GM12878 == 0) %>% 
  dplyr::select(SNP) %>% 
  unlist

hmec_hub_snp<-fromList(hub_snp_l) %>%
  mutate(SNP=unique(unlist(hub_snp_l))) %>% 
  filter(H1 ==0 & HMEC == 1 & GM12878 == 0) %>% 
  dplyr::select(SNP) %>% 
  unlist

gm12878_hub_snp<-fromList(hub_snp_l) %>%
  mutate(SNP=unique(unlist(hub_snp_l))) %>% 
  filter(H1 ==0 & HMEC == 0 & GM12878 == 1) %>% 
  dplyr::select(SNP) %>% 
  unlist

tot_snp_set<-unique(unlist(pheno_snp_set_l))
snp_enrich_tbl<-GO_set_enrich_fn(h1_hub_snp,tot_snp_set,pheno_snp_set_l)
print(snp_enrich_tbl %>% 
        filter(FDR<=0.01) %>% arrange(desc(OR)) %>% 
        left_join(ebi_gwas %>% distinct(MAPPED_TRAIT_URI,MAPPED_TRAIT),by=c("Gene.Set"="MAPPED_TRAIT_URI")),n=200)
