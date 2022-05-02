library(GenomicRanges)
library(tidyverse)
library(vroom)
library(furrr)
library(valr)
library(rtracklayer)
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
#-------------------------------------------------------------------------------------------------------
#Load EBI-GWAS data
ebi_gwas<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/EMBL_EBI/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv")

EFO_map<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/EMBL_EBI/gwas-efo-trait-mappings.tsv")
Cancer_efo_tbl<-EFO_map %>% 
  filter(`Parent term` == "Cancer") %>% 
  group_by(`EFO term`) %>% 
  dplyr::select(`EFO term`,`EFO URI`) %>% distinct

Cancer_efo_tbl<-Cancer_efo_tbl %>% 
  mutate(EFO=str_split(`EFO URI`,"/",simplify = T)[,5])

Cancer_ebi_gwas<-ebi_gwas %>% 
  filter(MAPPED_TRAIT_URI %in% Cancer_efo_tbl$`EFO URI`) %>% 
  filter(!(is.na(CHR_POS))) %>% 
  mutate(CHR_ID=paste0("chr",CHR_ID))

cancer_ebi_Grange<-   GRanges(seqnames=Cancer_ebi_gwas$CHR_ID,
                            ranges = IRanges(start=as.numeric(Cancer_ebi_gwas$CHR_POS),
                                             end=as.numeric(Cancer_ebi_gwas$CHR_POS)
                            ))
mcols(cancer_ebi_Grange)<-tibble(SNP=Cancer_ebi_gwas$SNPS)

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
cancer_ebi_Grange_hg19<-unlist(liftOver(cancer_ebi_Grange, ch))

ebi_Grange_hg19<-unlist(liftOver(ebi_Grange, ch))

hg19_coord <- read_delim("~/Documents/multires_bhicect/data/hg19.genome", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")
res_l<-vector('list',length(hub_GRanges_l))
names(res_l)<-names(hub_GRanges_l)

tmp_cl_tbl<-hub_GRanges_l[[2]] %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)
gap_bed<-read_delim("~/Documents/multires_bhicect/data/gap.bed", 
                    "\t", escape_double = FALSE, col_names = F, 
                    trim_ws = TRUE) %>% 
  dplyr::select(X2,X3,X4) %>% 
  dplyr::rename(chrom=X2,start=X3,end=X4)

res_l<-vector('list',length(hub_GRanges_l))
names(res_l)<-names(hub_GRanges_l)

for(i in names(hub_GRanges_l)){
  message(i)
  tmp_cl_tbl<-hub_GRanges_l[[i]] %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)
  plan(multisession,workers=5)
  
  rn_SNP_count<-future_map_int(1:1e3,function(x){
                    rn_pol<-bed_shuffle(tmp_cl_tbl,genome = hg19_coord,max_tries=1e8,within=T)
                    rn_GRange<-GRanges(seqnames=rn_pol$chrom,
                        ranges = IRanges(start=rn_pol$start,
                                         end=rn_pol$end))
    return(length(unique(subjectHits(findOverlaps(rn_GRange,ebi_Grange_hg19)))))
  })
  plan(sequential)

obs_count<-length(unique(subjectHits(findOverlaps(hub_GRanges_l[[i]],ebi_Grange_hg19))))
res_l[[i]]<-tibble(count=rn_SNP_count,set="boot") %>% 
         bind_rows(tibble(count=obs_count,set="obs")) %>% 
         mutate(line=i)
}
do.call(bind_rows,res_l) %>% 
  filter(set=="boot") %>% 
  ggplot(.,aes(line,count))+geom_violin()+
  geom_point(data=do.call(bind_rows,res_l) %>% 
               filter(set=="obs") )
ggsave("~/Documents/multires_bhicect/weeklies/weekly57/img/EMBL_GWAS_boot_enrich_all.png")

pnorm(obs_count,mean = mean(rn_SNP_count),sd = sd(rn_SNP_count),lower.tail = F)
#--------------------------
hub_snp<-ebi_Grange_hg19[unique(subjectHits(findOverlaps(hub_GRanges_l[[2]],ebi_Grange_hg19)))]

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno_hub_snp <- annotatePeak(hub_snp, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
peakAnno_snp <- annotatePeak(ebi_Grange_hg19, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
#--------------------------
# Loop through all categories of phenotypes
pheno_categ<-unique(EFO_map$`Parent term`)
ch = import.chain("~/Documents/multires_bhicect/data/epi_data/hg38ToHg19.over.chain")
gap_bed<-read_delim("~/Documents/multires_bhicect/data/gap.bed", 
                    "\t", escape_double = FALSE, col_names = F, 
                    trim_ws = TRUE) %>% 
  dplyr::select(X2,X3,X4) %>% 
  dplyr::rename(chrom=X2,start=X3,end=X4)

hg19_coord <- read_delim("~/Documents/multires_bhicect/data/hg19.genome", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")

res_l<-vector('list',length(hub_GRanges_l))
names(res_l)<-names(hub_GRanges_l)
for(i in names(hub_GRanges_l)){
  tmp_cl_tbl<-hub_GRanges_l[[i]] %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)
  
  enrich_pheno_l<-map(pheno_categ,function(x){
    
    tmp_efo_tbl<-EFO_map %>% 
      filter(`Parent term` == x) %>% 
      group_by(`EFO term`) %>% 
      dplyr::select(`EFO term`,`EFO URI`) %>% distinct
    
    tmp_efo_tbl<-tmp_efo_tbl %>% 
      mutate(EFO=str_split(`EFO URI`,"/",simplify = T)[,5])
    
    tmp_ebi_gwas<-ebi_gwas %>% 
      filter(MAPPED_TRAIT_URI %in% tmp_efo_tbl$`EFO URI`) %>% 
      filter(!(is.na(CHR_POS))) %>% 
      mutate(CHR_ID=paste0("chr",CHR_ID))
    
    tmp_ebi_Grange<-   GRanges(seqnames=tmp_ebi_gwas$CHR_ID,
                               ranges = IRanges(start=as.numeric(tmp_ebi_gwas$CHR_POS),
                                                end=as.numeric(tmp_ebi_gwas$CHR_POS)
                               ))
    mcols(tmp_ebi_Grange)<-tibble(SNP=tmp_ebi_gwas$SNPS)
    tmp_ebi_Grange_hg19<-unlist(liftOver(tmp_ebi_Grange, ch))
    
    message(i," Boot:", x)
    plan(multisession,workers=5)
    rn_SNP_count<-future_map_int(1:5e2,function(x){
      rn_pol<-bed_shuffle(tmp_cl_tbl,genome = hg19_coord,max_tries=1e8,within=T)
      rn_GRange<-GRanges(seqnames=rn_pol$chrom,
                         ranges = IRanges(start=rn_pol$start,
                                          end=rn_pol$end))
      return(length(unique(subjectHits(findOverlaps(rn_GRange,tmp_ebi_Grange_hg19)))))
    })
    plan(sequential)
    obs_count<-length(unique(subjectHits(findOverlaps(hub_GRanges_l[[i]],tmp_ebi_Grange_hg19))))
    return(tibble(count=rn_SNP_count,set="boot",pheno=x) %>% 
             bind_rows(tibble(count=obs_count,set="obs",pheno=x)))
    
  })
  res_l[[i]]<-do.call(bind_rows,enrich_pheno_l) %>% 
    mutate(line=i)
  
}

do.call(bind_rows,res_l) %>% 
  filter(set=="boot") %>% 
  ggplot(.,aes(pheno,count,fill=set))+
  geom_violin()+
  geom_point(data=do.call(bind_rows,res_l) %>% 
               filter(set=="obs"))+
  facet_grid(.~line)+
  coord_flip()
  
ggsave("~/Documents/multires_bhicect/weeklies/weekly57/img/EMBL_GWAS_boot_enrich_trait_class.png",width = 40,height = 23,units = "cm")



enrich_pheno_count_l<-map(pheno_categ,function(x){
    
    tmp_efo_tbl<-EFO_map %>% 
      filter(`Parent term` == x) %>% 
      group_by(`EFO term`) %>% 
      dplyr::select(`EFO term`,`EFO URI`) %>% distinct
    
    tmp_efo_tbl<-tmp_efo_tbl %>% 
      mutate(EFO=str_split(`EFO URI`,"/",simplify = T)[,5])
    
    tmp_tbl<-ebi_gwas %>% 
      filter(MAPPED_TRAIT_URI %in% tmp_efo_tbl$`EFO URI`) %>% 
      filter(!(is.na(CHR_POS))) %>% 
      distinct %>% summarise(n=n()) %>% mutate(set=x)
  })
do.call(bind_rows,enrich_pheno_count_l) %>% 
  mutate(set=fct_reorder(set,n)) %>% 
  ggplot(.,aes(n,set))+
  geom_bar(stat="identity")
ggsave("~/Documents/multires_bhicect/weeklies/weekly57/img/EMBL_GWAS_trait_class_count.png")
