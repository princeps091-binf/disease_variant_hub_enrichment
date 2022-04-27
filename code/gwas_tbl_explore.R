library(GenomicRanges)
library(tidyverse)
library(VariantAnnotation)
library(vroom)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
GWAS_file<-"~/Documents/multires_bhicect/data/epi_data/VCF/PLCO_GWAS_explorer/j_breast_cancer.tsv.gz"
gwas_tbl<-vroom(GWAS_file,col_select = c(1,2,13,15)) %>% 
  filter(CHR!=23)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)  <- paste0("chr",unique(gwas_tbl$CHR))

rn_gwas_ann<-function(txdb,gwas_tbl,n_samp){
  tmp_tbl<-gwas_tbl %>% 
    sample_n(n_samp)
  
  gwas_Grange<-   GRanges(seqnames=paste0("chr",tmp_tbl$CHR),
                          ranges = IRanges(start=as.numeric(tmp_tbl$POS),
                                           end=as.numeric(tmp_tbl$POS)
                          ))
  peakAnno_gwas <- annotatePeak(gwas_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
  return(peakAnno_gwas)
}
do.call(bind_rows,lapply(1:10,function(x)rn_gwas_ann(txdb,gwas_tbl,1e5))) %>% 
  group_by(Feature) %>% summarise(n=mean(Frequency))


gwas_tbl %>% 
  sample_n(1e6) %>% 
  ggplot(.,aes(BETA_European_all,-log10(P_European_all))) + 
  geom_point(alpha=0.1,size=0.1)+
  theme_minimal()
ggsave("~/Documents/multires_bhicect/weeklies/weekly56/img/PLCO_GWAS_volcano.png")

gg_dist<-gwas_tbl %>% 
  filter(CHR==22) %>% 
  mutate(eff=sign(BETA_European_all)* -log10(P_European_all),dir=sign(BETA_European_all)) %>% 
  ggplot(.,aes(POS,eff,color=dir))+
  geom_point(size=0.01)+
  facet_grid(CHR~.,scales="free")+
  theme(legend.position="none")
ggsave("~/Documents/multires_bhicect/weeklies/weekly56/img/PLCO_GWAS_gdist_chr22.png",width = 40,height = 23,units = "cm",gg_dist)


unique(summary_stats_R6_manifest$category)
summary_stats_R6_manifest %>% 
  filter(category %in% c("III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)", "Diseases marked as autimmune origin")) %>% 
  arrange(desc(n_cases))
