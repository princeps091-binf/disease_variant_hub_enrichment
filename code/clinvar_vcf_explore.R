library(GenomicRanges)
library(traseR)
library(tidyverse)
library(VariantAnnotation)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
data(taSNP)
vcf_file<-"~/Documents/multires_bhicect/data/epi_data/VCF/clinvar_2022_04_16.vcf"
vcf <- readVcf(vcf_file, "hg19")
test<-tibble(as.data.frame(info(vcf)))
test<-test %>% 
  dplyr::select(ALLELEID,CLNSIG,CLNDN,ORIGIN) %>% 
  mutate(ID=names(rowRanges(vcf)))

var_GRanges<-rowRanges(vcf)
var_GRanges<-renameSeqlevels(var_GRanges,mapSeqlevels(seqlevels(var_GRanges),"UCSC"))

breast_var_tbl<-test%>% 
  mutate(breast.var=map_dbl(CLNDN,function(x){
    length(grep('[B-b]reast|[M-m]ammary',x))
  })) %>% 
  filter(breast.var>0)
breast_var_GRange<-rowRanges(vcf)[which(names(rowRanges(vcf)) %in% as.character(breast_var_tbl$ID))]
breast_var_GRange<-renameSeqlevels(breast_var_GRange,mapSeqlevels(seqlevels(breast_var_GRange),"UCSC"))


immune_var_tbl<-test%>% 
  mutate(immune.var=map_dbl(CLNDN,function(x){
    length(grep('[I-i]mmun|[L-l]ymph',x))
  })) %>% 
  filter(immune.var>0)
immune_var_GRange<-rowRanges(vcf)[which(names(rowRanges(vcf)) %in% as.character(immune_var_tbl$ID))]
immune_var_GRange<-renameSeqlevels(immune_var_GRange,mapSeqlevels(seqlevels(immune_var_GRange),"UCSC"))

breast_cancer_gwas<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/j_breast_cancer.tsv.gz",n_max = 2e4)  
breast_cancer_gwas<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/PLCO_GWAS_explorer/j_breast_cancer.tsv.gz",col_select = 1:2)  

breast_gwas_Grange<-   GRanges(seqnames=paste0("chr",breast_cancer_gwas$CHR),
                               ranges = IRanges(start=as.numeric(breast_cancer_gwas$POS),
                                                end=as.numeric(breast_cancer_gwas$POS)
                               ))

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno_breast_clinvar <- annotatePeak(breast_var_GRange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
peakAnno_immune_clinvar <- annotatePeak(immune_var_GRange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
peakAnno_clinvar <- annotatePeak(var_GRanges[sample(1:length(var_GRanges),2e4)], tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
peakAnno_traser <- annotatePeak(taSNP, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
peakAnno_gwas <- annotatePeak(breast_gwas_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat

peakAnno_breast_clinvar %>% 
  as_tibble %>%
  mutate(set="Breast") %>% 
  bind_rows(.,
            peakAnno_clinvar %>% 
              as_tibble %>%
              mutate(set="All")) %>% 
  bind_rows(.,
            peakAnno_traser %>% 
              as_tibble %>%
              mutate(set="traseR")) %>%
  bind_rows(.,
            peakAnno_immune_clinvar %>% 
              as_tibble %>%
              mutate(set="Immune")) %>%
  bind_rows(.,
            peakAnno_gwas %>% 
              as_tibble %>%
              mutate(set="Breast_GWAS")) %>%
  
  ggplot(.,aes(set,Frequency,fill=Feature))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Paired")
#------------------------------------------
#Viz of genome distributions for variants

as.data.frame(ranges(immune_var_GRange)) %>% 
  mutate(chr=as.character(seqnames(immune_var_GRange))) %>%
  mutate(set="Immune") %>% 
  bind_rows(.,
            as.data.frame(ranges(breast_var_GRange)) %>% 
              mutate(chr=as.character(seqnames(breast_var_GRange))) %>%
              mutate(set="breast")) %>% 
  mutate(chr=fct_relevel(chr,paste0("chr",c(1:22,"X","Y","M")))) %>% 
  ggplot(.,aes(start,chr,color=set))+
  geom_point(size=0.5)+
  facet_wrap(set~.)

gg_tmp<-breast_cancer_gwas %>% 
mutate(CHR=paste0("chr",CHR)) %>%
  ggplot(.,aes(POS,CHR))+
  geom_point(size=0.5)
ggsave("./breast_gwas_coord.png",gg_tmp)
