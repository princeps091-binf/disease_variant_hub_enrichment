library(GenomicRanges)
library(traseR)
library(tidyverse)
library(VariantAnnotation)
library(scales)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
vcf_file<-"~/Documents/multires_bhicect/data/epi_data/VCF/clinvar2022_04_17.vcf"
vcf <- readVcf(vcf_file, "hg19")
test<-tibble(as.data.frame(info(vcf)))
test<-test %>% 
  dplyr::select(CLNSIG,CLNDN,ORIGIN) %>% 
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

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno_breast_clinvar <- annotatePeak(breast_var_GRange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
peakAnno_immune_clinvar <- annotatePeak(immune_var_GRange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat
peakAnno_clinvar <- annotatePeak(var_GRanges[sample(1:length(var_GRanges),2e4)], tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)@annoStat

peakAnno_breast_clinvar %>% 
  as_tibble %>%
  mutate(set="Breast") %>% 
  bind_rows(.,
            peakAnno_clinvar %>% 
              as_tibble %>%
              mutate(set="All")) %>% 
  bind_rows(.,
            peakAnno_immune_clinvar %>% 
              as_tibble %>%
              mutate(set="Immune")) %>%
  ggplot(.,aes(set,Frequency,fill=Feature))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Paired")
ggsave("~/Documents/multires_bhicect/weeklies/weekly56/img/clinvar_annotation.png")
#------------------------------------------
#Viz of genome distributions for variants

as.data.frame(ranges(immune_var_GRange)) %>% 
  mutate(chr=as.character(seqnames(immune_var_GRange))) %>%
  mutate(set="Immune") %>% 
  bind_rows(.,
            as.data.frame(ranges(breast_var_GRange)) %>% 
              mutate(chr=as.character(seqnames(breast_var_GRange))) %>%
              mutate(set="breast")) %>% 
  bind_rows(.,
            as.data.frame(ranges(var_GRanges)) %>% 
              mutate(chr=as.character(seqnames(var_GRanges))) %>%
              mutate(set="ALL")) %>% 
  mutate(chr=fct_relevel(chr,paste0("chr",c(1:22,"X","Y","M")))) %>% 
  ggplot(.,aes(start,chr,color=set))+
  geom_point(size=0.2)+
  scale_x_continuous(labels= label_number(scale = 1/1e6,suffix="Mb"))+
  facet_wrap(set~.)
ggsave("~/Documents/multires_bhicect/weeklies/weekly56/img/clinvar_genome_dist.png",width = 25,height=25,units = "cm",dpi = 500)
immune_var_tbl %>% 
  dplyr::select(ID,ORIGIN) %>% 
  unnest(cols=c(ORIGIN)) %>% 
  group_by(ORIGIN) %>% 
  summarise(n=n()) %>% 
  mutate(set="immune") %>% 
  bind_rows(.,
            breast_var_tbl %>% 
              dplyr::select(ID,ORIGIN) %>% 
              unnest(cols=c(ORIGIN)) %>% 
              group_by(ORIGIN) %>% 
              summarise(n=n()) %>%
              mutate(set="breast")) %>% 
  filter(!(grepl("^_",ORIGIN))) %>% 
  filter(n>0) %>% 
  ggplot(.,aes(set,n,fill=ORIGIN))+
  geom_bar(stat="identity",position="fill")+
  theme_minimal()
  scale_fill_brewer(palette="Set3")
