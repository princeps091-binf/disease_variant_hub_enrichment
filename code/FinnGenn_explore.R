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
GWAS_file<-"~/Documents/multires_bhicect/data/epi_data/VCF/FinnGen/summary_stats_finngen_R6_AUTOIMMUNE.gz"
#gwas_tbl<-vroom(GWAS_file)
gwas_tbl<-vroom(GWAS_file, col_select = c(1,2,7,9,10))
gwas_tbl<-vroom(GWAS_file, col_select = c(1,2))

gwas_tbl %>% 
  sample_n(1e6) %>% 
  ggplot(.,aes(beta,-log10(pval))) + 
  geom_point(alpha=0.1,size=0.1)+
  theme_minimal()
ggsave("~/Documents/multires_bhicect/weeklies/weekly56/img/Finngen_GWAS_volcano.png")

gg_dist<-gwas_tbl %>% 
  filter(`#chrom`==22) %>% 
  mutate(eff=sign(beta)* -log10(pval),dir=sign(beta)) %>% 
  ggplot(.,aes(pos,eff,color=dir))+
  geom_point(size=0.01)+
  theme(legend.position="none")
ggsave("~/Documents/multires_bhicect/weeklies/weekly56/img/Finngen_GWAS_gdist_chr22.png",width = 40,height = 23,units = "cm",gg_dist)
