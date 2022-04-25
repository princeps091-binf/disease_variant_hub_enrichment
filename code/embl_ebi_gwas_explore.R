library(tidyverse)
library(vroom)

ebi_gwas<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/EMBL_EBI/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv")

ebi_gwas %>% 
  filter(`OR or BETA` < 140000004000000) %>% 
  ggplot(.,aes(`OR or BETA`,-log10(`P-VALUE`))) + 
  geom_point(size=0.1)+
  scale_x_log10()

ebi_gwas %>% 
  group_by(MAPPED_TRAIT) %>% 
  summarise(n=n()) %>% 
#  filter(grepl("[B-b]reast ca",MAPPED_TRAIT)) %>% 
#  filter(grepl("[I-i]mmu|[L-l]ymph",MAPPED_TRAIT)) %>% 
  ggplot(.,aes(n)) + 
  geom_density()+
  scale_x_log10()


