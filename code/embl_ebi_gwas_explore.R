library(tidyverse)
library(vroom)
library(furrr)
ebi_gwas<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/EMBL_EBI/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv")
ebi_gwas_meta<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/EMBL_EBI/gwas_catalog_v1.0.2-studies_r2022-04-07.tsv")
ebi_gwas %>% 
#  filter(`OR or BETA` < 140000004000000) %>% 
  ggplot(.,aes(`OR or BETA`,-log10(`P-VALUE`))) + 
  geom_point(size=0.1)+
  scale_x_log10()

ebi_gwas %>% 
  group_by(`DISEASE/TRAIT`) %>% 
  summarise(n=n()) %>% 
  filter(grepl("[B-b]reast ca",MAPPED_TRAIT)) %>% 
#  filter(grepl("[I-i]mmu|[L-l]ymph|[L-l]euk",MAPPED_TRAIT)) %>% 
  ggplot(.,aes(n)) + 
  geom_density()+
  scale_x_log10()


EFO_map<-vroom("~/Documents/multires_bhicect/data/epi_data/VCF/EMBL_EBI/gwas-efo-trait-mappings.tsv")

EFO_map %>% 
  group_by(`Parent term`) %>% 
  distinct %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

efo_tbl<-EFO_map %>% 
  group_by(`EFO term`) %>% 
  dplyr::select(`EFO term`,`EFO URI`) %>% distinct

efo_tbl<-efo_tbl %>% 
  mutate(EFO=str_split(`EFO URI`,"/",simplify = T)[,5])


Cancer_efo_tbl<-EFO_map %>% 
  filter(`Parent term` == "Cancer") %>% 
  group_by(`EFO term`) %>% 
  dplyr::select(`EFO term`,`EFO URI`) %>% distinct

Cancer_efo_tbl<-Cancer_efo_tbl %>% 
  mutate(EFO=str_split(`EFO URI`,"/",simplify = T)[,5])

plan(multisession,workers=4)
Cancer_efo_tbl %>%
  mutate(n.variant=future_map_int(EFO,function(x){
  
  nrow(ebi_gwas %>% filter(grepl(x,MAPPED_TRAIT_URI)) %>% distinct)
  
}))
plan(sequential)
