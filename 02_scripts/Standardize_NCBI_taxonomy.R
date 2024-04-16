# Standardizing NCBI taxonomy:
library(taxonomizr)
library(tidyverse)

accessionTaxa.sql <- "./../../../99_div/accessionTaxa.sql"
ROD_v0.2 <- readRDS("ROD_with_seq_v0.2.rds")
ROD_v0.2$taxid
taxa_raw <- getRawTaxonomy(ROD_v0.2$taxid[1:100], accessionTaxa.sql)
# Trying to use the normalise function
desTax = c("superkingdom","kingdom","phylum", "subphylum", "class", "order", "family", "genus", "species")
taxa <- taxa <- getTaxonomy(ROD_v0.2$taxid, accessionTaxa.sql, desiredTaxa = desTax) %>% as_tibble()
unique(taxa$phylum)
unique(taxa$kingdom)
unique(ROD_v0.2$kingdom)
unique(ROD_v0.2$superkingdom)
unique(ROD_v0.2$phylum)
cbind(taxa$phylum,ROD_v0.2$phylum) %>% unique()
cbind(taxa$class,ROD_v0.2$class) %>% unique()
cbind(taxa$order,ROD_v0.2$order) %>% unique()
cbind(taxa$order,ROD_v0.2$order) %>% unique()
cbind(taxa$kingdom,ROD_v0.2$supergroup,ROD_v0.2$kingdom) %>% unique()

# Looks like it's ok. I have updated the taxonomy 
updated_taxonomy <- readxl::read_xlsx("./../assembly_summary_genbank.taxonomy_revised_V2.xlsx")
saveRDS(ROD_v0.2, "ROD_v0.2_with_blast.rds")


# add length: 
nchar(ROD_v0.2$sequence)

ROD_v0.3 <- ROD_v0.2 %>% select(seqid,assembly_id,size, taxid, org.copy.number ,sequence)
ROD_v0.3$length <- nchar(ROD_v0.3$sequence)
ROD_v0.3 <- ROD_v0.3 %>% relocate(seqid,assembly_id,taxid,size, length, org.copy.number ,sequence)
saveRDS(ROD_v0.2, "ROD_v0.3.rds")
le







