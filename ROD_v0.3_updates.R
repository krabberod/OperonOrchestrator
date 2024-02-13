# Standardizing NCBI taxonomy:
library(taxonomizr)
library(tidyverse)

accessionTaxa.sql <- "./../../../99_div/accessionTaxa.sql"
ROD_v0.2 <- readRDS("ROD_with_seq_v0.2.rds")
ROD_v0.2$taxid
# Looks like it's ok. I have updated the taxonomy 
updated_taxonomy <- readxl::read_xlsx("./../assembly_summary_genbank.taxonomy_revised_V2.xlsx")
# saveRDS(ROD_v0.2, "ROD_v0.2_with_blast.rds")


# add length: 
nchar(ROD_v0.2$sequence)

ROD_v0.3 <- ROD_v0.2 %>% select(seqid,assembly_id,size, taxid, sequence)
ROD_v0.3$length <- nchar(ROD_v0.3$sequence)
ROD_v0.3 <- ROD_v0.3 %>% relocate(seqid, size, assembly_id,taxid, length ,sequence)
ROD_v0.3 <- left_join(ROD_v0.3,updated_taxonomy)
ROD_v0.3$lineage <- paste(ROD_v0.3$Domain, ROD_v0.3$supergroup, ROD_v0.3$division,
                          ROD_v0.3$subdivision,
                          ROD_v0.3$class, ROD_v0.3$order, 
                          ROD_v0.3$family, ROD_v0.3$genus, 
                          ROD_v0.3$species, 
                         sep = ";") %>% str_replace_all(" ","_")
saveRDS(ROD_v0.3, "ROD_v0.3.rds")


### EXTRACT FASTA
df <- readRDS("ROD_v0.3")
df <- ROD_v0.3
header <- paste0(df$seqid,";","size=",df$size,"|",df$lineage , " taxid=",df$taxid, ";")
# header <- paste0(df$seqid,";","size=",df$size, " taxid=",df$taxid, ";")
Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
Xfasta[c(FALSE, TRUE)] <- df$sequence
writeLines(Xfasta, "ROD_v0.3.fasta")
getwd()

###
ROD_v0.3 <-readRDS("ROD_v0.3.rds")
ROD_v0.3 %>% filter(length >= 4000) %>% select(assembly_id) %>% unique()
sum(as.numeric(ROD_v0.3$size))




