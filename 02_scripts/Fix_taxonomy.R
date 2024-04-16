# Fix missing taxnomy in Goat: 
library(taxonomizr)
library(tidyverse)
library(readxl)

accessionTaxa.sql <- "./../../../99_div/accessionTaxa.sql"
assembly_info <- read.table("./../assembly_summary_genbank.taxid.txt", header = T, sep = "\t")
# taxaID <- read.table("./../kol7", header = F, stringsAsFactors = F)
# colnames(assembly_info) <- c("assembly_id","taxid")
assembly_info$taxid


# Add desired taxa: 
desTax = c("superkingdom","kingdom","phylum", "subphylum", "class", "order", "family", "genus", "species")
taxa <- getTaxonomy(assembly_info$taxid, accessionTaxa.sql, desiredTaxa = desTax) %>% 
  as_tibble()
# getTaxonomy(94885,accessionTaxa.sql)

tail(taxa)

assembly_info_taxa <- cbind(assembly_info,taxa)
# write.table(assembly_info_taxa, "./../assembly_summary_genbank.taxid_taxa.txt", col.names  = TRUE, 
#             row.names = FALSE, sep = "\t", quote = FALSE)

as_tibble(assembly_info_taxa)
assembly_bact <- assembly_info_taxa %>% filter(superkingdom == "Bacteria")
assembly_arch <- assembly_info_taxa %>% filter(superkingdom == "Archaea")
assembly_info_taxa <- assembly_info_taxa %>% filter(superkingdom != "Bacteria") %>% 
  filter(superkingdom != "Viruses") %>% filter(superkingdom != "Archaea")
assembly_info_taxa %>% select(kingdom) %>% unique()


# Find phylums without Kingdom in NCBI
assembly_info_taxa %>% 
  filter(is.na(kingdom)) %>% 
  select(genus) %>% unique()

# Some without phylum also lack kingdom
assembly_info_taxa %>% 
  filter(is.na(kingdom)) %>% 
  filter(is.na(phylum)) %>% 
  select(family) %>% unique()


# Etc.
# But I mainly need to fix the taxonomy for the groups that are in ROD. 
# So I filter only the needed genomes. Or I could fix it all... 

assembly_info_taxa %>% filter(phylum=="Apicomplexa") 

### 
# First the obvious "kingdoms"
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Apicomplexa", "Alveolata", kingdom)) 
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Fornicata", "Metamonada", kingdom)) 
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Heterolobosea", "Discoba", kingdom)) 
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Perkinsozoa", "Alveolata", kingdom)) 
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Oomycota", "Stramenopiles", kingdom)) 
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Rhodophyta", "Rhodophyta", kingdom)) 
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Ciliophora", "Alveolata", kingdom)) 
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Cercozoa", "Rhizaria", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Evosea", "Amoebozoa", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Euglenozoa", "Discoba", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Discosea", "Amoebozoa", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Foraminifera", "Rhizaria", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Bacillariophyta", "Stramenopiles", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Haptophyta", "Haptista", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Endomyxa", "Rhizaria", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Parabasalia", "Metamonada", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Preaxostyla", "Metamonada", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Imbricatea", "Rhizaria", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(phylum == "Tubulinea", "Amoebozoa", kingdom))

# And some "class"
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(class == "Bigyra", "Stramenopiles ", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(class == "Hyphochytriomycetes", "Stramenopiles", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(class == "Dinophyceae", "Alveolata", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(class == "Bolidophyceae", "Stramenopiles ", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(kingdom = ifelse(class == "Synurophyceae", "Stramenopiles ", kingdom))
assembly_info_taxa <- assembly_info_taxa %>% mutate(phylum = ifelse(class == "Dinophyceae", "Dinoflagellata", phylum))


assembly_info_taxa %>% filter(is.na(phylum))
assembly_info_taxa %>% filter(is.na(kingdom)) %>% select(class) %>% unique()
assembly_info_taxa %>% filter(kingdom=="Metazoa")
assembly_info_taxa %>% filter(class=="Dinophyceae")

write.csv(x = assembly_info_taxa, "./../assembly_summary_genbank.taxonomy.txt", row.names= FALSE, quote = FALSE)
# The work on this has been continued in Excel: 
# assembly_summary_genbank_taxonomy.xls


#### some are missing: 
missing_taxid <- read.table("./../list_of_missing", header=F)

desTax = c("superkingdom","kingdom","phylum", "subphylum", "class", "order", "family", "genus", "species")
taxa_missing <- getTaxonomy(missing_taxid$V2, accessionTaxa.sql, desiredTaxa = desTax) %>%as_tibble()
df_missing <- cbind(missing_taxid, taxa_missing)
write.table(df_missing, "./../list_of_missing_tax.csv", row.names = FALSE, quote = FALSE, sep="\t")
## updated in the file assembly_summary_genbank_taxonomy_revised
  