library(tidyverse)

# An error was discovered in 1.0. 
# Agaricomycotina was wrongly inserted in the subdivision column for 8 of cryptococcus sp. 
# Update to 1.1 involves correcting for this

ROD_v1.0 <- read_rds("./../04_github/ROD/ROD_v1.0_operon_variants.rds") %>% as_tibble()
ROD_v1.0_genome_stats <- read_rds("./../04_github/ROD/ROD_v1.0_genome_statistics.rds") %>% as_tibble()

# Locating the errors: 
colnames(ROD_v1.0)
ROD_v1.0 %>% filter(subdivision == "Agaricomycotina")
# Also in the lineage!
ROD_v1.0 %>% filter(subdivision == "Agaricomycotina") %>% select(lineage)



# And for the genome statistics: 
ROD_v1.0_genome_stats %>% filter(subdivision == "Agaricomycotina")
ROD_v1.0_genome_stats %>% filter(subdivision == "Agaricomycotina") %>% select(lineage) 

ROD_v1.0 %>% filter(if_any(where(is.character), ~str_detect(., "Agaricomycotina")))
ROD_v1.0_genome_stats %>% filter(if_any(where(is.character), ~str_detect(., "Agaricomycotina")))

# 8 Occurrences in two columns needs to be changed
# First the character columns (i.e "lineage")
ROD_v1.1 <- ROD_v1.0 %>%
  mutate(across(where(is.character), ~str_replace_all(., "Agaricomycotina", "Basidiomycota")))
ROD_v1.1_genome_stats <- ROD_v1.0_genome_stats %>%
  mutate(across(where(is.character), ~str_replace_all(., "Agaricomycotina", "Basidiomycota")))
# Then the factor columns: 
ROD_v1.1 <- ROD_v1.1 %>%
  mutate(subdivision = factor(str_replace_all(as.character(subdivision), "Agaricomycotina", "Basidiomycota")))
ROD_v1.1_genome_stats <- ROD_v1.1_genome_stats %>%
  mutate(subdivision = factor(str_replace_all(as.character(subdivision), "Agaricomycotina", "Basidiomycota")))



# Make sure there are no Agaricomycotina left
ROD_v1.1 %>% filter(if_any(where(is.character), ~str_detect(., "Agaricomycotina")))
ROD_v1.1 %>% filter(subdivision == "Agaricomycotina")
ROD_v1.1_genome_stats %>% filter(if_any(where(is.character), ~str_detect(., "Agaricomycotina")))
ROD_v1.1_genome_stats %>% filter(subdivision == "Agaricomycotina")


# Count Basidiomycota before and after (should increase with 8)
ROD_v1.0 %>% filter(subdivision == "Basidiomycota")
# 1,215
ROD_v1.0_genome_stats %>% filter(subdivision == "Basidiomycota")
# 764
ROD_v1.1 %>% filter(subdivision == "Basidiomycota")
# 1,215
ROD_v1.1_genome_stats %>% filter(subdivision == "Basidiomycota")
# 764

saveRDS(ROD_v1.1, "./../04_github/ROD/ROD_v1.1_operon_variants.rds")
write.table(ROD_v1.1, "./../04_github/ROD/ROD_v1.1_operon_variants.tab", quote = F, sep = "\t", row.names = F)
saveRDS(ROD_v1.1_genome_stats, "./../04_github/ROD/ROD_v1.1_genome_statistics.rds")
write.table(ROD_v1.1_genome_stats, "./../04_github/ROD/ROD_v1.1_genome_statistics.tab", quote = F, sep = "\t", row.names = F)

