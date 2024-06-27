# Script for updating the ROD from version 1.1 to 1.2
library(tidyverse)

# Read in the ROD_v1.1_operon_variants.rds and ROD_v1.1_genome_statistics.rds
# getwd()
# setwd("/Users/anderkkr/Library/CloudStorage/OneDrive-UniversitetetiOslo/00_Master_projects_OD/01_Active_Projects_OD/02_OMG/17_OperonOrchestrator/OperonOrchestrator/02_scripts")
ROD_v1.1 <- read_rds("./../04_github/ROD/ROD_v1.1_operon_variants.rds") %>% as_tibble()
ROD_v1.1_genome_stats <- read_rds("./../04_github/ROD/ROD_v1.1_genome_statistics.rds") %>% as_tibble()


# Update the ROD_v1.1 to ROD_v1.2
# Find which entries with wrongly spelled "Ochrophyta" in the subdivision column. It is spelled "Oochropyhta" in some instances!
# Filter on Stramenopiles in divison, then print unique values in subdivision 
ROD_v1.1_genome_stats %>% filter(division == "Stramenopiles") %>% select(subdivision) %>% unique()
# and check if "Oochrophyta" is present
# find any occurence if "Oochrophyta" in ROD_v1.1 main dataframe
ROD_v1.1 %>% filter(str_detect(subdivision, "Oochropyhta"))
# find any occurence if "Oochrophyta" in ROD_v1.1_genome_stats
ROD_v1.1_genome_stats %>% filter(str_detect(subdivision, "Oochropyhta"))
# find any occurence if "Oochropyhta" in any colum, idetnify  the column it is found in
ROD_v1.1_genome_stats %>% filter(if_any(where(is.character), ~str_detect(., "Oochropyhta")))

# Replace "Oochrophyta" with "Ochrophyta" in the ROD_v1.1 and ROD_v1.1_genome_stats
ROD_v1.2 <- ROD_v1.1 %>%
  mutate(across(where(is.factor), as.character), # Convert factors to characters first
         across(where(is.character), ~str_replace_all(., "Oochropyhta", "Ochrophyta")),
         across(where(is.character), as.factor)) # Optionally convert back to factors

ROD_v1.2_genome_stats <- ROD_v1.1_genome_stats %>%
  mutate(across(where(is.factor), as.character), # Convert factors to characters first
         across(where(is.character), ~str_replace_all(., "Oochropyhta", "Ochrophyta")),
         across(where(is.character), as.factor)) # Optionally convert back to factors
# Check if the replacement was successful
ROD_v1.2 %>% filter(division == "Stramenopiles") %>% select(subdivision) %>% unique()
ROD_v1.2 %>% filter(if_any(where(is.character), ~str_detect(., "Oochropyhta")))
ROD_v1.2_genome_stats %>% filter(if_any(where(is.character), ~str_detect(., "Oochropyhta")))


# Save the updated ROD_v1.2 and ROD_v1.2_genome_stats
ROD_v1.2$sequence <- as.character(ROD_v1.2$sequence)
ROD_v1.2$assembly_id <- as.character(ROD_v1.2$assembly_id)
ROD_v1.2$lineage <- as.character(ROD_v1.2$lineage)
ROD_v1.2_genome_stats$assembly_id <- as.character(ROD_v1.2_genome_stats$assembly_id)
ROD_v1.2_genome_stats$all_sizes <- as.character(ROD_v1.2_genome_stats$all_sizes)
ROD_v1.2_genome_stats$all_lengths <- as.character(ROD_v1.2_genome_stats$all_lengths)

str(ROD_v1.2_genome_stats)

saveRDS(ROD_v1.2, "./../04_github/ROD/ROD_v1.2_operon_variants.rds")
write.table(ROD_v1.2, gzfile("./../04_github/ROD/ROD_v1.2_operon_variants.tab.gz"), quote = F, sep = "\t", row.names = F)
saveRDS(ROD_v1.2_genome_stats, "./../04_github/ROD/ROD_v1.2_genome_statistics.rds")
write.table(ROD_v1.2_genome_stats, gzfile("./../04_github/ROD/ROD_v1.2_genome_statistics.tab.gz"), quote = F, sep = "\t", row.names = F)

# write assembly_id as character in tsv file
write.table(ROD_v1.2_genome_stats$assembly_id, "ROD_v1.2_genome_statistics_assembly_id.tsv", quote = F, sep = "\t", row.names = F)

read.table("ROD_v1.2.NCBI_metadata.tab", header = F, sep = "\t") %>% str()

