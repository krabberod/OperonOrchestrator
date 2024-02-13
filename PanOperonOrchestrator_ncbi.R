library(tidyverse)
library(Biostrings)
library(ape)
library(stringdist)
#library(DECIPHER)
library(ips)



#### IMPORT GOAT ####
# NB: Taxonomy updated from taxid in 
# "./../assembly_summary_genbank_taxonomy.xlsx"
Goat <- readxl::read_xlsx("Goat.xlsx")
head(Goat)
colnames(Goat)
table(Goat$phylum)

#### IMPORT FROM NCBI ####
fasta_file <- "./../01_results/02_ncbi_results/ribres_4000_cluster.fasta"
# fasta_file <- "./../01_results/02_ncbi_results/ribres_4000_precluster.fasta"

sequences <- readDNAStringSet(fasta_file)

df_ncbi <- NULL
df_ncbi$header <- as.character(names(sequences))
# colnames(df) <- "header"
df_ncbi$sequence <- unname(as.character(sequences))

# colnames(df) <- c("header", "sequence")
str(df_ncbi)
df_ncbi <- as_tibble(df_ncbi)

parts <- str_split(df_ncbi$header, "\\|", simplify = TRUE)[,1]
df_ncbi$assembly_id <- str_split(df_ncbi$header, "\\|", simplify = TRUE)[,1]


x<- str_split_fixed(df_ncbi$header, "\\|", 2)[,2] %>% 
  str_split_fixed(" ",2)
# IF these have been clustered, there is "size" info in the header. I.e. copynumber
x <- str_split_fixed(df_ncbi$header, "/", 2)[,2] %>% 
  str_split_fixed(";",2)
dim(x)
# df_ncbi$sub_sequence <- x[,1]


# copynumber <- as.data.frame(table(df_ncbi$assembly_id))
# If size info!
# df_ncbi$rDNA_CopyNumber <-  str_split_fixed(x[,2], "=", 2)[,2]
# df_ncbi$rDNA_CopyNumber <- as.numeric(df_ncbi$rDNA_CopyNumber )

#df_ncbi <- df_ncbi %>% select(-header)
df_ncbi <- df_ncbi %>% relocate(assembly_id,sub_sequence)
df_ncbi[493,]

# Use the metadata from GOAT, or Skip to Metadata from taxid only! 
# head(Goat)
# Goat$assembly_id
# assembly_id <- Goat$assembly_id %>% str_split_fixed("\\.",2)
# Goat$assembly_id <- assembly_id[,1]
# df_ncbi[50201,]$assembly_id
# df_ncbi <- left_join(df_ncbi, Goat, relationship = "many-to-many")
### ALTERNATIVE: 

# assembly_info_taxa <- read.csv("./../assembly_summary_genbank.taxonomy.txt")
assembly_info_taxa <- readxl::read_xlsx("./../assembly_summary_genbank.taxonomy_revised.xlsx")

#### SYNCHRONIZE METADATA ####
# full_join(mycocosm_genome_info,Goat)
colnames(assembly_info_taxa)

df_ncbi <- left_join(df_ncbi,assembly_info_taxa, relationship = "many-to-many")
df_ncbi[is.na(df_ncbi$taxid),]
missing_taxid <- df_ncbi[is.na(df_ncbi$taxid),]$assembly_id %>% unique()
missing_taxid

table(df_ncbi$kingdom)
table(df_ncbi$phylum)
table(df_ncbi$class)
df_ncbi$seq_length <- nchar(df_ncbi$sequence)
max(df_ncbi$seq_length)
min(df_ncbi$seq_length)
hist(df_ncbi$seq_length, breaks = 100)

#### STATS FOR GENOMES ####
# OBS Most of these scripts are from pre-cleaning of ROD. 
ROD_Genome_stats <- table(df_ncbi$assembly_id) %>% as.data.frame()
colnames(ROD_Genome_stats) <- c("assembly_id", "rDNA_copynumber") 
# genome_tax <- df_ncbi %>% select(!c(sub_sequence,sequence,seq_length)) %>% unique()
ROD_Genome_stats <- left_join(ROD_Genome_stats, assembly_info_taxa)

sum(is.na(ROD_Genome_stats$taxid))

write.csv(ROD_Genome_stats,"ROD_Genome_stats.csv")
saveRDS(ROD_Genome_stats, "ROD_Genome_stats.rds")


#### Make figures etc ####

table(ROD_Genome_stats$kingdom)
table(ROD_Genome_stats$kingdom) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity") +
  labs(x = "Tax", y = "Frequency") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("All genomes.pdf")

table(ROD_Genome_stats$phylum) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity") +
  labs(x = "Tax", y = "Frequency") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ROD_Genome_stats %>% filter(kingdom == "Fungi")
ROD_Genome_stats_Fungi <- ROD_Genome_stats %>% filter(kingdom == "Fungi") 
table(ROD_Genome_stats_Fungi$class) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity") +
  labs(x = "Order", y = "Frequency") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ROD_Genome_stats %>% filter(kingdom == "Viridiplantae")
ROD_Genome_stats_Viridiplantae <- ROD_Genome_stats %>% filter(kingdom == "Viridiplantae") 
table(ROD_Genome_stats_Viridiplantae$order) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity") +
  labs(x = "Viridiplantae - order", y = "Frequency") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ROD_Genome_stats %>% filter(kingdom == "Metazoa")
ROD_Genome_stats_Metazoa <- ROD_Genome_stats %>% filter(kingdom == "Metazoa") 
table(ROD_Genome_stats_Metazoa$class) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity") +
  labs(x = "Metazoa - class", y = "Frequency") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

table(ROD_Genome_stats$kingdom)
ROD_Genome_stats %>% filter(kingdom == "Stramenopiles")
ROD_Genome_stats_Strameopiles <- ROD_Genome_stats %>% filter(kingdom == "Stramenopiles") 
table(ROD_Genome_stats_Strameopiles$genus) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity") +
  labs(x = "Oomycota - genus", y = "Frequency") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("OOmycots.pdf")

ROD_Genome_stats %>% filter(is.na(kingdom))

hist(log10(ROD_Genome_stats$rDNA_copynumber))

ROD_Genome_stats %>% 
  ggplot(aes(x=rDNA_copynumber)) +
  geom_histogram(bins = 100) +
  # scale_x_log10(
  #   breaks = scales::trans_breaks("log10", function(x) 10^x),
  #   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_minimal() + annotation_logticks() 
ggsave("copynumber.pdf")


ROD_Genome_stats %>% 
  ggplot(aes(x=assembly_id, y=sort(rDNA_copynumber)))  +
  geom_point(na.rm = TRUE)+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) + axis.text.x=element_blank()



### Copy number calculated pr. Genome accession: ####
copynumer <- table(df_ncbi$assembly_id) %>% as.data.frame()
hist(copynumer$Freq)
selection <- copynumer %>% filter(Freq > 1000) 
hist(selection$Freq)

### The clustered data have =size info

# df_ncbi$rDNA_CopyNumber <- as.numeric(df_ncbi$rDNA_CopyNumber)
selection <- df_ncbi %>% filter(rDNA_CopyNumber < 60 ) %>% 
  filter(rDNA_CopyNumber > 2 )
hist(selection$rDNA_CopyNumber, breaks = 20)

selection <- df_ncbi %>% filter(rDNA_CopyNumber > 100 )
table(selection$phylum)

selection <- test %>% filter(rDNA_CopyNumber < 100 )
selection <- selection %>% filter(rDNA_CopyNumber > 2 )
table(selection$phylum)
hist(selection$rDNA_CopyNumber)


colnames(df_ncbi)
df_ncbi_NA <- df_ncbi[is.na(df_ncbi$taxon_id),]
table(df_ncbi_NA$assembly_id)
write_delim(df_ncbi_NA, "missing.taxid.csv", delim = "\t", quote = "none")
df_ncbi <- df_ncbi %>% filter(!is.na(taxon_id))
df_ncbi <- df_ncbi %>% filter(!is.na(species_taxid))
df_ncbi$lineage <- paste(df_ncbi$superkingdom, df_ncbi$kingdom, df_ncbi$phylum, df_ncbi$class, df_ncbi$order, 
                         df_ncbi$family, df_ncbi$genus, df_ncbi$species, 
                 sep = ";") %>% str_replace_all(" ","_")

# NCBI <- df_ncbi
# mycocosm <- df_mycocosm

# saveRDS(NCBI, "NCBI.rds")

## 

df_mycocosm$rDNA_CopyNumber <- as.numeric(df_mycocosm$rDNA_CopyNumber)
df_ncbi$rDNA_CopyNumber <- as.numeric(df_ncbi$rDNA_CopyNumber)

test <- full_join(df_mycocosm,df_ncbi)

test$rDNA_CopyNumber <- as.numeric(test$rDNA_CopyNumber)
hist(test$rDNA_CopyNumber) 
hist(test$seq_length, breaks = 100) 
median(test$seq_length)
table(test$phylum)


#### EXPORT FASTA ####
# Export fasta

#df <- test #%>% filter(seq_length > 8000)
df <- test #%>% filter(genus=="Neurospora")
# header <- paste(df$species,df$assembly_id,"|",df$sub_sequence," taxid=",df$taxon_id, ";")
df <- df_ncbi
header <- paste0(df$assembly_id,"|",df$sub_sequence,"|",df$lineage , " taxid=",df$taxon_id, ";")
header <- paste0(df$assembly_id,"|",df$sub_sequence,"|",df$lineage , " taxid=",df$species_taxid, ";")
Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
Xfasta[c(FALSE, TRUE)] <- df$sequence

writeLines(Xfasta, "all.fasta")
getwd()

