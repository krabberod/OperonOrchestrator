library(tidyverse)
library(Biostrings)
library(ape)
library(stringdist)
#library(DECIPHER)
library(ips)


#### IMPORT FROM MYCOCOSM ####
# Importing Fasta:
# The fasta file has "size" info that equals the number of copies for that sequence in the genome (at 99% clustering)

fasta_file <- "./../01_results/01_mycocosm_results/mycocosm_all.fasta"
sequences <- readDNAStringSet(fasta_file)

df_mycocosm <- NULL
df_mycocosm$header <- as.character(names(sequences))
df_mycocosm$sequence <- unname(as.character(sequences))
str(df_mycocosm)
df_mycocosm <- as_tibble(df_mycocosm)
parts <- str_split(df_mycocosm$header, "\\|", simplify = TRUE)[,1]

df_mycocosm$assembly_id <- str_split(df_mycocosm$header, "\\|", simplify = TRUE)[,1]
x <- str_split_fixed(df_mycocosm$header, "\\|", 2)[,2] %>% 
  str_split_fixed(" ",2)
# IF these have been clustered, there is "size" info in the header. I.e. copynumber
x <- str_split_fixed(df_mycocosm$header, "\\|", 2)[,2] %>% 
  str_split_fixed(";",2)
df_mycocosm$sub_sequence <- x[,1]
# df_mycocosm <- df_mycocosm %>% select(-header)
df_mycocosm$rDNA_CopyNumber <-  str_split_fixed(x[,2], "=", 2)[,2]
df_mycocosm$rDNA_CopyNumber <- as.numeric(df_mycocosm$rDNA_CopyNumber)

# Check if any rDNA copynumbers are missing
df_mycocosm[is.na(df_mycocosm$rDNA_CopyNumber), ]
df_mycocosm[443,]

df_mycocosm$data_origin <- "Mycocosm"

df_mycocosm$seq_length <- nchar(df_mycocosm$sequence)
max(df_mycocosm$seq_length)
min(df_mycocosm$seq_length)

df_mycocosm <- df_mycocosm %>% relocate(header, assembly_id, sub_sequence, rDNA_CopyNumber, data_origin)
df_mycocosm

#### IMPORT METADATA FROM MYCOCOSM ####
# For Goat see next section.
mycocosm_genome_info <- read.table("metadata.csv", header = TRUE, sep = "\t") %>% as_tibble()
head(mycocosm_genome_info)

mycocosm_genome_info <- mycocosm_genome_info %>% rename("Genome" = "assembly_id") %>%
  rename("NCBI_Taxid" = "taxon_id")


df_mycocosm <- left_join(df_mycocosm, mycocosm_genome_info)
table(df_mycocosm$phylum)
table(df_mycocosm$class)
max(df_mycocosm$rDNA_CopyNumber)
min(df_mycocosm$rDNA_CopyNumber)

colnames(df_mycocosm)
df_mycocosm[is.na(df_mycocosm$taxon_id), ]
df_mycocosm <- df_mycocosm %>% filter(!is.na(taxon_id))

df_mycocosm$lineage <- paste(df_mycocosm$superkingdom, df_mycocosm$phylum, df_mycocosm$class, 
                             df_mycocosm$order, df_mycocosm$family, 
                             df_mycocosm$genus, df_mycocosm$species, sep = ";") %>% 
  str_replace_all(" ","_")
# OBS TODO: REDO ALL TAXIDS BASED ON TAXONOMIZR

hist(as.numeric(df_mycocosm$rDNA_CopyNumber), breaks = 100)
hist(as.numeric(df_mycocosm$seq_length), breaks = 100)

# saveRDS(df_mycocosm, "mycocosm.rds")


#### IMPORT GOAT ####
Goat <- readxl::read_xlsx("Goat.xlsx")
head(Goat)
colnames(Goat)
table(Goat$phylum)





#### IMPORT FROM NCBI ####
# fasta_file <- "ncbi.single.copy.fasta"
# fasta_file <- "./../plants.4000.single.copy.fasta"
# fasta_file <- "./../test.fasta"
fasta_file <- "./../01_results/02_ncbi_results/ribres_4000_precluster.fasta"
# fasta_file <- "ribres_4000_cluster.fasta"

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
df_ncbi$sub_sequence <- x[,1]

# copynumber <- as.data.frame(table(df_ncbi$assembly_id))
# If size info!
# df_ncbi$rDNA_CopyNumber <-  str_split_fixed(x[,2], "=", 2)[,2]
# df_ncbi$rDNA_CopyNumber <- as.numeric(df_ncbi$rDNA_CopyNumber )

df_ncbi <- df_ncbi %>% select(-header)
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
assembly_info_taxa <- read.csv("./../assembly_summary_genbank.taxonomy.txt")


#### SYNCHRONIZE METADATA ####
# full_join(mycocosm_genome_info,Goat)
colnames(assembly_info_taxa)



df_ncbi <- left_join(df_ncbi,assembly_info_taxa, relationship = "many-to-many")
df_ncbi[is.na(df_ncbi$species_taxid),]
table(df_ncbi$kingdom)
table(df_ncbi$phylum)
table(df_ncbi$class)
df_ncbi$seq_length <- nchar(df_ncbi$sequence)
max(df_ncbi$seq_length)
min(df_ncbi$seq_length)
hist(df_ncbi$seq_length, breaks = 100)

#### STATS FOR GENOMES ####
ROD_Genome_stats <- table(df_ncbi$assembly_id) %>% as.data.frame()
colnames(ROD_Genome_stats) <- c("assembly_id", "rDNA_copynumber") 
# genome_tax <- df_ncbi %>% select(!c(sub_sequence,sequence,seq_length)) %>% unique()
ROD_Genome_stats <- left_join(ROD_Genome_stats, assembly_info_taxa)
write.csv(ROD_Genome_stats,"ROD_Genome_stats.csv")

table(ROD_Genome_stats$kingdom)

ROD_Genome_stats %>% filter(kingdom == "Fungi")
ROD_Genome_stats_Fungi <- ROD_Genome_stats %>% filter(kingdom == "Fungi") 
table(ROD_Genome_stats_Fungi$order) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity") +
  labs(x = "Order", y = "Frequency") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

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

ROD_Genome_stats %>% filter(is.na(kingdom))


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

