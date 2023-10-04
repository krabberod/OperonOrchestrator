library(tidyverse)
library(Biostrings)
library(ape)
library(stringdist)
#library(DECIPHER)
library(ips)



# install.packages("remotes")
# remotes::install_github("mikemc/tidyseq")

# Importing Fasta:
fasta_file <- "./../myco.rep.4000.fasta"
sequences <- readDNAStringSet(fasta_file)

df <- NULL
df$header <- as.character(names(sequences))
# colnames(df) <- "header"
df$sequence <- unname(as.character(sequences))

# colnames(df) <- c("header", "sequence")
str(df)
df <- as_tibble(df)

parts <- str_split(df$header, "\\|", simplify = TRUE)[,1]

df$Genome <- str_split(df$header, "\\|", simplify = TRUE)[,1]
x <- str_split_fixed(df$header, "\\|", 2)[,2] %>% 
  str_split_fixed(" ",2)
# IF these have been clustered, there is "size" info in the header. I.e. copynumber
x <- str_split_fixed(df$header, "\\|", 2)[,2] %>% 
  str_split_fixed(";",2)
df$sub_sequence <- x[,1]
df <- df %>% select(-header)
df <- df %>% relocate(Genome,sub_sequence)
df[493,]

# read metadata
mycocosm_genome_info <- read.table("./../metadata.csv", header = TRUE, sep = "\t")
df <- left_join(df, mycocosm_genome_info)
table(df$phylum)
table(df$class)
df$seq_length <- nchar(df$sequence)
max(df$seq_length)
min(df$seq_length)


copynumber <- as.data.frame(table(df$Genome))
# If size info!
colnames(copynumber) <- c("Genome", "rDNA_CopyNumber")
 str_split_fixed(x[,2], "=", 2)[,2]
copynumber$rDNA_CopyNumber <-  str_split_fixed(x[,2], "=", 2)[,2]

mycocosm_genome_info <- left_join(mycocosm_genome_info, copynumber)
max(na.omit(as.integer(mycocosm_genome_info$rDNA_CopyNumber)))
selection <- df %>% filter(Genome =="Abobie1") %>% select(sub_sequence, sequence)
colnames(df)
df[is.na(df$NCBI_Taxid), ]
df <- df %>% filter(!is.na(NCBI_Taxid))
lineage <- paste(df$superkingdom, df$phylum, df$class, df$order, df$family, df$genus, df$species, 
   sep = ";") %>% str_replace_all(" ","_")




# Export fasta
header <- paste(df$species,df$Genome,"|",df$sub_sequence," taxid=",df$NCBI_Taxid, ";")
header <- paste0(df$Genome,"|",df$sub_sequence,"|",lineage, " taxid=",df$NCBI_Taxid, ";")

Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
Xfasta[c(FALSE, TRUE)] <- df$sequence


writeLines(Xfasta, "mycocosm_ribo_testV1.fasta")
