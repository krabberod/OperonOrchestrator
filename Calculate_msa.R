library("tidyverse")
library("Biostrings")
library("msa")

ROD <- readRDS("ROD_with_seq_v0.2.rds")

# ROD %>% filter(rDNA_copy_number > 1) %>% select(assembly_id)

# test <- ROD %>% filter(assembly_id == "GCA_000001215") %>% select(seqid,sequence)
## Test
# msa_test <- msa(test$sequence, type = "dna")
# pairwiseAlignment(test$sequence[1],test$sequence[2])
## Not promising. 
# I rather want to take sequences out of ROD. pr. assembly id, write fastas, and manipulate on saga: 

### EXPORT FASTA
# genomes <- unique(ROD$assembly_id)
# genomes <- genomes[1:3]
# More than 1 copy: 
genomes <- table(ROD$assembly_id) %>% as.data.frame()
colnames(genomes) <- c("assembly_id", "Freq")
genomes$assembly_id<- as.character(genomes$assembly_id)

genomes <- genomes %>% filter(Freq > 1) %>% select(assembly_id)
genomes <- unique(genomes)
genomes$assembly_id[1]
df <-NULL
for (item in genomes$assembly_id) {
  print(item)
  df <- ROD %>% filter(assembly_id == item)
  header <- paste0(df$seqid, ";size=",df$size)
  Xfasta <- character(nrow(df) * 2)
  Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
  Xfasta[c(FALSE, TRUE)] <- df$sequence
  writeLines(Xfasta, paste0("pr_genome_fasta/",item, ".fasta"))
}
