library(tidyverse)
ROD_v0.4 <- read_rds("./../04_github/ROD/ROD_v0.4.rds") %>% as_tibble()
 
### EXPORT FASTA FULL ROD
df <- ROD_v0.4 #%>% filter(genus=="Neurospora")
header <- paste0(df$assembly_id,"|",df$seqid,";size=",df$size ,"|",df$lineage , " taxid=",df$taxid, ";")
Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
Xfasta[c(FALSE, TRUE)] <- df$sequence
writeLines(Xfasta, "./../04_github/ROD/ROD_v0.4.fasta")
getwd()


#### PR SELECTIO
df <- ROD_v0.4
colnames(ROD_v0.4)
selections <- df$family %>% unique()


for (selected in selections){
  df_sub <- df %>% filter(family==selected) 
  # header <- paste0(df_sub$seqid,";","size=",df_sub$size,"|",df_sub$lineage , " taxid=",df_sub$taxid, ";")
  header <- paste0(df_sub$seqid,"|",df_sub$lineage,"|","size=",df_sub$size)
  # header <- paste0(df_sub$seqid,";","size=",df_sub$size)
  Xfasta <- character(nrow(df_sub) * 2)
  Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
  Xfasta[c(FALSE, TRUE)] <- df_sub$sequence
  writeLines(Xfasta, paste0("./../05_goldenROD/pr_family/ROD_v0.4_", selected ,".fasta"))
  # getwd()
  #print(header)
  print(paste0("./../05_goldenROD/pr_family/ROD_v0.4_", selected ,".fasta"))
}
