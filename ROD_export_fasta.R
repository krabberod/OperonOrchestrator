# Export Fasta 

### EXTRACT FASTA

df <- readRDS("./../04_github/ROD/ROD_v0.4.rds")
header <- paste0(df$seqid,";","size=",df$size,"|",df$lineage , " taxid=",df$taxid, ";")
# header <- paste0(df$seqid,";","size=",df$size, " taxid=",df$taxid, ";")
Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
Xfasta[c(FALSE, TRUE)] <- df$sequence
writeLines(Xfasta, "./../04_github/ROD/ROD_v0.4.fasta")
getwd()


#### 
df <- ROD_v0.4
genomes <- df$assembly_id %>% unique()


for (genome in genomes){
  df_sub <- df %>% filter(assembly_id==genome) 
  # header <- paste0(df_subf$seqid,";","size=",df_sub$size,"|",df_sub$lineage , " taxid=",df_sub$taxid, ";")
  header <- paste0(df_sub$seqid,";","size=",df_sub$size)
  Xfasta <- character(nrow(df_sub) * 2)
  Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
  Xfasta[c(FALSE, TRUE)] <- df_sub$sequence
  writeLines(Xfasta, paste0("./../05_goldenROD/pr_genome/ROD_v0.4_", genome ,".fasta"))
  getwd()
  #print(header)
  print(paste0("./../05_goldenROD/pr_genome/ROD_v0.4_", genome ,".fasta"))
}
  
