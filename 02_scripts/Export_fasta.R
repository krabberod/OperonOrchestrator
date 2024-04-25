# Load the tidyverse package
library(tidyverse)

# Read the RDS file and convert it to a tibble
ROD <- read_rds("ROD_v1.0_operon_variants.rds") %>% as_tibble()

# Prepare the data for FASTA export
df <- ROD

# Create a header for the FASTA file using the seqid, lineage, and size columns. Values are separated by "|".
header <- paste0(df$seqid,"|",df$lineage,"|","size=",df$size)
# Create a character vector for the FASTA file
Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
Xfasta[c(FALSE, TRUE)] <- df$sequence

# Write the FASTA file
writeLines(Xfasta, "ROD_v1.0_reference_sequences.fasta")

# If you want to split the FASTA file by family, you can use the following code
# Select unique families from the data
selections <- df$family %>% unique()

# For each unique family, create a FASTA file
for (selected in selections){
  df_sub <- df %>% filter(family==selected) 
  header <- paste0(df_sub$seqid,"|",df_sub$lineage,"|","size=",df_sub$size)
  Xfasta <- character(nrow(df_sub) * 2)
  Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
  Xfasta[c(FALSE, TRUE)] <- df_sub$sequence

  # Write the FASTA file for the selected family
  writeLines(Xfasta, paste0("./../05_goldenROD/pr_family/ROD_v0.4_", selected ,".fasta"))

  # Print the path of the written FASTA file
  print(paste0("./../05_goldenROD/pr_family/ROD_v0.4_", selected ,".fasta"))
}