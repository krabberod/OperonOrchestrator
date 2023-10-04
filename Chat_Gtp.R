# Install and load the required packages
install.packages("DBI")
install.packages("RSQLite")
library(DBI)
library(RSQLite)

# Specify the name of your SQLite database
db_name <- "OperonOrchestrator.db"

# Create a connection to the database
con <- dbConnect(RSQLite::SQLite(), db_name)

# Create a table for sequence data
dbExecute(con, "CREATE TABLE sequences (id INTEGER PRIMARY KEY, header TEXT, sequence TEXT);")

# Create a table for metadata
dbExecute(con, "CREATE TABLE metadata (id INTEGER PRIMARY KEY, genome_name TEXT, 
          taxon TEXT, length INTEGER, superkingdom TEXT, phylum TEXT, class TEXT, 
          `order` TEXT, family TEXT, genus TEXT, species TEXT);")



# Load the Bioconductor Biostrings package for reading FASTA files
library(Biostrings)

# Read the FASTA file and insert data into the database
fasta_file <- "./../all.fasta"
sequences <- readDNAStringSet(fasta_file)

# Define a function to extract genome and sequence names from the FASTA header
extract_genome_and_sequence_names <- function(header) {
  parts <- unlist(strsplit(header, "\\|"))
  genome_name <- parts[1]
  sequence_name <- parts[2]
  return(list(genome_name = genome_name, sequence_name = sequence_name))
}

metadata_file <- "./../metadata.csv"
metadata <- read.table(metadata_file, header = TRUE, sep = "\t")





# Insert metadata into the database
for (i in 1:nrow(metadata)) {
  genome_name <- metadata$Genome[i]
  ncbi_taxid <- metadata$NCBI_Taxid[i]
  superkingdom <- metadata$superkingdom[i]
  phylum <- metadata$phylum[i]
  class <- metadata$class[i]
  order <- metadata$order[i]
  family <- metadata$family[i]
  genus <- metadata$genus[i]
  species <- metadata$species[i]
  
  dbExecute(con, "INSERT INTO metadata (genome_name, ncbi_taxid, superkingdom, phylum, 
            class, order, family, genus, species) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);", 
            genome_name, ncbi_taxid, superkingdom, phylum, class, order, family, genus, species)
}



# Read the FASTA file and insert data into the database
fasta_file <- "sequences.fasta"
sequences <- readDNAStringSet(fasta_file)


for (i in 1:length(sequences)) {
  # Match the genome name with the metadata and retrieve associated metadata
  metadata_query <- dbGetQuery(con, "SELECT * FROM metadata WHERE genome_name = ?;", genome_name)
  
  if (!is.null(metadata_query)) {
    # Insert metadata and sequence data into their respective tables
    genome_name <- metadata$Genome[i]
    ncbi_taxid <- metadata_query$taxon
    superkingdom <- metadata_query$superkingdom
    phylum <- metadata_query$phylum
    class <- metadata_query$class
    order <- metadata_query$order
    family <- metadata_query$family
    genus <- metadata_query$genus
    species <- metadata_query$species
    
    dbExecute(con, "INSERT INTO metadata (genome_name, ncbi_taxid, superkingdom, phylum, class, `order`, family, genus, species) 
              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);")
  }
}

# Close the database connection
dbDisconnect(con)


