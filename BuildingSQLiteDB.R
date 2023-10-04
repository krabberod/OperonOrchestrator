library(RSQLite)

con <- dbConnect(RSQLite::SQLite(), db_name)



dbListTables(con)
dbListObjects(con)
dbReadTable(con,"metadata")
dbReadTable(con,"sequences")

extract_genome_and_sequence_names <- function(header) {
  parts <- unlist(strsplit(header, "\\|"))
  genome_name <- parts[1]
  sequence_name <- parts[2]
  return(list(genome_name = genome_name, sequence_name = sequence_name))
}

sql_insert_sequence <- "INSERT INTO sequences (header, sequence) VALUES (?, ?);"

for (i in 1:length(sequences)) {
  header <- names(sequences)[i]
  sequence <- as.character(sequences[i])
  
  # Extract genome name and sequence name from the FASTA header
  header_info <- extract_genome_and_sequence_names(header)
  genome_name <- header_info$genome_name
  
  # Execute the SQL statement with placeholders
  dbExecute(con, sql_insert_sequence, list(header, sequence))
}

# Add the genome_name and sequence_name columns to the sequences table
dbExecute(con, "ALTER TABLE sequences ADD COLUMN genome_name TEXT;")
dbExecute(con, "ALTER TABLE sequences ADD COLUMN sequence_name TEXT;")

# Commit the changes to the database

dbExecute(con, "UPDATE sequences SET genome_name = SUBSTR(header, 1, INSTR(header, '|') - 1), 
          sequence_name = SUBSTR(header, INSTR(header, '|') + 1);")

# Query for unique genome names
query <- "SELECT DISTINCT genome_name FROM sequences;"
unique_genome_names <- dbGetQuery(con, query)

dbExecute(con, "ALTER TABLE metadata RENAME COLUMN taxon TO ncbi_taxid")



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



