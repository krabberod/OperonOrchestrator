
library(DBI)
library(RSQLite)

db_file <- "OperonOrchestrator.db"
con <- dbConnect(RSQLite::SQLite(), dbname = db_file)
dbWriteTable(con, "Mycocosm Sequence data", df, overwrite = TRUE) # Replace "table_name" with your desired table name
dbWriteTable(con, "Myocosm metadata", mycocosm_genome_info, overwrite = TRUE)
dbDisconnect(con)
