### TABLE TO NEWICK

library(data.tree)
library(tidyverse)

df <-paste0(ROD_v0.4_genome_stats$lineage,";",ROD_v0.4_genome_stats$assembly_id)

df <- df %>% as_data_frame()

addNode <- function(pathString, tree) {
  taxa <- unlist(strsplit(pathString, ";"))
  node <- tree
  for (taxon in taxa) {
    if (taxon != "") {
      if (is.null(node[[taxon]])) {
        node <- node$AddChild(name = taxon)
      } else {
        node <- node[[taxon]]
      }
    }
  }
  return(tree)  # Ensure the modified tree is returned
}

root <- Node$new("Root")
apply(df, 1, function(row) {
  # Create a path string, concatenating taxonomic levels separated by semicolons
  pathString <- paste(na.omit(row), collapse = ";")
  addNode(pathString, root)
})


for (i in 1:nrow(df)) {
  row <- df[i, ]
  pathString <- paste(na.omit(row), collapse = ";")
  root <- addNode(pathString, root)  # Ensure to update the root with returned tree
}


print(root)


treeToNewick <- function(node) {
  if (node$isLeaf) {
    return(node$name)
  } else {
    childrenNewick <- sapply(node$children, function(child) treeToNewick(child))
    return(paste0("(", paste(childrenNewick, collapse = ","), ")", node$name))
  }
}

newickString <- treeToNewick(root)
newickString <- paste0(newickString, ";")  # Newick format ends with a semicolon
cat(newickString)
write_lines(newickString,"/Users/anderkkr/Downloads/ROD.nwk")


