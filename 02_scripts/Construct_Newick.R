# Script: Construct_Newick.R
# Author: Anders K. Krabberød
# Date: March 2024  
# Description: This script takes a data frame with lineage and assembly_id information, 
# and constructs a Newick tree structure. The script uses the data.tree and tidyverse 
# libraries for creating and manipulating tree-like hierarchical data structures and 
# for data manipulation, respectively.

library(data.tree)
library(tidyverse)

df <-paste0(ROD_v0.4_genome_stats$lineage,";",ROD_v0.4_genome_stats$assembly_id)
colnames(ROD_v0.4_genome_stats)


stripped_list <- sapply(df, function(x) {
  parts <- strsplit(x, ";")[[1]]  # Split each entry by ';'
  paste(parts[1:4], collapse = ";")  # Recombine the first four elements
})

df <- stripped_list %>% as_data_frame() %>% unique()
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

# Create a new root node for the tree
root <- Node$new("Root")

# Apply the addNode function to each row of df
apply(df, 1, function(row) {
  # Create a path string by concatenating the non-missing values in the row, separated by semicolons
  pathString <- paste(na.omit(row), collapse = ";")
  # Add the path to the tree
  addNode(pathString, root)
})

# Loop over each row of df
for (i in 1:nrow(df)) {
  # Get the current row
  row <- df[i, ]
  # Create a path string by concatenating the non-missing values in the row, separated by semicolons
  pathString <- paste(na.omit(row), collapse = ";")
  # Add the path to the tree and update the root
  root <- addNode(pathString, root)
}

# Print the tree
print(root)

# Define a function to convert a tree to Newick format
treeToNewick <- function(node) {
  if (node$isLeaf) {
    # If the node is a leaf, return its name
    return(node$name)
  } else {
    # If the node is not a leaf, recursively convert its children to Newick format
    childrenNewick <- sapply(node$children, function(child) treeToNewick(child))
    # Return the Newick format of the node and its children
    return(paste0("(", paste(childrenNewick, collapse = ","), ")", node$name))
  }
}

# Convert the tree to Newick format
newickString <- treeToNewick(root)
# Add a semicolon at the end, as required by the Newick format
newickString <- paste0(newickString, ";")
# Print the Newick string
cat(newickString)
# Write the Newick string to a file
write_lines(newickString,"/Users/anderkkr/Downloads/ROD_class.nwk")