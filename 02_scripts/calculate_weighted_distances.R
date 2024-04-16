# A function that takes a distance matrix as input
# and returns a matrix weighted by the number of copies in the "size" information
# it calculates a weighted average distance between each pair of entries, using the sizes as weights. 
# Author: Anders K. Krabberød 
# March 2024

calculate_weighted_distances <- function(distance_matrix) {
  # Extract sizes from row names
  sizes <- as.numeric(sub(".*;size=", "", rownames(distance_matrix)))
  
  # Initialize a matrix to hold the weighted distances
  weighted_distance_matrix <- matrix(nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
  
  # Calculate the weighted distances
  for (i in 1:nrow(distance_matrix)) {
    for (j in 1:ncol(distance_matrix)) {
      # Geometric mean:
       weighted_distance_matrix[i, j] <- distance_matrix[i, j] / sqrt(sizes[i] * sizes[j]) 
      # Arithmetic mean: 
      # weighted_distance_matrix[i, j] <- distance_matrix[i, j] / ((sizes[i] + sizes[j]) / 2)
      # p=0.5
      #weighted_distance_matrix[i, j] <- distance_matrix[i, j] / ((sizes[i]^p + sizes[j]^p) / 2)^(1/p)
      
      
    }
  }
  
  # Assign row and column names
  rownames(weighted_distance_matrix) <- rownames(distance_matrix)
  colnames(weighted_distance_matrix) <- colnames(distance_matrix)
  
  # Return the weighted distance matrix
  return(weighted_distance_matrix)
}

