# Function: calculate_weighted_distances
# Author: Anders K. Krabberød
# Date: March 2024
# Description: This function takes a distance matrix and returns a matrix of 
# weighted distances. The weights are determined by the "size" information. 
# It calculates a weighted average distance between each pair of entries, 
# using the sizes as weights. The method of calculation can be specified.
# Usage: calculate_weighted_distances(distance_matrix, method = "geometric", p = 0.5)
# Arguments:
#   distance_matrix: A matrix of distances between entries. The row and column
#   names should contain the "size" information.
#   method: The method of calculating the weighted distances. Can be "geometric",
#   "arithmetic", or "power_mean". Default is "geometric".
#   p: The power for the "power_mean" method. Default is 0.5.
# Returns: A matrix of weighted distances between entries.
# Example:
#   distance_matrix <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3)
#   rownames(distance_matrix) <- c("A;size=1", "B;size=2", "C;size=3")
#   colnames(distance_matrix) <- c("A;size=1", "B;size=2", "C;size=3")
#   weighted_distance_matrix <- calculate_weighted_distances(distance_matrix, method = "power_mean", p = 1)
#   print(weighted_distance_matrix)
#   # Output:
#   #           A         B         C
#   # A 0.0000000 0.5000000 0.6666667

calculate_weighted_distances <- function(distance_matrix, method = "geometric", p = 0.5) {
  # Extract sizes from row names
  sizes <- as.numeric(sub(".*;size=", "", rownames(distance_matrix)))
  
  # Initialize a matrix to hold the weighted distances
  weighted_distance_matrix <- matrix(nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
  
  # Calculate the weighted distances
  for (i in 1:nrow(distance_matrix)) {
    for (j in 1:ncol(distance_matrix)) {
      if (method == "geometric") {
        # Geometric mean:
        weighted_distance_matrix[i, j] <- distance_matrix[i, j] / sqrt(sizes[i] * sizes[j]) 
      } else if (method == "arithmetic") {
        # Arithmetic mean: 
        weighted_distance_matrix[i, j] <- distance_matrix[i, j] / ((sizes[i] + sizes[j]) / 2)
      } else if (method == "power_mean") {
        # Power mean with p:
        weighted_distance_matrix[i, j] <- distance_matrix[i, j] / ((sizes[i]^p + sizes[j]^p) / 2)^(1/p)
      }
    }
  }
  
  # Set the row and column names of the weighted distance matrix
  rownames(weighted_distance_matrix) <- rownames(distance_matrix)
  colnames(weighted_distance_matrix) <- colnames(distance_matrix)
  
  return(weighted_distance_matrix)
}

