data <- summarized_data
data_long <- data %>%
  pivot_longer(cols = c("class", "order"), names_to = "Category", values_to = "Value") %>%
  filter(!is.na(Value))  # Remove NA values that may result from pivoting

subdivisions <- unique(data$subdivision)

plot_list <- list()
for (subdivision in subdivisions) {
  # Filter the data for the current subdivision
  data_sub <- data_long %>%
    filter(subdivision == !!subdivision)
  
  # Create the plot
  p <- ggplot(data_sub, aes(x = Category, y = total_size, fill = Value)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = paste("Subdivision:", subdivision), x = "Class/Order", y = "Total Size") 
    #scale_fill_brewer(palette = "Paired")  # Use a color palette that distinguishes the stacks well
  
  # Save the plot in the list
  plot_list[[subdivision]] <- p
}

print(plot_list[[3]])
###
unique_subdivisions <- unique(data$subdivision)
plots <- list()

for (i in seq_along(unique_subdivisions)) {
  subdivision_data <- data %>% 
    filter(subdivision == unique_subdivisions[i]) %>%
    group_by(class, order) %>%
    summarise(total_size = sum(total_size, na.rm = TRUE))  # Summing up in case of multiple entries for a class-order combination
  
  # Generate the plot
  p <- ggplot(subdivision_data, aes(x = class, y = total_size, fill = order)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    labs(title = paste("Subdivision: ", unique_subdivisions[i]), x = "Class", y = "Total Size") 
    #scale_fill_brewer(palette = "Paired")  # Choose a palette that offers clear distinction between "orders"
  
  # Store the plot in the list
  plots[[unique_subdivisions[i]]] <- p
}
print(plots[[4]])


#####

for (subdivision in subdivisions) {
  # Filter the data for the current subdivision
  data_sub <- data %>%
    filter(subdivision == !!subdivision)

p <- ggplot(data_sub, aes(x = order, y = total_size, fill = order)) +  # x-axis uses 'order' to ensure each stack is distinguishable
  geom_bar(stat = "identity", position = "stack") +
  #facet_wrap(~class, scales = "free", nrow = 2) +  # Facets by 'class', adjust 'nrow' as needed
  theme_minimal() +
  labs(x = "Order", y = "Total Size") +
  #scale_fill_brewer(palette = "Paired") +  # Use a clear and distinct color palette
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angle for readability

plot_list[[subdivision]] <- p
}

plot_list[[1]]

#### 
divisions <- unique(ROD_v0.4_genome_stats$division)
divisions <- "Viridiplantae"

for (division in divisions) {
  # Filter the data for the current subdivision
  data_sub <- ROD_v0.4_genome_stats %>%
    filter(division == !!division)
  
  p <- ggplot(data_sub, aes(x = class, y = sum.size, fill = order)) +  # x-axis uses 'order' to ensure each stack is distinguishable
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    #facet_wrap(~class, scales = "free", nrow = 2) +  # Facets by 'class', adjust 'nrow' as needed
    theme_minimal() +
    labs(x = "Order", y = "Total Size") +
    #scale_fill_brewer(palette = "Paired") +  # Use a clear and distinct color palette
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angle for readability
  p
  plot_list[[subdivision]] <- p
}

plot_list[[4]]

