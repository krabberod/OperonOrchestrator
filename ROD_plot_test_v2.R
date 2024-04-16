library(ggplot2)
library(dplyr)
library(patchwork)



summarized_data <- ROD_v0.4 %>% #filter(division == "Metazoa") %>% 
  group_by(subdivision, class, order) %>%
  summarise(total_size = length(unique(assembly_id)), .groups = 'drop') 

ggplot(summarized_data, aes(x = class, y = total_size, fill = order)) +
  geom_bar(stat = "identity") +
  facet_wrap(subdivision ~ ., scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Class", y = "Total Size", fill = "Order") +
  ggtitle("Stacked Barplot of Sizes by Class and Order") +
  theme(legend.position = "none")


ggplot(summarized_data, aes(x = class, y = total_size, fill = order)) +
  geom_col(position = position_dodge(width = 0.7)) +  # Use geom_col with position_dodge
  facet_wrap(subdivision ~ ., scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Class", y = "Total Size", fill = "Order") +
  ggtitle("Stacked Barplot of Sizes by Class and Order") +
  theme(legend.position = "none")


#complete_data <- df %>%
complete_data <- ROD_v0.4 %>% filter(division == "Metazoa") %>% 
  #complete(subdivision, class, order, fill = list(size = 0)) %>%
  group_by(subdivision, class, order) %>%
  summarise(total_size = sum(size), .groups = 'drop')

# Plot with consistent bar widths across facets
ggplot(complete_data, aes(x = class, y = total_size, fill = order)) +
  geom_col(width = 0.7) +  # Set a fixed width for bars
  facet_wrap(subdivision ~ ., scales = "free_x") +  # Free scales only for x
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Class", y = "Total Size", fill = "Order") +
  ggtitle("Stacked Barplot of Sizes by Class and Order") +
  theme(legend.position = "none")

 #####

df <-  ROD_v0.4 %>% filter(division == "Metazoa")
plot_list <- list()
# Loop through each unique subdivision
for(subdivision in unique(df$subdivision)) {
  # Filter data for the current subdivision
  sub_data <- df %>% filter(subdivision == !!subdivision)
  
  # Summarize data
  summarized_data <- sub_data %>%
    group_by(class, order) %>%
    summarise(total_size = sum(size), .groups = 'drop')
  
  # Create plot
  p <- ggplot(summarized_data, aes(x = class, y = total_size, fill = order)) +
    geom_bar(stat = "identity") +
    labs(x = "Class", y = "Total Size", fill = "Order") +
    ggtitle(paste(subdivision)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  # Add plot to the list
  plot_list[[subdivision]] <- p
}

# Remove NULL elements if any
plot_list <- plot_list[!sapply(plot_list, is.null)]

combined_plot <- wrap_plots(plot_list)

# Print the combined plot
combined_plot

#### 



library(ggplot2)
library(dplyr)

# Assuming 'df' is your data frame
# Summarize data first
summarized_data <- ROD_v0.4 %>% filter(division == "Metazoa") %>%
  group_by(assembly_id, subdivision, class, order) %>%
  summarise(total_size = length(unique(assembly_id)), .groups = 'drop')

# Create a new factor for x-axis combining subdivision and class
summarized_data <- summarized_data %>%
  mutate(x_group = factor(paste(subdivision, class, sep = " - "), 
                          levels = unique(paste(subdivision, class, sep = " - "))))

# Plot
ggplot(summarized_data, aes(x = x_group, y = total_size, fill = order)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Subdivision - Class", y = "Total Size", fill = "Order") +
  ggtitle("Stacked Barplot of Sizes by Subdivision and Class") +
#  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")


#### 
# Assuming 'df' is your data frame
# Summarize data first
# summarized_data <-  ROD_v0.3 %>% filter(division == "Metazoa") %>%
#   group_by(subdivision, class, order) %>%
#   summarise(total_size = sum(size), .groups = 'drop')
summarized_data <- ROD_v0.4 %>% filter(division == "Metazoa") %>%
  group_by(subdivision, class, order) %>%
  summarise(total_size = sum(length(unique(assembly_id))), .groups = 'drop')

# Create a new factor for x-axis combining subdivision and class
summarized_data <- summarized_data %>%
  mutate(x_group = paste(subdivision, class, sep = " - ")) %>%
  mutate(x_group = as.factor(x_group))

# Calculate total size for each x_group to order them
ordering_data <- summarized_data %>%
  group_by(x_group) %>%
  summarise(total_size = sum(total_size)) %>%
  arrange(desc(total_size)) %>%
  mutate(x_group = fct_inorder(x_group))

# Apply the ordering to the original data
summarized_data$x_group <- factor(summarized_data$x_group, levels = ordering_data$x_group)

# Plot
ggplot(summarized_data %>% filter(total_size>3), aes(x = x_group, y = total_size, fill = order)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Subdivision - Class", y = "Total Size", fill = "Order") +
  ggtitle("Stacked Barplot of Genomes by Subdivision and Class, Ordered by Size") +
#  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")

#####
summarized_data <- ROD_v0.4 %>% filter(supergroup == "Opisthokonta") %>%
  group_by(division, subdivision, class) %>%
  summarise(total_size = sum(length(unique(assembly_id))), .groups = 'drop')

# Create a new factor for x-axis combining division and subdivision
summarized_data <- summarized_data %>%
  mutate(x_group = paste(division, subdivision, sep = " - ")) %>%
  mutate(x_group = as.factor(x_group))

# Calculate total size for each x_group to class them
classing_data <- summarized_data %>%
  group_by(x_group) %>%
  summarise(total_size = sum(total_size)) %>%
  arrange(desc(total_size)) %>%
  mutate(x_group = fct_inorder(x_group))

# Apply the classing to the original data
summarized_data$x_group <- factor(summarized_data$x_group, levels = classing_data$x_group)

# Plot
ggplot(summarized_data %>% filter(total_size > 10), aes(x = x_group, y = total_size, fill = class)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Subdivision - Class", y = "Total Size", fill = "Order") +
  ggtitle("Stacked Barplot of Genomes by Subdivision and Class, Ordered by Size") +
  #  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "right")


  