library(tidyverse)
library(wesanderson)
library(ggplot2)

ROD_v1.2 <- read_rds("./../04_github/ROD/ROD_v1.2.rds") %>% as_tibble()
ROD_v1.2_genome_stats <- read_rds("./../04_github/ROD/ROD_v1.2_genome_stats.rds") %>% as_tibble()

ROD_v1.2_genome_stats$median_distance[is.na(ROD_v1.2_genome_stats$median_distance)] <- 0
ROD_v1.2_genome_stats$mean_distance[is.na(ROD_v1.2_genome_stats$mean_distance)] <- 0

df <- ROD_v1.2_genome_stats

df_data <- ROD_v1.2_genome_stats %>%
  group_by(supergroup,division,subdivision,class,order,family,genus,species) %>%
  summarise(genomes = n(),  
            max.rDNA_copies = max(rDNA_copies), 
            max.rDNA_variants = max(rDNA_variants),
            avg.distance = mean(mean_distance),
            min.length = min(length_min),
            max.length = max(length_max), mean.length = mean(length_mean)
            ) %>%
  ungroup()
write.table(df_data, "01_ROD_results/Stat_pr_species_ROD_v1.2.tab", quote = F, sep = "\t", row.names = F)



df_data <- ROD_v1.2_genome_stats %>%
  group_by(supergroup,division,subdivision) %>%
  summarise(genomes = n(),  
            max.rDNA_copies = max(rDNA_copies), 
            max.rDNA_variants = max(rDNA_variants),
            avg.distance = mean(mean_distance),
            min.length = min(length_min),
            max.length = max(length_max), mean.length = mean(length_mean)
  ) %>%
  ungroup()
write.table(df_data, "01_ROD_results/Stat_pr_subdivision_ROD_v1.2.tab", quote = F, sep = "\t", row.names = F)

# Lenght pr. other groups
group_of_interest <- "subdivision"
df <- ROD_v1.2%>% 
  #arrange(.data[[group_of_interest]], desc(size)) %>%
  group_by(.data[[group_of_interest]]) %>%
  reframe(rDNA_variants = n(),
          genomes = length(unique(assembly_id)),
          total_size = sum(size),
          rDNA_copies_mean = sum(size) /length(unique(assembly_id)),
          rDNA_variants_mean = n() /length(unique(assembly_id)),
          size_max = max(size),
          size_min = min(size), 
          size_second_max = ifelse(n_distinct(size) > 1, sort(size, decreasing = TRUE)[2], NA), 
          size_max_second_diff = size_max - ifelse(n_distinct(size) > 1, sort(size, decreasing = TRUE)[2], NA),
          size_max_prop =  size_max / total_size,
          size_second_prop =  size_second_max / total_size, 
          remaining_size = total_size - size_max - size_second_max, 
          all_sizes = paste(size, collapse = ";"),
          all_lengths = paste(length, collapse = ";"),
          supergroup = paste(unique(supergroup),collapse = ";"),
          division = paste(unique(division),collapse = ";"),
          subdivision = paste(unique(subdivision),collapse = ";"),
          class = paste(unique(class),collapse = ";"),
          order = paste(unique(order),collapse = ";"),
          family = paste(unique(family),collapse = ";"),
          genus = paste(unique(genus),collapse = ";"),
          species = paste(unique(species),collapse = ";")
  ) %>% ungroup() 

# df %>% select(subdivision,rDNA_copies_mean, rDNA_variants_mean)
write.table(df, "01_ROD_results/Stat_pr_subdivision_ROD_v1.2.tab", quote = F, sep = "\t", row.names = F)


# Create a summary data frame with counts for each division
division_counts <- ROD_v1.2_genome_stats %>%
  group_by(division) %>%
  summarise(count = n()) %>%
  ungroup()

# Merge the counts back with the original data frame
df_with_counts <- merge(ROD_v1.2_genome_stats, division_counts, by = "division") %>% as_tibble()

# Create a new label combining division name and count
df_with_counts$division_label <- paste(df_with_counts$division, "\n(n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(x = reorder(division_label, -sum.size), y = sum.size, fill = division)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ supergroup, scales = "free", nrow = 2) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  #scale_fill_manual(values="black") + # Use a color palette for different divisions
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Adjusting x labels for readability
        #strip.background = element_blank(), # Optional: remove background of facet labels
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")

ggsave("01_ROD_results/ROD_operons_pr.genomes_overview.pdf")

####

split.factor="assembly_id"
summarized_data <- ROD_v1.2%>% 
  group_by(.data[[split.factor]]) %>%
  summarise(sum.size = sum(size), .groups = 'drop')


#### Opisthokonta
#### FUNGI ####

fungi <- ROD_v1.2 %>% filter(division == "Fungi")
fungi_stats <- ROD_v1.2_genome_stats %>% filter(division == "Fungi")
sum(fungi$size)

fungi %>% select(supergroup, assembly_id)
df <- left_join(fungi,summarized_data) %>% select(division,subdivision,class,sum.size) %>% unique()

# Create a summary data frame with counts for each class
class_counts <- ROD_v1.2 %>% filter(division == "Fungi")  %>% 
  group_by(class) %>%
  summarise(count = n()) %>%
  ungroup() #%>% View()

# Merge the counts back with the original data frame
df_with_counts <- merge(df, class_counts, by = "class")

# Create a new label combining division name and count
df_with_counts$class_label <- paste(df_with_counts$class, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(x = reorder(class_label, -count), y = sum.size, fill = class)) + 
  geom_boxplot(fill="#F1BB7B", notch = FALSE) + 
  # geom_jitter(width = 0.1, alpha = 0.1) +  
  facet_wrap(~ subdivision, scales = "free", nrow = 2) + # Group by supergroup in separate panels
  labs(x = "", y = "") +
  #scale_fill_manual(values="black") + # Use a color palette for different divisions
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), # Adjusting x labels for readability
        #strip.background = element_blank(), # Optional: remove background of facet labels
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_fungus.pdf")

#### METAZOA ####
metazoa <- ROD_v1.2 %>% filter(division == "Metazoa")
metazoa_stats <- ROD_v1.2_genome_stats %>% filter(division == "Metazoa")
sum(metazoa$size)

metazoa %>% select(supergroup,.data[[split.factor]]) 
metazoa$seqid %>% unique() %>% length()

df <- left_join(metazoa,summarized_data) %>% select(division,subdivision,class,.data[[split.factor]],sum.size) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each class
class_counts <- ROD_v1.2_genome_stats %>%
  group_by(class) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, class_counts, by = "class")

# Create a new label combining division name and count
df_with_counts$class_label <- paste(df_with_counts$class, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(x = reorder(class_label, -count), y = sum.size, fill = class)) + 
  geom_boxplot(fill="#F1BB7B") +  
  facet_wrap(~ subdivision, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_metazoa.pdf")


#### VIRIDIPLANTAE Division ####
viridiplantae <- ROD_v1.2 %>% filter(division == "Viridiplantae")
viridiplantae_stats <- ROD_v1.2_genome_stats %>% filter(division == "Viridiplantae")
sum(viridiplantae$size)

viridiplantae %>% select(supergroup,.data[[split.factor]]) 
viridiplantae$seqid %>% unique() %>% length()

df <- left_join(viridiplantae,summarized_data) %>% select(division,subdivision,class,.data[[split.factor]],sum.size) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each class
class_counts <- ROD_v1.2_genome_stats %>%
  group_by(class) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, class_counts, by = "class")

# Create a new label combining division name and count
df_with_counts$class_label <- paste(df_with_counts$class, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(x = reorder(class_label, -count), y = sum.size, fill = class)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ subdivision, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_viridiplantae.pdf")

#### Magnoliopsida SUBDIVISION ####

magnoliopsida <- ROD_v1.2 %>% filter(class == "Magnoliopsida")
magnoliopsida_stats <- ROD_v1.2_genome_stats %>% filter(class == "Magnoliopsida")
sum(magnoliopsida$size)

magnoliopsida %>% select(supergroup,.data[[split.factor]]) 
magnoliopsida$seqid %>% unique() %>% length()

df <- left_join(magnoliopsida,summarized_data) %>% select(division,subdivision,order,class,.data[[split.factor]],sum.size) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$order <- factor(df$order, levels = unique(df$order))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each order
order_counts <- magnoliopsida_stats %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts %>% filter(count > 1), aes(y = reorder(order_label, count), x = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_magnoliopsida.pdf")


#### Ascomycota - Order  ####

ascomycota <- ROD_v1.2 %>% filter(subdivision == "Ascomycota")
ascomycota_stats <- ROD_v1.2_genome_stats %>% filter(subdivision == "Ascomycota")
sum(ascomycota$size)

ascomycota$seqid %>% unique() %>% length()

df <- left_join(ascomycota,summarized_data)  %>% select(division,subdivision,order,class,.data[[split.factor]],sum.size) %>% unique()

# Create a summary data frame with counts for each order
order_counts <- ascomycota_stats %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(y = reorder(order_label, count), x = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(y = "", x = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_ascomycota.pdf")


#### Basidiomycota - Order  ####

basidiomycota <- ROD_v1.2 %>% filter(subdivision == "Basidiomycota")
basidiomycota_stats <- ROD_v1.2_genome_stats %>% filter(subdivision == "Basidiomycota")
sum(basidiomycota$size)

basidiomycota$seqid %>% unique() %>% length()

df <- left_join(basidiomycota,summarized_data)  %>% select(division,subdivision,order,class,.data[[split.factor]],sum.size) %>% unique()

# Create a summary data frame with counts for each order
order_counts <- basidiomycota_stats %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts %>% filter(count > 2), aes(y = reorder(order_label, count), x = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  #xfacet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(y = "", x = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_basidiomycota.pdf")


#### ARHTHROPODA #####

arthropoda <- ROD_v1.2 %>% filter(subdivision == "Arthropoda")
arthropoda_stats <- ROD_v1.2_genome_stats %>% filter(class == "Arthropoda")
sum(arthropoda$size)

arthropoda %>% select(supergroup,.data[[split.factor]]) 
arthropoda$seqid %>% unique() %>% length()

df <- left_join(arthropoda,summarized_data) %>% select(division,subdivision,order,class,.data[[split.factor]],sum.size) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$order <- factor(df$order, levels = unique(df$order))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each order
order_counts <- arthropoda %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts , aes(y = reorder(order_label, count), x = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  # facet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_arthropoda.pdf")

#### CHORDATA #####

chordata <- ROD_v1.2 %>% filter(subdivision == "Chordata")
chordata_stats <- ROD_v1.2_genome_stats %>% filter(subdivision == "Chordata")
sum(arthropoda$size)

# chordata %>% select(supergroup,.data[[split.factor]]) 
chordata$seqid %>% unique() %>% length()

df <- left_join(chordata,summarized_data) 

# Create a summary data frame with counts for each order
order_counts <- chordata %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(y = reorder(order_label, count), x = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  # facet_wrap(~ class, scales = "free", nrow = 2) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=4), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_chordata.pdf")



# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(x = reorder(subdivision, -count), y = sum.size, fill = class)) + 
  geom_bar(stat = "identity", position = "stack") +
  #  facet_wrap(~ class, scales = "free", nrow = 3) #+ # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=4), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) #+
#theme(legend.position = "none")

#### INSECTA #####

insecta <- ROD_v1.2 %>% filter(class == "Insecta")
insecta_stats <- ROD_v1.2_genome_stats %>% filter(class == "Insecta")
sum(insecta$size)

# chordata %>% select(supergroup,.data[[split.factor]]) 
insecta$seqid %>% unique() %>% length()

df <- left_join(insecta,summarized_data) 

# Create a summary data frame with counts for each order
order_counts <- insecta %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(y = reorder(order_label, count), x = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 2) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_insecta.pdf")


#### INSECTA #####

insecta <- ROD_v1.2 %>% filter(class == "Insecta")
insecta_stats <- ROD_v1.2_genome_stats %>% filter(class == "Insecta")
sum(insecta$size)
insecta$seqid %>% unique() %>% length()

df <- left_join(insecta,summarized_data) 

# Create a summary data frame with counts for each order
order_counts <- insecta %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(y = reorder(order_label, count), x = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 2) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_insecta.pdf")

###
summarized_data <- ROD_v1.2 %>% filter(division == "Metazoa") %>% 
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
