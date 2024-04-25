### STATS FOR ROD
library(tidyverse)
library(wesanderson)
library(ggplot2)
library(metacoder)
library(ape)
library(ggridges)

# Read the current version:
ROD_v0.4 <- read_rds("./../04_github/ROD/ROD_v0.4.rds") %>% as_tibble()
ROD_v0.4_genome_stats <- read_rds("./../04_github/ROD/ROD_v0.4_genome_stats.rds") %>% as_tibble()


#### General Stats ####
# Total number of rDNA copies
sum(ROD_v0.4_genome_stats$rDNA_copies)
sum(ROD_v0.4$size)


# Counting copies: 
# Pr.genome: 
# copy pr. genome, 1
ROD_v0.4_genome_stats %>% filter(rDNA_copies == 1) %>% 
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))
  
# copy pr. genome, 2-10
ROD_v0.4_genome_stats %>% filter(rDNA_copies >= 2) %>% filter(rDNA_copies <= 10) %>% 
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))

# copy pr. genome, 11-100
ROD_v0.4_genome_stats %>%  filter(rDNA_copies >= 11) %>% filter(rDNA_copies <= 100) %>%
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))

# copy pr. genome, 101-1000
ROD_v0.4_genome_stats %>%  filter(rDNA_copies >= 101) %>% filter(rDNA_copies <= 1000) %>%
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))

ROD_v0.4_genome_stats %>%  filter(rDNA_copies >= 1001) %>% #filter(rDNA_copies <= 1000) %>%
  summarize(total_rDNA_copies = sum(rDNA_copies), average_rDNA_copies = mean(rDNA_copies))

# Counting VARIANTS
# copy pr. genome, 1
ROD_v0.4_genome_stats %>% filter(rDNA_variants == 1) %>% 
  summarize(total_rDNA_copies = sum(rDNA_variants),average_rDNA_copies = mean(rDNA_variants))

# copy pr. genome, 2-10
ROD_v0.4_genome_stats %>% filter(rDNA_variants >= 2) %>% filter(rDNA_variants <= 10) %>% 
  summarize(total_rDNA_copies = sum(rDNA_variants),average_rDNA_copies = mean(rDNA_variants))

# copy pr. genome, 11-100
ROD_v0.4_genome_stats %>%  filter(rDNA_variants >= 11) %>% filter(rDNA_variants <= 100) #%>%
  summarize(total_rDNA_copies = sum(rDNA_variants),average_rDNA_copies = mean(rDNA_variants))

# copy pr. genome, 101-1000
ROD_v0.4_genome_stats %>%  filter(rDNA_variants >= 100) %>% filter(rDNA_variants <= 1000) %>%
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))

ROD_v0.4_genome_stats %>%  filter(rDNA_variants >= 1000) %>% #filter(rDNA_copies <= 1000) %>%
  summarize(total_rDNA_copies = sum(rDNA_copies), average_rDNA_copies = mean(rDNA_copies))


# Max copynumber in genome
(ROD_v0.4_genome_stats %>% filter(total_size == max(total_size)))$assembly_id
max(ROD_v0.4_genome_stats$length_max)
min(ROD_v0.4_genome_stats$length_min)
ROD_v0.4 %>% filter(length ==min(ROD_v0.4$length))
ROD_v0.4 %>% filter(length ==max(ROD_v0.4$length))
ROD_v0.4_genome_stats %>% filter(length_min == min(ROD_v0.4_genome_stats$length_min)) %>% .$mean_distance

ROD_v0.4_genome_stats %>% filter(ROD_v0.4_genome_stats$division == "Fungi") %>% 
  filter(rDNA_variants == 81)

# Fungi percentage
ROD_v0.4_genome_stats %>% filter(ROD_v0.4_genome_stats$division == "Fungi")
dim(ROD_v0.4_genome_stats %>% filter(ROD_v0.4_genome_stats$division == "Fungi"))[1]/dim(ROD_v0.4_genome_stats)[1]
dim(ROD_v0.4_genome_stats %>% filter(ROD_v0.4_genome_stats$division == "Metazoa"))[1]/dim(ROD_v0.4_genome_stats)[1]
dim(ROD_v0.4_genome_stats %>% filter(ROD_v0.4_genome_stats$division == "Viridiplantae"))[1]/dim(ROD_v0.4_genome_stats)[1]

ROD_v0.4_genome_stats %>% filter(size_max_second_diff ==0)
max(na.omit(ROD_v0.4_genome_stats$max_distance))
ROD_v0.4_genome_stats %>% filter( max_distance == max(na.omit(ROD_v0.4_genome_stats$max_distance))) %>% .$mean_distance
ROD_v0.4 %>% filter(assembly_id =="GCA_009804295")

ROD_v0.4_genome_stats %>% filter( mean_distance == max(na.omit(ROD_v0.4_genome_stats$mean_distance)))
ROD_v0.4_genome_stats %>% filter( median_distance == max(na.omit(ROD_v0.4_genome_stats$median_distance)))
ROD_v0.4_genome_stats %>% filter(subdivision == "Echinamoebida")  %>% .$max_distance

boxplot(ROD_v0.4$length)
plot(density(ROD_v0.4$length))
median(ROD_v0.4$length)
ROD_v0.4 %>% filter(length == min(length))
ROD_v0.4 %>% filter(length == max(length))

ROD_v0.4 %>% filter(length > 5000) %>% filter(length < 7000)

# 5000-7000
61739/69480

ROD_v0.4 %>% filter(length > 5000) %>% filter(length < 6500)
59781/69480

# Lenght pr. other groups
group_of_interest <- "subdivision"
df <- ROD_v0.4 %>% 
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

df %>% select(subdivision,rDNA_copies_mean, rDNA_variants_mean)

View(df[842,])
head(df$species)
length(df$species)
sum(df$genomes)

ROD_v0.4_genome_stats %>% filter(species =="Brassica juncea")

# Order data according to some variable
ordered_data <- df %>% 
#  filter(rDNA_variants ==1 ) %>% 
  arrange(desc(total_size))

ordered_data
#df$group_of_interest <- factor(df$group_of_interest, levels = df$group_of_interest)
df[[group_of_interest]] <- factor(df[[group_of_interest]], levels = unique(df[[group_of_interest]]))


# Make data long, for plotting
long_data <- ordered_data %>%
  pivot_longer(cols = c(size_max, size_second_max, remaining_size), names_to = "category", values_to = "value")

# Plotting (only makes sense if group_of_interest = "assembly_id")
# ggplot(long_data, aes(x = .data[[group_of_interest]], y = value, fill = category)) +
#   geom_bar(stat = "identity") +  # Stacked bar plot
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot of Sizes by Assembly ID",
#        x = "Assembly ID",
#        y = "Total Size",
#        fill = "Category") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjusting x labels for readability
# 
# long_data <- ordered_data %>%
#   pivot_longer(cols = c(size_max, size_min), names_to = "variable", values_to = "value")
# 
# # Plotting
# ggplot(long_data, aes(x = assembly_id, y = value, color = variable, group = variable)) +
#   geom_line() +  # Using geom_line() to connect points; you can change to geom_point() if you prefer dots
#   theme_minimal() +
#   labs(title = "Comparison of Maximum Size and Difference to Second Max",
#        x = "Assembly ID",
#        y = "Value",
#        color = "Variable") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjusting x labels for readability


