library(wesanderson)
library(ggplot2)

# ROD - statistics
# ROD_v0.0_unclean <- readRDS("ROD_v0.0_unclean.rds")
# colnames(ROD_v0.0_unclean)
# ROD_v0.0_unclean <- ROD_v0.0_unclean %>% select(sub_sequence, rDNA_CopyNumber, assembly_id, taxid,seq_length, sequence)
updated_taxonomy <- readxl::read_xlsx("./../assembly_summary_genbank.taxonomy_revised_V2.xlsx")
ROD_v0.0_unclean <- left_join(df_ncbi, updated_taxonomy, relationship = "many-to-many")
ROD_v0.0_unclean
ROD_v0.3 <- readRDS("ROD_v0.3.rds")
ROD_v0.3$size <- as.numeric(ROD_v0.3$size)

ROD_v0.3$division <- factor(ROD_v0.3$division, levels = unique(ROD_v0.3$division))
ROD_v0.3$subdivision <- factor(ROD_v0.3$subdivision, levels = unique(ROD_v0.3$subdivision))
ROD_v0.3$class <- factor(ROD_v0.3$class, levels = unique(ROD_v0.3$class))


contaminants <- anti_join(ROD_v0.0_unclean %>% select(assembly_id,seqid,size,taxid,length),ROD_v0.3)
sum(contaminants$size)
length(unique(contaminants$assembly_id))
sum(ROD_v0.3$size)
length(unique(ROD_v0.3$assembly_id))

ROD_v0.3_Genome_stats <- ROD_v0.3 %>% select(-c(seqid, length, size, sequence)) %>% unique()

rDNAcopies <- ROD_v0.3 %>%                  
  group_by(assembly_id) %>%        
  summarise(copynumbers = sum(size))
ROD_v0.3_Genome_stats$size <- rDNAcopies$copynumbers
ROD_v0.3_Genome_stats %>% filter(size > 100)
max(ROD_v0.3_Genome_stats$size)
View(ROD_v0.3_Genome_stats %>% filter(size == 6534))
max(ROD_v0.3$length)
ROD_v0.3 %>% filter(length == 16463) %>% View()
ROD_v0.3_Genome_stats %>% filter(assembly_id =="GCA_932526495")


split.factor="subdivision"
summarized_data <- ROD_v0.3_Genome_stats %>% 
  group_by(.data[[split.factor]]) %>%
  summarise(sum.size = sum(size), .groups = 'drop')

ROD_v0.3_Genome_stats %>% select(supergroup,.data[[split.factor]]) 
df <- left_join(ROD_v0.3_Genome_stats,summarized_data) %>% select(supergroup, division,.data[[split.factor]],sum.size) %>% unique()

ggplot(df, aes(x = .data[[split.factor]], y = sum.size)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  facet_wrap(supergroup ~ .) +
  labs(x = split.factor, y = "Sum of Size", 
       title = paste("Bar Chart of", split.factor, "vs. Sum of Size"))

df %>% filter(division == "Metazoa") %>% 
ggplot(aes(x = division, y = sum.size, fill = .data[[split.factor]])) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set3", name = "Supergroup") +  # Customizing the fill scale and legend title
  theme_minimal() +
  labs( y = "Sum of Size", 
       title = "Sum of Size by Subdivision and Supergroup") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") 


df %>% filter(division == "Metazoa") %>%
ggplot( aes(x = subdivision, y = sum.size, fill = subdivision)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Supergroup", y = "Sum of Size", fill = "Subdivision", 
       title = "Sum of Size by Supergroup with Subdivision Fill") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df %>% filter(division == "Metazoa") %>%
  ggplot( aes(x = reorder(class, -sum.size), y = sum.size)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Class", y = "Sum of operons", fill = "Subdivision", 
       title = "Operons in Metazoa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position = "none")

df %>% filter(division == "Viridiplantae") %>%
  ggplot( aes(x = reorder(class, -sum.size), y = sum.size)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Class", y = "Sum of operons", fill = "Subdivision", 
       title = "Operons in Viridiplantae") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position = "none")


### COUNTING PER GENOME
split.factor="assembly_id"
summarized_data <- ROD_v0.3_Genome_stats %>% 
  group_by(.data[[split.factor]]) %>%
  summarise(sum.size = sum(size), .groups = 'drop')

df <- left_join(ROD_v0.3_Genome_stats,summarized_data) %>% select(supergroup, division,.data[[split.factor]],sum.size) %>% unique()

df %>% filter(division == "Fungi") %>%
  ggplot(aes(x = reorder(assembly_id, -sum.size), y = sum.size)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Class", y = "Sum of operons", fill = "Subdivision", 
       title = "Operons in fungi") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position = "none")

# VIOLIN PLOT

df$supergroup <- factor(df$supergroup, levels = unique(df$supergroup))
df$division <- factor(df$division, levels = unique(df$division))

ggplot(df, aes(x = division, y = log(sum.size), fill = division)) + 
  geom_violin(trim = FALSE) + 
  facet_wrap(~ supergroup, scales = "free_x", nrow = 1) + # Group by supergroup in separate panels
  labs(title = "Violin plots of sum.size by division, grouped by supergroup",
       x = "Division",
       y = "Sum Size") +
  #scale_fill_brewer(palette = "Set1") + # Use a color palette for different divisions
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Adjusting x labels for readability
        strip.background = element_blank(), # Optional: remove background of facet labels
        strip.text.x = element_text(face = "bold")) # Optional: bold facet labels



### BOXPLOT WITH n=operons.pr.genome ####
df$supergroup <- factor(df$supergroup, levels = unique(df$supergroup))
df$division <- factor(df$division, levels = unique(df$division))

# Create a summary data frame with counts for each division
division_counts <- ROD_v0.3_Genome_stats %>%
  group_by(division) %>%
  summarise(count = n()) %>%
  ungroup()

# Merge the counts back with the original data frame
df_with_counts <- merge(df, division_counts, by = "division")

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


#### Opisthokonta
#### FUNGI ####

fungi <- ROD_v0.3 %>% filter(division == "Fungi")
fungi_stats <- ROD_v0.3_Genome_stats %>% filter(division == "Fungi")
sum(fungi$size)

fungi %>% select(supergroup,.data[[split.factor]]) 
df <- left_join(fungi,summarized_data) %>% select(division,subdivision,class,.data[[split.factor]],sum.size) %>% unique()
df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each class
class_counts <- ROD_v0.3_Genome_stats %>%
  group_by(class) %>%
  summarise(count = n()) %>%
  ungroup()

# Merge the counts back with the original data frame
df_with_counts <- merge(df, class_counts, by = "class")

# Create a new label combining division name and count
df_with_counts$class_label <- paste(df_with_counts$class, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(x = reorder(class_label, -count), y = sum.size, fill = class)) + 
  geom_boxplot(fill="goldenrod") +  
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
metazoa <- ROD_v0.3 %>% filter(division == "Metazoa")
metazoa_stats <- ROD_v0.3_Genome_stats %>% filter(division == "Metazoa")
sum(metazoa$size)

metazoa %>% select(supergroup,.data[[split.factor]]) 
metazoa$seqid %>% unique() %>% length()

df <- left_join(metazoa,summarized_data) %>% select(division,subdivision,class,.data[[split.factor]],sum.size) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each class
class_counts <- ROD_v0.3_Genome_stats %>%
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
ggsave("01_ROD_results/ROD_operons_pr.genomes_metazoa.pdf")


#### VIRIDIPLANTAE Division ####
viridiplantae <- ROD_v0.3 %>% filter(division == "Viridiplantae")
viridiplantae_stats <- ROD_v0.3_Genome_stats %>% filter(division == "Viridiplantae")
sum(viridiplantae$size)

viridiplantae %>% select(supergroup,.data[[split.factor]]) 
viridiplantae$seqid %>% unique() %>% length()

df <- left_join(viridiplantae,summarized_data) %>% select(division,subdivision,class,.data[[split.factor]],sum.size) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each class
class_counts <- ROD_v0.3_Genome_stats %>%
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

magnoliopsida <- ROD_v0.3 %>% filter(class == "Magnoliopsida")
magnoliopsida_stats <- ROD_v0.3_Genome_stats %>% filter(class == "Magnoliopsida")
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
ggplot(df_with_counts, aes(y = reorder(order_label, count), x = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_magnoliopsida.pdf")


#### ASCOMYCOTA #####

ascomycota <- ROD_v0.3 %>% filter(subdivision == "Ascomycota")
ascomycota_stats <- ROD_v0.3_Genome_stats %>% filter(class == "ascomycota")
sum(ascomycota$size)

ascomycota %>% select(supergroup,.data[[split.factor]]) 
ascomycota$seqid %>% unique() %>% length()

df <- left_join(ascomycota,summarized_data) %>% select(division,subdivision,order,class,.data[[split.factor]],sum.size) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$order <- factor(df$order, levels = unique(df$order))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each order
order_counts <- ascomycota %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(x = reorder(order_label, -count), y = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_ascomycota.pdf")




#### BASIDIOMYCOTA #####

basidiomycota <- ROD_v0.3 %>% filter(subdivision == "Basidiomycota")
basidiomycota_stats <- ROD_v0.3_Genome_stats %>% filter(class == "Basidiomycota")
sum(basidiomycota$size)

basidiomycota %>% select(supergroup,.data[[split.factor]]) 
basidiomycota$seqid %>% unique() %>% length()

df <- left_join(basidiomycota,summarized_data) %>% select(division,subdivision,order,class,.data[[split.factor]],sum.size) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$order <- factor(df$order, levels = unique(df$order))
df$class <- factor(df$class, levels = unique(df$class))

# Create a summary data frame with counts for each order
order_counts <- basidiomycota %>%
  group_by(order) %>%
  summarise(count = n(),.groups = 'drop') %>% 
  arrange(desc(count))

# Merge the counts back with the original data frame
df_with_counts <- merge(df, order_counts, by = "order")

# Create a new label combining division name and count
df_with_counts$order_label <- paste(df_with_counts$order, " (n=", df_with_counts$count, ")", sep="")

# my_colors <- wes_palette("GrandBudapest1", n = length(unique(df$division)), type = "continuous")
ggplot(df_with_counts, aes(x = reorder(order_label, -count), y = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_basidiomycota.pdf")

#### ARHTHROPODA #####

arthropoda <- ROD_v0.3 %>% filter(subdivision == "Arthropoda")
arthropoda_stats <- ROD_v0.3_Genome_stats %>% filter(class == "Arthropoda")
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
ggplot(df_with_counts, aes(x = reorder(order_label, -count), y = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
  labs(x = "", y = "Operon copies pr. genome") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), # Adjusting x labels for readability
        strip.text.x = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("01_ROD_results/ROD_operons_pr.genomes_arthropoda.pdf")

#### CHORDATA #####

chordata <- ROD_v0.3 %>% filter(subdivision == "Chordata")
chordata_stats <- ROD_v0.3_Genome_stats %>% filter(class == "Chordata")
sum(arthropoda$size)

chordata %>% select(supergroup,.data[[split.factor]]) 
chordata$seqid %>% unique() %>% length()

df <- left_join(chordata,summarized_data) %>% select(assembly_id,size,subdivision,class,order,) %>% unique()

df$division <- factor(df$division, levels = unique(df$division))
df$subdivision <- factor(df$subdivision, levels = unique(df$subdivision))
df$order <- factor(df$order, levels = unique(df$order))
df$class <- factor(df$class, levels = unique(df$class))

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
ggplot(df_with_counts, aes(x = reorder(order_label, -count), y = sum.size, fill = order)) + 
  geom_boxplot(fill="goldenrod") +  
  facet_wrap(~ class, scales = "free", nrow = 3) + # Group by supergroup in separate panels
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
