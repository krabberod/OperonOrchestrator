library(wesanderson)
library(ggplot2)
library(tidyverse)
library("ggvenn") 

ROD_v0.3 <- readRDS("ROD_v0.3.rds")
seqid <- read.table('./../05_goldenROD/golden_rod_18S_AND_28S.seqid', header =FALSE)



colnames(seqid) <- "seqid"
colnames(ROD_v0.3)
ROD_v0.4 <- left_join(seqid, ROD_v0.3)

saveRDS(ROD_v0.4, './../04_github/ROD/ROD_v0.4.rds')
ROD_v0.4 <- readRDS("./../04_github/ROD/ROD_v0.4.rds") %>% as_tibble()
length(ROD_v0.4$assembly_id %>% unique())

cutoff <- 4300
shorty <- ROD_v0.4 %>% filter(length <= cutoff) %>% .$assembly_id %>% unique()
longy <- ROD_v0.4 %>% filter(length > cutoff) %>% .$assembly_id %>% unique()

B <- list("shorty"=shorty, "longy"=longy)
# create venn diagram and display all the sets 
ggvenn(B)
intersect(shorty,longy)

ROD_v0.4 %>% filter(assembly_id == "GCA_932526495") %>% View()
ROD_v0.4 %>% filter(species =="Alces alces") %>% select(assembly_id) %>% unique()


length(ROD_v0.3$supergroup %>% unique())
sum(ROD_v0.4$size)

max(ROD_v0.4$size)
length(ROD_v0.4$species %>% unique())
species <- table(ROD_v0.4$species) %>% as.data.frame()
arrange(species,desc(Freq))

# Calculate pr. genome copynuber
rDNA_copynumber <- ROD_v0.4 %>% 
  group_by(assembly_id) %>%
  summarise(rDNA_copynumber = sum(size), .groups = 'drop')
max(rDNA_copynumber$rDNA_copynumber)

rDNA_copynumber %>% filter(rDNA_copynumber > 100)

min(ROD_v0.4$length)

ggplot(ROD_v0.4, aes(x=length, y=size)) +
  geom_point(size=2, shape=23)


### FASTA #### 
df <- ROD_v0.4 #%>% filter(genus=="Neurospora")
header <- paste0(df$assembly_id,"|",df$seqid,";size=",df$size ,"|",df$lineage , " taxid=",df$taxid, ";")
Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
Xfasta[c(FALSE, TRUE)] <- df$sequence
writeLines(Xfasta, "./../04_github/ROD/ROD_v0.4.fasta")
getwd()

ROD_Genome_stats <- table(ROD_v0.4$assembly_id)  %>% as.data.frame()
colnames(ROD_Genome_stats) <- c("assembly_id","rep99")
left_join(ROD_Genome_stats, ROD_v0.4) %>% as_tibble()




