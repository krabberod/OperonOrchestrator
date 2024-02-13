library(tidyverse)
library(Biostrings)
# Fix missing taxnomy in Goat: 
library(taxonomizr)
library(tidyverse)
library(readxl)

accessionTaxa.sql <- "./../../../99_div/accessionTaxa.sql"

### Look for contaminants and wrongly annotated ribosomes
# All copies have been blasted agains NCBI_nr, PR2 and Silva. 
# Simplest test is to match the taxid for the blast hits 1) intragenomical, 2) against the taxonomic assignment

# Goat <- readxl::read_xlsx("Goat.xlsx")
# assembly_info_taxa <- readxl::read_xlsx("./../assembly_summary_genbank.taxonomy_revised.xlsx")
ROD_Genome_stats <- readRDS("ROD_Genome_stats.rds") %>% as_tibble()

# Read blasthits with taxonomy and taxid:  
# blastResults_tax <- read_rds("./../01_results/02_ncbi_results/ribres_4000_cluster.blastn.tsv.besthit.tax.rds") %>% as_tibble()
blastResults <- read.table("./../01_results/02_ncbi_results/ribres_4000_cluster.blastn.tsv.besthit", header = F) %>% as_tibble()

colnames(blastResults) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", 
                             "sstart", "send", "evalue", "bitscore", "slen", "staxid")
 
desTax = c("superkingdom","kingdom","phylum", "subphylum", "class", "order", "family", "genus", "species")
taxa <- getTaxonomy(blastResults$staxid %>% str_split(";", simplify = T) %>% .[,1], accessionTaxa.sql, desiredTaxa = desTax) %>% 
   as_tibble()
# # Following the nomenclautre of blast: add s- for subject to the result columns
colnames(taxa) <- c("ssuperkingdom","skingdom","sphylum", "ssubphylum", "sclass", "sorder", "sfamily", "sgenus", "sspecies")
blastResults_tax <- cbind(blastResults, taxa) %>% as_tibble()
blastResults_tax$assembly_id <- blastResults_tax$qseqid %>% str_split("\\|", simplify = TRUE) %>% .[,1]

# Match the taxid from ROD to the t blast results to verify the correspondence between balst and original assignment. 

blastResults_tax <- left_join(blastResults_tax, ROD_Genome_stats)
blastResults_tax <- blastResults_tax %>% relocate("assembly_id","qseqid","taxid","staxid","sseqid")

# write.table(blastResults_tax, "./../01_results/02_ncbi_results/ribres_4000_cluster.blastn.tsv.besthit.tax", 
#             col.names = TRUE, sep = "\t", row.names = FALSE)


saveRDS(blastResults_tax, "./../01_results/02_ncbi_results/ribres_4000_cluster.blastn.tsv.besthit.tax.rds")

# check if any lack taxid:             
blastResults_tax[is.na(blastResults_tax$staxid),]
blastResults_tax[is.na(blastResults_tax$taxid),]

subset(blastResults_tax, blastResults_tax$taxid == blastResults_tax$staxid) 
subset(blastResults_tax, blastResults_tax$genus == blastResults_tax$sgenus)
subset(blastResults_tax, blastResults_tax$family == blastResults_tax$sfamily)


# Some have no "sperkingdom"
dim(blastResults_tax[is.na(blastResults_tax$superkingdom),])

# Remove the strangest hits: 
blastResults_tax <- blastResults_tax[!is.na(blastResults_tax$ssuperkingdom),]

taxids.counts <- blastResults_tax %>%                  
  # filter(!is.na(aa)) %>%    
  group_by(assembly_id) %>%         
  summarise(Unique_Elements = n_distinct(staxid))

# Find those with only one taxis (might still have multiple copies)
# And they might still be wrong (i.e. a contaminant)
taxids.unique <- taxids.counts %>% filter(Unique_Elements == 1)

# Look for identical taxid in Genome annotation and blast from the unique taxids:
# Subset the blast.results
blastResults_tax.uniq.taxids <- blastResults_tax %>% filter(assembly_id %in% taxids.unique$assembly_id)
# Check if taxid (i.e. genome annot) is identitcal to blast search (i.e. staxis)
subset(blastResults_tax.uniq.taxids, blastResults_tax.uniq.taxids$taxid == blastResults_tax.uniq.taxids$staxid)

# 15530 has identical taxids 

# More than one taxid
taxids.multiple <- taxids.counts %>% filter(Unique_Elements > 1)
blastResults_tax.mutiple <- blastResults_tax %>% filter(assembly_id %in% taxids.multiple$assembly_id)

# Some of those that have more than one taxids can still be the same genus: 
blastResults_tax %>% filter(assembly_id == "GCA_000001765")

# Or same family
blastResults_tax %>% filter(assembly_id == "GCA_000003195")

# Or widely different 
blastResults_tax %>% filter(assembly_id == "GCA_000004655")


na.genus <- blastResults_tax %>%                  
  filter(is.na(genus))

uniq.genus <- blastResults_tax %>%                  
  filter(!is.na(genus)) %>%    
  group_by(assembly_id) %>%        
  summarise(Unique_Elements = n_distinct(genus))
uniq.genus %>% filter(Unique_Elements > 1)
  

# Find those with same genus
same_taxid <- subset(blastResults_tax, blastResults_tax$taxid == blastResults_tax$staxid)
sum(is.na(same_taxid$assembly_id))
dim(same_taxid)
length(unique(same_taxid$assembly_id))

same_genus <- subset(blastResults_tax, blastResults_tax$genus == blastResults_tax$sgenus)
same_genus <- blastResults_tax[blastResults_tax$genus == blastResults_tax$sgenus, ]
sum(is.na(same_genus$assembly_id))
same_genus <- same_genus[!is.na(same_genus$assembly_id),]
dim(same_genus)
length(unique(same_genus$assembly_id))

different_genus <- blastResults_tax[is.na(blastResults_tax$genus) | is.na(blastResults_tax$sgenus) | 
                   blastResults_tax$genus != blastResults_tax$sgenus, ]

same_family <- subset(blastResults_tax, blastResults_tax$family == blastResults_tax$sfamily)
dim(same_family)
length(unique(same_family$assembly_id))


# different_family <- subset(blastResults_tax, blastResults_tax$family != blastResults_tax$sfamily)
different_family <- blastResults_tax[is.na(blastResults_tax$family) | is.na(blastResults_tax$sfamily) | 
                                       blastResults_tax$family != blastResults_tax$sfamily, ]
# same_class <- subset(blastResults_tax, blastResults_tax$class == blastResults_tax$sclass)

different_family %>%                  
  filter(!is.na(staxid)) %>%    
  group_by(assembly_id) %>%        
  summarise(Unique_Elements = n_distinct(staxid))


same_ribo <- blastResults_tax[blastResults_tax$taxid == blastResults_tax$staxid | 
                                blastResults_tax$species == blastResults_tax$sspecies |
                                blastResults_tax$genus == blastResults_tax$sgenus |
                                blastResults_tax$family == blastResults_tax$sfamily,]
same_ribo <- same_ribo[!is.na(same_ribo$assembly_id),]

true_ribo <- same_ribo %>%                  
  filter(!is.na(staxid)) %>%    
  group_by(assembly_id) %>%        
  summarise(Unique_Elements = n_distinct(staxid))

ROD <- same_ribo

# Extract the contaminant sequences (included here are genomes that also have the "true" sequence)
blastResults_tax.contam_seq <- blastResults_tax %>% filter(!(qseqid %in% same_ribo$qseqid))
# anti_join(blastResults_tax,same_ribo)
blastResults_tax.contam_seq[is.na(blastResults_tax.contam_seq$qseqid),]
blastResults_tax.contam_seq %>%                  
  group_by(assembly_id) %>%        
  summarise(Unique_Elements = n_distinct(staxid))

# 4192 Genomes have contaminants

# The Genomes
blastResults_tax.contam_ass <- blastResults_tax %>% filter(!(assembly_id %in% same_ribo$assembly_id))
blastResults_tax.contam_ass %>%                  
  # filter(!is.na(staxid)) %>%    
  group_by(assembly_id) %>%        
  summarise(Unique_Elements = n_distinct(staxid))


api <- blastResults_tax.contam_ass %>% filter(phylum =="Apicomplexa" & sphylum =="Apicomplexa")
SAR <- blastResults_tax.contam_ass %>% filter(supergroup == "SAR") %>% filter(class == sclass)
SAR_2 <- blastResults_tax.contam_ass %>% filter(supergroup == "SAR") %>% filter(phylum == sphylum)
metazoa <- blastResults_tax.contam_ass %>% filter(skingdom == "Metazoa") %>% filter(class == sclass)
metazoa_2 <-blastResults_tax.contam_ass %>% filter(skingdom == "Metazoa") %>% filter(order == sorder)
fungi <- blastResults_tax.contam_ass %>% filter(skingdom == "Fungi") %>% filter(class == sclass)
plants <- blastResults_tax.contam_ass %>% filter(skingdom == "Viridiplantae") %>% filter(order == sorder)
plants_2 <- blastResults_tax.contam_ass %>% filter(skingdom == "Viridiplantae") %>% filter(class == sclass)
rhodophyta <- blastResults_tax.contam_ass %>% filter(phylum == "Rhodophyta" & sphylum == "Rhodophyta")
chlorophyta <- blastResults_tax.contam_ass %>% filter(phylum == "Chlorophyta" & sphylum == "Chlorophyta")
haptophyta <- blastResults_tax.contam_ass %>% filter(phylum == "Haptophyta" & sphylum == "Haptophyta")
evosea <- blastResults_tax.contam_ass %>% filter(phylum == "Evosea" & sphylum == "Evosea")
Prasinodermales <- blastResults_tax.contam_ass %>% filter(order == "Prasinococcales" & sorder == "Prasinodermales")
keep <- blastResults_tax.contam_ass %>% filter(assembly_id %in% c("GCA_001179505", "GCA_004369235", "GCA_021439945", " GCA_900128395"))
Taphrinomycotina <- blastResults_tax.contam_ass %>% filter(assembly_id %in% c("GCA_000312925", "GCA_000836195", "GCA_001929475", "GCA_003717165", 
                                                                              "GCA_005281575", "GCA_005281585", "GCA_005281805", "GCA_008802775"))



ROD_1 <- rbind(ROD, api, SAR, SAR_2, metazoa, metazoa_2, fungi, plants,
               plants_2, rhodophyta, chlorophyta, haptophyta, evosea, 
               Prasinodermales, keep, Taphrinomycotina)
ROD_1 <- unique(ROD_1)
to_be_sorted <- anti_join(blastResults_tax.contam_ass, ROD_1)
contam <- to_be_sorted %>% filter(kingdom == "Fungi" & skingdom!="Fungi" | 
                                    skingdom == "Fungi" & kingdom!="Fungi" |
                                    kingdom == "Metazoa" & skingdom!="Metazoa"| 
                                    skingdom == "Metazoa" & kingdom!="Metazoa"|
                                    kingdom == "Viridiplantae" & skingdom!="Viridiplantae" |
                                    skingdom == "Viridiplantae" & kingdom!="Viridiplantae" |
                                    supergroup == "SAR" & skingdom =="Viridiplantae" | 
                                    supergroup == "SAR" & skingdom =="Fungi" | 
                                    supergroup == "SAR" & skingdom =="Metazoa" |
                                    sphylum =="Nematoda" & phylum != "Nematode" |
                                    phylum == "Evosea" & sphylum != "Evosea" |
                                    supergroup == "Eukaryota_XXXX"
                                    )
contam_1 <- to_be_sorted %>% filter(kingdom =="Metazoa") %>% filter(sphylum != phylum)
contam_2 <- to_be_sorted %>% filter(kingdom =="Fungi") %>% filter(sphylum != phylum)
contam_3 <- unique(rbind(contam_1, contam, contam_2))
contam_3 <- contam_3[!is.na(contam_3$assembly_id),]



to_be_sorted_1 <- anti_join(to_be_sorted, contam_3)
# For those alreay in ROD, it is not  with another copy
to_be_sorted_1 <- to_be_sorted_1 %>% filter(!(assembly_id %in% ROD_1$assembly_id))
write.table(to_be_sorted_1, "to_be_sorted_1.csv", row.names = FALSE, quote = FALSE, sep = "\t")


ROD_1 %>% filter(assembly_id =="GCA_018136815")
sorted_taxid <- to_be_sorted_1 %>% select(assembly_id,taxid,staxid)

taxid <- getTaxonomy(sorted_taxid$taxid, accessionTaxa.sql, desiredTaxa = desTax)
staxid <- getTaxonomy(sorted_taxid$staxid, accessionTaxa.sql, desiredTaxa = desTax)
x <- cbind(sorted_taxid, staxid, taxid)
colnames(x) <- c("assembly_id","taxid","staxid","ssuperkingdom","skingdom","sphylum","ssubphylum","sclass", "sorder","sfamily", "sgenus", "sspecies",
                 "superkingdom","kingdom","phylum","subphylum","class", "order","family", "genus", "species")



desTax
x %>% select(class, sclass) 

x %>% filter(kingdom == "Fungi" & skingdom!="Fungi" | 
                     skingdom == "Fungi" & kingdom!="Fungi" |
                     kingdom == "Metazoa" & skingdom!="Metazoa"| 
                     skingdom == "Metazoa" & kingdom!="Metazoa"|
                     kingdom == "Viridiplantae" & skingdom!="Viridiplantae" |
                     skingdom == "Viridiplantae" & kingdom!="Viridiplantae" )

contam_4 <- x %>% filter(class != sclass | 
               order != sorder |
               sphylum != phylum ) 

contam_4 %>% select(assembly_id)
siste_rest <- anti_join(x, contam_4) %>% as_tibble()
siste_rest %>% filter(ssubphylum == "Taphrinomycotina") %>% select(assembly_id)

      
      

# write.table(contam, "contam.csv", row.names = FALSE, quote = FALSE, sep = "\t")
# length(unique(ROD_1$assembly_id))

saveRDS(ROD_1,"ROD_no_seq_v0.1.rds")

# Add sequence: 
library("Biostrings")
fasta_file <- "./../01_results/02_ncbi_results/ribres_4000_cluster.fasta"
sequences <- readDNAStringSet(fasta_file)
seq_db <- NULL
seq_db$header <- as.character(names(sequences))
seq_db$sequence <- unname(as.character(sequences))
str(seq_db)
seq_db <- as_tibble(seq_db)
ROD_seq <- left_join(ROD_1, seq_db, join_by(qseqid == header))

# x<- str_split_fixed(seq_db$header, "\\|", 2)[,2] %>% 
#   str_split_fixed(" ",2)
# IF these have been clustered, there is "size" info in the header. I.e. original copynumber
ROD_seq$seqid <- ROD_seq$qseqid %>% str_split_fixed(";",2) %>% .[,1]
ROD_seq$size <- ROD_seq$qseqid %>% str_split_fixed(";",2) %>% .[,2]  %>% str_split_fixed("=",2) %>% .[,2]
ROD_seq <- ROD_seq %>% dplyr::rename(org.copy.number = rDNA_copynumber)

# Need to adjust the copycount. Summ size pr. genome: 
# rDNA_copy_count <- table(ROD_seq$assembly_id) %>% as.data.frame()
# colnames(rDNA_copy_count) <- c("assembly_id", "rDNA_copy_number")
# ROD_seq <- left_join(ROD_seq, rDNA_copy_count)
# ROD_seq$contaminants <- ROD_seq$org.copy.number - ROD_seq$rDNA_copy_number

ROD_seq$lineage <- paste(ROD_seq$superkingdom, ROD_seq$supergroup, ROD_seq$kingdom, 
                         ROD_seq$phylum, ROD_seq$class, ROD_seq$order, 
                         ROD_seq$family, ROD_seq$genus, ROD_seq$species, 
                         sep = ";") %>% str_replace_all(" ","_")



ROD_seq <- relocate(ROD_seq, seqid, assembly_id, sequence, size,
                    # rDNA_copy_number,
                    org.copy.number, 
                    #contaminants,taxid, 
                    lineage ,superkingdom, supergroup, 
                    kingdom, phylum, class, order, family, genus, species, qseqid)

# ROD_seq %>% filter(contaminants > 0) 
# max(ROD_seq$contaminants)
# ROD_seq %>% filter(contaminants == 368) 
# max(ROD_seq$rDNA_copy_number)
# ROD_seq %>% filter(rDNA_copy_number == 3058) 


##### SIZE SIZE SIZE SIZE !!!!
saveRDS(ROD_seq,"ROD_with_seq_v0.2.rds")

### EXPORT FASTA
df <- readRDS("ROD_with_seq_v0.2.rds")

# header <- paste0(df$seqid,";","size=",df$size,"|",df$lineage , " taxid=",df$taxid, ";")
header <- paste0(df$seqid,";","size=",df$size, " taxid=",df$taxid, ";")
Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
Xfasta[c(FALSE, TRUE)] <- df$sequence
writeLines(Xfasta, "ROD_v0.2.fasta")
getwd()




