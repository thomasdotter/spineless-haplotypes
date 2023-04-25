#      /^\    /^\
#     {  O}  {  O}
#      \ /    \ /
#      //     //       _------_
#     //     //     ./~        ~-_
#    / ~----~/     /              \
#   /         :   ./       _---_    ~-
#   |  \________) :       /~     ~\   |
#   |        /    |      |  :~~\  |   |
#   |       |     |      |  \___-~    |
#   |        \ __/`^\______\.        ./
#    \                     ~-______-~\.
#    .|                                ~-_
#     /_____________________________________~~____
#
# 

############################################################################
##        Anna Thomasdotter                                               ##
##        IMBRSea Spring 2022 - UNIBO                                     ##
##                                                                        ##
##        Thesis script 1                                                 ##
##        Bioinformatics pipeline                                         ##
##        Haplotype generation from metabarcoding data                    ##
##        Raw data from SeaMoBB project                                   ##
############################################################################

# Pipeline components:
#   1. Denoising (unoise) + clustering (swarm)
#   2. Decontamination (decontam)
#   3. Co-occurrence filtering (lulu)
#   4. Abundance filtering
#   5. Remove stop codons
#   6. Taxonomic assignment part 1 (RDP classifier)
#   7. Taxonomic assignment part 2 (BOLDigger)

# Input files: quality filtered fasta files, parsed by sample 
# Obitools steps have already been taken prior to this pipeline; see manuscript  


library(Biostrings)
library(decontam)
library(phyloseq)
library(JAMP)
library(tidyverse)
library(lulu)
library(seqinr)


########################## 1. DENOISE + CLUSTER #########################

wd <- "~/Documents/test_again" # working directory; jamp files will output here
input_folder <- "/store/usr7eco1/Anna_store/Jamp_inputfiles/ITL02" # location of input files
script_location <- "Denoise_SWARM13_280322.R" # location of modified JAMP::Denoise() script, currently not public

# Create a new folder with JAMP and copy input files into it
setwd(wd)
Empty_folder() 
new_folder <- "./A_Empty_Folder/_data"
list_of_files <- list.files(input_folder) 
file.copy(file.path(input_folder,list_of_files), new_folder)

ee <- list.files("A_Empty_Folder/_data", full.names = TRUE) 
ee # check that all the names of your files are here

# Denoise using modified script
source(script_location)
Denoise(files = ee, minsize = 2, minrelsize = 1e-13, OTUmin = 1e-13, 
        minHaploPresence = 1, poolsamples = F, unoise_alpha = 5, 
        withinOTU = 1, minOTUPresence = 1, threads = 4)

# Unoise has already removed chimeras, but SWARM does it once more; 
# so just rename "Chimera" to "OTU" in the output file
# In terminal:
cd Anna_R/real_ITL01
sed 's/Chimera/OTU/g' E_haplo_SWARM_alpha4.txt > 1_denoised_haplotable.txt


########################## 2. DECONTAMINATION #########################

setwd("~/Anna_R/real_ITL01")
df <- read.table("1_denoised_haplotable.txt", sep = ";")
df <- as_tibble(df) %>%
  rename_with(~gsub("_derep.fasta_size_2_denoised.txt", "", .x))

# Format for phyloseq
# Phyloseq works based on samples, OTUs, and (possibly) taxonomic assignment
# In this case, OTU are haplotypes for my data (confusing, but yes)
# It should look like this and become otu_table, where the numbers are abundances:
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       96      50      36      35      59      80      83      63
## OTU2       52      67      39      39      37      57      20      15
## OTU3       94      18      15      11      14      75       1      12
## OTU4       27      88      98     100      59      27      30      30
# Note: when shum does this, the sample names are on the rows. This is adjusted
# when importing by putting OTU = otu_table(as.matrix(count_tab), taxa_are_rows = FALSE)

# Sample information (such as the example in decontam, with control column)
# goes into sample_data() and should look like this for my data:
## sample_ID        sample_or_control     batch
## CRO_SAN_03_2     true_sample           ITL01
## Tneg_1_3         neg_control           ITL02

# Create OTU (haplotype) table:
colnames(df) # Check which cols contain sample names; its 4:579
otu_counts <- df %>% select(c(4:579))
otu_counts_df <- as.data.frame(otu_counts)
row.names(otu_counts_df) <- df$haplotype

# Create sample info table:
sample_info_table <- tibble(sample_ID = names(df[4:579]))
sample_info_table <- sample_info_table %>% 
  separate(sample_ID, into = c("sample_area", "sample_location", 
                               "sample_unit", "sample_plate", "sample_replicate"),
           remove = FALSE) %>%
  mutate(sample_or_control = case_when(sample_area == "Tneg" ~ "negative",
                                       TRUE ~ "real_sample")) 
sample_info_table
# Add batch information
batch = read_csv("batch_index.csv") # external file: contains sample_ID and corresponding batch
sample_info_table = left_join(sample_info_table, batch)
sample_info_df <- as.data.frame(sample_info_table)
row.names(sample_info_df) <- sample_info_df$sample_ID

# Create phyloseq object
OTU = otu_table(otu_counts_df, taxa_are_rows = TRUE)
SAM = sample_data(sample_info_df)
haplodf <- merge_phyloseq(phyloseq(OTU), SAM)

## Start of decontam
# Inspect library sizes
head(sample_data(haplodf))
df1 <- as.data.frame(sample_data(haplodf)) # Put sample_data into a ggplot-friendly data.frame
df1$LibrarySize <- sample_sums(haplodf)
df1 <- df1[order(df1$LibrarySize),]
df1$Index <- seq(nrow(df1))
ggplot(data=df1, aes(x=Index, y=LibrarySize, color=sample_or_control)) + 
  geom_jitter(alpha = 0.5, size = 1, height = 2000) +
  #geom_smooth()+
  ggtitle("Library sizes")

# Remove negative controls with library sizes over 2000
# This is done due to cross-contamination, skip if samples are clean
discard = subset(df1, sample_or_control == "negative" & LibrarySize > 2000)
discard = discard$sample_ID
sample_data(haplodf) = subset(sample_data(haplodf), !(sample_ID %in% discard)) 

# Identify contaminants (prevalence method)
sample_data(haplodf)$is.neg <- sample_data(haplodf)$sample_or_control == "negative"
contamdf.prev <- isContaminant(haplodf, method="prevalence", neg="is.neg", batch = "batch")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# Identify contaminants at stricter threshold (prevalence 0.5)
contamdf.prev05 <- isContaminant(haplodf, method="prevalence", neg="is.neg", 
                                 batch = "batch", threshold=0.5)
table(contamdf.prev05$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
haplodf.pa <- transform_sample_counts(haplodf, function(abund) 1*(abund>0))
haplodf.pa.neg <- prune_samples(sample_data(haplodf.pa)$sample_or_control == "negative", haplodf.pa)
haplodf.pa.pos <- prune_samples(sample_data(haplodf.pa)$sample_or_control == "real_sample", haplodf.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(haplodf.pa.pos), pa.neg=taxa_sums(haplodf.pa.neg), contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant, label = row.names(df.pa))) + 
  #geom_text(hjust=1, vjust=0, size = 3) +
  geom_jitter(alpha = 0.5) + 
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") +
  ggtitle("Prevalence in controls vs real, threshold 0.5, alpha 5")
# Ideally, samples should separate (see decontam tutorial), but again, cross-contamination

# Write prevalence data at 0.5% to file
write.csv(contamdf.prev05, file="decontam_stats.prev05.csv")

# Filter the original dataframe (df) to keep only non-contaminant haplos
not_contamination <- row.names(contamdf.prev05)[contamdf.prev05$contaminant == FALSE]
df_not_contamination <- df %>% filter(haplotype %in% not_contamination)
write_csv(df_not_contamination, "2_decontaminated_prev05.csv") # save as csv



########################## 3. LULU FILTERING #########################


setwd("~/Anna_R/")

# LULU needs:
# 1. R (or R-studio) - LULU was developed in R version 3.3.2.
# 2. LULU r-package 
# 3. An OTU table - produced with any algorithm. In this case, ESV table
# 4. OTU sequences - a file with representavei sequence of each OTU 
# 5. Access to a tool for making a match list (see below, e.g. VSEARCH)

# Create the OTU table (which in this case is an ESV table) from the 
# decontaminated haplotype table without sequences
esv_table1 <- read_csv("2_decontaminated_prev05.csv")
colnames(esv_table1)
esv_table <- select(esv_table1, c(2,4:579)) # Haplotype ID and sample columns
colnames(esv_table)
write_csv(esv_table, "lulu_ESV_table.csv")

# Make a fasta file with haplotype sequences
write.fasta(sequences = as.list(esv_table1$sequences), 
            names = esv_table1$haplotype, "lulu_ESV_sequences.fasta", 
            nbchar = 400, as.string = TRUE)

## In terminal: 
# Make a matchlist using vsearch
vsearch --usearch_global lulu_ESV_sequences.fasta --db lulu_ESV_sequences.fasta --self --id .84 --iddef 1 --userout lulu_match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

## Back in RStudio
# Import files and run lulu
esv_table <- read.csv("lulu_ESV_table.csv", header=TRUE,as.is=TRUE, 
                      row.names = 1)
matchlist <- read.table("lulu_match_list.txt", header=FALSE,as.is=TRUE, 
                        stringsAsFactors=FALSE)
curated_result <- lulu(esv_table, matchlist)

# Inspect results
head(curated_result$curated_table) # Access curated table
head(curated_result$original_table) # Access original table
head(curated_result$curated_otus) # Access list of remaining ESVs
head(curated_result$otu_map) # Information about fate of haplotypes (ESVs)
curated_result$discarded_count
curated_result$curated_count

# Save results in RStudio format
saveRDS(curated_result, file = "lulu_curated_result.rds")

# Save curated haplotype table
curated_table <- curated_result$curated_table
curated_table <- curated_table %>% mutate(haplotype = row.names(curated_table)) %>%
  select(haplotype, everything())
write_csv(curated_table, "3_lulu_curated_table.csv")



########################## 4. ABUNDANCE FILTERING #########################

# Create a dataframe in long format with OTU information and reads

# Load files and merge
curated_table <- read_csv("3_lulu_curated_table.csv") # Reads table
# Add sequences and OTUs
temp = read_csv("2_decontaminated_prev05.csv")
temp = select(temp, OTU, haplotype, sequences)
master = left_join(curated_table, temp)

# Pivot longer
colnames(master) # Columns 2:577 contain samples
lulu_long <- master %>% 
  pivot_longer(c(2:577), names_to = "sample_ID", values_to = "count")

write_csv(lulu_long, "lulu_curated_long.csv")
# lulu_long <- read_csv("lulu_curated_long.csv")

# Filter 1: Relative haplotype abundance in replicate 
# Set haplotype observations in a replicate to 0 (i.e. remove haplotype
# observations) if they constitute less than 0.01% of the total reads in that replicate
lulu_long1 <- lulu_long %>% filter(count > 0) %>% group_by(sample_ID) %>%
  mutate(replicate_total_reads = sum(count)) %>% ungroup() %>% 
  mutate(haplo_rel_replicate = count/replicate_total_reads)
lulu_long1 <- filter(lulu_long1, haplo_rel_replicate > 0.0001)

# Filter 2: Absolute haplotype abundance in replicate
# Remove haplotype observations if they have less than 5 reads in a replicate
lulu_long1 <- lulu_long1 %>% ungroup %>% filter(count > 4)

# Filter 3: Occurrence filter within sample (group of 3 technical replicates)
# Remove haplotype observations if the haplotype does not occur in at least two out
# of three technical replicates for a given sample
lulu_long1 <- lulu_long1 %>% 
  mutate(sample = substr(sample_ID,1,nchar(sample_ID)-2)) %>% 
  group_by(OTU, haplotype, sample) %>% mutate(occurrence = n())
lulu_long1 <- filter(lulu_long1, occurrence > 1)

# Pivot the final dataframe and format
lulu_long2 <- lulu_long1 %>% ungroup() %>% select(haplotype, sample, count)
lulu_long2 <- lulu_long2 %>% group_by(haplotype, sample) %>%
  summarise(count = sum(count)) %>% 
  pivot_wider(names_from = sample, values_from = count, values_fill = 0)
OTU_seq = select(master, OTU, haplotype, sequences)
lulu_filtered <- left_join(lulu_long2, OTU_seq, by = "haplotype")
colnames(lulu_filtered) # columns containing samples: 2:193
lulu_filtered <- lulu_filtered %>% select(OTU, haplotype, 
                                          c(2:193), sequences)

write_csv(lulu_filtered, "4_filtered_table.csv")


########################## 5. REMOVE STOP CODONS #########################
## Remove sequences containing stop codons
# Write fasta to inspect in MEGA for stop codons
write.fasta(sequences = as.list(lulu_filtered$sequences), 
            names = lulu_filtered$haplotype, 
            "mega_check.fasta", nbchar = 400, as.string = TRUE)

# Translate the DNA sequences in MEGA into amino acid seqs
# Here, used invertebrate code but double-checked which haplotypes were removed
# Export result in fasta format and import here again
df_align = readAAMultipleAlignment("mega_checked.fas", format = "fasta")
seq <- as.character(df_align)
df3 = data.frame(haplotype = names(seq), sequences = seq)
df3 = df3 %>% mutate(haplotype = gsub(" ", "_", haplotype))
# Remove haplotypes with * in sequence (i.e. termination or stop)
haplo_stop = df3$haplotype[grepl("\\*", df3$sequences)]
lulu_filtered = filter(lulu_filtered, !(haplotype %in% haplo_stop))

write_csv(lulu_filtered, "5_stops_removed.csv")

########################## 6. TAX. ASSIGNMENT 1 #########################
## Taxonomic assignment using RDP classifier

# Create fasta file for input
temp <- read_csv("5_stops_removed.csv")
write.fasta(sequences = as.list(temp$sequences), names = temp$haplotype, 
            "tax_input.fasta", nbchar = 400, as.string = TRUE)

## In terminal (double-check file locations so that they are correct):
java -Xmx128g -jar /store/usr7eco1/Fra_store/tools/rdp_classifier_2.12/dist/classifier.jar classify -t /store/usr7eco1/Fra_store/tools/COI_mydata_trained/rRNAClassifier.properties -o /home/usr7eco1/Anna_R/final_taxo_output_nofilter.txt /home/usr7eco1/Anna_R/tax_input.fasta  

## Back in RStudio console: Import taxonomy 
tax <- read_tsv("final_taxo_output_nofilter.txt", col_names = FALSE)
tax[1,]

# This file is in a weird format; again, wrangle to get it right
tax1 <- select(tax, -c(2,3,))
new_colnames <- c("haplotype", "x", "p_x", "superkingdom", "p_superkingdom", "kingdom", "p_kingdom",
                  "phylum", "p_phylum", "class", "p_class", "order", "p_order", "family", "p_family",
                  "genus", "p_genus", "species", "p_species")
tax1 <- select(tax, c(1, 4:6, 8,9,11,12,14,15,17,18,20,21,23,24,26,27,29))
names(tax1) <- new_colnames
write_csv(tax1, "RDP_taxonomy_nofilters.csv")

# Filter to only keep metazoans with >0.80 bootstrap at phylum level
unique(tax1$kingdom)
tax1 = tax1 %>% filter(kingdom == "Metazoa" & p_phylum >= 0.80)

write_csv(tax1, "6_RDP_taxonomy.csv")

########################## 7. TAX. ASSIGNMENT 2 #########################
## BOLDigger
# Format for input: fasta file with unique names and sequence in one line (not wrapped)
lulu_filtered <- read_csv("5_stops_removed.csv")

write.fasta(sequences = as.list(lulu_filtered$sequences), 
            names = lulu_filtered$haplotype, "lulu_tax_input.fasta", 
            nbchar = 400, as.string = TRUE)

# Run BOLDigger with jamp hits and flags
# Manually curate haplotypes with questionable assignments at species level
# Export the flagged tab from excel to csv

# Import flagged tab
bold <- read_csv("BOLDResults_lulu_tax.csv")
bold <- bold %>% dplyr::rename("haplotype" = "ID") %>% mutate(haplotype = gsub(">", "", haplotype))

# Add both taxonomies to the haplotable
rdp <- read_csv("6_RDP_taxonomy.csv")
tax_master = select(lulu_filtered, OTU, haplotype)
tax_master = left_join(tax_master, rdp)
tax_master = left_join(tax_master, bold)

## Set all the conditions
tax_master = tax_master %>% 
  mutate(phylum1 = case_when(Similarity >= 98 & (Similarity/100 > p_species | is.na(species) == TRUE) ~ Phylum,
                             p_class >= 0.85 ~ phylum,
                             TRUE ~ ""),
         class1 = case_when(Similarity >= 98 & (Similarity/100 > p_species | is.na(species) == TRUE) ~ Class,
                            p_class >= 0.85 ~ class,
                             TRUE ~ ""),
         order1 = case_when(Similarity >= 98 & (Similarity/100 > p_species | is.na(species) == TRUE) ~ Order,
                            p_class >= 0.85 ~ order,
                             TRUE ~ ""),
         family1 = case_when(Similarity >= 98 & (Similarity/100 > p_species | is.na(species) == TRUE) ~ Family,
                             p_species >= 0.90 ~ family, 
                             TRUE ~ ""),
         genus1 = case_when(Similarity >= 98 & (Similarity/100 > p_species | is.na(species) == TRUE) & is.na(Species) == FALSE ~ Genus,
                            p_species >= 0.90 ~ genus, 
                            TRUE ~ ""),
         species1 = case_when(Similarity >= 98 & (Similarity/100 > p_species | is.na(species) == TRUE) & is.na(Species) == FALSE ~ paste(Genus, Species, sep =" "),
                              p_species >= 0.90 ~ species, 
                              TRUE ~ ""))
# Remove unassigned at class
tax_master = filter(tax_master, class1 != "")
tax_master = select(tax_master, OTU, haplotype, contains("1"))
colnames(tax_master) = gsub("1", "", colnames(tax_master))
tax_master = mutate(tax_master, species = gsub("_", " ", species))

# Filter unwanted classes
unique(tax_master$phylum)
phylum_drop = c("Rhodophyta", "Ochrophyta")
tax_master = tax_master %>% filter(!(phylum %in% phylum_drop))
unique(tax_master$class)
class_drop = c("Insecta", "Clitellata", "Chilopoda")
tax_master = tax_master %>% filter(!(class %in% class_drop))
unique(tax_master$order)
tax_master = tax_master %>% filter(order != "Opiliones")
unique(tax_master$family)


write_csv(tax_master, "tax_RDP_BOLD.csv")

# Bind the curated taxonomy to the final haplotable
tax_master = read_csv("tax_RDP_BOLD.csv")
temp = tax_master %>% select(-OTU)
haplotable = left_join(temp, lulu_filtered)
haplotable = haplotable %>% mutate(species = gsub("_", " ", species)) %>%
  relocate(OTU, .before = "haplotype") %>% 
  arrange(phylum, class, order, genus, species)

## Final clean-up
# Remove OTUs contain control sequences
pos <- read_delim("../positive-ctrls-SEAMoBB.txt", delim = "\t")
colnames(pos)[4] <- "sequences"
test <- left_join(haplotable, pos, by = "sequences")
test <- test %>% select(OTU, mock_type, taxon) %>% filter(!is.na(taxon))
haplotable <- filter(haplotable, !(OTU %in% test$OTU))

# Remove OTUs with conflicting class assignments
test <- haplotable %>% group_by(OTU) %>% mutate(n_class = n_distinct(class)) %>%
  select(OTU, haplotype, class, n_class) %>% filter(n_class > 1)
# 7 OTUs with 17 haplos total contain conflicting class-level assignment 
# and will be excluded from analysis
haplotable <- haplotable %>% filter(!(OTU %in% test$OTU))

# Exclude unwanted classes
unique(haplotable$class)
haplotable <- haplotable %>% filter(class != "Actinopteri" & class != "Actinopterygii")

# Export haplotable
write_csv(haplotable,"7_final_haplotable_reads.csv")

# Convert haplotable to presence absence per sample and export
colnames(haplotable) # 9:200 contain sample names
presence_absence <- function(x, na.rm = FALSE) (ifelse(x>0,1,0))
haplotable_pa <- haplotable %>% mutate_at(c(9:200), .funs = presence_absence)
write_csv(haplotable_pa, "7_final_haplotable_pa.csv")


## Done!

