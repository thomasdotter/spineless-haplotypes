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
##        Thesis script 2                                                 ##
##        Statistical analysis                                            ##
##        Haplotype generation from metabarcoding data and some plots     ##
##        Raw data from SeaMoBB project                                   ##
############################################################################

# AIM: investigate intraspecific patterns (by class and by species),
# general outputs, regional differences, comparison with vagile and sessile
# records

# INPUT FILES: haplotype table with presence/absence occurrences
# In this case, "7_final_haplotable_pa.csv"
# For comparison with morphologial records: "vagilefauna.csv" and 
# "sessilefauna.csv"


library(pegas)
library(Biostrings)
library(msa)
library(tidyverse)
library(seqinr)
library(vegan)
setwd("~/Documents/stats_again")


########################## Data import  #########################

df <- read_csv("7_final_haplotable_pa.csv")
df_reads <- read_csv("7_final_haplotable_reads.csv")

# Filter to retain real samples
df <- df %>% select(-c(contains("Tneg"), contains("Tpos"))) %>%
  filter(rowSums(across(c(where(is.numeric))))>0)

df_reads <- df_reads %>% 
  select(-c(contains("Tneg"), contains("Tpos"))) %>%
  filter(rowSums(across(c(where(is.numeric))))>0) 

########################## Class-level diversity ############################

# Make a summary by class: number of species, number of OTUs, number of 
# haplotypes, number of occurrences for each region

summary <- df %>% group_by(class) %>% 
  summarise(n_haplo = n_distinct(haplotype),  
            n_OTU = n_distinct(OTU),
            n_species = n_distinct(species, na.rm = TRUE),
            LIV = rowSums(across(contains("LIV"))),
            CRO = rowSums(across(contains("CRO"))),
            PAL = rowSums(across(contains("PAL"))))
summary <- summary %>% group_by(class, n_OTU, n_species, n_haplo) %>% 
  summarise(LIV = sum(LIV), CRO = sum(CRO), PAL = sum(PAL), all = LIV+PAL+CRO)
write_csv(summary, "summary_class.csv")


### Calculate diversity by class using species as units
df_c <- filter(df, !is.na(species)) # remove haplos without species assigned

# Loop to calculate pi, not by region (pooled data)
df_pi <- tibble(class = character(), order = character(), OTU = character(), 
               species = character(),
               n_haplotypes = numeric(), n_OTUs = numeric(),
               occur = numeric(), pi = numeric())
spp <- (unique(df_c$species))
spp <- na.omit(spp)
for (i in 1:length(spp)){
  temp_df <- filter(df_c, species == spp[i]) %>% 
    mutate(occur = rowSums(across(c(where(is.numeric))))) %>% 
    select(class, OTU, haplotype, species, occur, sequences)
  seqs <- temp_df %>% uncount(occur) 
  seqs <- as.DNAbin(DNAStringSet(seqs$sequences))
  bind <- tibble(class = toString(unique(temp_df$class)),
                OTU = toString(unique(temp_df$OTU)), 
                species = toString(unique(temp_df$species)),
                n_haplotypes = n_distinct(temp_df$haplotype),
                n_OTUs = n_distinct(temp_df$OTU),
                occur = sum(temp_df$occur), pi = nuc.div(seqs))
  df_pi <- rbind(df_pi, bind)
} 

df_pi$pi[df_pi$pi == "NaN"]= 0
write_csv(df_pi, "class_species_pi.csv")


# Loop to calculate pi by region
df_pi <- tibble(class = character(), order = character(), OTU = character(), 
               species = character(), n_haplotypes = numeric(),
               reads = numeric(), pi = numeric())
spp <- (unique(df_c$species))
spp <- na.omit(spp)
regions <- c("CRO", "LIV", "PAL")

for (i in 1:length(spp)){
  temp_df <- filter(df_c, species == spp[i])
  for (j in 1:3){
    temp_df1 <- temp_df %>% select(class, order, OTU, haplotype, species, 
                                  contains(regions[j]), sequences) %>%
      mutate(reads = rowSums(across(c(where(is.numeric))))) %>% 
      select(class, order, OTU, haplotype, species, reads, sequences)
    if (sum(temp_df1$reads)>0) {
      seqs <- temp_df1 %>% uncount(reads) 
      seqs <- select(seqs, haplotype, sequences)
      seqs <- as.DNAbin(DNAStringSet(seqs$sequences))
      bind <- tibble(class = toString(unique(temp_df1$class)),
                    order = toString(unique(temp_df1$order)),
                    OTU = toString(unique(temp_df1$OTU)), 
                    species = toString(unique(temp_df1$species)),
                    n_haplotypes = n_distinct(temp_df1$haplotype),
                    reads = sum(temp_df1$reads), pi = nuc.div(seqs),
                    region = regions[j])
      df_pi <- rbind(df_pi, bind)
  } else {
  } #end if else
 } # end regions loop
} # end species or otu loop

df_pi$pi[df_pi$pi == "NaN"]<- 0
write_csv(df_pi, "class_species_pi_region.csv")



########################## AMOVA using ade4 ############################
library(ade4)
# Run loop with AMOVA for each species of interest

# Create empty dataframe to populate in loop
df_var <- tibble(Species = character(), "Variations Between region" = numeric(),
                "Variations Between samples Within region" = numeric(), 
                "Variations Within samples" = numeric(),
                "Total variations" = numeric(), p_withinsample = numeric(),
                p_betweensamples = numeric(),
                p_betweenregion = numeric())
# Select the species of interest, here occurring with 5 or more in each region
# First, reformat dataframe
dat <- df %>% filter(!is.na(species)) %>%
  pivot_longer(cols = c(9:168), names_to = "sample_ID", values_to = "reads")
dat <- dat %>% separate(sample_ID, into = c("region", "site", "unit", "plate"),
                       remove = FALSE) %>% 
  select(species, haplotype, sample_ID, region, site, unit, plate, sequences, reads) %>%
  #filter(!is.na(plate)) %>% # to remove natural substrate 
  mutate(site = paste(region, site, sep = "_"), unit = paste(site, unit, sep = "_"))

# Second, select species of interest
spp <- dat %>% group_by(species, region) %>% filter(sum(reads) > 4) %>%
  group_by(species) %>% filter(n_distinct(region) == 3)
spp <- unique(spp$species)

# Set structures (i.e. group levels in the nested AMOVA)
## "Structures: a data frame containing, in the jth row and the kth column, 
## the name of the group of level k to which the jth population belongs"
structures <- dat %>% select(site, region) %>% distinct() # to site
structures <- as.data.frame(structures, drop = FALSE)
structures$region <- as.factor(structures$region)
row.names(structures) <- structures[,1]
structures <- structures[,-1, drop = FALSE]
structures

# Loop to run AMOVA for each species
for (i in 1:length(spp)){
## Distances: create distance matrix
# convert sequences to dnabin for distance matrix
temp <- dat %>% filter(species == spp[i]) %>% select(haplotype, sequences) %>% distinct()
seqs <- temp$sequences
names(seqs) <- paste0(temp$haplotype)
seqs <- as.DNAbin(DNAStringSet(seqs))
# calculate distance matrix
distances <- cailliez(dist.dna(seqs))

## Samples: a data frame with haplotypes (or genotypes) as rows, 
## populations as columns and abundance as entries
temp <- dat %>% filter(species == spp[i])
samples <- temp %>% select(haplotype, site, reads) %>% group_by(haplotype, site) %>%
  summarise(abundance = sum(reads))
samples <- pivot_wider(samples, values_from = abundance, names_from = site)
samples <- as.data.frame(samples)
row.names(samples) <- samples[,1]
samples <- samples[,-1]

# Run AMOVA
(am <- amova(samples, distances, structures))

# Extract the variation percentages 
var <- am$componentsofcovariance[,2]
test <- randtest(am, nrepet = 999)
species <- spp[i]
bind <- tibble(Species = species, "Variations Between region" = var[1], 
              "Variations Between samples Within region" = var[2], 
              "Variations Within samples" = var[3], 
              "Total variations" = var[4],
              p_withinsample = test$pvalue[1],
              p_betweensamples = test$pvalue[2],
              p_betweenregion = test$pvalue[3])
# Bind to empty df
df_var <- bind_rows(df_var, bind)

} # end loop



###################### LOOP: Haplotype networks ##########################
df_s <- df %>% filter(!is.na(species)) # select haplos w/ species assigned

# Select the species of interest, here any spp with more than 4 haplotypes
multispp <- df_s %>% group_by(species) %>% 
  summarise(n_haplo = n_distinct(haplotype)) %>%
  ungroup() %>% filter(n_haplo > 4)

# Set pie colors and haplonet options; see pegas manual for details
pie_colors <- c("#ffca8b","#8fbaba","#995ebc")
setHaploNetOptions(link.type.alt = 0, haplotype.outer.color = 0,
                   show.mutation = 2, pie.colors.function = pie_colors)
#dir.create("plots")
setwd("./plots") # Get into the folder you want to save the plots in

# Loop to plot and save for each species of interest
for (i in 1:nrow(multispp)){
  occur <- df_s %>% filter(species == multispp$species[i]) # FIlter for the current spp
  occur <- pivot_longer(occur, 9:168, names_to = "sample_location", values_to = "reads")
  occur <- occur %>% filter(reads > 0) %>% 
    mutate(substrate = case_when(str_length(sample_location)>10 ~ "arms", TRUE ~ "natural"),
           location_id = substr(sample_location, 1,7),
           sample_location = paste0(location_id, "_",substrate)) %>% 
    select(-c(substrate, location_id))
  # summarise the dataframe
  pop1 <- occur %>% select(haplotype, sample_location, reads) %>% 
    group_by(haplotype, sample_location) %>% summarise(reads = sum(reads))
  pop1 <- pop1 %>% pivot_wider(names_from = sample_location, values_from = reads,
                              values_fill = 0)
  # Add sequences
  temp <- select(df_s, haplotype, sequences)
  pop2 <- left_join(pop1, temp, by = "haplotype")
  
  # Get occurrences per location
  pop2 <- pop2 %>% mutate(CRO = rowSums(across(contains("CRO_"))),
                         LIV = rowSums(across(contains("LIV_"))),
                         PAL = rowSums(across(contains("PAL_")))) %>%
    select(haplotype, LIV, PAL, CRO, sequences)
  
  # Format for pegas
  seq <- pop2$sequences
  names(seq) <- pop2$haplotype
  alignSeq1 <- as.DNAbin(DNAStringSet(seq)) # get into DNAbin format
  haploSeq <- haplotype(alignSeq1, labels = names(alignSeq1)) # create haplotypes
  haplonet <- haploNet(haploSeq) # calculate network links
  
  # Create matrix for location to set pie segments later on
  # Row names: haplotype label, column names: population (here, region)
  location <- pop2 %>% ungroup %>% select(haplotype, contains("LIV"), 
                                         contains("PAL"), contains("CRO"))
  location <- as.data.frame(location)
  row.names(location) <- location[,1]
  location <- location[,-1]
  location <- as.matrix(location, drop = FALSE)
  
  # Create matrix or vector for "abundances" to set circle sizes later
  # Row names: haplotype label
  abundance <- pop2 %>% ungroup() %>% 
    mutate(reads <- rowSums(across(c(where(is.numeric)))))
  abundance <- as.data.frame(select(abundance, haplotype, reads), drop = FALSE)
  abundance <- as.data.frame(abundance)
  row.names(abundance) <- abundance[,1]
  abundance <- abundance[,-1, drop = FALSE]
  abundance <- as.matrix(abundance, drop = FALSE)
  
  # Plot the haplotype network and export as svg and png 
  svg(paste0(gsub(" ", "_", multispp$species[i]), ".svg"))
  plot(haplonet, pie = location, size = abundance, 
       labels = FALSE,  fast = TRUE, threshold = 0, legend = F)
  mtext(text = multispp$species[i], side = 3)
  dev.off()
  png(paste0(gsub(" ", "_", multispp$species[i]), ".png"))
  plot(haplonet, pie = location, size = abundance, 
       labels = FALSE,  fast = TRUE, threshold = 0, legend = F)
  mtext(text = multispp$species[i], side = 3)
  dev.off()
  
}

setwd("..") # Return to original wd


#################### Haplotype networks, not looped ######################
# Load data
occurrence <- read_csv("7_final_haplotable_pa.csv") # presence_absence
unique(occurrence$species)
occur <- occurrence %>% filter(inferred_spp == "Halisarca dujardini") # Select species here
# Pivot longer
colnames(occur)
occur <- occur %>% pivot_longer(9:200, names_to = "sample_location", values_to = "reads")
occur <- occur %>% filter(!grepl("Tpos|Tneg", sample_location) & reads > 0) %>% 
  mutate(substrate = case_when(str_length(sample_location)>10 ~ "arms", TRUE ~ "natural"),
         location_id = substr(sample_location, 1,7),
         sample_location = paste0(location_id, "_",substrate)) %>% select(-c(substrate, location_id))
# Summarise for haplotype
pop1 <- occur %>% select(haplotype, sample_location, reads) %>% 
  group_by(haplotype, sample_location) %>% summarise(reads = sum(reads))
pop1 <- pop1 %>% pivot_wider(names_from = sample_location, values_from = reads,
                            values_fill = 0)
# Add sequences
temp <- select(occurrence, haplotype, sequences)
pop2 <- left_join(pop1, temp, by = "haplotype")

# Summarise for location
pop2 <- pop2 %>% mutate(CRO = rowSums(across(contains("CRO_"))),
                       LIV = rowSums(across(contains("LIV_"))),
                       PAL = rowSums(across(contains("PAL_")))) %>%
  select(haplotype, LIV, PAL, CRO, sequences)
#Format for haplonet
seq <- pop2$sequences
names(seq) <- pop2$haplotype
alignSeq1 <- as.DNAbin(DNAStringSet(seq))
haploSeq <- haplotype(alignSeq1, labels = names(alignSeq1)) # create haplotypes
haploSeq # inspect
haplonet <- haploNet(haploSeq) # create haplonet
haplonet # inspect

# Location matrix for pies
colnames(pop2)
location <- pop2 %>% ungroup %>% select(haplotype, contains("LIV"), contains("PAL"), contains("CRO"))
location <- as.data.frame(location)
row.names(location) <- location[,1]
location <- location[,-1]
location <- as.matrix(location, drop = FALSE)

# Abundance matrix (or vector?) for size
colnames(pop2)
abundance <- pop2 %>% ungroup() %>% mutate(reads = rowSums(across(c(where(is.numeric)))))
abundance <- as.data.frame(select(abundance, haplotype, reads), drop = FALSE)
abundance <- as.data.frame(abundance)
row.names(abundance) <- abundance[,1]
abundance <- abundance[,-1, drop = FALSE]
abundance <- as.matrix(abundance, drop = FALSE)

# Pie colors function
pie_colors <- c("#ffca8b",  "#8fbaba", "#995ebc") # LIV PAL CRO
attr(haplonet, "labels") # check order; should match for all three
row.names(location)
row.names(abundance)
# Set options and plot
setHaploNetOptions(link.type.alt = 0, haplotype.outer.color = 0,
                   show.mutation = 2)
plot(haplonet, pie = location, size = abundance, 
     labels = F, col = pie_colors, fast = TRUE, threshold = 0, legend = TRUE) 
xy <- replot()
replot(xy)


############### Rarefaction and accumulation curves using vegan ###############
# Dataframe with haplotypes and reads is df_reads
# Species (or haplotypes, or MOTUs) should be columns, and samples (or plates, or replicates)
# should be rows. Transpose for correct input format.
str(df_reads)
names(df_reads)

#### By species 
temp <- df_reads %>% select(8:168) %>% filter(!is.na(species)) %>% 
  group_by(species) %>% summarise(across(2:160, sum))
  
temp <- as_data_frame(temp)
names <- temp$species
temp1 <- as.data.frame(t(temp[,-1]))
names(temp1) <- names # this worked out
# are there na's
test <- temp1[rowSums(is.na(temp1)) > 0, ] 
nrow(test) # zero so this is good


df1 <- temp1
# Plot rarefaction
rarecurve(df1, xlab = "Number of reads", ylab = "Number of species", 
          label = TRUE, log = "x")

# Make separate rarefaction curves per region
regions <- c("PAL", "LIV", "CRO")
titles = c("Palinuro", "Livorno", "Rovinj")

for (i in 1:length(regions)){
  region = regions[i]
  temp <- df_reads %>% select(species, contains(region)) %>% 
    filter(!is.na(species)) %>% group_by(species) %>%
    summarise(across(-1, sum))
  
  temp <- as_data_frame(temp)
  names <- temp$species
  temp1 <- as.data.frame(t(temp[,-1]))
  names(temp1) <- names # this worked out
  
  df1 <- temp1
  # Plot rarefaction
  rarecurve(df1, xlab = "Number of reads", ylab = "Number of species", 
            label = FALSE, log = "x", main = titles[i])
}


### By haplotypes
names(df_reads)
temp <- df_reads %>% select(haplotype, 9:168) %>% group_by(haplotype) %>%
  summarise(across(2:160, sum))

temp <- as_data_frame(temp)
names <- temp$haplotype
temp1 <- as.data.frame(t(temp[,-1]))
names(temp1) <- names # this worked out

df1 <- temp1
# Plot rarefaction
# Make separate rarefaction curves per region
regions <- c("PAL", "LIV", "CRO")
titles = c("Palinuro", "Livorno", "Rovinj")

for (i in 1:length(regions)){
  region = regions[i]
  temp <- df_reads %>% 
    filter(!is.na(species)) %>% # filter for haplotypes with species assigned
    select(haplotype, contains(region)) %>% 
    group_by(haplotype) %>% summarise(across(-1, sum))
  
  temp <- as_data_frame(temp)
  names <- temp$haplotype
  temp1 <- as.data.frame(t(temp[,-1]))
  names(temp1) <- names # this worked out
  
  df1 <- temp1
  # Plot rarefaction
  rarecurve(df1, xlab = "Number of reads", ylab = "Number of haplotypes with species assigned", 
            label = FALSE, log = "x", main = titles[i])
}



############### Species accumulation curves using vegan ##################
# Input here is a table with row names corresponding to "site" (in this case, ARMS or sample)
# Columns are "species" units (species, haplotypes, or MOTUs)

# Start with species and samples
names(df_reads)
temp <- df_reads %>% 
  filter(!is.na(species)) %>%
  select(species, 9:168) %>% group_by(species) %>%
  summarise(across(2:160, sum))

temp <- as_data_frame(temp)
names <- temp$species
temp1 <- as.data.frame(t(temp[,-1]))
names(temp1) <- names # this worked out


sac <- specaccum(temp1)
plot(sac, ci.type = "polygon", ci.col = "#b6dbff", ylab = "Number of species", xlab = "Number of samples")

# Haplotypes and samples
names(df_reads)
temp <- df_reads %>% 
  filter(!is.na(species)) %>% # remove haplotypes without species assignment
  select(haplotype, 9:168) %>% group_by(haplotype) %>%
  summarise(across(2:160, sum))

temp <- as_data_frame(temp)
names <- temp$haplotype
temp1 <- as.data.frame(t(temp[,-1]))
names(temp1) <- names # this worked out


sac <- specaccum(temp1)
plot(sac, ci.type = "polygon", ci.col = "#b6dbff", ylab = "Number of haplotypes", xlab = "Number of samples")


