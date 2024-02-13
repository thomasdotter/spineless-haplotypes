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
##        Thesis script 3                                                 ##
##        Plots                                                           ##
##        Raw data from SeaMoBB project                                   ##
############################################################################
setwd("~/Documents/stats_again")

library(tidyverse)
library(ggplot2)
library(pegas)
library(BioStrings)

############## FIGURE 4: Plot pi per class
df_pi <- read_csv("class_species_pi.csv")

palette10 = c("#f5a543","#ec8035","#ff6d7c","#ff6db6","#d85ac9","#bf74ea",
              "#af89f8","#7491ea","#5aaaf6","#58b9d7","#5dc3a1","#6bc878",
              "#94ce5b","#924900")
# If subsetting for classes with more than 4 OTUs
#df_pi = df_pi %>% group_by(class) %>% filter(n_distinct(OTU) > 4)
# If subsetting for haplotypes with more than 5 occurrences
#df_pi = df_pi %>% filter(occur > 4)
# If subsetting for classes with more than 4 haplotype occurrences
df_pi = df_pi %>% group_by(class) %>% mutate(class_occur = sum(occur)) %>%
  filter(class_occur > 4)
ggplot(df_pi, aes(x=pi)) +
  geom_histogram(aes(fill=class),color="grey80", bins=30) + # change bins here
  facet_grid(class~., switch = "y") + # add class~region if plotting by region
  scale_fill_manual(values = palette10)+
  theme_bw() +
  theme(strip.text.y.left = element_text(angle = 0), legend.position = "none") +
  ylab("Number of species")+
  xlab(expression("Nucleotide diversity "*pi)) +
  #coord_cartesian(xlim = c(-0.001,0.03)) + # comment out to see the whole plot
  scale_y_continuous(breaks = c(2,8,16), limits = c(0, NA), trans = "sqrt")
  # edit last line to set breaks and limits, maybe change transformation (trans)

# Save plot
ggsave("./plots/class_MOTU.png")


############## FIGURE 6: Haplotype networks
# Load data
occurrence <- read_csv("7_final_haplotable_pa.csv") # presence_absence
unique(occurrence$species)
occur <- occurrence %>% filter(inferred_spp == "Halisarca dujardini") # Select species here

## Format data
## Needed: dataframe with columns for haplo id, each location (population) containing frequency, and sequences.
## Each row is a unique haplotype and contains info on its frequency
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

## Create matrices describing each haplotype's location (population) for colors and 
## abundance (frequency) for pegas plotting
## See pegas vignette for info: https://emmanuelparadis.github.io/pegas/PlotHaploNet.pdf
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
pie_colors <- c("#ffca8b",  "#8fbaba", "#995ebc") # LIV PAL CRO; edit here for other colors & add more if more populations
attr(haplonet, "labels") # check order here,output from next two lines should match
row.names(location) 
row.names(abundance) 

# Set options and plot
setHaploNetOptions(link.type.alt = 0, haplotype.outer.color = 0,
                   show.mutation = 2) # see vignette for options
plot(haplonet, pie = location, size = abundance, 
     labels = F, col = pie_colors, fast = TRUE, threshold = 0, legend = TRUE) 
xy <- replot() # if replotting
replot(xy)



############## SUPPLEMENTARY MATERIAL

library(ggpubr) # necessary to get equations for linear reg. (R and p)
# Make a dataframe that combines all values for pi
reg <- read_csv("class_species_pi_region.csv")
reg <- select(reg, -order)
all <- read_csv("class_species_pi.csv")
all <- all %>% mutate(region = "All", reads = occur) %>% select(-c(occur, n_OTUs))
df_pi <- rbind(reg, all)
df_pi <- df_pi %>% mutate(region = factor(region, levels = c("All", "LIV", "PAL", "CRO"))) %>%
  mutate(region = gsub("CRO", "ROV", region))

# Frequency distribution haplotypes per species
ggplot(df_pi, aes(x=n_haplotypes)) + 
  facet_grid(region~., switch = "y") +
  geom_bar(stat = "count") +
  labs(title = "Frequency of haplotypes per species") +
  xlab("Number of haplotypes in species") +
  ylab("Frequency")+
  theme_bw()

# Frequency distribution nuclear diversity per species
ggplot(df_pi, aes(x=pi)) + 
  facet_grid(region~., switch = "y") +
  geom_histogram(bins = 40) +
  labs(title = "Frequency of nucleotide diversity per species") +
  xlab("Nucleotide diversity in species") +
  ylab("Frequency")+
  theme_bw()

# Scatterplot no haplos vs sum sample occurrence
ggplot(df_pi, aes(x=reads, y = n_haplotypes, col = region)) +
  geom_point() +
  stat_smooth(method="lm", se = FALSE) +
  facet_wrap(region~.) +
  labs(title = "Number of haplotypes vs haplotype occurrences") +
  xlab("Haplotype occurrences in species") +
  ylab("Number of haplotypes in species")+
  theme_bw()+
  theme(legend.position = "none")+
  stat_cor(label.x= 100, label.y = 20) +
  stat_regline_equation(label.x = 101, label.y = 23)

# Scatterplot nuc div vs sum sample occurrence
ggplot(df_pi, aes(x=reads, y = pi, col = region)) +
  geom_point() +
  stat_smooth(method="lm", se = FALSE) +
  facet_wrap(region~.) +
  labs(title = "Nucleotide diveristy vs haplotype occurrences") +
  xlab("Haplotype occurrences in species") +
  ylab("Nucleotide diversity in species") +
  theme_bw()+
  theme(legend.position = "none") +
  stat_cor(label.x= 100, label.y = 0.15) +
  stat_regline_equation(label.x = 101, label.y = 0.18)


