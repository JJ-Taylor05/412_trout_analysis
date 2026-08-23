rm(list = ls()) #wipe data

# Load packages
install.packages(c("betapart", "ade4", "labdsv", "ape", "ggplot2", "vegan"))

library(betapart)
library(ade4)
library(labdsv)
library(ape)
library(ggplot2)
library(vegan)
library(readr)
library(tidyr)
library(dplyr)

# Read data all(data set too big for calculations)
records = read.csv("/Users/theomissen/Desktop/PostGrad/GENE412/Module 2/records.csv")

# Read SI rivers
records = read.csv("/Users/theomissen/Desktop/PostGrad/GENE412/Module 2/WilderlabPublicData SI Rivers Fish/records.csv")

# Read NI rivers
records = read.csv("/Users/theomissen/Desktop/PostGrad/GENE412/Module 2/WilderlabPublicData NI Rivers Fish/records.csv")

# Read All still bodies (ponds, wetlands, lakes)
records = read.csv("/Users/theomissen/Desktop/PostGrad/GENE412/Module 2/WilderlabPublicData still bodies Fish/records.csv")




# rename commonname values where it is blank
records$CommonName[is.na(records$CommonName) | records$CommonName == ""] <- records$Name[is.na(records$CommonName) | records$CommonName == ""] #claude

# Reorders the dataframe taken from dylan's code
commonNamesList <- as.list(unique(records$CommonName))
uidsList <- as.list(unique(records$UID))
data_matrix = matrix(NA, nrow = length(uidsList), ncol = length(commonNamesList)) #taken from Dylan's code
data = data.frame(data_matrix)

rownames(data) <- uidsList
colnames(data) <- commonNamesList

# Reorders the dataframe, no filtereing necessary because no names missing
data <- records %>%
  pivot_wider(
    id_cols = UID,
    names_from = CommonName,
    values_from = Count,
    values_fn = sum,
    values_fill = 0
  )


View(data)

# Marine vector all
# Marine vector N.I.
Presence_threshold = 7
Marine_filter = c("Blue cod","Yellowtail kingfish","Bluefin Gurnard","Trevally","Tarakihi","Banded Parrotfish","Blue warehou","Swordfish","Red cod","Orange roughy","Indian mackerel","Mahi-Mahi","Sunfish","Yelloweye mullet")
data$Marine_species_presence = as.integer(rowSums(data[, Marine_filter] >= Presence_threshold) > 0)

# Marine vector S.I Tarakihi`, `Swordfish`, `Red cod`, `Indian mackerel`, `Mahi-Mahi`,`Sunfish` aren't present for S.I river samples
Presence_threshold = 7
Marine_filter = c("Blue cod","Yellowtail kingfish","Bluefin Gurnard","Trevally","Banded Parrotfish","Blue warehou","Orange roughy","Yelloweye mullet")
data$Marine_species_presence = as.integer(rowSums(data[, Marine_filter] >= Presence_threshold) > 0)

# Marine vector still bodies
Presence_threshold = 7
Marine_filter = c("Blue cod","Yellowtail kingfish","Bluefin Gurnard","Trevally","Banded Parrotfish","Orange roughy","Yelloweye mullet")
data$Marine_species_presence = as.integer(rowSums(data[, Marine_filter] >= Presence_threshold) > 0)

# Filter out marine species
data_no_marine = data[data$Marine_species_presence == 0, ]

# Remove Marine variable
data_no_marine$Marine_species_presence <- NULL

# Create a new uid list to name the rows and remove UID column
nm_uidsList = data_no_marine$UID
data_no_marine <- data_no_marine[, -1]

# Create binomial version of dnm
dnmb = ifelse(data_no_marine <7,0,1)
dnmb = as.data.frame(dnmb)
rownames(dnmb) = nm_uidsList

#Trout presence vector
Presence_threshold = 1
Trout_filter = c("Rainbow trout", "Brown trout", "Salmon/Trout", "Trout")
dnmb$Trout_filter_presence = as.integer(rowSums(dnmb[, Trout_filter] >= Presence_threshold) > 0)

# Remove extra trout variables
dnmb$"Rainbow trout" <- NULL
dnmb$"Brown trout" <- NULL
dnmb$"Salmon/Trout" <- NULL
dnmb$"Trout" <- NULL

# Identify columns and rows that are entirely 0 across all remaining samples
empty_cols <- colSums(dnmb) == 0
sum(empty_cols) # How many
colnames(dnmb)[empty_cols] # which ones
empty_rows <- rowSums(dnmb) == 0
sum(empty_rows) # How many
rownames(dnmb)[empty_rows] # which ones

# Remove empty columns and rows
dnmbc <- dnmb[, !empty_cols]
dnmbc <- dnmbc[!empty_rows, ]
View(dnmbc)

# Convert to matrix for beta diversity functions
dnmb_matrix = as.matrix(dnmbc)
View(dnmb_matrix)

# Calculate multi-site beta diversity metrics(for all samples no workie, too big, see claude sampling alternative)
# works for S.I rivers
# works for N.I Rivers
#works for still bodies
beta_multi <- beta.multi(dnmb_matrix)
print(beta_multi)

# Calculate pairwise Sorensen dissimilarity matrix (graham's code)
#works for S.I rivers
# works for N.I Rivers
# works for still bodies
soren <- beta.pair(
  dnmb_matrix,
  index.family = "sorensen"
)
Bsor <- soren$beta.sor



library(betapart)

set.seed(123)          # makes the whole sequence of random draws reproducible
n_iterations = 100    # how many times to repeat the process
sample_size = 100     # how many sites (rows) to randomly draw each time

# Pre-allocate a data frame to hold one row of results per iteration
results <- data.frame(
  turnover = numeric(n_iterations),
  nestedness = numeric(n_iterations),
  total = numeric(n_iterations)
)
for (i in 1:n_iterations) {
  # Randomly pick 'sample_size' row indices, no repeats, from the full matrix
  sample_rows <- sample(1:nrow(dnmb_matrix), size = sample_size)
  # Subset the matrix down to just those randomly chosen rows
  subset_matrix <- dnmb_matrix[sample_rows, ]
  # Run the multi-site beta diversity calculation on just this subset
  beta_result <- beta.multi(subset_matrix, index.family = "sorensen")
  # Store this iteration's three values in row i of the results data frame
  results$turnover[i] <- beta_result$beta.SIM
  results$nestedness[i] <- beta_result$beta.SNE
  results$total[i] <- beta_result$beta.SOR
}
# Average values across all 100 iterations
colMeans(results)
# Variability (standard deviation) across iterations
apply(results, 2, sd)

#Claude's alternative which should work with the vegdist() command we used instead of graham's code on big samples
#this will take a while to run
Bsor <- vegdist(dnmb_matrix, method = "bray", binary = TRUE)




# Histogram of dissimilarities
hist(Bsor,breaks = nrow(dnmb_matrix),main = "Histogram of Sorensen distances",xlab = "Sorensen dissimilarity")
#image manually saved

# Principal Coordinates Analysis (PCoA)
pcoa_result <- pcoa(Bsor)

# Create data frame for plotting PCoA results
dnmbpcoa <- data.frame(
  p1 = pcoa_result$vectors[, 1],
  p2 = pcoa_result$vectors[, 2],
  p3 = pcoa_result$vectors[, 3],
  p4 = pcoa_result$vectors[, 4],
  Trout_present = factor(dnmbc$Trout_filter_presence),
  UID = rownames(dnmbc)
)


# PCoA plot with labels
ggplot(dnmbpcoa, aes(x = p1, y = p2)) +
 
  geom_point(
    aes(fill = Trout_present),
    shape = 21,
    size = 2
  ) +
 
  geom_text(
    aes(label = UID),
    vjust = -1,
    size = 3
  ) +
 
  theme_classic() +
 
  labs(
    x = "PCoA 1",
    y = "PCoA 2",
    fill = ""
  ) +
 
  scale_fill_manual(
    values = c(
      "0" = "forestgreen",
      "1" = "goldenrod"
    )
  ) +
 
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black")
  )


# PCoA plot without labels
ggplot(dnmbpcoa, aes(x = p1, y = p2)) +
 # instruct PCoA groups
  geom_point(
    aes(fill = Trout_present),
    shape = 21,
    size = 2
  ) +
  theme_classic() +
 #axis labels
  labs(
    x = "PCoA 1",
    y = "PCoA 2",
    fill = ""
  ) +
 #Colouring the nodes
  scale_fill_manual(
    values = c(
      "0" = "forestgreen",
      "1" = "goldenrod"
    )
  ) +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black")
  )


# PERMANOVA test for differences between trout presence groups
permanova_result <- adonis2(
  Bsor ~ Trout_present,
  data = dnmbpcoa,
  permutations = 9999
)
print(permanova_result)



