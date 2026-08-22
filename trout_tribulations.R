library(leaflet)
library(htmlwidgets)
library(ggplot2)
library(ggVennDiagram)
library(vegan)

## Trout dataframe ##
# Get trout time
records <- read.csv("records.csv")
samples <- read.csv("samples.csv")

# Make df
records <- records[!is.na(records$CommonName) & records$CommonName != "", ]

commonNamesList <- as.list(unique(records$CommonName))
uidsList <- as.list(unique(records$UID))

orgDF <- data.frame(matrix(NA, nrow = length(uidsList), ncol = length(commonNamesList)))
rownames(orgDF) <- uidsList
colnames(orgDF) <- commonNamesList

for (i in 1:nrow(records)) {
  id <- records$UID[i]
  name <- records$CommonName[i]
  count <- records$Count[i]
  orgDF[as.character(id), as.character(name)] <- count
}

# Flag each site (UID) as trout-present if trout was recorded there
orgDF$UID <- as.character(rownames(orgDF))
orgDF$Has_Trout <- !is.na(orgDF$"Rainbow trout") | !is.na(orgDF$"Brown trout" | !is.na(orgDF$"Trout"))

# Exclude marine species - likely mislabelled
marineSpecies <- c(
  "Blue cod", "Yellowtail kingfish", "Bluefin Gurnard", "Trevally",
  "Tarakihi", "Banded Parrotfish", "Blue warehou", "Swordfish", "Red cod",
  "Orange roughy", "Indian mackerel", "Mahi-Mahi", "Sunfish", "Yelloweye mullet"
)
marineCols <- intersect(colnames(orgDF), marineSpecies)

hasMarine <- rowSums(!is.na(orgDF[, marineCols, drop = FALSE])) > 0 
cat(sprintf(
  "%d sites with marine species removed from analysis\n",
  length(marineCols)
))

orgDF <- orgDF[!hasMarine, ]

## Trout map ##
# Bring in site coordinates from samples.csv
samples$UID <- as.character(samples$UID)
mapData <- merge(orgDF[, c("UID", "Has_Trout")], samples, by = "UID")

# Colour palette: blue = trout present, grey = trout absent
pal <- colorFactor(c("grey15", "mediumblue"), domain = c(FALSE, TRUE))

troutMap <- leaflet(mapData) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~Longitude, lat = ~Latitude,
    color = ~pal(Has_Trout),
    radius = 5, stroke = FALSE, fillOpacity = 0.75,
    popup = ~paste0(
      "UID: ", UID,
      "<br>Trout present: ", Has_Trout,
      "<br>Site: ", ClientSampleID
    )
  ) %>%
  addLegend(
    "bottomright", pal = pal, values = ~Has_Trout,
    title = "Trout Present", labels = c("No", "Yes")
  )

troutMap  

htmlwidgets::saveWidget(troutMap, "trout_terrors.html", selfcontained = TRUE)

## Overlap between trout and predator/prey species ##
# Presence flag for species
orgDF$Has_Eel <- !is.na(orgDF$"Longfin eel")
orgDF$Has_Koaro <- !is.na(orgDF$"Koaro")
orgDF$Has_Inanga <- !is.na(orgDF$"Inanga")

# Site ID list per species
troutsites <- orgDF$UID[orgDF$Has_Trout]
eelsites <- orgDF$UID[orgDF$Has_Eel]
koarosites <- orgDF$UID[orgDF$Has_Koaro]
inangasites <- orgDF$UID[orgDF$Has_Inanga]

# Trout vs Predator
trout_v_eel <- ggVennDiagram(
  list(Trout = troutsites, "Longfin eel" = eelsites), label = "count"
) +
  scale_fill_gradient(low = "white", high = "white") +
  labs(title = "Trout vs Longfin Eel Site Overlap") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

trout_v_eel
ggsave("trout_v_eel_venn.png", plot = trout_v_eel, width = 6, height = 6, dpi = 300, bg = "white")

# Trout vs Inanga
trout_v_inanga <- ggVennDiagram(
  list(Trout = troutsites, "Inanga" = inangasites), label = "count"
) +
  scale_fill_gradient(low = "white", high = "white") +
  labs(title = "Trout vs Inanga Site Overlap") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

trout_v_inanga
ggsave("trout_v_inanga_venn.png", plot = trout_v_inanga, width = 6, height = 6, dpi = 300, bg = "white")


# Trout vs Koaro
trout_v_koaro <- ggVennDiagram(
  list(Trout = troutsites, "Koaro" = koarosites), label = "count"
) +
  scale_fill_gradient(low = "white", high = "white") +
  labs(title = "Trout vs Koaro Site Overlap") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

trout_v_koaro
ggsave("trout_v_koaro_venn.png", plot = trout_v_koaro, width = 6, height = 6, dpi = 300, bg = "white")

## Diversity calculations ##
# Site x species community matrix
nonSpeciesCols <- c("UID", "Has_Trout", "Has_Eel", "Has_Koaro", "Has_Inanga")
speciesCols <- setdiff(colnames(orgDF), nonSpeciesCols)
commMatrix <- orgDF[, speciesCols]
commMatrix[is.na(commMatrix)] <- 0

# Shannon Diversity Index per site
orgDF$Shannon <- diversity(commMatrix, index = "shannon")

# Boxplot template
plotDiversity <- function(data, groupVar, groupLabels, title, colours) {
  data$Group <- factor(data[[groupVar]], labels = groupLabels)
  ggplot(data, aes(x = Group, y = Shannon, fill = Group)) +
    geom_boxplot(outlier.alpha = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
    scale_fill_manual(values = colours) +
    labs(title = title, x = NULL, y = "Shannon Diversity Index") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Diversity of all sites w and without trout (uncorrected)
allsites_diversity <- plotDiversity(
  orgDF, "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - All Sites",
  c("grey15", "mediumblue")
)

allsites_diversity
ggsave("allsites_diversity.png", plot = allsites_diversity, width = 6, height = 5, dpi = 300)

# Diversity of Koaro sites w and without trout (uncorrected)
koarosites_diversity <- plotDiversity(
  orgDF[orgDF$Has_Koaro, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Koaro Sites",
  c("grey15", "lightblue4")
)

koarosites_diversity
ggsave("koarosites_diversity.png", plot = koarosites_diversity, width = 6, height = 5, dpi = 300)

# Diversity of Inanga sites w and without trout (uncorrected)
inangasites_diversity = plotDiversity(
  orgDF[orgDF$Has_Inanga, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Inanga Sites",
  c("grey15", "skyblue")
)

inangasites_diversity
ggsave("inangasites_diversity.png", plot = inangasites_diversity, width = 6, height = 5, dpi = 300)

# Diversity of Trout sites w and without eel (uncorrected)
troutsites_diversity <- plotDiversity(
  orgDF[orgDF$Has_Trout, ], "Has_Eel", c("No Eel", "Eel"),
  "Shannon Diversity - Trout Sites",
  c("grey15", "forestgreen")
)

troutsites_diversity
ggsave("troutsites_diversity.png", plot = troutsites_diversity, width = 6, height = 5, dpi = 300)

## Check for variation in eDNA sampling ##
# Add SeqDepth to df
orgDF$SeqDepth <- rowSums(commMatrix)

# Variation in seq depth
summary(orgDF$SeqDepth)
seqdepth_var <- ggplot(orgDF, aes(x = SeqDepth)) +
  geom_histogram(bins = 50, fill = "grey") +
  labs(title = "eDNA sequencing depth across sites",
       x = "Number of sequences", y = "Number of sites") +
  theme_minimal()

seqdepth_var
ggsave("seqdepth_var.png", plot = seqdepth_var, width = 6, height = 5, dpi = 300)

# Check if Shannon diversity is confounded by seq depth
depthTest <- cor.test(orgDF$SeqDepth, orgDF$Shannon, method = "spearman")
depthTest

# Plot corr
corr_depth_shannon <- ggplot(orgDF, aes(x = SeqDepth, y = Shannon)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", colour = "grey") +
  annotate(
    "label",
    x = Inf, y = Inf,
    label = paste0(
      "Spearman rho = ", round(depthTest$estimate, 2),
      "\np = ", signif(depthTest$p.value, 3)
    ),
    hjust = 1, vjust = 1,
    fill = "white", label.size = 0.3
  ) +
  labs(
    title = "Shannon diversity vs sequencing depth",
    x = "Number of sequences (SeqDepth)", y = "Shannon Diversity Index"
  ) +
  theme_minimal()

corr_depth_shannon
ggsave("corr_depth_shannon.png", plot = corr_depth_shannon, width = 6, height = 5, dpi = 300)

## Rarefy (normalisation transformation) community matrix ##
## Recompute Shannon with common seq depth ##
# Remove sites with low sequencing depth
minDepth <- round(quantile(orgDF$SeqDepth, 0.30, na.rm = TRUE))
keep <- !is.na(orgDF$SeqDepth) & orgDF$SeqDepth >= minDepth
cat(sprintf(
  "Rarefying to %d sequences; dropping %d of %d sites below this depth\n",
  minDepth, sum(!keep), nrow(orgDF)
))

# Rarefaction curve check
set.seed(1)
rarecurve(commMatrix[sample(nrow(commMatrix), 30), ], step = 50, label = FALSE)
abline(v = minDepth, col = "firebrick", lty = 2)

commMatrixRarefy <- commMatrix[keep, ]
orgDFRarefy <- orgDF[keep, ]

set.seed(2)
rarefyMatrix <- rrarefy(commMatrixRarefy, sample = minDepth)
orgDFRarefy$ShannonRarefied <- diversity(rarefyMatrix, index = "shannon")

# Boxplot template
plotDiversityRarefied <- function(data, groupVar, groupLabels, title, colours) {
  data$Group <- factor(data[[groupVar]], labels = groupLabels)
  ggplot(data, aes(x = Group, y = ShannonRarefied, fill = Group)) +
    geom_boxplot(outlier.alpha = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
    scale_fill_manual(values = colours) +
    labs(title = title, x = NULL, y = "Shannon Diversity Index (rarefied)") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Diversity of all sites w and without trout (rarefied)
allsites_diversity_rarefied <- plotDiversityRarefied(
  orgDFRarefy, "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - All Sites (with correction)",
  c("grey15", "mediumblue")
)

allsites_diversity_rarefied
ggsave("allsites_diversity_rarefied.png", plot = allsites_diversity_rarefied, width = 6, height = 5, dpi = 300)

# Diversity of Koaro sites w and without trout (rarefied)
koarosites_diversity_rarefied <- plotDiversity(
  orgDFRarefy[orgDFRarefy$Has_Koaro, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Koaro Sites (with correction)",
  c("grey15", "lightblue4")
)

koarosites_diversity_rarefied
ggsave("koarosites_diversity_rarefied.png", plot = koarosites_diversity_rarefied, width = 6, height = 5, dpi = 300)

# Diversity of Inanga sites w and without trout (rarefied)
inangasites_diversity_rarefied = plotDiversity(
  orgDFRarefy[orgDFRarefy$Has_Inanga, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Inanga Sites (with correction)",
  c("grey15", "skyblue")
)

inangasites_diversity_rarefied
ggsave("inangasites_diversity_rarefied.png", plot = inangasites_diversity_rarefied, width = 6, height = 5, dpi = 300)

# Diversity of Trout sites w and without eel (rarefied)
troutsites_diversity_rarefied <- plotDiversity(
  orgDFRarefy[orgDFRarefy$Has_Trout, ], "Has_Eel", c("No Eel", "Eel"),
  "Shannon Diversity - Trout Sites (with correction)",
  c("grey15", "forestgreen")
)

troutsites_diversity_rarefied
ggsave("troutsites_diversity_rarefied.png", plot = troutsites_diversity_rarefied, width = 6, height = 5, dpi = 300)

## Histogram of # of species per site ##
# Number of species per site
orgDF$Richness <- rowSums(commMatrix > 0)

# Standardise axes
richnessXLim <- c(min(orgDF$Richness, na.rm = TRUE) -1, max(orgDF$Richness, na.rm = TRUE) + 1)
richnessYMax <- max(
  table(orgDF$Richness[orgDF$Has_Trout]),
  table(orgDF$Richness[!orgDF$Has_Trout])
)

# Plot
plotSpeciesNo <- function(data, title, colour, xlim, ymax) {
  meanRichness <- mean(data$Richness, na.rm = TRUE)
  ggplot(data, aes(x = Richness)) +
    geom_histogram(binwidth = 1, boundary = -0.5, fill = colour, colour = "white") +
    scale_x_continuous(limits = xlim) +
    coord_cartesian(ylim = c(0, ymax + 1.05)) +
    annotate(
      "label",
      x = xlim[2], y = ymax + 1.05,
      label = paste0("Mean species per site = ", round(meanRichness, 2)),
      hjust = 1, vjust = 1,
      fill = "white", label.size = 0.3
    ) +
    labs(title = title, x = "Number of species", y = "Number of sites") +
    theme_minimal()
}

richness_trout <- plotSpeciesNo(
  orgDF[orgDF$Has_Trout, ], "Species Richness - Trout Sites", "mediumblue",
  xlim = richnessXLim, ymax = richnessYMax
)

richness_trout
ggsave("richness_hist_trout.png", plot = richness_trout, width = 6, height = 5, dpi = 300)

richness_notrout <- plotSpeciesNo(
  orgDF[!orgDF$Has_Trout, ], "Species Richness - No Trout Sites", "grey15",
  xlim = richnessXLim, ymax = richnessYMax
)

richness_notrout
ggsave("richness_hist_notrout.png", plot = richness_notrout, width = 5, height = 5, dpi = 300)



