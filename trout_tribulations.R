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
orgDF$Has_Trout <- !is.na(orgDF$"Rainbow trout") | !is.na(orgDF$"Brown trout") | !is.na(orgDF$"Trout")

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

## Native species-specific richness ##
# Identify native species in the dataset (based on NZFFD dataset)
nativeSpecies <- c(
  # Eels 
  "Eel", "Eels", "Longfin eel", "Shortfin eel", "Speckled longfin eel",
  
  # Galaxiids, kokopu, and torrentfish 
  "Galaxiids", "Inanga", "Koaro",
  "Banded kokopu", "Giant kokopu", "Giant or shortjaw kokopu", "Shortjaw kokopu",
  "Alpine galaxias", "Bignose galaxias", "Canterbury galaxias",
  "Clutha flathead galaxias", "Clutha or Teviot flathead galaxias",
  "Clutha, Pomahaka or Taieri flathead galaxias", "Dune lakes galaxias",
  "Dusky galaxias", "Dusky or roundhead galaxias", "Dwarf galaxias",
  "Dwarf or Alpine galaxias", "Eldons galaxias", "Gollum galaxias",
  "Longjawed galaxias", "Lowland longjaw galaxias", "Nevis gollum galaxias",
  "Northern flathead galaxias", "Pomahaka galaxias", "Roundhead galaxias",
  "Southern flathead galaxias", "Southern, Clutha, Pomahaka or Taieri flathead galaxias",
  "Taieri flathead galaxias", "Taieri or Clutha flathead galaxias",
  "Taieri or Teviot flathead galaxias", "Taieri or southern flathead galaxias",
  "Waitaki lowland longjaw galaxias", "Waitaki upland longjaw galaxias",
  "Torrentfish", "Torrentfishes",
  
  # Mudfish 
  "Black mudfish", "Brown mudfish", "Canterbury mudfish", "Northland mudfish",
  "Mudfish", "mudfish",
  
  # Bullies 
  "Bluegilled bully", "Common bully", "Common/Cran/Dinahs bully", "Crans bully",
  "Giant bully", "Kaharore bully", "Redfin bully", "Upland bully",
  "Upland or kaharore bully", "Bullies",
  
  # Smelt 
  "Common smelt", "Stokell's smelt", "Smelt",
  
  # Lamprey 
  "Lamprey", "Pouched lamprey",
  
  # Flounder 
  "Black Flounder, freshwater flounder"
)

nativeCols <- intersect(colnames(commMatrix), nativeSpecies)

# Native species richness
orgDF$NativeRichness <- rowSums(commMatrix[, nativeCols, drop = FALSE] > 0)

# Shared axis limits
nativeXLim <- c(min(orgDF$NativeRichness, na.rm = TRUE) - 1, max(orgDF$NativeRichness, na.rm = TRUE) + 1)
nativeYMax <- max(
  table(orgDF$NativeRichness[orgDF$Has_Trout]),
  table(orgDF$NativeRichness[!orgDF$Has_Trout])
)

# Plot template
plotNativeSpeciesNo <- function(data, title, colour, xlim, ymax) {
  meanNative <- mean(data$NativeRichness, na.rm = TRUE)
  
  ggplot(data, aes(x = NativeRichness)) +
    geom_histogram(binwidth = 1, boundary = -0.5, fill = colour, colour = "white") +
    scale_x_continuous(limits = xlim) +
    coord_cartesian(ylim = c(0, ymax + 1.05)) +
    annotate(
      "label",
      x = xlim[2], y = ymax + 1.05,
      label = paste0("Mean native species/site = ", round(meanNative, 2)),
      hjust = 1, vjust = 1,
      fill = "white", label.size = 0.3
    ) +
    labs(title = title, x = "Number of native species", y = "Number of sites") +
    theme_minimal()
}

nativerichness_trout <- plotNativeSpeciesNo(
  orgDF[orgDF$Has_Trout, ], "Native Species Richness - Trout Sites", "mediumblue",
  xlim = nativeXLim, ymax = nativeYMax
)

nativerichness_trout
ggsave("native_richness_hist_trout.png", plot = nativerichness_trout, width = 6, height = 5, dpi = 300)

nativerichness_notrout <- plotNativeSpeciesNo(
  orgDF[!orgDF$Has_Trout, ], "Native Species Richness - No Trout Sites", "grey15",
  xlim = nativeXLim, ymax = nativeYMax
)

nativerichness_notrout
ggsave("nativerichness_hist_notrout.png", plot = nativerichness_notrout, width = 6, height = 5, dpi = 300)

## Habitat covariates ##
# Pull habitat data
habitatData <- samples[, c("UID", "EnvironmentType", "Region")]
habitatData$UID <- as.character(habitatData$UID)
modelDF <- merge(orgDF, habitatData, by = "UID")

# Shannon diversity ~ trout presence + habitat
shannonModel <- lm(Shannon ~ Has_Trout + EnvironmentType + Region, data = modelDF)
summary(shannonModel)

# Overall richness ~ trout presence + habitat 
richnessModel <- glm(Richness ~ Has_Trout + EnvironmentType + Region, data = modelDF, family = poisson)
summary(richnessModel)

# Native richness ~ trout presence + habitat
nativeRichnessModel <- glm(NativeRichness ~ Has_Trout + EnvironmentType + Region, data = modelDF, family = poisson)
summary(nativeRichnessModel)


## Stratified comparison of trout/no trout in each habitat
# Plot 
plotDiversityStratified <- ggplot(modelDF, aes(x = Has_Trout, y = Shannon, fill = Has_Trout)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.6) +
  scale_x_discrete(labels = c("No Trout", "Trout")) +
  scale_fill_manual(values = c("grey15", "mediumblue")) +
  facet_wrap(~ EnvironmentType, scales = "free_y") +
  labs(
    title = "Shannon Diversity - Trout vs No Trout in each habitat type",
    x = NULL, y = "Shannon Diversity Index"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
    
plotDiversityStratified
ggsave("diversity_stratified_by_habitat.png", plot = plotDiversityStratified, width = 8, height = 6, dpi = 300)

## Upstream vs downstream of waterfalls effect ##
## Using a small subset of the data that has manually assigned location ##
aboveWaterfallUIDs <- c(
  "406581", "406582", "406583", "406584", "406585", "406586",  # Moutere River U/S Moutere Weir
  "511455",                                                     # Immediately upstream St Ronans Ave weir
  "519349", "519350", "519351", "519352", "519353", "519354",  # Pukatea Stream Site 2 - above ornamental weirs
  "519358",                                                     # Pukatea Stream Site 3 - above barrier 6
  "524825",                                                     # Waitawhara Stream above weir
  "527455", "527459",                                           # silverstream - US of weir
  "531234", "531235", "531236", "531239", "531240", "531242",  # Wharerangi Stream above waterfall
  "538756", "538760",                                           # Upstream River Rd Barrier
  "558166", "558167", "558168", "558169", "558170", "558171",  # Mangatu Stream US of waterfall
  "703861", "703862", "703863", "703864", "703865", "703866",  # Te Arai River Bush Intake Above Weir
  "714002", "714006",                                           # U/S Weir
  "716092", "716094",                                           # Wharekopae above falls
  "727781", "727782", "727783", "727784", "727785", "727786",  # Little Huia Stream above weir
  "727841", "727842", "727843", "727844", "727845", "727846",  # Big Huia Stream above weir
  "731731", "731732", "731733", "731734", "731735", "731736",  # Mangahere u/s waterfall
  "736781", "736782", "736783", "736784", "736785", "736786",  # Orongorongo River above weir
  "739911", "739912", "739913", "739914", "739915", "739916"   # Waterfall Rd upstream of culvert
)

belowWaterfallUIDs <- c(
  "519343", "519344", "519345", "519346", "519347", "519348",  # Pipitea Stream Site 1 - below all barriers
  "520155", "520275", "520278",                                 # Mangaiwi Stream downstream of Weir
  "528206",                                                     # Falls Ck Downstream
  "539593", "539638",                                           # Silverstream - DS of weir
  "703911", "703912", "703913", "703914", "703915", "703916",  # Te Arai Trib Below Intake Weir
  "731701", "731702", "731703", "731704", "731705", "731706",  # Boundary stream below Shine falls
  "736684"                                                      # Orongorongo River below weirs
)

waterfallDF <- orgDF[orgDF$UID %in% c(aboveWaterfallUIDs, belowWaterfallUIDs), ]
waterfallDF$Barrier <- ifelse(waterfallDF$UID %in% aboveWaterfallUIDs, "Above Waterfall", "Below Waterfall")

if (nrow(waterfallDF) == 0) {
  cat("No UIDs entered yet - fill in aboveWaterfallUIDs / belowWaterfallUIDs to run this comparison\n")
} else {
  waterfallPlot <- ggplot(waterfallDF, aes(x = Barrier, y = Shannon, fill = Barrier)) +
    geom_boxplot(outlier.alpha = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    scale_fill_manual(values = c("Above Waterfall" = "steelblue", "Below Waterfall" = "grey15")) +
    labs(
      title = "Shannon Diversity - Above vs Below Waterfall",
      x = NULL, y = "Shannon Diversity Index"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  waterfallPlot
  ggsave("diversity_waterfall_barrier.png", plot = waterfallPlot, width = 6, height = 5, dpi = 300)
  
  nAbove <- sum(waterfallDF$Barrier == "Above Waterfall")
  nBelow <- sum(waterfallDF$Barrier == "Below Waterfall")
  
  if (nAbove >= 3 && nBelow >= 3) {
    waterfallTest <- wilcox.test(Shannon ~ Barrier, data = waterfallDF)
    cat(sprintf(
      "Above vs below waterfall: n(above)=%d, n(below)=%d, Wilcoxon p = %.4g\n",
      nAbove, nBelow, waterfallTest$p.value
    ))
  } else {
    cat(sprintf(
      "Too few sites for a test yet (n(above)=%d, n(below)=%d) - add more UIDs\n",
      nAbove, nBelow
    ))
  }
}
