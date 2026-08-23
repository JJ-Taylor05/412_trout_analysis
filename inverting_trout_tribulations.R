library(leaflet)
library(htmlwidgets)
library(ggplot2)
library(ggVennDiagram)
library(vegan)
library(Matrix)

## Trout presence (from the fish dataset) ##
records <- read.csv("records.csv")
samples <- read.csv("samples.csv")

records <- records[!is.na(records$CommonName) & records$CommonName != "", ]
troutUIDs <- unique(records$UID[records$CommonName %in% c("Rainbow trout", "Brown trout", "Trout")])

## Invertebrate community dataframe ##
# Get invert time
invert_records <- read.csv("invert_records.csv")
invert_samples <- read.csv("invert_samples.csv")

# Unlike the fish CommonName column, ~19% of invert rows have no CommonName
# at all (many hits are only identified to genus/species level, e.g. a
# specific Deleatidium sequence variant). Dropping every blank-CommonName
# row the way the original script does would throw away a fifth of the
# invertebrate data, so instead we fall back to the scientific Name when
# CommonName is missing. TaxonKey is then our "species column" for the rest
# of the script, in place of CommonName.
invert_records$TaxonKey <- ifelse(
  !is.na(invert_records$CommonName) & invert_records$CommonName != "",
  invert_records$CommonName,
  invert_records$Name
)
invert_records <- invert_records[!is.na(invert_records$TaxonKey) & invert_records$TaxonKey != "", ]

# A small lookup table, one row per taxon, carrying its Order/Family/Genus.
# We'll use this a few times below to classify *columns* of the community
# matrix (e.g. "is this taxon a mayfly?") without re-scanning all ~814k raw
# rows every time.
taxonLookup <- unique(invert_records[, c("TaxonKey", "Order", "Family", "Genus")])
taxonLookup <- taxonLookup[!duplicated(taxonLookup$TaxonKey), ]

# Make df
# The original fish script fills the site x species matrix with a for-loop,
# assigning one cell of a data.frame at a time. That works fine for ~130k
# fish rows, but the invert file has ~814k rows, and single-cell data.frame
# assignment gets much slower as the table grows (each assignment re-scans
# the whole frame), so a loop like that would take a very long time here.
# Instead we build the same site x taxon count matrix with a sparse matrix:
# every (site, taxon, count) triplet is placed directly by row/column index,
# with no cell-by-cell searching. It produces an identical table, just much
# faster (this build takes well under a second here).
uidFactor <- factor(invert_records$UID)
taxonFactor <- factor(invert_records$TaxonKey)

sparseCounts <- sparseMatrix(
  i = as.integer(uidFactor),
  j = as.integer(taxonFactor),
  x = invert_records$Count,
  dims = c(nlevels(uidFactor), nlevels(taxonFactor)),
  dimnames = list(levels(uidFactor), levels(taxonFactor))
)

orgDF <- as.data.frame(as.matrix(sparseCounts))
orgDF$UID <- rownames(orgDF)

# Flag each site (UID) as trout-present if trout was recorded there in the
# fish dataset
orgDF$Has_Trout <- orgDF$UID %in% as.character(troutUIDs)

## ----------------------------------------------------------------------
## Note on "native invertebrate species", used throughout this script:
## eDNA samples pick up everything the water passed over or near - alongside
## genuine stream fauna, invert_records.csv also contains terrestrial
## incidentals (aphids, moths, garden slugs), marine/aquaculture species
## (Pacific oysters, tiger shrimp), and other non-target hits. There's no
## reliable column to separate native from exotic taxon-by-taxon, so per
## your steer we're using a simplified proxy: the classic EPTOM insect
## orders (Ephemeroptera/mayflies, Plecoptera/stoneflies,
## Trichoptera/caddisflies, Odonata/dragon-&-damselflies,
## Megaloptera/dobsonflies) plus three native non-insect groups (koura
## freshwater crayfish, Paratya freshwater shrimp, and kakahi freshwater
## mussel/native mud snail). This is a coarse simplification, not a full
## taxonomic native/exotic classification - it will miss some genuinely
## native taxa (e.g. native worms, native beetles) and could in principle
## include a self-introduced non-native within one of these orders. Treat
## "NativeRichness" below as a proxy for "core native stream macroinvertebrate
## richness", and refine the aquaticOrders/aquaticFamilies/aquaticGenera
## lists below if you have a better taxon reference (e.g. the MCI taxa list).
## ----------------------------------------------------------------------
aquaticOrders <- c("Ephemeroptera", "Plecoptera", "Trichoptera", "Odonata", "Megaloptera")
aquaticFamilies <- c("Parastacidae", "Atyidae", "Hyriidae")
aquaticGenera <- c("Potamopyrgus")

isAquaticTaxon <- taxonLookup$Order %in% aquaticOrders |
  taxonLookup$Family %in% aquaticFamilies |
  taxonLookup$Genus %in% aquaticGenera
aquaticTaxonKeys <- taxonLookup$TaxonKey[isAquaticTaxon]

cat(sprintf(
  "%d of %d detected taxa classed as core native aquatic macroinvertebrates\n",
  length(aquaticTaxonKeys), nrow(taxonLookup)
))

## Invertebrate map ##
# Bring in site coordinates from invert_samples.csv
invert_samples$UID <- as.character(invert_samples$UID)
mapData <- merge(orgDF[, c("UID", "Has_Trout")], invert_samples, by = "UID")

# Colour palette: blue = trout present, grey = trout absent
pal <- colorFactor(c("grey15", "mediumblue"), domain = c(FALSE, TRUE))

invertMap <- leaflet(mapData) %>%
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

invertMap

htmlwidgets::saveWidget(invertMap, "invert_terrors.html", selfcontained = TRUE)

## Overlap between trout and native invertebrates threatened by trout ##
# These three groups are well documented in NZ freshwater ecology as
# affected by trout:
#  - Koura / freshwater crayfish (Paranephrops spp.): brown trout are a
#    documented predator of small/juvenile koura and suppress their
#    foraging activity (Usio & Townsend 2000, NZ J Mar Freshw Res 34:557-567,
#    https://doi.org/10.1080/00288330.2000.9516956; Shave et al. 1994, NZ J
#    Ecol 18:1-10).
#  - Deleatidium mayfly: the best-studied invertebrate prey item in NZ
#    trout-stream research - trout presence changes its drift behaviour and
#    grazing activity (McIntosh & Townsend 1996, Oecologia 108:174-181,
#    https://doi.org/10.1007/BF00333229; Flecker & Townsend 1994, Ecological
#    Applications 4:798-807, https://doi.org/10.2307/1942009).
#  - Paratya freshwater shrimp: NZ's only native freshwater shrimp, with a
#    range now reduced and patchy largely attributed to trout introduction
#    (De Grave et al., IUCN Red List account for Paratya curvirostris).
# We use the Genus column here rather than TaxonKey/CommonName, because
# CommonName is inconsistent for these taxa (e.g. many Deleatidium hits are
# only labelled with a sequence-variant code, not "Deleatidium"), while
# Genus reliably groups all of a taxon's records together regardless of
# how finely it was identified.
koaraUIDs <- unique(invert_records$UID[invert_records$Genus == "Paranephrops"])
deleatidiumUIDs <- unique(invert_records$UID[invert_records$Genus == "Deleatidium"])
paratyaUIDs <- unique(invert_records$UID[invert_records$Genus == "Paratya"])

orgDF$Has_Koura <- orgDF$UID %in% as.character(koaraUIDs)
orgDF$Has_Deleatidium <- orgDF$UID %in% as.character(deleatidiumUIDs)
orgDF$Has_Paratya <- orgDF$UID %in% as.character(paratyaUIDs)

# Site ID list per species
troutsites <- orgDF$UID[orgDF$Has_Trout]
kourasites <- orgDF$UID[orgDF$Has_Koura]
deleatidiumsites <- orgDF$UID[orgDF$Has_Deleatidium]
paratyasites <- orgDF$UID[orgDF$Has_Paratya]

# Trout vs Koura
trout_v_koura <- ggVennDiagram(
  list(Trout = troutsites, "Koura" = kourasites), label = "count"
) +
  scale_fill_gradient(low = "white", high = "white") +
  labs(title = "Trout vs Koura Site Overlap") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

trout_v_koura
ggsave("trout_v_koura_venn.png", plot = trout_v_koura, width = 6, height = 6, dpi = 300, bg = "white")

# Trout vs Deleatidium
trout_v_deleatidium <- ggVennDiagram(
  list(Trout = troutsites, "Deleatidium" = deleatidiumsites), label = "count"
) +
  scale_fill_gradient(low = "white", high = "white") +
  labs(title = "Trout vs Deleatidium Mayfly Site Overlap") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

trout_v_deleatidium
ggsave("trout_v_deleatidium_venn.png", plot = trout_v_deleatidium, width = 6, height = 6, dpi = 300, bg = "white")

# Trout vs Paratya
trout_v_paratya <- ggVennDiagram(
  list(Trout = troutsites, "Paratya" = paratyasites), label = "count"
) +
  scale_fill_gradient(low = "white", high = "white") +
  labs(title = "Trout vs Paratya Shrimp Site Overlap") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

trout_v_paratya
ggsave("trout_v_paratya_venn.png", plot = trout_v_paratya, width = 6, height = 6, dpi = 300, bg = "white")

## Diversity calculations ##
# Two community matrices, used for different things below:
#  - commMatrixFull: every detected taxon, including terrestrial/marine
#    incidentals. Useful as a QC signal (e.g. sequencing depth) but not a
#    meaningful measure of *stream* diversity, since it's dominated by
#    whatever eDNA happened to wash through.
#  - commMatrixAquatic: restricted to the core native aquatic macroinvertebrate
#    taxa defined above. This is what we use for Shannon diversity, since
#    that's the ecologically meaningful "is trout affecting the stream
#    invertebrate community" question you're asking - including hundreds of
#    aphid/slug/moth incidentals in that number would just add noise
#    unrelated to trout's effect on the stream.
nonSpeciesCols <- c("UID", "Has_Trout", "Has_Koura", "Has_Deleatidium", "Has_Paratya")
speciesCols <- setdiff(colnames(orgDF), nonSpeciesCols)
commMatrixFull <- orgDF[, speciesCols]
commMatrixAquatic <- orgDF[, intersect(speciesCols, aquaticTaxonKeys)]

# Shannon Diversity Index per site (core aquatic taxa only, see note above)
orgDF$Shannon <- diversity(commMatrixAquatic, index = "shannon")

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
ggsave("allsites_diversity_invert.png", plot = allsites_diversity, width = 6, height = 5, dpi = 300)

# Diversity of Koura sites w and without trout (uncorrected)
kourasites_diversity <- plotDiversity(
  orgDF[orgDF$Has_Koura, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Koura Sites",
  c("grey15", "lightblue4")
)

kourasites_diversity
ggsave("kourasites_diversity.png", plot = kourasites_diversity, width = 6, height = 5, dpi = 300)

# Diversity of Deleatidium sites w and without trout (uncorrected)
deleatidiumsites_diversity <- plotDiversity(
  orgDF[orgDF$Has_Deleatidium, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Deleatidium Sites",
  c("grey15", "skyblue")
)

deleatidiumsites_diversity
ggsave("deleatidiumsites_diversity.png", plot = deleatidiumsites_diversity, width = 6, height = 5, dpi = 300)

# Diversity of Paratya sites w and without trout (uncorrected)
paratyasites_diversity <- plotDiversity(
  orgDF[orgDF$Has_Paratya, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Paratya Sites",
  c("grey15", "forestgreen")
)

paratyasites_diversity
ggsave("paratyasites_diversity.png", plot = paratyasites_diversity, width = 6, height = 5, dpi = 300)

## Check for variation in eDNA sampling ##
# Add SeqDepth to df (full matrix here, since this is a QC check on total
# sequencing effort at the site, not on the aquatic-only subset)
orgDF$SeqDepth <- rowSums(commMatrixFull)

# Variation in seq depth
summary(orgDF$SeqDepth)
seqdepth_var <- ggplot(orgDF, aes(x = SeqDepth)) +
  geom_histogram(bins = 50, fill = "grey") +
  labs(title = "eDNA sequencing depth across sites",
       x = "Number of sequences", y = "Number of sites") +
  theme_minimal()

seqdepth_var
ggsave("seqdepth_var_invert.png", plot = seqdepth_var, width = 6, height = 5, dpi = 300)

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
ggsave("corr_depth_shannon_invert.png", plot = corr_depth_shannon, width = 6, height = 5, dpi = 300)

## Rarefy (normalisation transformation) community matrix ##
## Recompute Shannon with common seq depth ##
# Rarefy on the aquatic-only matrix, since that's the matrix Shannon above
# is calculated from - rarefying commMatrixFull instead would normalise for
# depth of *all* eDNA (including incidentals), not depth of aquatic taxa.
orgDF$SeqDepthAquatic <- rowSums(commMatrixAquatic)

# Remove sites with low sequencing depth
minDepth <- round(quantile(orgDF$SeqDepthAquatic, 0.30, na.rm = TRUE))
keep <- !is.na(orgDF$SeqDepthAquatic) & orgDF$SeqDepthAquatic >= minDepth
cat(sprintf(
  "Rarefying to %d sequences; dropping %d of %d sites below this depth\n",
  minDepth, sum(!keep), nrow(orgDF)
))

# Rarefaction curve check
set.seed(1)
rarecurve(commMatrixAquatic[sample(nrow(commMatrixAquatic), 30), ], step = 50, label = FALSE)
abline(v = minDepth, col = "firebrick", lty = 2)

commMatrixRarefy <- commMatrixAquatic[keep, ]
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
ggsave("allsites_diversity_rarefied_invert.png", plot = allsites_diversity_rarefied, width = 6, height = 5, dpi = 300)

# Diversity of Koura sites w and without trout (rarefied)
kourasites_diversity_rarefied <- plotDiversityRarefied(
  orgDFRarefy[orgDFRarefy$Has_Koura, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Koura Sites (with correction)",
  c("grey15", "lightblue4")
)

kourasites_diversity_rarefied
ggsave("kourasites_diversity_rarefied.png", plot = kourasites_diversity_rarefied, width = 6, height = 5, dpi = 300)

# Diversity of Deleatidium sites w and without trout (rarefied)
deleatidiumsites_diversity_rarefied <- plotDiversityRarefied(
  orgDFRarefy[orgDFRarefy$Has_Deleatidium, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Deleatidium Sites (with correction)",
  c("grey15", "skyblue")
)

deleatidiumsites_diversity_rarefied
ggsave("deleatidiumsites_diversity_rarefied.png", plot = deleatidiumsites_diversity_rarefied, width = 6, height = 5, dpi = 300)

# Diversity of Paratya sites w and without trout (rarefied)
paratyasites_diversity_rarefied <- plotDiversityRarefied(
  orgDFRarefy[orgDFRarefy$Has_Paratya, ], "Has_Trout", c("No Trout", "Trout"),
  "Shannon Diversity - Paratya Sites (with correction)",
  c("grey15", "forestgreen")
)

paratyasites_diversity_rarefied
ggsave("paratyasites_diversity_rarefied.png", plot = paratyasites_diversity_rarefied, width = 6, height = 5, dpi = 300)

## Histogram of # of species per site ##
# Number of species per site (full matrix - total detected taxa at a site,
# same "all taxa" definition the original fish script used)
orgDF$Richness <- rowSums(commMatrixFull > 0)

# Standardise axes
richnessXLim <- c(min(orgDF$Richness, na.rm = TRUE) - 1, max(orgDF$Richness, na.rm = TRUE) + 1)
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
ggsave("richness_hist_trout_invert.png", plot = richness_trout, width = 6, height = 5, dpi = 300)

richness_notrout <- plotSpeciesNo(
  orgDF[!orgDF$Has_Trout, ], "Species Richness - No Trout Sites", "grey15",
  xlim = richnessXLim, ymax = richnessYMax
)

richness_notrout
ggsave("richness_hist_notrout_invert.png", plot = richness_notrout, width = 5, height = 5, dpi = 300)

## Native invertebrate species-specific richness ##
# Native richness here uses the same aquaticTaxonKeys defined near the top
# of the script (see the note there on how "native" is being approximated)
orgDF$NativeRichness <- rowSums(commMatrixAquatic > 0)

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
ggsave("native_richness_hist_trout_invert.png", plot = nativerichness_trout, width = 6, height = 5, dpi = 300)

nativerichness_notrout <- plotNativeSpeciesNo(
  orgDF[!orgDF$Has_Trout, ], "Native Species Richness - No Trout Sites", "grey15",
  xlim = nativeXLim, ymax = nativeYMax
)

nativerichness_notrout
ggsave("nativerichness_hist_notrout_invert.png", plot = nativerichness_notrout, width = 6, height = 5, dpi = 300)

## Habitat covariates ##
# Pull habitat data
habitatData <- invert_samples[, c("UID", "EnvironmentType", "Region")]
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
ggsave("diversity_stratified_by_habitat_invert.png", plot = plotDiversityStratified, width = 8, height = 6, dpi = 300)

## Upstream vs downstream of waterfalls effect ##
## Using a small subset of the data that has manually assigned location ##
## (site list unchanged from the fish script - these are the same physical
## sites, just joined here to the invertebrate community data instead) ##
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
      title = "Shannon Diversity - All Sites, Above vs Below Waterfall",
      x = NULL, y = "Shannon Diversity Index"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  waterfallPlot
  ggsave("diversity_waterfall_barrier_invert.png", plot = waterfallPlot, width = 6, height = 5, dpi = 300)
  
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
