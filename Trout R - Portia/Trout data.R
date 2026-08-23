#Load packages
library(pacman)
p_load("tidyverse", "data.table")

#Load data
records <- read.csv("~/GENE 412 Trout/records.csv") 
samples <- read.csv("~/GENE 412 Trout/samples.csv")

#Set as data table
records <- data.table(records)
samples <- data.table(samples)

#Remove first column (HID)
records <- records[,-1]


#Table with count of species (eDNA) present per UID
UIDdiversity <- table(records$UID) # create a count table of UID
UIDdiversity <- data.table(UIDdiversity) #change to a data frame
colnames(UIDdiversity) <- c("UID", "Species Count") #rename columns
UIDdiversity <- cbind(UIDdiversity, samples[,c(3, 4)]) #include lattitude and longitude



#Maps
p_load(leaflet, sf, rnaturalearth, rnaturalearthdata)
leaflet(data = UIDdiversity)

nz_map <- ne_countries(scale = "medium", country = "New Zealand", returnclass = "sf")


#Heat map over area, not really relevant
{
#ggplot() +
#  # Base Map Layer
#  geom_sf(data = nz_map, fill = "gray95", color = "gray60") +
#  # Density/Heat Map Layer
#  stat_density_2d(data = UIDdiversity, aes(x = Longitude, y = Latitude, fill = after_stat(level)), 
#                  geom = "polygon", alpha = 0.4) +
#  # Color palette for the heat map
#  scale_fill_gradientn(colors = c("blue", "green", "yellow", "red"), name = "Density")  +
#  # Optional: Overlay individual data points
#  geom_point(data = UIDdiversity, aes(x = Longitude, y = Latitude), size = 0.8, color = "maroon1", alpha = 0.3) +
#  # Zoom the map view tightly around your data points
#  coord_sf(xlim = c(165, 180), ylim = c(-34, -48), expand = FALSE) +
#  labs(title = "NZ Coordinate Heat Map Displaying Count of Fish Density", x = "Longitude", y = "Latitude") +
#  theme_minimal()
}


#Density per location
ggplot() +
  # Base Map Layer
  geom_sf(data = nz_map, fill = "gray95", color = "gray60") +
  # Maps the point color directly to your count variable
  geom_point(data = UIDdiversity, aes(x = Longitude, y = Latitude, color = `Species Count`), 
             size = 2, alpha = 0.7) +
  # Color palette for the heat map
  scale_color_gradientn(colors = c("blue", "yellow", "red"), name = UIDdiversity$`Species Count`)  +
  # Zoom the map view tightly around your data points
  coord_sf(xlim = c(165, 180), ylim = c(-34, -48), expand = FALSE) +
  labs(title = "NZ Coordinate Heat Map Displaying Count of Fish Density", x = "Longitude", y = "Latitude") +
  theme_minimal()





#Longitude 
longitude <- UIDdiversity$Longitude
longitude <- data.table(longitude)

#Latitude
latitude <- UIDdiversity$Latitude
latitude <- data.table(latitude)

#Number of times species appear
species <- table(records$Name)
species <- data.table(species)
colnames(species) <- c("Name", "Species Count") #rename columns

#Oncorhynchus mykiss = Rainbow Trout
#Salmo trutta = Brown Trout


#tibble <- records[,c("UID", "Name", "CommonName", "Count")]



#need something like if uid and species match then add count and combine CommonName

tibble <- data.frame(tibble)
tibble <- data.table(tibble)


browntrout <- records[records$Name == "Salmo trutta", ]
rainbowtrout <- records[records$Name == "Oncorhynchus mykiss", ]

combinedtrout <- rbind(browntrout, rainbowtrout)

records$presence <- ifelse(records$Count <10, 0, 1)


#want to have a look at how many species are present in trout vs non trout water bodies


Marine_filter = c("Blue.cod","Yellowtail.kingfish","Bluefin.Gurnard","Trevally","Tarakihi","Banded.Parrotfish","Blue.warehou","Swordfish","Red.cod","Orange.roughy","Indian.mackerel","Mahi.Mahi","Sunfish","Yelloweye.mullet")
