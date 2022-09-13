# Mapping species ranges and overlap

#### Setup ####
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(tidyr)
library(dplyr)
library(robis)
library(leaflet)
library(spocc)
library(rgbif)
library(scrubr)
library(qlcMatrix)
library(CoordinateCleaner)
library(terra)
library(geodata)
library(sdmpredictors)
library(raster)


#### Occurrence data ####

#my_species <- c("Penaeus vannamei", "Ruditapes philippinarum")   #test

my_species <- c("Penaeus vannamei", "Ruditapes philippinarum", "Salmo salar", "Chanos chanos", "Oncorhynchus mykiss", "Penaeus monodon", "Magallana gigas", "Anadara granosa", "Mytilus chilensis", "Apostichopus japonicus")


spocc_data <- occ(my_species, from = c("gbif", "obis"), limit = 10000, has_coords = TRUE)
spocc_data

names(spocc_data$gbif$data[[1]])

# columns we actually want:
interest <- c('species', 'longitude', 'latitude', 'prov', 'eventDate', 'occurrenceStatus', 'countryCode', 'basisOfRecord')
final <- data.frame(matrix(ncol = length(interest), nrow = 0))
colnames(final) <- interest
# limit to databases of interest
db <- c(spocc_data[['gbif']], spocc_data[['obis']])
db <- db[c(2,4)]

# combines data into one data frame with actual non-ambiguous species name
for (j in db) {
  for (i in j) {
    
    df <- as.data.frame(i)
    df <- df[,names(df) %in% interest]
    final <- rbind(df, final)
  
  }
}
spocc_df <- final
summary(spocc_df)

# spocc_df <- occ2df(spocc_data, what = 'all')
# spocc_df

# names(occ)
# spocc <- occ[,c("scientificName","sst","date_year","decimalLatitude", "decimalLongitude")]
# spocc <- spocc[ which(spocc$date_year >= 1995 & spocc$date_year <= 2020), ]

# leaflet() %>%
#   addTiles() %>%
#   addCircles(data = subset(spocc_df, name == "Penaeus vannamei"), lng = ~ longitude, lat = ~ latitude, col = "blue") %>%
#   addCircles(data = subset(spocc_df, name == "Ruditapes philippinarum"), lng = ~ longitude, lat = ~ latitude, col = "red") %>%
#   addCircles(data = subset(spocc_df, name == "Salmo salar"), lng = ~ longitude, lat = ~ latitude, col = "green") %>%
#   addCircles(data = subset(spocc_df, name == "Chanos chanos"), lng = ~ longitude, lat = ~ latitude, col = "yellow") %>%
#   addCircles(data = subset(spocc_df, name == "Oncorhynchus mykiss"), lng = ~ longitude, lat = ~ latitude, col = "pink") %>%
#   addCircles(data = subset(spocc_df, name == "Penaeus monodon"), lng = ~ longitude, lat = ~ latitude, col = "orange") %>%
#   addCircles(data = subset(spocc_df, name == "Magallana gigas"), lng = ~ longitude, lat = ~ latitude, col = "purple") %>%
#   addCircles(data = subset(spocc_df, name == "Anadara granosa"), lng = ~ longitude, lat = ~ latitude, col = "brown") %>%
#   addCircles(data = subset(spocc_df, name == "Mytilus chilensis"), lng = ~ longitude, lat = ~ latitude, col = "black") %>%
#   addCircles(data = subset(spocc_df, name == "Apostichopus japonicus"), lng = ~ longitude, lat = ~ latitude, col = "grey")
# 

# Map occurrence data to ensure every species is there ----
# USE THIS FOR FUTURE PLOTS (colors automatically)

# bright magenta = NA
pal <- colorFactor(
  palette = "Set3",
  na.color = 'magenta',
  domain = my_species
)

leaflet() %>%
  addTiles() %>%
    addCircles(data = spocc_df, lng = ~ longitude, lat = ~ latitude, col = ~pal(species), popup = paste(spocc_df$species)) %>%
  addLegend(position = 'bottomleft', pal = pal, values = my_species, title = 'Species', na.label = 'N/A')

  

#### Cleaning data ####

#spocc <- occ[ which(spocc$date_year >= 1995 & spocc$date_year <= 2020), ]
nrow(spocc_df)

# remove other records
unique(spocc_df$basisOfRecord)
occ_clean <- subset(spocc_df, !(basisOfRecord %in% c("PreservedSpecimen", "D", "DerivedFromLiterature", "MaterialSample", "NomenclaturalChecklist" , "FossilSpecimen", "MATERIAL_SAMPLE", "PRESERVED_SPECIMEN", "MATERIAL_CITATION", "FOSSIL_SPECIMEN"))) 

# keep only presence data
occ_clean <- subset(spocc_df, !(occurrenceStatus %in% c("ABSENT", NA))) 

# remove NA coordinates
occ_clean <- subset(occ_clean, !is.na(latitude) & !is.na(longitude))
nrow(occ_clean)

# remove duplicates
dups <- which(duplicated(occ_clean[ , c("longitude", "latitude")]))
length(dups)
if (length(dups) > 0)  occ_clean <- occ_clean[-dups, ]
nrow(occ_clean)

# remove lat and long that are zero
occ_clean <- subset(occ_clean, !(latitude == 0 & longitude == 0))
nrow(occ_clean)

# do some automatic data cleaning with functions of the 'scrubr' package:
occ_clean <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(occ_clean))))
nrow(occ_clean)

pal <- colorFactor(
  palette = "Set3",
  na.color = 'magenta',
  domain = my_species
)

leaflet() %>%
  addTiles() %>%
  addCircles(data = occ_clean, lng = ~ longitude, lat = ~ latitude, col = ~pal(species), popup = paste(occ_clean$species)) %>%
  addLegend(position = 'bottomleft', pal = pal, values = my_species, title = 'Species', na.label = 'N/A')

# more automatic data cleaning
occ_coordclean <- clean_coordinates(occ_clean, lon = "longitude", lat = "latitude", countries = "countryCode", tests = c("centroids", "equal", "institutions", "outliers", "zeros"), inst_rad = 100, zeros_rad = 1, outliers_method = "quantile")

head(occ_coordclean)

leaflet() %>%
  addTiles() %>%
  addCircles(data = occ_coordclean, lng = ~ longitude, lat = ~ latitude, col = ~pal(species), popup = paste(occ_coordclean$species)) %>%
  addLegend(position = 'bottomleft', pal = pal, values = my_species, title = 'Species', na.label = 'N/A')

nrow(occ_clean)
occ_clean <- occ_clean[occ_coordclean$.summary == TRUE, ]
nrow(occ_clean)

# remove coordinates outside region of interest
library(rgdal)
library(raster)
library(ggplot2)
library(rgeos)

download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip", 
              destfile = 'oceans.zip')

unzip(zipfile = "oceans.zip", 
      exdir = 'ne-ocean-10m')
oceans <- readOGR("ne-ocean-10m/ne_10m_ocean.shp")
class(oceans)
crs(oceans)

#Map again
pal <- colorFactor(
  palette = "Set3",
  na.color = 'magenta',
  domain = my_species
)

leaflet() %>%
  addTiles() %>%
  addCircles(data = occ_clean, lng = ~ longitude, lat = ~ latitude, col = ~pal(species), popup = paste(occ_clean$species)) %>%
  addLegend(position = 'bottomleft', pal = pal, values = my_species, title = 'Species', na.label = 'N/A')

## Next Steps:
# clean data - remove inland points
# map according to country or separate species to minimise overlap - random so no overlap or by classification?
# add legend 
# do the same for algae?

