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
install.packages("qlcMatrix")
library(qlcMatrix)
library(CoordinateCleaner)
library(terra)
library(geodata)

#### Occurrence data ####

#my_species <- c("Penaeus vannamei", "Ruditapes philippinarum")   #test

my_species <- c("Penaeus vannamei", "Ruditapes philippinarum", "Salmo salar", "Chanos chanos", "Oncorhynchus mykiss", "Penaeus monodon", "Magallana gigas", "Anadara granosa", "Mytilus chilensis", "Apostichopus japonicus")


spocc_data <- occ(my_species, from = c("gbif", "obis"), has_coords = TRUE)
spocc_data

#combines data into one data frame
spocc_df <- occ2df(spocc_data)
spocc_df

#names(occ)
#spocc <- occ[,c("scientificName","sst","date_year","decimalLatitude", "decimalLongitude")]
#spocc <- spocc[ which(spocc$date_year >= 1995 & spocc$date_year <= 2020), ]

leaflet() %>%
  addTiles() %>%
  addCircles(data = subset(spocc_df, name == "Penaeus vannamei"), lng = ~ longitude, lat = ~ latitude, col = "blue") %>%
  addCircles(data = subset(spocc_df, name == "Ruditapes philippinarum"), lng = ~ longitude, lat = ~ latitude, col = "red") %>%
  addCircles(data = subset(spocc_df, name == "Salmo salar"), lng = ~ longitude, lat = ~ latitude, col = "green") %>%
  addCircles(data = subset(spocc_df, name == "Chanos chanos"), lng = ~ longitude, lat = ~ latitude, col = "yellow") %>%
  addCircles(data = subset(spocc_df, name == "Oncorhynchus mykiss"), lng = ~ longitude, lat = ~ latitude, col = "pink") %>%
  addCircles(data = subset(spocc_df, name == "Penaeus monodon"), lng = ~ longitude, lat = ~ latitude, col = "orange") %>%
  addCircles(data = subset(spocc_df, name == "Magallana gigas"), lng = ~ longitude, lat = ~ latitude, col = "purple") %>%
  addCircles(data = subset(spocc_df, name == "Anadara granosa"), lng = ~ longitude, lat = ~ latitude, col = "brown") %>%
  addCircles(data = subset(spocc_df, name == "Mytilus chilensis"), lng = ~ longitude, lat = ~ latitude, col = "black") %>%
  addCircles(data = subset(spocc_df, name == "Apostichopus japonicus"), lng = ~ longitude, lat = ~ latitude, col = "grey")
  
# Map occurrence data to ensure every species is there

leaflet() %>%
  addTiles() %>%
  addCircles(data = subset(spocc, scientificName == "Penaeus vannamei"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "blue") %>%
  addCircles(data = subset(spocc, scientificName == "Ruditapes philippinarum"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "red") %>%
  addCircles(data = subset(spocc, scientificName == "Salmo salar"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "green") %>%
  addCircles(data = subset(spocc, scientificName == "Chanos chanos"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "yellow") %>%
  addCircles(data = subset(spocc, scientificName == "Oncorhynchus mykiss"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "pink") %>%
  addCircles(data = subset(spocc, scientificName == "Penaeus monodon"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "orange") %>%
  addCircles(data = subset(spocc, scientificName == "Magallana gigas"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "purple") %>%
  addCircles(data = subset(spocc, scientificName == "Anadara granosa"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "brown") %>%
  addCircles(data = subset(spocc, scientificName == "Mytilus chilensis"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "black") %>%
  addCircles(data = subset(spocc, scientificName == "Apostichopus japonicus"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "grey")

#### Cleaning data ####

spocc <- occ[ which(spocc$date_year >= 1995 & spocc$date_year <= 2020), ]

names(spocc)
unique(spocc$occurrenceStatus) 

# keep only presence data
occ_clean <- subset(spocc, !(occurrenceStatus %in% NA)) 

# remove NA coordinates
occ_clean <- subset(occ_clean, !is.na(decimalLatitude) & !is.na(decimalLongitude))
nrow(occ_clean)

# remove duplicates
dups <- which(duplicated(occ_clean[ , c("decimalLongitude", "decimalLatitude")]))
length(dups)
if (length(dups) > 0)  occ_clean <- occ_clean[-dups, ]
nrow(occ_clean)

# remove lat and long that are zero
occ_clean <- subset(occ_clean, !(decimalLatitude == 0 & decimalLongitude == 0))
nrow(occ_clean)

# do some automatic data cleaning with functions of the 'scrubr' package:

?coord_incomplete

occ_clean <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(occ_clean))))
nrow(occ_clean)

#remove 
unique(spocc$basisOfRecord) 
occ_clean2 <- subset(occ_clean, !(basisOfRecord %in% c("PreservedSpecimen", "DerivedFromLiterature", "MaterialSample", "D", "FossilSpecimen", "NomenclaturalChecklist")))

leaflet() %>%
  addTiles() %>%
  addCircles(data = subset(occ_clean2, scientificName == "Penaeus vannamei"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "blue") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Ruditapes philippinarum"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "red") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Salmo salar"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "green") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Chanos chanos"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "yellow") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Oncorhynchus mykiss"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "pink") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Penaeus monodon"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "orange") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Magallana gigas"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "purple") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Anadara granosa"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "brown") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Mytilus chilensis"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "black") %>%
  addCircles(data = subset(occ_clean2, scientificName == "Apostichopus japonicus"), lng = ~ decimalLongitude, lat = ~ decimalLatitude, col = "grey")

names(occ_clean)

## Next Steps:
# clean data - remove inland and repeat occurrences
# clean data for each species rather than as a group??????
# add gbif data?
# map according to country or separate species to minimise overlap - random so no overlap or by classification?
# add legend 
# do the same for algae?
# fuzzy occurrences per species - too many figures?

