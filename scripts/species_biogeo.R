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

# will want to rerun with a full my_species, remove code referencing seaweed specifically
# my_species <- c("Saccharina japonica", "Saccharina latissima")   #test

my_species <- c("Salmo salar", "Oncorhynchus mykiss", "Chanos chanos", "Penaeus vannamei","Penaeus monodon", "Ruditapes philippinarum",     "Magallana gigas", "Saccharina japonica", "Undaria pinnatifida", "Kappaphycus alvarezii", "Porphyra tenera") # use Saccharina japonica next time


spocc_data <- occ(my_species, from = c("gbif", "obis"), limit = 10000, has_coords = TRUE)
spocc_data

tst_jap <- occ("Laminaria japonica", from = c("gbif", "obis"), limit = 10000, has_coords = TRUE)

names(spocc_data$gbif$data[[1]])

saveRDS(spocc_data, paste0("../output niches/spocc_", Sys.Date(), "my_species_raw.rds"))  # we name it "raw" because the data haven't been cleaned yet

# spocc_data_restricted <- readRDS(paste0("./outputs/spocc_", my_species, "_raw.rds"))


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
occ_clean <- subset(occ_clean, !(occurrenceStatus %in% c("ABSENT", NA))) 

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

nrow(occ_clean)
occ_clean <- occ_clean[occ_coordclean$.summary == TRUE, ]
nrow(occ_clean)

# remove coordinates outside region of interest
library(rgdal)
library(raster)
library(ggplot2)
library(rgeos)
library(sf)
sf_use_s2(FALSE)


# download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip", 
#               destfile = 'oceans.zip')
# 
# unzip(zipfile = "oceans.zip", 
#       exdir = 'ne-ocean-10m')
oceans <- readOGR("ne-ocean-10m/ne_10m_ocean.shp")
class(oceans)
crs(oceans)

oceans <- st_as_sf(oceans)
oceans_clean <- oceans[1,"min_zoom"]
points <- subset(occ_clean, select = c("latitude", "longitude", "species"))
points <- st_as_sf(points, coords = c("longitude", "latitude"), crs = 4326)
# 1 = in ocean, 0 = on land
points_plot <- points %>% mutate(in_ocean = lengths(st_within(points, oceans_clean)))


st_write(points_plot, paste0('points_plot',strftime(Sys.time(), "%m-%d-%y_%H%M"),'.shp'))

# STOP HERE AND APPEND #
algae <- c("Laminaria japonica", "Undaria pinnatifida", "Kappaphycus alvarezii", "Porphyra tenera")
noAlgae <- c("Salmo salar", "Oncorhynchus mykiss", "Chanos chanos", "Penaeus vannamei","Penaeus monodon",  "Ruditapes philippinarum",     "Magallana gigas")
points_ocean <- points_plot[points_plot$in_ocean==1,]
points_ocean$type <-NA
points_ocean$type[points_ocean$species %in% algae] <- "algae"
points_ocean$type[!(points_ocean$species %in% algae)] <- "other"

# PLOTTING ----
library(viridis)
levels(my_species) <- my_species
my_species <- factor(my_species)
plot <- st_transform(points_ocean,4326)
plot_other <- plot[plot$type=="other",]
plot_algae <- plot[plot$type=="algae",]

pal <- colorFactor(
  palette = viridis(11),
  na.color = 'magenta',
  levels = my_species,
  ordered = TRUE
)
pal_my <- colorFactor(
  palette = my_colors,
  na.color = 'magenta',
  domain = my_species
)
# pal_other <- colorFactor(
#   palette = heat.colors(length(noAlgae)),
#   na.color = 'magenta',
#   domain = noAlgae,
#   
# )
# pal_algae <- colorFactor(
#   palette = viridis(length(algae)),
#   na.color = 'magenta',
#   domain = algae
# )

leaflet(plot) %>%
  addTiles() %>%
  addCircleMarkers(col = ~pal(species), stroke = FALSE, radius = 2, fillOpacity = 1, popup = paste(plot$species)) %>%
  addLegend(position = 'topright', pal = pal, values = my_species, title = 'Species', na.label = 'N/A', labels = my_species)

# try with ggplot
# goal is getting one plot grouped by species type in legend w/colors
library(RColorBrewer)
world_coordinates <- map_data("world")
plot$species <- factor(plot$species, levels = my_species)
plot_other$species <- factor(plot_other$species, levels = noAlgae)
plot_algae$species <- factor(plot_algae$species, levels = algae)

my_colors <- c(inferno(length(noAlgae), begin = 0.35, end = .9, direction = -1), colors()[c(614,518,610,494)])
# brewer.pal(length(noAlgae), name = 'YlOrRd')

library(sdmpredictors)
library(raster)
library(rgdal)

# make raster plot from BioORACLE layers
SST_2050_8.5 <- read.asciigrid(fname = file_tst)
SST_2050_8.5 <- load_layers("BO22_RCP85_2050_tempmean_ss")
SST_2100_8.5 <- load_layers("BO22_RCP85_2100_tempmean_ss")
SST_max2050_8.5 <- load_layers("BO22_RCP85_2050_tempmax_ss")
SST_max2100_8.5 <- load_layers("BO22_RCP85_2100_tempmax_ss")
SST_present <- load_layers("BO2_tempmean_ss")
SST_maxpresent <- load_layers("BO2_tempmax_ss")
SST_diff <- SST_2100_8.5 - SST_present
SST_maxdiff <- SST_max2100_8.5 - SST_maxpresent
diff_df <- as.data.frame(SST_diff, xy=TRUE)
diff_maxdf <- as.data.frame(SST_maxdiff, xy=TRUE)
names(diff_df) <- c('x', 'y', 'temp. diff (mean)')
names(diff_maxdf) <- c('x', 'y', 'temp. diff (max)')

# first pass plotting
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
plot(SST_present, col = my.colors(1000), axes = FALSE, box = FALSE)
plot(SST_2050_8.5, col = my.colors(1000), axes = FALSE, box = FALSE)
plot(SST_diff, col = my.colors(1000), axes = FALSE, box = FALSE)

# outline black, bit smaller, like 10-20% transparent, fill as is
# plotting occurrences
ggplot() +
  geom_map(
    data = world_coordinates, map = world_coordinates,
    aes(long, lat, map_id = region),  color = '#616161', fill = 'white', alpha = .5)  +
  # geom_raster(data = diff_df, mapping=aes(x=x, y=y, fill=`temp. diff (mean)`), alpha = .5)+
  # # scale_fill_distiller(palette = 'Spectral',name='temp diff (deg C)', limits = c(-2,6)) +
  geom_sf(data = plot, aes(fill = species), shape = 21, size =2.25, colour = '#9E9E9E', alpha = .8, inherit.aes = FALSE) + 
            scale_fill_manual(values = my_colors, labels = my_species) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
# raster mean temp with labels
ggplot() +
  geom_raster(data = diff_df, mapping=aes(x=x, y=y, fill=`temp. diff (mean)`), alpha = .5)+
  scale_fill_distiller(palette = 'Spectral',name='temp diff (deg C)', limits = c(-2,6)) +
  geom_point(data = occ_clean, aes(x=longitude, y=latitude, colour=species)) +
  scale_colour_manual(values = my_colors)

# raster mean temp without labels or background
ggplot() +
  geom_raster(data = diff_df, mapping=aes(x=x, y=y, fill=`temp. diff (mean)`), alpha = 1)+
  scale_fill_distiller(palette = 'Spectral',name='temp diff (deg C)', limits = c(-.5,7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# raster max temp
ggplot() +
  geom_raster(data = diff_maxdf, mapping=aes(x=x, y=y, fill=`temp. diff (max)`), alpha = 1)+
  scale_fill_distiller(palette = 'Spectral',name='max temp diff (deg C)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## Next Steps:
# map according to country or separate species to minimise overlap - random so no overlap or by classification?
# generalize approach to deal with synonymous names - check unique species matches length of my_species



