# Wrangling relevant temperature data for marine species of interest ----

## Setup ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
library(tidyr)
library(dplyr)
library(robis)

checks <- read.delim("../input data/OBIS test marine species data.txt")
# Thermal guilds from overlapping temp. ranges ---------------------------------
dat <- read.delim("Marine species data.txt")
df <- dat[,! colnames(dat) %in% c("Common.name","Farming.Environment","Temp.Ref")]
# df <- na.omit(df)
overlap <- function(A, B) {
  shared <- pmax(0, min(A[2], B[2]) - max(A[1], B[1]))
  max(shared / c(diff(A), diff(B)))
}

eg <- expand.grid(a = seq_len(nrow(df)), b = seq_len(nrow(df)))
eg <- eg[eg$a < eg$b,]

together <- cbind(
  setNames(df[eg$a,], paste0(names(df), "1")),
  setNames(df[eg$b,], paste0(names(df), "2"))
)
together <- within(together, {
  shared  = pmax(0, pmin(Max..Temp1, Max..Temp2) - pmax(Min..Temp1, Min..Temp2))
  overlap = pmax(shared / (Max..Temp1 - Min..Temp1), shared / (Max..Temp2 - Min..Temp2))
})[, c("Species1", "Species2", "overlap")]

bigenough <- together[together$overlap >= 0.75,]
groups <- split(bigenough$Species2, bigenough$Species1)

for (ltr in df$Species) {
  ind <- (ltr == names(groups)) | sapply(groups, `%in%`, x = ltr)
  groups <- c(
    setNames(list(unique(c(ltr, names(groups[ind]), unlist(groups[ind])))), ltr),
    groups[!ind]
  )
}

groups <- data.frame(
  ID = rep(seq_along(groups), lengths(groups)),
  Species = unlist(groups)
)

guilds <- merge(df, groups, by = "Species")

isUnique <- function(vector){
  return(!any(duplicated(vector)))
}

guild_groups <- guilds %>% group_by(ID) %>% summarise(
  group_mean = mean(Mean.Temp),
  group_min = mean(Min..Temp),
  group_max = mean(Max..Temp)
)

# Next steps
# by classification
# directionality of overlap?
# stats by group
 
 
# Thermal niches from obis -----------------------------------------------------

## testing prior to implementation
#tst <- occurrence("Tilapia guineensis")
#tst_sst <- tst[which(tst$date_year >= 1995 & tst$date_year <= 2020 & tst$marine == TRUE),
#                      c("scientificName", "sst")]
#tst_sum <- as.data.frame(summary(tst_sst))
#species = "Tilapia guineensis"
#target <- tst_sum[8:13,"Freq"]
# row <- c(species)
# as.numeric(sub(".*:","",target))
# for (stat in target) {
#   add <- as.numeric(sub(".*:","",stat))
#   row <- append(row, add)
# }

# function to make thermal niches from a data frame with one column called "Species"
# Returns "final_df" and writes to file as .csv
df <- head(df)
dat <- df

make_niches <- function(df, n = 20) {
  ## prepare empty final data frame
  ## establish columns
  file_name <- strftime(Sys.time(), "%m-%d-%y_%H%M")
  log <- file(paste0("../output niches/logs/",file_name,".log"), open = 'a')
  cat(file_name, " run messages","\n", file = log)
  columns <- c("Species", "Min", "1stQ", "Median", "Mean", "3rdQ", "Max")
  final_df <- data.frame(matrix(ncol = length(columns)))
  err <- c()
  colnames(final_df) <- columns
  prog <- 0
  db <- df$Species
  ## For-loop processing temperature summary
  
  for (species in db) { # outer for loop: obtain occurrence, summarize, extract stats and add to final
    tryCatch({ # catch errors without breaking the loop and print the error and species that caused it
      if (species %in% final_df$Species) { # skip duplicates
        cat(paste0("\n","Duplicate of ", species, " found and omitted from final dataset"), 
            file = log, sep="\n")
        next
      } else if (is.na(species) | species=="") {
        next
      }
      cat("\n", species, "\n")
      occ <- occurrence(species)
      # filter to date range matching BIO-ORACLE and only presence records
      occ_sum <- occ[which(occ$date_year >= 1995 & occ$date_year <= 2020 & occ$absence == FALSE & occ$marine == TRUE),
                     c("scientificName", "sst")]
      if (nrow(occ_sum) < n) { # this should skip the species if it has less than n occurrences
        cat(paste("\n", "Number of occurences less than threshold of", n, ",", species, "omitted from final dataset."), 
            file = log, sep="\n")
        next
      }
      occ_sum <- as.data.frame(summary(occ_sum))
      row <- c(species)
      target <- occ_sum[8:13,"Freq"] # prepare mini dataframe with predictably placed summary stats
      for (stat in target) { # inner for-loop: take summary stats and extract only the values
        add <- as.numeric(sub(".*:","",stat))
        row <- append(row, add)
      }
      if (any(is.na(row))) {
        cat(paste("\n", species, "summary contained NAs and was omitted from the final dataset"), 
            file = log, sep="\n")
        next
      } else {
        final_df <- rbind(row, final_df)
      }
      prog <- prog + 1
      # print status
      cat(paste("\n", round((prog/length(db))*100),"% of species processed in final dataset"),
          file = log, sep="\n")
    }, error=function(e){cat("ERROR in:", species, "\n", conditionMessage(e), "\n", 
                             file = log, sep="\n")
      err <- c(err,species)})
  }
  write.csv(final_df, file = paste0("../output niches/",file_name,"_temp_niche.csv"))
  closeAllConnections()
}
# RUN ON CHECKS BEFORE PUSHING
system.time({make_niches(checks)})
make_niches(dat)


# Next steps: 
# wrap in function that takes first dataframe and the name of the
# species containing column and does the rest possibility of including options
# to query other databases - could adapt to work with terrestrial data as well
# or multi-database via spocc
# Refine error handling

# # make_niches <- function(df) {
# ## prepare empty final data frame
# final_df <- rbind(row, final_df)
# ## establish columns
# columns <- c("Species", "Min", "1stQ", "Median", "Mean", "3rdQ", "Max")
# final_df <- data.frame(matrix(ncol = length(columns)))
# colnames(final_df) <- columns
# prog <- 0
# db <- df$Species
# startTime <- Sys.time()
# ## For-loop processing temperature summary
# 
# for (species in db) { # outer for loop: obtain occurrence, summarize, extract stats and add to final
#   tryCatch({ # catch errors without breaking the loop and print the error and species that caused it
#     occ <- occurrence(species)
#     occ_sum <- occ[,c("scientificName", "sst")]
#     occ_sum <- as.data.frame(summary(occ_sum))
#     row <- c(species)
#     target <- occ_sum[8:13,"Freq"] # prepare mini dataframe with predictably placed summary stats
#     for (stat in target) { # inner for-loop: take summary stats and extract only the values
#       add <- as.numeric(sub(".*:","",stat))
#       row <- append(row, add)
#     }
#     final_df <- rbind(row, final_df)
#     prog <- prog + 1
#     # print status
#     cat(paste("\n", round((prog/length(db))*100),"% of species processed", 
#               "\n Species just processed: ", species, "\n"))
#   }, error=function(e){cat("ERROR in:", species, "\n", conditionMessage(e), "\n")})
# }
# endTime <- Sys.time()
# cat(paste0("Total run time: ", endTime - startTime))


# Visualization to identify guilds ----
library(ggplot2)
library(cowplot)
niches <- read.delim2("fish_temp_niche.csv", header = TRUE, sep = ",")
names <- c("Index","Species", "Min", "1stQ", "Median", "Mean", "3rdQ", "Max")
names(niches) <- names
niches <- niches[order(niches$Median),]
top <- head(niches, 20)


p <- ggplot(niches, aes(x = reorder(Species, as.numeric(Median)), y = as.numeric(Median))) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggtitle("Aquaculture Fish") +
  ylab("Median Temperature (C)") +
  xlab("Species") +
  geom_errorbar(aes(ymax = as.numeric(`3rdQ`), ymin = as.numeric(`1stQ`)),
                position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_point(mapping = aes(y = as.numeric(`Min`)), color = "blue") +
  geom_point(mapping = aes(y = as.numeric(`Max`)), color = "red") +
  coord_flip()
p


amy <- read.delim2("2022-07-15_temp_niche.csv", header = TRUE, sep = "")
amy_nums <- read.delim2("Temp niche data_Amy.txt", header = TRUE, sep = "", row.names = NULL)
amy_nums$Species <- paste(amy_nums$row.names, amy_nums$Species)
amy_nums[296,] <-c("Anguilla", "Anguilla bicolor bicolor",	25.72,	27.59,	28.34,	27.97,	28.61,	29.6)
amy_nums[202,] <- c("Holothuria","Holothuria (Metriatyla) scabra",	22.57,	26.29,	27.78,	27.43,	28.9,	30.17)

amy_hot <- amy_nums[as.numeric(amy_nums$Median)>=20,]
amy_cold <- amy_nums[as.numeric(amy_nums$Median)<20,]
a <- ggplot(amy_hot, aes(x = reorder(Species, as.numeric(Median)), y = as.numeric(Median))) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggtitle("All Aquaculture >= 20C median T") +
  ylab("Median Temperature (C)") +
  xlab("Species") +
  geom_errorbar(aes(ymax = as.numeric(`X3rdQ`), ymin = as.numeric(`X1stQ`)),
                position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_point(mapping = aes(y = as.numeric(`Min`)), color = "blue") +
  geom_point(mapping = aes(y = as.numeric(`Max`)), color = "red") +
  coord_flip()

b <- ggplot(amy_cold, aes(x = reorder(Species, as.numeric(Median)), y = as.numeric(Median))) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggtitle("All Aquaculture < 20C median T") +
  ylab("Median Temperature (C)") +
  xlab("Species") +
  geom_errorbar(aes(ymax = as.numeric(`X3rdQ`), ymin = as.numeric(`X1stQ`)),
                position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_point(mapping = aes(y = as.numeric(`Min`)), color = "blue") +
  geom_point(mapping = aes(y = as.numeric(`Max`)), color = "red") +
  coord_flip()
c <- ggplot(amy, aes(x = reorder(Species, as.numeric(Median)), y = as.numeric(Median))) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggtitle("All Aquaculture") +
  ylab("Median Temperature (C)") +
  xlab("Species") +
  geom_errorbar(aes(ymax = as.numeric(`X3rdQ`), ymin = as.numeric(`X1stQ`)),
                position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_point(mapping = aes(y = as.numeric(`Min`)), color = "blue") +
  geom_point(mapping = aes(y = as.numeric(`Max`)), color = "red") +
  coord_flip()
c
plot_grid(a,b)


# problems:
# n obs
# plotting from summary statistics more difficult than generating summary statistics and then plotting
# revise original function to filter
# query different databases
# downstream applications - practical
# other database tools: rfishbase, raquamaps, spocc multi query
# temperature data: robis pulls from BIO-Oracle which is long term averaged 1995-2020


