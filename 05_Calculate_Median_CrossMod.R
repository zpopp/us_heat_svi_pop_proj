# Created by: Zach Popp
# Date created: 01.05.2023
# Date modified: 04.16.2024
# Version no.: v1
#
# Overview: We want to calculate median metrics by SSP and across models.
# This script calculates the median, and transforms some of the large outputs
# at the population scale into metrics for 1000s or 100,000s. There is also
# code provided to evaluate which models were most represented in the final
# dataset
#

library(biscale)
library(cowplot)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(data.table)

# Set directories
#
proc_popgrid_dir <- "" # Where all model aggregates linked with population is stored
svi_dir <- ""          # Where SVI data is stored 

# Derivation of cross-model median measures for key hazard variables
# Read in processed dataset with exposure and vulnerability
#
us_temp_svi <- readRDS(paste0(proc_popgrid_dir, "all_mods_pop_20_50_v5.rds"))

# Convert to data.table
#
us_temp_svi_dt <- as.data.table(us_temp_svi)

# Calculate the median for each variable stratified by ssp and across models
#
us_temp_svi_med_proj <- us_temp_svi_dt[, lapply(.SD, function(uu) (median(uu, na.rm=TRUE))), 
                            .SDcols = names(us_temp_svi_dt)[grepl("proj|change", names(us_temp_svi_dt))], 
                            by = c("FIPS", "ssp")]

names(us_temp_svi_med_proj)[3:30] <- paste0("med_", names(us_temp_svi_med_proj)[3:30])
names(us_temp_svi_med_proj)[3:30]  <- gsub("_annual|_annual|_annual|_20|_50", "" , names(us_temp_svi_med_proj)[3:30])

# Calculate the median for historic variables based on subset to SSP585, where
# all models are available
#
us_temp_svi_dt_585 <- filter(us_temp_svi_dt, ssp == "ssp585")

us_temp_svi_med_hist <- us_temp_svi_dt_585[, lapply(.SD, function(uu) (median(uu, na.rm=TRUE))), 
                                       .SDcols = names(us_temp_svi_dt_585)[grepl("hist", names(us_temp_svi_dt_585))], 
                                       by = c("FIPS", "ssp")]

names(us_temp_svi_med_hist)[3:16] <- paste0("med_", names(us_temp_svi_med_hist)[3:16])
names(us_temp_svi_med_hist)[3:16]  <- gsub("_annual|_annual|_annual|_20|_50", "" , names(us_temp_svi_med_hist)[3:16])

us_temp_svi_med_hist$ssp <- NULL

# Join historic and projected data
#
us_temp_svi_med <- merge(us_temp_svi_med_proj, us_temp_svi_med_hist, by = "FIPS")

# Create exposure variables as by 1000 for easier reading
# Determination of hundred thousand or thousand denominator based on evaluating
# what will lead to max being between 100 - 1000
#
med_PDD_hundthou <- us_temp_svi_med[, lapply(.SD, function(uu) (uu / 100000)), 
                                    .SDcols = names(us_temp_svi_med)[grepl("med_PDD", names(us_temp_svi_med))], 
                                    by = c("FIPS", "ssp")]

names(med_PDD_hundthou)[3:8] <- paste0(names(med_PDD_hundthou)[3:8], "_hundthou")

# Process to be in 1000s
#
med_popexp_thou <- us_temp_svi_med[, lapply(.SD, function(uu) (uu / 1000)), 
                                    .SDcols = names(us_temp_svi_med)[grepl("med_PHD|med_PD95|med_PD_HI95", names(us_temp_svi_med))], 
                                    by = c("FIPS", "ssp")]

names(med_popexp_thou)[3:17] <- paste0(names(med_popexp_thou)[3:17], "_thou")

# Join together thousand and hundred thousands
#
us_temp_svi_med <- merge(us_temp_svi_med, med_PDD_hundthou, by = c("FIPS", "ssp"))
us_temp_svi_med <- merge(us_temp_svi_med, med_popexp_thou, by = c("FIPS", "ssp"))

# Save median output for hazards
#
saveRDS(us_temp_svi_med, paste0(proc_popgrid_dir, "all_mods_pop_20_50_med_v5.rds"))

# Read in SVI data for 2020
#
svi1 <- readRDS(paste0(svi_dir, "svi1.rds"))
svi1 <- svi1[c("FIPS", "SVI")]

# Merge median dataset with SVI
#
us_temp_svi <- left_join(us_temp_svi_med, svi1, by = "FIPS")
us_temp_svi$SVI <- round(us_temp_svi$SVI, digits=2)

# Create SVI quintile variable
#
us_temp_svi$SVI_quint <- cut(us_temp_svi$SVI, breaks=c(quantile(us_temp_svi$SVI, probs = seq(0, 1, by = 0.20))), 
                             labels=c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"))

# Read in regions
us_temp_svi$ST <- substr(us_temp_svi$FIPS, 1, 2)

#Read in processed dataset with exposure and vulnerability, including County labels
new_eng <- list("09","23","25","33","44","50")
mid_atl <- list("34","36","42")
east_nc <- list("18","17","26","39","55")
west_nc <- list("19","20","27","29","31","38","46")
south_atl <- list("10","11","12","13","24","37","45","51","54")
east_sc <- list("01","21","28","47")
west_sc <- list("05","22","40","48")
mountain <- list("04","08","16","35","30","49","32","56")
pacific <- list("06","41","53")
all_states <- unique(us_temp_svi$ST)

west <- c(mountain, pacific)
south <- c(south_atl, east_sc, west_sc)
midwest <- c(east_nc, west_nc)
northeast <- c(new_eng, mid_atl)

# Classify tracts
us_temp_svi$region <- ifelse(us_temp_svi$ST %in% west, "West",
                             ifelse(us_temp_svi$ST %in% south, "South",
                                    ifelse(us_temp_svi$ST %in% midwest, "Midwest",
                                           ifelse(us_temp_svi$ST %in% northeast, "Northeast", NA))))

# Classify tracts subregion
us_temp_svi$subregion <- ifelse(us_temp_svi$ST %in% new_eng, "New England",
                                ifelse(us_temp_svi$ST %in% mid_atl, "Mid Atlantic",
                                       ifelse(us_temp_svi$ST %in% east_nc, "East North Central",
                                              ifelse(us_temp_svi$ST %in% west_nc, "West North Central", 
                                                     ifelse(us_temp_svi$ST %in% south_atl, "South Atlantic", 
                                                            ifelse(us_temp_svi$ST %in% east_sc, "East South Central", 
                                                                   ifelse(us_temp_svi$ST %in% west_sc, "West South Central", 
                                                                          ifelse(us_temp_svi$ST %in% mountain, "Mountain", 
                                                                                 ifelse(us_temp_svi$ST %in% pacific, "Pacific", NA)))))))))


# Read in population data
#
final_pop_20_50 <- readRDS(paste0(proc_popgrid_dir, "final_pop_20_50.rds"))

# Merge population and hazard data
#
us_temp_svi_pop <- left_join(us_temp_svi, final_pop_20_50, by = "FIPS")

# Update population based on ssp
#
us_temp_svi_pop$Pop50 <- ifelse(us_temp_svi_pop$ssp == "ssp245", us_temp_svi_pop$Pop50_ssp2,
                                ifelse(us_temp_svi_pop$ssp == "ssp370", us_temp_svi_pop$Pop50_ssp3,
                                       ifelse(us_temp_svi_pop$ssp == "ssp585", us_temp_svi_pop$Pop50_ssp5, NA)))

# Save relevant columns
#
us_temp_svi_pop <- as.data.frame(us_temp_svi_pop)

# Save final dataset
#
saveRDS(us_temp_svi_pop, paste0(proc_popgrid_dir, "all_mods_pop_20_50_med_update_v5.rds"))


