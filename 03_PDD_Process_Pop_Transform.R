# Created by: Zach Popp
# Date created: 12.07.2023
# Date modified: 
# Version no: v1
#
# Overview: We want to calculate population hazard. For the contemporary period,
#           we will use 2020 population data from TIGRIS. 
#           For the projected period, we will use decadal projections generated 
#           by Gao et al. which are SSP-specific to estimate the change rate 
#           from 2050 to 2020. We will use this change rate to 
#           calculate estimated 2050 population.
#           
#           Population exposure metrics are then calculated by multiplying 
#           heat hazard metrics by the tract 2020 and 2050 populations
#
require(terra)
require(sf)
require(data.table)
require(ggplot2)
library(dplyr)
library(tidycensus)

# Set directories
#
svi_dir <- "" # Directory where SVI data are stored
census_dir <- "" # Directory where state FIPS list is stored (you can build yourself with tigris::states())
proj_dir <- "" # Directory where 2020 population will be stored
popgrid_dir <- "" # Directory where population projections were stored and where outputs will go ( downloaded from SEDAC https://doi.org/10.7927/q7z9-9r6, accessed 1-25-24)
models_dir <- "" # Directory where model aggregated (all mod/ssp) from script 1 and 2a are stored

# Read in original SVI file
#
svi1 <- readRDS(paste0(svi_dir, "svi1.rds"))

# Save relevant columns
#
svi1 <- svi1[c("FIPS", "Shape")]

# Step 1: Read in tract population in 2020 for census tracts by state.

################# Read in tract population data ################################
# Read in state FIPS
states <- read.csv(paste0(census_dir, "US_States_FIPS_Codes.csv"), stringsAsFactors = FALSE)
states <- states %>% filter(StAbbrev != "AK" & StAbbrev != "HI")
states$StFIPS <- formatC(as.numeric(states$StFIPS), width = 2, format = "d", flag = "0")
StFIPS <- states$StFIPS

for (i in c(StFIPS)) {
  
  # Read in tract population for 2020 from tidycensus
  #
  tractpop <- get_decennial(geography = "tract",
                            variables = "P1_001N",
                            state = i,
                            year = 2020)
  
  # Rename and select relevant variables
  tractpop <- tractpop %>%
    rename(Pop20 = value) %>%
    select(GEOID, Pop20)
  
  # Merge tract population data
  #
  if (i == StFIPS[1]) {
    final_pop <- tractpop
  } else {
    final_pop <- rbind(final_pop, tractpop)
  }
  
}

# Save output
#
saveRDS(final_pop, paste0(proj_dir, "pop2020_tidycensus.rds"))

################################################################################

# Step 2: Assess weighted mean population projections in 2020 and 2050,
# using a for loop for the projections across SSPs. Calculate fractional change
# because projections started in 2000 and differ from 2020 actual populations
# 
################# Processing Gau Proj. Population data #########################

# Determine scenario
#
scenarios <- c("ssp245", "ssp370", "ssp585")

# Run loop through scenarios
#
for (scen in c(scenarios)) {
  
  # Output step in process for bash QC
  #
  cat("Processing pop ", scen, "\n")

  # Create indicator of scenario in naming convention of population grid
  #
  pop_scenario <- substr(scen, 1, 4)
  
  # Read in 1km gridded population
  #
  pop_grid_20 <- rast(paste0(popgrid_dir, toupper(pop_scenario),"_1km/",pop_scenario,"_total_2020.tif"))
  pop_grid_50 <- rast(paste0(popgrid_dir,toupper(pop_scenario),"_1km/",pop_scenario,"_total_2050.tif"))
  
  # Convert SVI tracts to vect object in order to extract raster values
  #
  svi1_v <- vect(svi1)
  
  # Project raster data to NAD 83
  #
  pop_grid_20_proj <- project(pop_grid_20, crs(svi1_v))
  pop_grid_50_proj <- project(pop_grid_50, crs(svi1_v))
  
  # Crop raster to extent of US
  #
  pop_grid_US_20 <- crop(pop_grid_20_proj, ext(svi1_v))
  pop_grid_US_50 <- crop(pop_grid_50_proj, ext(svi1_v))
  
  # Extract sum population across census tract
  #
  poly_zon_20 <- terra::extract(pop_grid_US_20, svi1_v, fun=sum, na.rm=TRUE, weights=TRUE)
  poly_zon_50 <- terra::extract(pop_grid_US_50, svi1_v, fun=sum, na.rm=TRUE, weights=TRUE)
  
  # Bind zonal output to dataset
  svi_pop <-cbind(svi1, poly_zon_20)
  svi_pop <-cbind(svi_pop, poly_zon_50)
  
  # Convert dataset to df object
  #
  svi_pop <- as.data.frame(svi_pop)
  svi_pop$Shape <- NULL
  
  # Rename population input
  #
  newvarname_20 <- paste0('Pop_2020_sum_', pop_scenario)
  newvarname_50 <- paste0('Pop_2050_sum_', pop_scenario)
  setnames(svi_pop, paste0(pop_scenario, '_total_2020'),newvarname_20)
  setnames(svi_pop, paste0(pop_scenario, '_total_2050'),newvarname_50)
  
  # Calculate fractional change from 2020 to 2050
  #
  fract_var_name <- paste0("change_50_20_", pop_scenario)
  svi_pop[[fract_var_name]] <- svi_pop[[newvarname_50]] / svi_pop[[newvarname_20]]

  # Rebuild all model/scenario dataset
  #
  if (scen == scenarios[1]) {
    final_gau <- svi_pop
  } else {
    final_gau <- left_join(final_gau, svi_pop, by="FIPS")
  }
  
}

# Output population change measures
# Only needs to be saved once
#
saveRDS(final_gau, paste0(popgrid_dir, "Processed/pop_gau_fractional_change_20_50.rds"))

################################################################################

# Step 3: There are some tracts where the fractional change is registered as
# infinite/missing. Where these arise, we can impute the change rate for the 
# county. After cleaning, apply fractional change to 2020 population to get
# 2050 population projections by SSP
#
# Read in products from above if changing only step 3 below
tract_pop <- readRDS(paste0(proj_dir, "pop2020_tidycensus.rds"))
final_gau <- readRDS(paste0(popgrid_dir, "Processed/pop_gau_fractional_change_20_50.rds"))

################# Save 2050 to 2020 fractional change ##########################
# Prep tract_pop for join
tract_pop$FIPS <- tract_pop$GEOID
tract_pop <- tract_pop[c("FIPS", "Pop20")]

# Prep gau for join
#
final_gau_fract <- final_gau[c("FIPS", "change_50_20_ssp2", "change_50_20_ssp3", "change_50_20_ssp5" )]

# Merge tract_pop to gau
#
final_gau_fract_tract <- left_join(final_gau_fract, tract_pop, by = "FIPS")

# Create county FIPS in dataset
#
final_gau_fract_tract$CoFIPS <- substr(final_gau_fract_tract$FIPS, 1, 5)

# Check for errant values
#
summary(final_gau_fract_tract$change_50_20_ssp2)
summary(final_gau_fract_tract$change_50_20_ssp3)
summary(final_gau_fract_tract$change_50_20_ssp5)

# Update infinity change to NA, as well as tract where change is outlier and product of 
# extremely small base population
# 
final_gau_fract_tract$change_50_20_ssp5 <- ifelse(final_gau_fract_tract$change_50_20_ssp5 == Inf | final_gau_fract_tract$change_50_20_ssp5 > 20, NA, final_gau_fract_tract$change_50_20_ssp5)

# Group by county and calculate average change measures
#
final_gau_fract_tract_county <- final_gau_fract_tract %>%
  group_by(CoFIPS) %>%
  summarise(county_change_50_20_ssp2 = mean(change_50_20_ssp2, na.rm=T),
            county_change_50_20_ssp3 = mean(change_50_20_ssp3, na.rm=T),
            county_change_50_20_ssp5 = mean(change_50_20_ssp5, na.rm=T))

# Subset to tracts with 2020 population but NA change data
#
final_gau_missing <- final_gau_fract_tract %>%
  filter(is.na(change_50_20_ssp5))

# Subset to those with complete data
#
final_gau_complete <- final_gau_fract_tract %>%
  filter(!is.na(change_50_20_ssp5))

# Bring tract-level measures in for missing tracts
#
final_gau_missing <- left_join(final_gau_missing, final_gau_fract_tract_county, by = "CoFIPS")

# Reassign county changes as FIPS fractional change for missing data
#
final_gau_missing <- final_gau_missing %>%
  mutate(change_50_20_ssp2 = county_change_50_20_ssp2,
         change_50_20_ssp3 = county_change_50_20_ssp3,
         change_50_20_ssp5 = county_change_50_20_ssp5,)

# Save relevant columns
#
final_gau_missing <- final_gau_missing[c("FIPS", "change_50_20_ssp2", "change_50_20_ssp3", "change_50_20_ssp5", "Pop20", "CoFIPS")]

# Merge complete and missing data
#
final_gau_fract_tract <- rbind(final_gau_complete, final_gau_missing)

# Calculate 2050 pop based on fractional change measures from Gau
#
final_gau_fract_tract$Pop50_ssp2 <- final_gau_fract_tract$Pop20 * final_gau_fract_tract$change_50_20_ssp2 
final_gau_fract_tract$Pop50_ssp3 <- final_gau_fract_tract$Pop20 * final_gau_fract_tract$change_50_20_ssp3
final_gau_fract_tract$Pop50_ssp5 <- final_gau_fract_tract$Pop20 * final_gau_fract_tract$change_50_20_ssp5

# Limit data to only Pop 20 and Pop 50
#
final_pop_20_50 <- final_gau_fract_tract[c("FIPS", "Pop20", "Pop50_ssp2", "Pop50_ssp3", "Pop50_ssp5")]

# Save population estimates output
#
#saveRDS(final_pop_20_50, paste0(popgrid_dir, "Processed/final_pop_20_50.rds"))

################################################################################

# Step 4: Final step! Read in the 2050 population, read in all models
# merged together, join the data and calculate population exposure measures
#
################# Process population exposure metrics ##########################

# Read in population change data
#
final_pop_20_50 <- readRDS(paste0(popgrid_dir, "Processed/final_pop_20_50.rds"))

# Join final population measures to SVI data
#
svi5 <- readRDS(paste0(models_dir, "all_models.rds"))

# Join model outputs to tract measures
#
svi5 <- left_join(svi5, final_pop_20_50, by = "FIPS")

# Determine scenario
scenarios <- c("ssp245", "ssp370", "ssp585")

for (scen in c(scenarios)) {
  
  # Output step in process for bash QC
  #
  cat("Processing ", scen, "\n")
  
  # Subset data to scenario
  #
  svi_scen <- svi5 %>% filter(ssp == scen)
  
  # Create indicator of scenario in naming convention of population vars
  #
  pop_scenario <- substr(scen, 1, 4)
  
  # Establish new variable names
  #
  Pop50 <- paste0("Pop50_", pop_scenario)
  
  # Remove ID columns
  #
  ID_cols <- names(svi_scen)[grepl("ID", names(svi_scen))]
  keep_names <- names(svi_scen)[!names(svi_scen) %in% ID_cols]
  
  # Keep only relevant names
  #
  svi_scen <- svi_scen[c(keep_names)]
  
  # Calculate population change for SSP
  #
  svi_scen$Pop_Change <-  svi_scen[[Pop50]] - svi_scen$Pop20
  
  # Set variable inputs
  #
  haz <- c("CDD_18", "CDD_24", "hotday_37", "tasmax95", "heatindex95", "hotday_hi", "hotday_hi_ex")
  popexp <- c("PDD_18", "PDD_24", "PHD_37", "PD95", "PD_HI95", "PHD_HI40", "PHD_HI52")
  
  for (j in 1:length(haz)) {
    
    cat("Processing ", haz[j], "\n")
      
    hist_hazvar <- paste0(haz[j], "_hist")
    proj_hazvar <- paste0(haz[j], "_proj")
    change_hazvar <- paste0(haz[j], "_change")
    
    hist_popvar <- paste0(popexp[j], "_hist")
    proj_popvar <- paste0(popexp[j], "_proj")
    change_popvar <- paste0(popexp[j], "_change")
    
    # Run computation
    #
    svi_scen[[hist_popvar]] <- svi_scen[[hist_hazvar]] * svi_scen$Pop20
    svi_scen[[proj_popvar]] <- svi_scen[[proj_hazvar]] * svi_scen[[Pop50]]
    svi_scen[[change_popvar]] <- (svi_scen[[hist_hazvar]] * svi_scen$Pop_Change) + (svi_scen[[change_hazvar]] * svi_scen$Pop20) + (svi_scen[[change_hazvar]] * svi_scen$Pop_Change)
    
  }

  # Rebuild all model/scenario dataset
  #
  if (scen == scenarios[1]) {
    final_all_mod_pop <- svi_scen
  } else {
    final_all_mod_pop <- rbind(final_all_mod_pop, svi_scen)
  }
  
}

# Save output to model and scenario specific file including all measures
#
saveRDS(final_all_mod_pop, paste0(popgrid_dir, "Processed/all_mods_pop_20_50_v5.rds"))
