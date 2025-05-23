# Created by: Zach Popp, Ian Sue Wing
# Date created: 12.07.2023
# Date modified: 04.01.2025
# Version no.:
# Updates since last version:   
# Overview: This script will take a set of projected climate netCDF files and 
#           historic climate netCDF, calculate the 95%ile max temperature and
#           heat index, the average annual number of cooling degree days, and
#           multiple metrics of extreme heat (days over 37.5 dry temp, days
#           with heat index > 40F) for both the contemporary (1995 - 2014) and
#           projected (2041 - 2060) period.
#           
#           For each period, then calculate the mean of these variables across
#           US census tracts. 
#           
#           The result will be a dataset with social vulnerability
#           index data, 95 percentile maximum temperature + cooling degree day
#           data for 1995-2014, 2041-2060, and their difference at the tract level
#           Separate files will be generated for each model and ssp scenario
#           input through bash scripting.
#
#           This script assume that the NEX-CMIP6 daily average (tas), daily max
#           (tasmax) and daily specific humidity (huss)are downloaded and stored 
#           on your system (download from https://nex-gddp-cmip6.s3.us-west-2.amazonaws.com/index.html#NEX-GDDP-CMIP6/)
#           with a file structure storing models by SSP and measure. 
#           In this application all models for a given SSP (ie: ssp585) and 
#           measure (ie: tmax) are in a single directory. Modifications may 
#           be required to the file structure based on your local file system.
#
#           The SVI data and geometry are used for analyses. These data can be 
#           accessed at: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/MZF7XQ  

require(terra)
require(sf)
require(data.table)
require(ggplot2)
require(dplyr)
require(weathermetrics)
require(HeatStress)

# The below inputs will come from a bash script. This will allow processing 
# across distinct models and ssp scenarios to be conducted across different
# nodes on the SCC
# 
# The input [b] is an index from 1:32 which will be run against a list of all
# available CMIP6 models in order to select which model to process. The input
# scenario is either "ssp245", "ssp370", or "ssp585".
#
# This code will need to be modified based on the storage system for exposure-
# level CMIP6 model outputs. For this project, models were stored by exposure,
# by SSP, and with labels of model and year.
#
# ie:
# directory/tasmax/historical/tasmax_day_<MODEL>_<SSP>_r1i1p1f1_gn_<YEAR>.nc
#
args <- commandArgs(trailingOnly = TRUE)
b <- as.numeric(args[1]) # index for model
scenario <- as.character(args[2]) # SSP scenario

# Set base directory for climate model data
#
base_cmip_dir <- "" # where your cmip6 data are stored
out_mod_dir <- "" # where your outputs should go (script assumes that mod-specific dirs are created)
svi_dir <- "" # where your SVI data are stored

# The below process is used to set the list of models against which to index
# in bash.
# Extract applicable model names to be selected via bash. 
#
check <- list.files(paste0(base_cmip_dir, 'tasmax/historical'),
                    pattern=paste0('.*2000.nc'),
                    full.names = F)

# Establish naming scheme to extract models
#
pattern <- "day_(.*?)_historical"
match <- regexpr(pattern, check, perl=TRUE)

# Remove the prefix 'day_' and suffix '_historical' from the result
#
models <- regmatches(check, match)

models <- sub("day_", "", models)
models <- sub("_historical", "", models)

models <- models[models %in% c("HadGEM3-GC31-MM", "GFDL-ESM4", "EC-Earth3",
                   "CNRM-ESM2-1", "CNRM-CM6-1", "CMCC-ESM2",
                   "ACCESS-CM2", "ACCESS-ESM1-5")]

# Subset to index from bash
#
model_spec <- models[b]

# Add text for index 12 or 13 which are both GFDL-CM4
#
if ( b == 12) {
  model_spec <- "GFDL-CM4_.*_r1i1p1f1_gr1"
} else if (b == 13) {
  model_spec <- "GFDL-CM4_.*_r1i1p1f1_gr2"
}

############# Preparing SVI Data ##############
# Read in SVI data for 2020. We use the SVI native shapefile to run the process-
# ing. 2020 census tract geographies are used.
#
svi1<-readRDS(paste0(svi_dir, "svi1.rds"))

# Capture extent of SVI object, which will represent the United States
#
e0 <- st_bbox(svi1)

# Modify extent to be comparable to climate models, which use an x axis running
# from 0 to 360 instead of -180 to +180
#
e1 <- ext(c(e0['xmin']+360,e0['xmax']+360,e0['ymin'],e0['ymax']))

# Save a new SVI file to tag hazard metrics too
#
svi2 <- svi1

# Update the geometry to be consistent with projected climate models
#
svi2$Shape <- st_geometry(svi2) + c(360, 0)

# Convert sf object to SpatVector in order to use terra functions
#
svi2 <- vect(svi2)

# Read in example raster to get CRS
#
example_file <- list.files(paste0(base_cmip_dir, 'tasmax/historical'),
                                        pattern=paste0('.*',model_spec,'.*199[5-9].nc'),
                                        full.names = T)[1]
example_rast <- rast(example_file)

# Project SVI to CRS of Rasters
#
crs(svi2) <- crs(example_rast)
svi2 <- project(svi2, crs(example_rast))

# Save SF version of object to bind results to
#
svi3 <- st_as_sf(svi2)

################################################

############# Processing Max Temperature TASMAX95 ##############

# The function below will take a list of netCDF files and an extent, then generate
# the 95 percentile maximum temperature for the range of time represented by the files
# and range of locations represented by the extent
#
prep_tasmax <- function(file_list, extent) {

  # Create rasters from listed files
  #
  s0 <- rast(file_list)
  
  # Crop raster to extent of US
  #
  s1 <- crop(s0,e1)
  
  # Apply calculation for 95th percentile Temperature, and convert from K to C
  #
  s2 <- app(s1,quantile,0.95,na.rm=T) - 273.15
   
  # Return output for further use
  #
  return(s2)
}

# This list includes historical data for 1995 to 1999
#
tasmax_file_list_hist_90s <- list.files(base_cmip_dir, 'tasmax/historical',
                                        pattern=paste0('.*',model_spec,'.*199[5-9].nc'),
                                        full.names = T)

# This list includes historical data for 2000 to 2014
#
tasmax_file_list_hist_00s <- list.files(base_cmip_dir, 'tasmax/historical',
                                        pattern=paste0('.*',model_spec,'.*20[0-1].*.nc'),
                                        full.names = T)

# We combine the lists above for our historical data files
#
tasmax_file_list_hist <- c(tasmax_file_list_hist_90s, tasmax_file_list_hist_00s)

# This list includes projected data for 2040 to 2060
#
tasmax_file_list_proj <- list.files(paste0(base_cmip_dir, 'tasmax/',scenario),
                                    pattern=paste0('.*',model_spec,'.*20[4-6].*.nc'),
                                    full.names = T)

# This will subset data for 2041 to 2060
#
tasmax_file_list_proj <- tasmax_file_list_proj[2:21]

# Output step in process for bash QC
#
cat("Processing tasmax", scenario, "of", model_spec, "\n")

# Calculate the tasmax95 for historical data and projected data for the US.
#
tmax_hist <- prep_tasmax(tasmax_file_list_hist, e1)
tmax_proj <- prep_tasmax(tasmax_file_list_proj, e1)

# Calculate the change across epochs in tasmax95
#
tmax_change <- tmax_proj - tmax_hist

# Use terra extract function to calculate mean across census tract for historic,
# projected, and change in tasmax over time.
#
x0_max <- terra::extract(tmax_change, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_max_hist <- terra::extract(tmax_hist, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_max_proj <- terra::extract(tmax_proj, svi2, fun=mean, na.rm=TRUE, exact=TRUE)

# Bind zonal outputs with SVI dataset, and rename variables as they are added
# to reflect the correct measure
#
svi3<-cbind(svi3,x0_max)
setnames(svi3,'lyr.1','tasmax95_change')
svi3<-cbind(svi3,x0_max_hist)
setnames(svi3,'lyr.1','tasmax95_hist')
svi3<-cbind(svi3,x0_max_proj)
setnames(svi3,'lyr.1','tasmax95_proj')

################################################

############# Processing Max Temperature HOT DAYS (37.5) ##############

# The function below will take a list of netCDF files and an extent, then generate
# the number of days per year > 37.5 C for the range of time represented by the files
# and range of locations represented by the extent
#
prep_hotdays <- function(file_list, extent) {
  
  # Create rasters from file list
  #
  d0 <- rast(file_list)
  
  # Crop raster to extent of US
  #
  d1 <- crop(d0,e1)
  
  # Modify values to Celsius
  #
  d2 <- d1 - 273.15
  
  # Create a new SpatRaster with hot day value (yes/no for over 37.5C)
  #
  result <- ifel(d2 > 37.5, 1, 0)
  
  # Sum hot days across study period, then divide by 20 to get annual metric
  #
  d3 <- app(result, sum, na.rm = T)/20
  
  # Return output for further use
  #
  return(d3)
}

# Calculate the tasmax for historical data and projected data for the US.
#
hotday_hist <- prep_hotdays(tasmax_file_list_hist, e1)
hotday_proj <- prep_hotdays(tasmax_file_list_proj, e1)

# Calculate the change over time in tasmax
#
hotday_change <- hotday_proj - hotday_hist

# Use terra zonal function to calculate mean across census tract for historic,
# projected, and change in hot days over time.
#
x0_hds <- terra::extract(hotday_change, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_hds_hist <- terra::extract(hotday_hist, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_hds_proj <- terra::extract(hotday_proj, svi2, fun=mean, na.rm=TRUE, exact=TRUE)

# Bind zonal outputs with SVI dataset, and rename variables as they are added
# to reflect the correct measure
#
svi3<-cbind(svi3,x0_hds)
setnames(svi3,'sum','hotday_37_change_annual')
svi3<-cbind(svi3,x0_hds_hist)
setnames(svi3,'sum','hotday_37_hist_annual')
svi3<-cbind(svi3,x0_hds_proj)
setnames(svi3,'sum','hotday_37_proj_annual')

################################################

############# Processing Average Temperature CDD24s ##############

# The function below will take a list of netCDF files and an extent, then generate
# the number of cooling degree days (cumulative number of average daily temperature
# degrees over a threshold across some time span, here 24 C and per year (after division by 20))
#
prep_cdd_24 <- function(file_list, extent) {
  
  # Create rasters from file list
  #
  r0<- rast(file_list)
  
  #C rop raster to extent of US
  #
  r1 <- crop(r0, e1)
  
  # Modify values to Celsius
  #
  r2 <- r1 - 273.15
  
  # Create a new SpatRaster with CDD value
  #
  result <- ifel(r2 > 24, r2 - 24, 0)
  
  # Sum CDD across study period, then divide by 20 to get annual metric
  #
  r3 <- app(result,sum,na.rm=T)/20
  
  return(r3)
}

# This list includes historical data for 1995 to 1999, note use of tas instead of tasmax
#
tas_file_list_hist_90s <- list.files(base_cmip_dir, 'tas/historical',
                                     pattern=paste0('.*',model_spec,'.*199[5-9].nc'),
                                     full.names = T)

# This list includes historical data for 2000 to 2014
#
tas_file_list_hist_00s <- list.files(base_cmip_dir, 'tas/historical',
                                     pattern=paste0('.*',model_spec,'.*20[0-1]..nc'),
                                     full.names = T)

# We combine the lists above for our historical data files
#
tas_file_list_hist <- c(tas_file_list_hist_90s, tas_file_list_hist_00s)

# This list includes projected fata for 2040 to 2060
#
tas_file_list_proj <- list.files(paste0(base_cmip_dir, 'tas/',scenario),
                                 pattern=paste0('.*',model_spec,'.*20[4-6]..nc'),
                                 full.names = T)

# This will subset data for 2041 to 2060
#
tas_file_list_proj<- tas_file_list_proj[2:21]

# Output step in process for bash QC
#
cat("Processing tas", scenario, "of", model_spec, "\n")

# Calculate the annual CDDs for historical data and projected data for the US.
#
cdd_hist_24 <- prep_cdd_24(tas_file_list_hist, e1)
cdd_proj_24 <- prep_cdd_24(tas_file_list_proj, e1)

# Calculate the change over time in annual CDDs
#
cdd_change_24 <- cdd_proj_24 - cdd_hist_24

# Use terra zonal function to calculate mean across census tract for historic,
# projected, and change in CDD over time.
#
x0_cdd_24 <- terra::extract(cdd_change_24, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_cdd_hist_24 <- terra::extract(cdd_hist_24, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_cdd_proj_24 <- terra::extract(cdd_proj_24, svi2, fun=mean, na.rm=TRUE, exact=TRUE)

# Bind zonal outputs with SVI dataset, and rename variables as they are added
# to reflect the correct measure
#
svi3<-cbind(svi3,x0_cdd_24)
setnames(svi3,'sum','CDD_24_change_annual')
svi3<-cbind(svi3,x0_cdd_hist_24)
setnames(svi3,'sum','CDD_24_hist_annual')
svi3<-cbind(svi3,x0_cdd_proj_24)
setnames(svi3,'sum','CDD_24_proj_annual')

################################################

############# Processing Average Temperature CDD18s ##############

# The function below will take a list of netCDF files and an extent, then generate
# the number of cooling degree days (cumulative number of average daily temperature
# degrees over a threshold across some time span, here 24 C and per year (after division by 20))
#
prep_cdd_18 <- function(file_list, extent) {
  
  # Create rasters from file list
  #
  r0<- rast(file_list)
  
  #C rop raster to extent of US
  #
  r1 <- crop(r0, e1)
  
  # Modify values to Celsius
  #
  r2 <- r1 - 273.15
  
  # Create a new SpatRaster with CDD value
  #
  result <- ifel(r2 > 18, r2 - 18, 0)
  
  # Sum CDD across study period, then divide by 20 to get annual metric
  #
  r3 <- app(result,sum,na.rm=T)/20
  
  return(r3)
}


# Output step in process for bash QC
#
cat("Processing tas", scenario, "of", model_spec, "\n")

# Calculate the annual CDDs for historical data and projected data for the US.
#
cdd_hist_18 <- prep_cdd_18(tas_file_list_hist, e1)
cdd_proj_18 <- prep_cdd_18(tas_file_list_proj, e1)

# Calculate the change over time in annual CDDs
#
cdd_change_18 <- cdd_proj_18 - cdd_hist_18

# Use terra zonal function to calculate mean across census tract for historic,
# projected, and change in CDD over time.
#
x0_cdd_18 <- terra::extract(cdd_change_18, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_cdd_hist_18 <- terra::extract(cdd_hist_18, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_cdd_proj_18 <- terra::extract(cdd_proj_18, svi2, fun=mean, na.rm=TRUE, exact=TRUE)


# Bind zonal outputs with SVI dataset, and rename variables as they are added
# to reflect the correct measure
#
svi3<-cbind(svi3,x0_cdd_18)
setnames(svi3,'sum','CDD_18_change_annual')
svi3<-cbind(svi3,x0_cdd_hist_18)
setnames(svi3,'sum','CDD_18_hist_annual')
svi3<-cbind(svi3,x0_cdd_proj_18)
setnames(svi3,'sum','CDD_18_proj_annual')

################################################

############ Save result for non - Humidity ####################################

# There are several models that do not have humidity available. To avoid 
# having our script crash for these models, we output them here and then end
# the script. The rest of the models continue to heat index processing below
#
if (model_spec == "MIROC6" |
    model_spec == "NESM3" |
    model_spec == "IPSL-CM6A-LR") {
  
  # Add scenario and model to mark output
  #
  svi3$ssp <- scenario
  svi3$mod <- model_spec
  
  # Add NA for humidity vars to assist with join later
  #
  svi3$heatindex95_change <- svi3$heatindex95_hist <- svi3$heatindex95_proj <-
    svi3$hotday_hi_change <- svi3$hotday_hi_hist <- svi3$hotday_hi_proj <- 
    svi3$hotday_hi_ex_change <- svi3$hotday_hi_ex_hist <- svi3$hotday_hi_ex_proj <- NA
  
  # Save output to model and scenario specific file including all measures
  #
  saveRDS(svi3, paste0(out_mod_dir, model_spec,"/us_climate_svi_",model_spec,"_",scenario,"_v5.rds"))
  
  # Break to move to next scenario
  #
  break
}

################################################

############# Processing Heat Index (HIMAX95) #####################################

# The function below will take lists of netCDF files for max temperature and 
# for specific humidity, an extent, and then generate heat index and 
# the 95th percentile heat index over a threshold across some time span, 
# here 24 C and per year (after division by 20))
# 
# Reference for equation:
#
# Set function based on reference: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
#
# Set constants
#
p_base <- 101325

rh <- function(temp_in, huss_in) {
  rh <- (0.263 * huss_in * p_base) / exp((17.67*(temp_in - 273.15)/(temp_in - 29.65)))
  
  return(rh)
}

prep_hi_max <- function(file_list_huss_in, file_list_tmax_in, extent) {
  
  # Create rasters from file lists
  #
  r0_huss <- rast(file_list_huss_in)
  r0_tmax <- rast(file_list_tmax_in)
  
  # Crop rasters to extent of US
  #
  r1_huss <- crop(r0_huss, e1)
  r1_tmax <- crop(r0_tmax, e1)
  
  # Calculate relative humidity from specific humidity and temperature 
  # assuming constant pressure
  # Source: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
  #
  r1_rh <- xapp(r1_tmax, r1_huss, fun = function(temp, huss) {
    rh(temp_in = temp, huss_in = huss)
  })
  
  # If relative humidity is greater than 100, reset to 100. RH > 100 is 
  # representative of rain, for the purposes of heat index calculation, this
  # supersaturation should not be accounted for in further adjusting the heat
  # index
  #
  r1_rh <- ifel(r1_rh > 100, 100, r1_rh)
  
  # Modify values to Celsius
  #
  r1_tmax <- r1_tmax - 273.15
  
  # Compute heat index from rh and tmax
  #
  r1_hi <- xapp(r1_tmax, r1_rh, fun = function(temp, relhum) {
    heat.index(t = temp, rh = relhum,
               temperature.metric = "celsius",
               output.metric = "celsius", round = 6)
  })
  
  # Apply calculation for 95th percentile Heat Index
  #
  hi_95 <- app(r1_hi,quantile,0.95,na.rm=T)
  
  # Return output for further use
  #
  return(hi_95)
}


# This list includes historical data for 1995 to 1999, note use of huss 
#
huss_file_list_hist_90s <- list.files(base_cmip_dir, 'huss/historical',
                                     pattern=paste0('.*',model_spec,'.*199[5-9].nc'),
                                     full.names = T)

# This list includes historical data for 2000 to 2014
#
huss_file_list_hist_00s <- list.files(base_cmip_dir, 'huss/historical',
                                     pattern=paste0('.*',model_spec,'.*20[0-1]..nc'),
                                     full.names = T)

# We combine the lists above for our historical data files
#
huss_file_list_hist <- c(huss_file_list_hist_90s, huss_file_list_hist_00s)

# This list includes projected fata for 2040 to 2060
#
huss_file_list_proj <- list.files(paste0(base_cmip_dir, 'huss/',scenario),
                                 pattern=paste0('.*',model_spec,'.*20[4-6]..nc'),
                                 full.names = T)

# This will subset data for 2041 to 2060
#
huss_file_list_proj<- huss_file_list_proj[2:21]

# Output step in process for bash QC
#
cat("Processing huss", scenario, "of", model_spec, "\n")

# Calculate the annual CDDs for historical data and projected data for the US.
#
hi_hist<- prep_hi_max(file_list_huss_in = huss_file_list_hist, 
                      file_list_tmax_in = tasmax_file_list_hist, 
                      e1)
hi_proj<- prep_hi_max(file_list_huss_in = huss_file_list_proj, 
                      file_list_tmax_in = tasmax_file_list_proj, 
                      e1)

# Calculate the change over time in annual CDDs
#
hi_change <- hi_proj - hi_hist

# Use terra zonal function to calculate mean across census tract for historic,
# projected, and change in 95th percentile heat index over time.
#
x0_hi <- terra::extract(hi_change, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_hi_hist <- terra::extract(hi_hist, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_hi_proj <- terra::extract(hi_proj, svi2, fun=mean, na.rm=TRUE, exact=TRUE)

# Bind zonal outputs with SVI dataset, and rename variables as they are added
# to reflect the correct measure
#
svi3<-cbind(svi3,x0_hi)
setnames(svi3,'lyr.1','heatindex95_change')
svi3<-cbind(svi3,x0_hi_hist)
setnames(svi3,'lyr.1','heatindex95_hist')
svi3<-cbind(svi3,x0_hi_proj)
setnames(svi3,'lyr.1','heatindex95_proj')

################################################################################

############# Processing Max Heat Index HOT DAYS (HI40) ###############################

# The function below will take a list of netCDF files and an extent, then generate
# the number of days per year > 37.5 C for the range of time represented by the files
# and range of locations represented by the extent
#
prep_hotdays_hi <- function(file_list_huss_in, file_list_tmax_in, extent) {
  
  # Create rasters from file lists
  #
  r0_huss <- rast(file_list_huss_in)
  r0_tmax <- rast(file_list_tmax_in)
  
  # Crop rasters to extent of US
  #
  r1_huss <- crop(r0_huss, e1)
  r1_tmax <- crop(r0_tmax, e1)
  
  # Calculate relative humidity from specific humidity and temperature 
  # assuming constant pressure
  # Source: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
  #
  r1_rh <- xapp(r1_tmax, r1_huss, fun = function(temp, huss) {
    rh(temp_in = temp, huss_in = huss)
  })
  
  # If relative humidity is greater than 100, reset to 100. RH > 100 is 
  # representative of rain, for the purposes of heat index calculation, this
  # supersaturation should not be accounted for in further adjusting the heat
  # index
  #
  r1_rh <- ifel(r1_rh > 100, 100, r1_rh)
  
  # Modify values to Celsius
  #
  r1_tmax <- r1_tmax - 273.15
  
  # Compute heat index from rh and tmax
  #
  r1_hi <- xapp(r1_tmax, r1_rh, fun = function(temp, relhum) {
    heat.index(t = temp, rh = relhum,
               temperature.metric = "celsius",
               output.metric = "celsius", round = 6)
  })
  
  # Create a new SpatRaster with hot day value (yes/no for over 37.5C)
  #
  result <- ifel(r1_hi > 40, 1, 0)
  
  # Sum hot days across study period, then divide by 20 to get annual metric
  #
  d3 <- app(result, sum, na.rm = T)/20
  
  # Return output for further use
  #
  return(d3)
}

# Calculate the tasmax for historical data and projected data for the US.
#
hotday_hi_hist <- prep_hotdays_hi(huss_file_list_hist, tasmax_file_list_hist, e1)
hotday_hi_proj <- prep_hotdays_hi(huss_file_list_proj, tasmax_file_list_proj, e1)

# Calculate the change over time in tasmax
#
hotday_hi_change <- hotday_hi_proj  - hotday_hi_hist

# Use terra zonal function to calculate mean across census tract for historic,
# projected, and change in 95th percentile heat index over time.
#
x0_hotday_hi <- terra::extract(hotday_hi_change, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_hotday_hi_hist <- terra::extract(hotday_hi_hist, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_hotday_hi_proj <- terra::extract(hotday_hi_proj, svi2, fun=mean, na.rm=TRUE, exact=TRUE)


# Bind zonal outputs with SVI dataset, and rename variables as they are added
# to reflect the correct measure
#
svi3<-cbind(svi3,x0_hotday_hi)
setnames(svi3,'sum','hotday_hi_change')
svi3<-cbind(svi3,x0_hotday_hi_hist)
setnames(svi3,'sum','hotday_hi_hist')
svi3<-cbind(svi3,x0_hotday_hi_proj)
setnames(svi3,'sum','hotday_hi_proj')

################################################

############# Processing Max Heat Index EXTREME HOT DAYS ###############################

# The function below will take a list of netCDF files and an extent, then generate
# the number of days per year > 37.5 C for the range of time represented by the files
# and range of locations represented by the extent
#
prep_hotdays_hi_ex <- function(file_list_huss_in, file_list_tmax_in, extent) {
  
  # Create rasters from file lists
  #
  r0_huss <- rast(file_list_huss_in)
  r0_tmax <- rast(file_list_tmax_in)
  
  # Crop rasters to extent of US
  #
  r1_huss <- crop(r0_huss, e1)
  r1_tmax <- crop(r0_tmax, e1)
  
  # Calculate relative humidity from specific humidity and temperature 
  # assuming constant pressure
  # Source: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
  #
  r1_rh <- xapp(r1_tmax, r1_huss, fun = function(temp, huss) {
    rh(temp_in = temp, huss_in = huss)
  })
  
  # If relative humidity is greater than 100, reset to 100. RH > 100 is 
  # representative of rain, for the purposes of heat index calculation, this
  # supersaturation should not be accounted for in further adjusting the heat
  # index
  #
  r1_rh <- ifel(r1_rh > 100, 100, r1_rh)
  
  # Modify values to Celsius
  #
  r1_tmax <- r1_tmax - 273.15
  
  # Compute heat index from rh and tmax
  #
  r1_hi <- xapp(r1_tmax, r1_rh, fun = function(temp, relhum) {
    heat.index(t = temp, rh = relhum,
               temperature.metric = "celsius",
               output.metric = "celsius", round = 6)
  })
  
  # Create a new SpatRaster with hot day value (yes/no for over 37.5C)
  #
  result <- ifel(r1_hi > 40, 1, 0)
  
  # Sum hot days across study period, then divide by 20 to get annual metric
  #
  d3 <- app(result, sum, na.rm = T)/20
  
  # Return output for further use
  #
  return(d3)
}

# Calculate the tasmax for historical data and projected data for the US.
#
hotday_hi_ex_hist <- prep_hotdays_hi_ex(huss_file_list_hist, tasmax_file_list_hist, e1)
hotday_hi_ex_proj <- prep_hotdays_hi_ex(huss_file_list_proj, tasmax_file_list_proj, e1)

# Calculate the change over time in tasmax
#
hotday_hi_ex_change <- hotday_hi_ex_proj  - hotday_hi_ex_hist

# Use terra zonal function to calculate mean across census tract for historic,
# projected, and change in 95th percentile heat index over time.
#
x0_hotday_hi_ex <- terra::extract(hotday_hi_ex_change, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_hotday_hi_ex_hist <- terra::extract(hotday_hi_ex_hist, svi2, fun=mean, na.rm=TRUE, exact=TRUE)
x0_hotday_hi_ex_proj <- terra::extract(hotday_hi_ex_proj, svi2, fun=mean, na.rm=TRUE, exact=TRUE)

# Bind zonal outputs with SVI dataset, and rename variables as they are added
# to reflect the correct measure
#
svi3<-cbind(svi3,x0_hotday_hi_ex)
setnames(svi3,'sum','hotday_hi_ex_change')
svi3<-cbind(svi3,x0_hotday_hi_ex_hist)
setnames(svi3,'sum','hotday_hi_ex_hist')
svi3<-cbind(svi3,x0_hotday_hi_ex_proj)
setnames(svi3,'sum','hotday_hi_ex_proj')

################################################

############# Save Outputs #####################################################

# Add text for index 12 or 13 which are both GFDL-CM4
#
if ( b == 12) {
  model_spec <- "GFDL-CM4_gr1"
} else if (b == 13) {
  model_spec <- "GFDL-CM4_gr2"
}

# Add scenario and model to mark output
#
svi3$ssp <- scenario
svi3$mod <- model_spec

# Save output to model and scenario specific file including all measures
#
saveRDS(svi3, paste0(output_dir, model_spec, "/us_climate_svi_", model_spec, "_",scenario, ".rds"))
