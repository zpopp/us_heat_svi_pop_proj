# Created by: Zach Popp
# Date created: 10.02.2024
# Date modified: 
# Version no.: v1
#     
# Overview: Generate box plot that stratifies by SSP and by SVI quintile, to 
# demonstrate how across models and SSPs, the outputs compare. Assess
# additional stratification by region
#

library(biscale)
library(cowplot)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(data.table)
library(egg)
library(ggpubr)

# Set directories
#
proc_popgrid_dir <- "" # Where population linked model data stored
boxplot_out <- ""   # Where to save figures

# Derivation of cross-model median measures for key hazard variables
# Read in processed dataset with exposure and vulnerability
#
us_temp_svi <- readRDS(paste0(proc_popgrid_dir, "all_mods_pop_20_50_v5.rds"))

# Convert to data.table
#
us_temp_svi_dt <- as.data.table(us_temp_svi)

# Assess available names
#
us_temp_svi_dt$SVI_cat <- ifelse(us_temp_svi_dt$SVI >= 0.8, "SVI Q5 (High Vuln)",
                                 ifelse(us_temp_svi_dt$SVI > 0.2 & us_temp_svi_dt$SVI <0.8, "SVI Q2-4",
                                        ifelse(us_temp_svi_dt$SVI <= 0.2, "SVI Q1 (Low Vuln)", NA)))

# Remove missing SVI
#
us_temp_svi_dt <- us_temp_svi_dt[!is.na(us_temp_svi_dt$SVI), ]

# Set factor variables
#
us_temp_svi_dt$SVI_cat <- factor(us_temp_svi_dt$SVI_cat, levels = c("SVI Q1 (Low Vuln)",  "SVI Q2-4", "SVI Q5 (High Vuln)"))
us_temp_svi_dt$ssp <- factor(us_temp_svi_dt$ssp, levels = c("ssp245", "ssp370", "ssp585"), labels = c("SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))


################ ALL TRACT + MODELS ############################################

boxplot_ssp <- function(data_in, var_in, var_name) {

  ggplot(data = us_temp_svi_dt[us_temp_svi_dt$SVI_cat %in% c("SVI Q5 (High Vuln)", "SVI Q1 (Low Vuln)"), ], aes(y = ssp, x = var_in, fill = ssp, col = SVI_cat)) +
    geom_boxplot(outliers = FALSE, show.legend = TRUE) +
    theme_classic(base_size = 25) +
    scale_fill_manual(values = c("#fdffe6", "#feb39b", "#b94250"))+
    scale_color_manual(values = c("#24BAFB", "black"))+
    labs(title = paste0(var_name, " Distribution"),
         x = paste0("Change in ", var_name),
         y = "Shared Socioeconomic Pathway",
         fill = "SSP",
         col = "Tract SVI Category") 
  
}

# TRACT LEVEL PLOTS - ALL MODELS and SSPs
# Create plots
#
haz_plot_ssp_cdd <- boxplot_ssp(var_in = us_temp_svi_dt[us_temp_svi_dt$SVI_cat %in% c("SVI Q5 (High Vuln)", "SVI Q1 (Low Vuln)"), ]$CDD_24_change_annual, var_name = "Cooling Degree Days (24C)")
haz_plot_ssp_hd_40 <- boxplot_ssp(var_in = us_temp_svi_dt[us_temp_svi_dt$SVI_cat %in% c("SVI Q5 (High Vuln)", "SVI Q1 (Low Vuln)"), ]$hotday_hi_change, var_name = "Hot Days (HI > 40C)")
haz_plot_ssp_tmax <- boxplot_ssp(var_in = us_temp_svi_dt[us_temp_svi_dt$SVI_cat %in% c("SVI Q5 (High Vuln)", "SVI Q1 (Low Vuln)"), ]$tasmax95_change, var_name = "95th Percentile Max Temp")

pop_haz_plot_ssp_cdd <- boxplot_ssp(var_in = us_temp_svi_dt[us_temp_svi_dt$SVI_cat %in% c("SVI Q5 (High Vuln)", "SVI Q1 (Low Vuln)"), ]$PDD_24_change_annual, var_name = "Person Degree Days")
pop_haz_plot_ssp_hd_40 <- boxplot_ssp(var_in = us_temp_svi_dt[us_temp_svi_dt$SVI_cat %in% c("SVI Q5 (High Vuln)", "SVI Q1 (Low Vuln)"), ]$PHD_HI40_change, var_name = "Person Hot Days (HI40)")
pop_haz_plot_ssp_tmax <- boxplot_ssp(var_in = us_temp_svi_dt[us_temp_svi_dt$SVI_cat %in% c("SVI Q5 (High Vuln)", "SVI Q1 (Low Vuln)"), ]$PD95_change, var_name = "Person Degrees")

# Combine plots
#
plot_all <- ggpubr::ggarrange(haz_plot_ssp_cdd, haz_plot_ssp_hd_40, haz_plot_ssp_tmax,
                      pop_haz_plot_ssp_cdd, pop_haz_plot_ssp_hd_40, pop_haz_plot_ssp_tmax, 
                      nrow = 2, ncol = 3, common.legend = TRUE, legend = "right")

ggsave(paste0(boxplot_out, "boxplot_crossmod_alltract.png"),
       plot_all, width = 12000, height = 6000,
       dpi = 350, units = "px", limitsize = FALSE)


#################### CROSS-MODEL SUMMARY #######################################
# POPULATION HEAT HAZARD PLOTS - SUM ACROSS MODEL/SSP
#
us_temp_svi_dt_sum <- us_temp_svi_dt %>%
  filter(SVI_cat %in% c("SVI Q5 (High Vuln)", "SVI Q1 (Low Vuln)")) %>%
  dplyr::group_by(SVI_cat, mod, ssp) %>%
  dplyr::mutate(pop_prod_CDD_annual_change = CDD_change_annual * Pop20 / sum(Pop20, na.rm=T),
         pop_prod_hotday_change = hotdays_change_annual * Pop20 / sum(Pop20, na.rm=T),
         pop_prod_tmax95_change = tasmax95_change * Pop20 / sum(Pop20, na.rm=T))%>%
  dplyr::summarise(pop_mean_CDD_annual_change = sum(pop_prod_CDD_annual_change, na.rm=T),
            pop_mean_hotday_annual_change = sum(pop_prod_hotday_change, na.rm=T),
            pop_mean_tmax95_change = sum(pop_prod_tmax95_change, na.rm=T),
            sum_PDD_change = sum(PDD_change_annual, na.rm = T),
            sum_PHD_change = sum(PHD_change_annual, na.rm = T),
            sum_PD95_change = sum(PD95_change_annual, na.rm = T))


boxplot_ssp_nation <- function(var_in, var_name, var_abr) {
  
  ggplot(data = us_temp_svi_dt_sum, aes(y = ssp, x = var_in, fill = ssp, col = SVI_cat)) +
    geom_boxplot(outliers = FALSE, show.legend = FALSE) +
    theme_classic(base_size = 15) +
    scale_fill_manual(values = c("#fdffe6", "#feb39b", "#b94250"))+
    scale_color_manual(values = c("darkred", "#493657"))+
    labs(title = paste0("National Mean ", var_name, " Distribution \nAll CMIP6 Models"),
         x = paste0("Population-Weighted Mean Change in ", var_name),
         y = "Shared Socioeconomic Pathway",
         fill = "SSP",
         col = "Tract SVI Category") 
  
}

# NATIONWIDE PLOTS - ALL MODELS and SSPs
# Create plots
#
haz_plot_ssp_cdd_sum <- boxplot_ssp_nation(var_in = us_temp_svi_dt_sum$pop_mean_CDD_annual_change, var_name = "Cooling Degree Days", var_abr = "CDDs")
haz_plot_ssp_hd_sum <- boxplot_ssp_nation(var_in = us_temp_svi_dt_sum$pop_mean_hotday_annual_change, var_name = "Hot Days", var_abr = "Hot Days")
haz_plot_ssp_tmax_sum <- boxplot_ssp_nation(var_in = us_temp_svi_dt_sum$pop_mean_tmax95_change, var_name = "95th Percentile Temp", var_abr = "95th Percentile Max Temp")

pophaz_plot_ssp_cdd_sum <- boxplot_ssp_nation(var_in = us_temp_svi_dt_sum$sum_PDD_change, var_name = "Person Degree Days", var_abr = "Person Degree Days")
pophaz_plot_ssp_hd_sum <- boxplot_ssp_nation(var_in = us_temp_svi_dt_sum$sum_PHD_change, var_name = "Person Hot Days", var_abr = "Person Hot Days")
pophaz_plot_ssp_tmax_sum <- boxplot_ssp_nation(var_in = us_temp_svi_dt_sum$sum_PD95_change, var_name = "Person Degrees", var_abr = "Person Degrees")

# Combine plots
#
plot_all_sum <- ggarrange(haz_plot_ssp_cdd_sum, haz_plot_ssp_hd_sum, haz_plot_ssp_tmax_sum,
                          pophaz_plot_ssp_cdd_sum, pophaz_plot_ssp_hd_sum, pophaz_plot_ssp_tmax_sum, 
                      nrow = 2, ncol = 3)

ggsave(paste0(boxplot_out, "boxplot_crossmod_summeas.png"),
       plot_all_sum, width = 8000, height = 3000,
       dpi = 350, units = "px", limitsize = FALSE)
