# Copyright (c) 2026 SODE Akoeugnigan Idelphonse
# Licensed under the MIT License (see LICENSE file for details)
 

# Required packages
library(sf)
library(terra)
library(ggplot2)
library(dplyr)  
library(INLA)
library(fmesher)
library(patchwork)
library(GGally)
library(RColorBrewer)

#-----------------------------------------
#--- 1) The observations and vector data
#-----------------------------------------

# Color palette across maps
col_grad <- colorRampPalette(c('blue', 'cyan', 'yellow', 'red'))(200)

# Projection system: WGS84, Zone 31N and kilometer units
proj <- "+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

#--- Vector data
ben_utm <- read_sf("data/shapefile/BEN_adm0.shp") |> 
  st_transform(crs = proj)

zc_utm <- read_sf("data/shapefile/Clim_zone_merge.shp") |>
  st_transform(crs = proj)

#-- Clean presence-only data
point_df <- read.csv("data/Adansonia_occurrence.csv")
sum(duplicated(point_df))

# Project data into UTM for distance calculation
point_utm <- st_as_sf(x = point_df, 
                      coords = c("long", "lat"), 
                      crs = 4326) |> 
  st_transform(crs = proj)
nrow(point_utm)

#-- Counts data
abund_df <- read.csv("data/Adansonia_abundance.csv")
sum(duplicated(abund_df))

abund_utm <- abund_df |> 
  dplyr::distinct() |> 
  st_as_sf(coords = c("long", "lat"), crs = 4326) |> 
  st_transform(crs = proj)
dim(abund_utm)


#------------------------------
#--- 2.) Environmental data 
#------------------------------

# Import the selected covariates as raster files
r_tmp <- rast("data/covariates/bio1_wc30s.tif")
r_files <- list.files("data/covariates", pattern ='tif', full.names = TRUE)

# Align the list of covariates
r_sc <- lapply(r_files, function(f) {
  r <- rast(f)
  resample(r, r_tmp, method = "bilinear")
})

# Project and standardise outputs
r_sc <- rast(r_sc) |> 
  project(proj) |> 
  scale()
names(r_sc)

#---------------------------------------------
#--- 3) Data preparation for model fitting
#---------------------------------------------

# Datasets 
dataset <- list(Presence = point_utm, Count = abund_utm)

# We keep a standardised order for the selected covariates 
vars_pc <- c("bio1_wc30s", "bio14_wc30s", "srtm_slope", "SLTPPT_d2", "CLYPPT_d6")  
covariates_pc <- r_sc[[vars_pc]]

# For the scenario of dataset-specific covariates effects
vars_pp <- paste0(vars_pc, "_pp") 
covariates <- c(covariates_pc, covariates_pc)
names(covariates) <- c(vars_pc, vars_pp)
rm(r_files, r_tmp)

#--- Build the mesh for the 2D SPDE -----

# Search the optimal mesh parameters based on the Belmont (2022) tutorial 
xrange <- diff(range(st_coordinates(point_utm)[,"X"]))
yrange <- diff(range(st_coordinates(point_utm)[,"Y"]))

(max_edge <- min(xrange, yrange)/(3*5))
(bnd_outer <- min(xrange, yrange)/3)  # offset (domain extension)
ceiling(max_edge)*c(1, 3)  # max.edge (the outer layer has lower triangles density)
ceiling(max_edge)*1/2      # cutoff (avoid too many triangles around clustered points)

# The final mesh validated for model fitting (see script 02)
bndr <- fm_as_segm(ben_utm)
mesh <- fm_mesh_2d(
  boundary = bndr,
  max.edge = c(22, 64),   # alternative tested: max.edge = c(20, 60), c(20, 40) 
  offset = c(1e-3, 100),  # The 1e-3 has no effect on the triangles density
  cutoff = 10,            # alternative tested: cutoff = ceiling(max_edge/5)
  crs = proj
)
mesh$n

# Figure A3
plot(mesh)
png(file = "figures/fig_A3_mesh.jpeg", width = 500, height = 800)
plot(mesh)
dev.off()

#--- 1D SPDE for nonlinear component ----
x_bio14 <- seq(-1, 6, length = 15)     
mesh1D  <- fm_mesh_1d(x_bio14, boundary = "free", degree = 2) 
spde_bio14 <- inla.spde2.pcmatern(mesh = mesh1D,
                                  alpha = 2,
                                  constr = TRUE,
                                  prior.range = c(1, 0.01),
                                  prior.sigma = c(1, 0.01)
                                 )

# Figure A4 
p_m1 <- ggplot() + geom_fm(data = mesh1D)
ggsave("figures/fig_A4_mesh_bio14.jpeg",  p_m1, width = 5, height = 4, dpi = 300)
rm(p_m1)

#--- Figure A5 (map of selected covariates)
covar_final <- covariates_pc
names(covar_final) <- c("Annual temperature", "Rainfall of driest month", "SRTM slope", 
                        "Soil texture silt fraction", "Soil texture clay fraction")

png(file = "figures/fig_A5_covar_map.jpeg", width = 600, height = 600)
plot(covar_final, fun = function() lines(vect(ben_utm)))
dev.off()
rm(covar_final)

#--- Prepared count data with the covariates for EDA
data_abund <- extract(r_sc[[vars_pc]], abund_utm) |>  
  cbind(abund =  abund_utm$abund, counts = abund_utm$counts, area = abund_utm$area, 
        st_coordinates(abund_utm)[, c("X", "Y")]
  ) 
table(complete.cases(data_abund))
head(data_abund)

#--- Correlation between the final covariates
p_cor <- data_abund |> 
  dplyr::select(bio1_wc30s, bio14_wc30s, srtm_slope, CLYPPT_d6, SLTPPT_d2) |>
  ggpairs(upper = list(continuous ='points'), 
                  lower = list(continuous ='cor')
          )

# Figure 2 (correlation plot)
ggsave("figures/fig2_covar_corr.jpeg",  p_cor, width = 7, height = 6, dpi = 300)
rm(p_cor)
