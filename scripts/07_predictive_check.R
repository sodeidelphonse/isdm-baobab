# Copyright (c) 2026 SODE Akoeugnigan Idelphonse
# Licensed under the MIT License (see LICENSE file for details)
 

library(inlabru)
library(DHARMa)
library(isdmtools)

bru_options_set(control.compute = list(dic = TRUE, config = TRUE)) 

#----------------------------------------------------------
#--- Posterior predictive check for integrated models  
#----------------------------------------------------------

# Base model components 
jcmp0 <- ~ -1 + area(log(abund_utm$area), model = "offset")+
  bio1_wc30s(covariates_pc$bio1_wc30s, model = "linear") + 
  srtm_slope(covariates_pc$srtm_slope, model = "linear") + 
  SLTPPT_d2(covariates_pc$SLTPPT_d2, model = "linear") + 
  CLYPPT_d6(covariates_pc$CLYPPT_d6, model = "linear") + 
  bio14_wc30s(covariates_pc$bio14_wc30s, model = spde_bio14)+
  Presence_inter(1) + Count_inter(1)

pcmatern <- inla.spde2.pcmatern(mesh, 
                                prior.range = c(10, 0.01), 
                                prior.sigma = c(0.1, 0.01)  # final prior 
                               )

#--- Fit the shared spatial model 
jcmp1 <- update.formula(jcmp0, ~ . + spde(main = geometry, model = pcmatern))

# LGCP component
obs_pp <- bru_obs(
  formula = geometry ~ bio1_wc30s + bio14_wc30s + srtm_slope + SLTPPT_d2 +
                       CLYPPT_d6 + Presence_inter + spde, 
  family = "cp", 
  data = point_utm,
  domain = list(geometry = mesh),
  samplers = ben_utm 
)

# Poisson likelihood 
obs_pois <- bru_obs(
  formula = counts ~ bio1_wc30s + bio14_wc30s + srtm_slope + SLTPPT_d2 + CLYPPT_d6 + 
                      area + Count_inter + spde, 
  family = "poisson", 
  data = abund_utm
)

# Integrated model with Poisson likelihood 
system.time(jfit_pois <- bru(jcmp1, obs_pois, obs_pp, 
                              options = list(control.inla = list(int.strategy = "eb"), 
                                             bru_max_iter = 20))
            )
summary(jfit_pois) 

# Negative binomial (NB) likelihood 
obs_nb <- bru_obs(
  formula = counts ~ bio1_wc30s + bio14_wc30s +srtm_slope +SLTPPT_d2 + CLYPPT_d6 +
                    Count_inter +area +spde, 
  family = "nbinomial", 
  control.family = list(variant = 1),
  data = abund_utm
)

# Integrated model with NB likelihood 
system.time(jfit_nb <- bru(jcmp1, obs_nb, obs_pp, 
                            options = list(control.inla = list(int.strategy = "eb"),
                                           bru_max_iter = 20
                                           ))
            )
summary(jfit_nb) 

# Difference in DIC: The Poisson family fits the data well.
jfit_pois$dic$dic - jfit_nb$dic$dic


#--- Residuals analysis for the fitted models -----

# Generate 1000 posterior samples for each family of count observation part
set.seed(234)
samples <- lapply(list(pois = jfit_pois, nb = jfit_nb),  
                  function(model) {
                    generate(model, newdata = abund_utm, 
                             formula = ~ exp(Count_inter +spde +bio1_wc30s +bio14_wc30s +
                                            srtm_slope +SLTPPT_d2 +CLYPPT_d6)*abund_utm$area[1],
                             n.samples = 1000)
                  })

#--- a) Check the Poisson family with replicated data
set.seed(234)
ppc_pois <- simulate_replicates(samples$pois, family = "poisson")

res_dharm_pois <- createDHARMa(simulatedResponse = ppc_pois,    
                               observedResponse = abund_utm$counts,        
                               fittedPredictedResponse = apply(samples$pois, 1, median), # or mean
                               integerResponse = TRUE)

# Figure A2
jpeg(file = "figures/fig_A2_resid_shared_pois.jpeg", width = 800, height = 400)
par(mfrow = c(1,2))
plotQQunif(res_dharm_pois)    
testDispersion(res_dharm_pois)
dev.off()

# No spatial autocorrelation present
testSpatialAutocorrelation(res_dharm_pois, 
                           x = st_coordinates(abund_utm)[,"X"], 
                           y = st_coordinates(abund_utm)[,"Y"], 
                           plot = FALSE)

#--- b) Check the negative binomial family with replicated data
set.seed(234)
ppc_nb <- simulate_replicates(samples$nb, family = "nbinomial", 
                              jfit_nb$summary.hyperpar[1,"mean"]
                              )

res_dharm_nb <- createDHARMa(simulatedResponse = ppc_nb, 
                             observedResponse = abund_utm$counts,    
                             fittedPredictedResponse = apply(samples$nb, 1, median), # or mean 
                             integerResponse = TRUE)

# Residual plots for the NB likelihood
par(mfrow = c(1,2))
plotQQunif(res_dharm_nb)    
testDispersion(res_dharm_nb)

# Spatial autocorrelation remains in the residuals
testSpatialAutocorrelation(res_dharm_nb, 
                           x = st_coordinates(abund_utm)[,"X"], 
                           y = st_coordinates(abund_utm)[,"Y"], 
                           plot = FALSE)

#------------------------------ END -------------------------------------------
