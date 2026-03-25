# Copyright (c) 2026 SODE Akoeugnigan idelphonse
# Licensed under the MIT License (see LICENSE file for details)
 

library(inlabru)
library(PointedSDMs)

bru_options_set(control.compute = list(dic = TRUE, config = TRUE)) 

#-----------------------------------
#--- 1) The shared latent model
#-----------------------------------

# Instantiate the shared component model
model_shared <- startISDM(data = dataset, 
                          Mesh = mesh,  
                          Projection = proj, 
                          responseCounts = 'counts',
                          spatialCovariates = covariates_pc, 
                          pointsSpatial = "shared",  
                          Offset = "area",
                          Boundary = ben_utm
                          )

# Define a PC prior for the shared latent (SPDE)
model_shared$specifySpatial(sharedSpatial = TRUE, 
                            PC = TRUE,
                            prior.range = c(10, 0.01),
                            prior.sigma = c(0.1, 0.01)
                           )

# Specify a nonlinear effect 
model_shared$changeComponents("bio14_wc30s(main = bio14_wc30s, model = spde_bio14)")

# Fit the model
spat_shared <- fitISDM(model_shared, 
                       options = list(control.inla = list(int.strategy ='eb'))
                      )
summary(spat_shared)


#------------------------------------------------
#--- 2) Bayesian LGCP model for presence-only 
#------------------------------------------------

model_pp <- startISDM(data = list(Presence = point_utm), 
                      Mesh = mesh,  
                      Projection = proj, 
                      spatialCovariates = covariates_pc,  
                      pointsSpatial = "shared", 
                      Boundary = ben_utm
                      )

model_pp$specifySpatial(sharedSpatial = TRUE, 
                        PC = TRUE,
                        prior.range = c(10, 0.01),
                        prior.sigma = c(0.1, 0.01)
                       )

model_pp$changeComponents("bio14_wc30s(main = bio14_wc30s, model = spde_bio14)")

spat_pp <- fitISDM(model_pp, 
                   options = list(control.inla = list(int.strategy ='eb', diagonal=0.1))  
                  ) 
summary(spat_pp)


#---------------------------------------------------
#--- 3) Bayesian Poisson model for abundance data
#---------------------------------------------------

model_pc <- startISDM(data = list(Count = abund_utm), 
                      Mesh = mesh,  
                      Projection = proj, 
                      responseCounts = 'counts',
                      spatialCovariates = covariates_pc, 
                      pointsSpatial = "shared", 
                      Offset = "area",
                      Boundary = ben_utm
                      )

model_pc$specifySpatial(sharedSpatial = TRUE, 
                        PC = TRUE,
                        prior.range = c(10, 0.01),
                        prior.sigma = c(0.2, 0.01)  
                        )
model_pc$changeComponents("bio14_wc30s(main = bio14_wc30s, model = spde_bio14)")

spat_pc <- fitISDM(model_pc, 
                   options = list(control.inla = list(int.strategy ='eb'))
                  ) 
summary(spat_pc)



