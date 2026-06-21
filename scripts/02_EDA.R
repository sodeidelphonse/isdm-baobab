# Copyright (c) 2026 SODE Akoeugnigan Idelphonse
# Licensed under the MIT License (see LICENSE file for details)
 

#--------------------------------------
#--- Exploratory data analysis (EDA)
#--------------------------------------

library(geoR)
library(spatstat.explore) 
library(spatstat.model)
library(ggspatial)
library(INLA)
library(inlabru)
library(isdmtools)  # Install the version v0.4.0 from the binary/source files or from GitHub
library(MASS)
library(spaMM)
library(DHARMa)

source("scripts/08_utils.r")

#--------------------------------
#--- 1) Plot of raw data sets 
#--------------------------------

# Bind datasets
xy_all <- bind_datasets(list(Presence = point_utm, Abundance = abund_utm)) |> 
  dplyr::select(datasetName, geometry) |>
  dplyr::mutate(Dataset = datasetName)

# Map of all point locations 
po <- ggplot() +
  geom_sf(data = ben_utm, fill = "NA") +
  geom_sf(data = xy_all, aes(colour = Dataset), size = 1.5) +  
  labs(x ="Longitude", y ="Latitude") +
  scale_x_continuous(breaks = seq(0, 4, 1)) +
  scale_y_continuous(breaks = seq(6, 13, 2)) +
  theme_bw(base_size = 14) +  
  annotate("text",  x = 580, y = 1400, label = "(a)", size = 5, fontface ="bold") +
  ggspatial::annotation_north_arrow(location = "tl", height = unit(1, "cm"), width = unit(0.7, "cm")) +
  ggspatial::annotation_scale(location = "br", bar_cols = c("grey60", "white")) +
  ggplot2::geom_sf(data = zc_utm, fill = NA, color = "grey20", linewidth = 0.3)

# Map of the abundance
pc <- ggplot() +
  geom_sf(data = ben_utm, fill = "NA") +
  geom_sf(data = abund_utm, aes(colour = abund), size = 2) + 
  scale_colour_gradientn(colours = col_grad, name = expression("# trees/km"^2)) +
  labs(x ="Longitude", y ="Latitude") +
  scale_x_continuous(breaks = seq(0, 4, 1)) +
  scale_y_continuous(breaks = seq(6, 13, 2)) +
  theme_bw(base_size = 14) +   
  annotate("text", x = 580, y = 1400, label = "(b)", size = 5, fontface ="bold") +
  ggspatial::annotation_north_arrow(location = "tl", height = unit(0.8, "cm"), width = unit(0.6, "cm")) +
  ggspatial::annotation_scale(location = "br", bar_cols = c("grey60", "white"))+
  ggplot2::geom_sf(data = zc_utm, fill = NA, color = "grey20", linewidth = 0.3)

# Save the combined plots (Figure 1)
(p_data <- po + pc)
ggsave("figures/fig1_map_datasets.jpeg", plot = p_data, width = 8, height = 6, dpi = 300)
rm(p_data, po, pc)


#------------------------------------------
#--- 2) Variogram of relative abundance
#------------------------------------------

# Compute the empirical variogram 
v_emp <- variog(coords = st_coordinates(abund_utm), 
             data = abund_utm$abund, option = "bin", 
             lam = 0.5,    # Box-Cox transformation
             uvec = seq(0, 200, by = 10)
             )

# Monte Carlo envelope
set.seed(231)
mc_env <- variog.mc.env(coords = st_coordinates(abund_utm), 
                     data = abund_utm$abund, 
                     obj = v_emp, 
                     nsim = 999, 
                     save.sim = TRUE
                     )

vario_df <- data.frame(dist = v_emp$u,       # Distance
                       semi = v_emp$v,       # semivariance
                       lo = mc_env$v.lower,  # Lower bound
                       hi = mc_env$v.upper   # Upper bound 
                       )  

# Fit variogram models to the abundance data
ini_vals <- expand.grid(seq(0, 2, l = 5), seq(0, 2, l = 5))

v_mat <- variofit(v_emp, ini = ini_vals, fix.nug = TRUE, wei = "equal", 
                  cov.model = "matern", kappa = 1)   # Matern covariance

v_exp <- variofit(v_emp, ini = ini_vals, fix.nug = TRUE, wei = "equal", 
                  cov.model = "matern", kappa = 0.5) # Exponential

# Create data frame to customize the variogram chart
dist_vals <- v_emp$u  
vm_mat <- (v_mat$cov.pars[1] + v_mat$nugget) - 
       cov.spatial(dist_vals, cov.pars = v_mat$cov.pars, 
              kappa = v_mat$kappa, cov.model = v_mat$cov.model)

vm_exp <- (v_exp$cov.pars[1] + v_exp$nugget) - 
       cov.spatial(dist_vals, cov.pars = v_exp$cov.pars, 
              kappa = v_exp$kappa, cov.model = v_exp$cov.model)

vm_df <- data.frame(dist = rep(dist_vals, 2), 
                    semivariance = c(vm_mat, vm_exp),
                    model = rep(c("Mat", "Exp"), each = length(dist_vals))
                    )

# Visualisation
plot_v <- ggplot(vario_df, aes(x = dist)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.5) +  
  geom_point(aes(y = semi), color = "black", shape = 1, size = 1.8) +  
  geom_line(data = vm_df, aes(x = dist, y = semivariance, color = model, linetype = model), 
            linewidth = 1) +
  annotate("text", x = min(vario_df$dist), y = max(vario_df$hi), label = "(a)", 
           hjust = 0, fontface = "bold", size = 5) +  
  labs(x = "Distance [km]", y = "Semivariance", color = "Model", linetype = "Model") +
  scale_color_manual(values = c("Mat" = "black", "Exp" = "blue"), 
                     labels = c("Matérn", "Exponential"), name = NULL) +
  scale_linetype_manual(values = c("Mat" = "solid", "Exp" = "dashed"), 
                        labels = c("Matérn", "Exponential"), name = NULL) +
  theme_bw(base_size = 14) + 
  theme(legend.position = c(0.8, 0.90),  
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent"))
plot_v


#------------------------------------------------------
#--- 3) K-Function and K envelope for point pattern
#------------------------------------------------------

# Create a point pattern object
win <- owin(xrange = range(st_coordinates(point_utm)[,"X"]), 
            yrange = range(st_coordinates(point_utm)[,"Y"]))
pp <- as.ppp(X = st_coordinates(point_utm), W = win)

# Simulate the K envelope
set.seed(123)
E <- envelope(pp, Kinhom, nsim = 999,  r = seq(0, 100, length.out = 50)) 

E_df <- data.frame(
  r = E$r,           
  theo = E$theo,     # Theoretical K-function
  obs = E$obs,       # Observed K-function
  lo = E$lo,         # Lower envelope
  hi = E$hi          # Upper envelope
)

plot_k <- ggplot(E_df, aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey60", alpha = 0.5) +            
  geom_line(aes(y = theo, color = "theo"), linetype = "dashed", linewidth = 0.8) +  
  geom_line(aes(y = obs, color = "obs"),  linetype = "solid", linewidth = 0.8) +    
  annotate("text", x = min(E_df$r), y = max(E_df$hi), label = "(b)", 
           hjust = 0, fontface = "bold", size = 5) +  
  scale_color_manual(
    values = c("theo" = "red", "obs" = "black"),  
    labels = c("theo" = bquote(K[inhom]^{theo} * (r)), 
               "obs"  = bquote(K[inhom]^{obs} * (r))), name = NULL) +
  labs(x = "Distance [km]", y = expression(K[inhom](r))) +
  theme_bw(base_size = 14) + 
  theme(legend.position = c(0.72, 0.90),  
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent"))
plot_k

# Combined plot of K-function and variogram (Figure 3)
(plot_vk <- plot_v + plot_k)
ggsave("figures/fig4_Kinhom_variog.jpeg", plot = plot_vk, width = 9, height = 5, dpi = 300)
rm(plot_v, plot_k, plot_vk)


#----------------------------------------------------------
#--- 4) Relation between LGCP covariance and K-function
#----------------------------------------------------------

# Note: Until now, we have been identifying a "baseline" spatial scale without  
# making any environmental assumptions. This also facilitates the comparison of the
# estimated spatial range with the distance from origin to the crossing point of
# the K_{inhom} function with its theoretical counterpart.

#--- A) Bayesian LGCP -----
bru_options_set(control.compute = list(dic = TRUE, config = TRUE)) 

#--- Some initial mesh configurations
bndr <- fm_as_segm(ben_utm)

mesh2d <- fm_mesh_2d(
  boundary = bndr,
  max.edge = ceiling(max_edge*c(1,3)),  # tested max_edge*c(1,2)  
  offset = c(1e-3, 100),  
  cutoff = ceiling(max_edge/4),         # tested max_edge/5, max_edge/2
  crs = proj
)
mesh2d$n

# The Matern model for spatial latent 
pcmatern0 <- inla.spde2.pcmatern(mesh2d, 
                                prior.range = c(10, 0.01), 
                                prior.sigma = c(0.1, 0.01)  # tried sigma_0 = 0.3, 0.5, 1  
                               )

# Model components and observation model
cmp_po <- ~ -1 + Presence_inter(1) + spde(geometry, model = pcmatern0)

obs_po <- bru_obs(
  family = "cp", 
  formula = geometry ~ Presence_inter + spde, 
  data = point_utm,
  domain = list(geometry = mesh2d),
  samplers =  ben_utm
)

# Fit the initial LGCP model with the empirical Bayes integration  
fit_lgcp_spde0 <- bru(cmp_po, obs_po, 
                     options = list(control.inla = list(int.strategy = "eb"))
                     )
summary(fit_lgcp_spde0)


#--- The final mesh configuration 
mesh <- fm_mesh_2d(
  boundary = bndr,
  max.edge = c(22, 64),   
  offset = c(1e-3, 100),  
  cutoff = 10,            
  crs = proj
)
mesh$n

pcmatern <- inla.spde2.pcmatern(mesh, 
                                prior.range = c(10, 0.01), 
                                prior.sigma = c(0.1, 0.01)   
                                )

cmp_po <- ~ -1 + Presence_inter(1) + spde(geometry, model = pcmatern)

obs_po <- bru_obs(formula = geometry ~ Presence_inter + spde, 
                  family = "cp", 
                  data = point_utm,
                  domain = list(geometry = mesh), 
                  samplers = ben_utm
                 )

fit_lgcp_spde <- bru(cmp_po, obs_po, 
                options = list(control.inla = list(int.strategy = "eb"))
                )
summary(fit_lgcp_spde)

#--- Difference of DIC between mesh configurations -----
# A configuration with the lower DIC is preferable
fit_lgcp_spde0$dic$dic - fit_lgcp_spde$dic$dic


#--- Credible intervals for Bayesian LGCP covariance parameters ------
fit_lgcp_spde$summary.hyperpar

# StDev and variance of spatial signal (checked)
marg_sd <- fit_lgcp_spde$marginals.hyperpar$`Stdev for spde`
inla.qmarginal(c(0.025, 0.5, 0.975), marg_sd)   

marg_var <- inla.tmarginal(function(x) x^2, marg_sd)
inla.qmarginal(c(0.025, 0.5, 0.975), marg_var)  

# Kappa
sqrt(8)/fit_lgcp_spde$summary.hyperpar["Range for spde", "mean"]
sqrt(8)/fit_lgcp_spde$summary.hyperpar["Range for spde", "0.975quant"]
sqrt(8)/fit_lgcp_spde$summary.hyperpar["Range for spde", "0.025quant"]

# 10% Practical range 
solve_practical_range(param_val = fit_lgcp_spde$summary.hyperpar["Range for spde", "mean"], 
                      nu = 1, engine = "inla", thresh = 0.10)

# 95% CI for the 10% practical range
solve_practical_range(param_val = fit_lgcp_spde$summary.hyperpar["Range for spde", "0.025quant"], 
                      nu = 1, engine = "inla", thresh = 0.10)

solve_practical_range(param_val = fit_lgcp_spde$summary.hyperpar["Range for spde", "0.975quant"], 
                      nu = 1, engine = "inla", thresh = 0.10)


#--- B) Relation between the LGCP covariance and K-function ---------

# Fit the classic LGCP to the point pattern
fit_lgcp_spat <- kppm(pp, cluster = "LGCP", model = "matern", nu = 1)
fit_lgcp_spat

dist_vals  <- seq(0, 100, length.out = 50)

# 1.) Normalized covariance from the frequentist LGCP 
cor_freq <- std_matern_corr(fit_lgcp_spat, "kppm", dist_vals)
summary(cor_freq$pair_cor_sc)

# 2.) Normalized covariance function from the Bayesian LGCP 
cor_bayes <- std_matern_corr(fit_lgcp_spde, "inla", r = dist_vals) 
summary(cor_bayes$pair_cor_sc)

# 3.) Compute the normalized empirical inhomogeneous K-function 
Ki <- Kinhom(pp,  r = dist_vals)  
Ki$iso_sc  <- (Ki$iso - min(Ki$theo, na.rm = T))/diff(range(Ki$theo, na.rm =T))
Ki$theo_sc <- (Ki$theo - min(Ki$theo, na.rm = T))/diff(range(Ki$theo, na.rm =T))
summary(Ki)

# 4.) Compute the normalized empirical pairwise correlation: K'(r)/(2r*pi)
pcf_k <- pcf(Ki, spar = 1, method = "b", r = dist_vals) 
pcf_k$g_sc <- (pcf_k$pcf - min(pcf_k$pcf, na.rm = T))/diff(range(pcf_k$pcf, na.rm =T, finite =T))
pcf_k$theo_sc <- (pcf_k$theo - min(pcf_k$pcf, na.rm = T))/diff(range(pcf_k$pcf, na.rm =T, finite =T))
summary(pcf_k)

# 5.) Visualization
range_k  <- 45 
plot_cor <- ggplot() +
  geom_line(aes(x = cor_freq$dist_vals, y = cor_freq$pair_cor_sc, color = "C(r) lgcp"), linewidth = 1) +
  geom_line(aes(x = cor_bayes$dist_vals, y = cor_bayes$pair_cor_sc, color = "C(r) Bayes"), linewidth = 0.8) +
  geom_line(aes(x = pcf_k$r, y = pcf_k$g_sc, color = "g(r)"), linewidth = 0.8, linetype ="dashed") + 
  geom_line(aes(x = Ki$r, y = Ki$iso_sc, color ="K(r)"), linewidth = 0.8) +
  geom_line(aes(x = Ki$r, y = Ki$theo_sc, color ="Kth"), linewidth = 0.8, linetype ="dashed") +
  geom_vline(xintercept = range_k, linetype = "dotted", linewidth = 0.5, color = "grey40") +  
  geom_vline(xintercept = cor_bayes$rho, linetype = "dotted", linewidth = 0.5, color = "grey40") + 
  
  labs(x = "Distance r (km)", y = "Normalized value", color = NULL) +
  scale_color_manual(values = c("C(r) Bayes"= "black", "C(r) lgcp" = "green", "g(r)" = "black", "g(r) lgcp" = "orange", "K(r)" = "blue", "Kth" = "red"),
                     breaks = c("C(r) Bayes", "C(r) lgcp", "g(r)", "g(r) lgcp", "K(r)", "Kth"),
                     labels = c(bquote(rho[Bayes-LGCP]*(r)), bquote(rho[Freq-LGCP]*(r)), bquote(g[inhom]^{obs}*(r)), bquote(g[inhom]^{lgcp}*(r)), 
                                bquote(K[inhom]^{obs}*(r)), bquote(K[inhom]^{theo}*(r)))
  ) +
  theme_grey(base_size = 14) +
  theme(legend.position = c(0.67, 0.80), legend.background = element_rect(fill ="transparent"))
plot_cor

# Figure 4
ggsave("figures/fig5_pair_corr.jpeg", plot_cor, width = 4.5, height = 4, dpi = 200)


#-----------------------------------------------------------------
#--- 5) Bootstrapping for estimating confidence intervals 
#-----------------------------------------------------------------

# The outputs from the analyses performed here have served to fill in the Table 2, 
# specifically the confidence intervals of classic/frequentist approaches.

#--- A) Simulation of covariance parameters for the variogram model

# Required data structure
dt <- as.geodata(data_abund, 
                 coords.col = c("X", "Y"),
                 data.col = which(colnames(data_abund) %in% c("abund", "counts")),
                 covar.col = vars_pc,
                 units.m.col = "area",
                 na.action = "ifdata"
                )

# Retrieve the estimated covariance parameters 
est_sigmasq <- v_mat$cov.pars[1]
est_phi <- v_mat$cov.pars[2]
est_tausq <- 0  

## Step 1: Generate 1000 spatial simulations based on the variogram model
set.seed(231)
nsim <- 1000
sim_grf <- grf(nrow(dt$coords), grid = dt$coords, nsim = nsim, 
               cov.model = "matern", kappa = 1,
               cov.pars = c(est_sigmasq, est_phi), 
               nugget = est_tausq, 
               lambda = 0.5)  # the same Box-Cox transformation 


## Step 2: Re-estimate covariance parameters for each simulation

# Pre-load replicates of parameters if they exist
if (file.exists("results/simulations_variofit.csv")) {
  message("Loading pre-calculated replicates from outputs folder...")
  sim_geor <- read.csv("results/simulations_variofit.csv") |> as.matrix() 
} else {
  message("No pre-calculated files found. Starting 1,000 simulations (this may take ~1 hour)...")
  sim_geor <- matrix(NA, nrow = nsim, ncol = 2)  
  colnames(sim_geor) <- c("sigmasq", "phi")

  system.time(
    for(i in 1:nsim) {
      temp_data <- dt
      temp_data$data <- sim_grf$data[,i]
      
      fit_sim <- likfit(temp_data, 
                        ini = c(est_sigmasq, est_phi), 
                        fix.kappa = TRUE, 
                        kappa = 1, 
                        lambda = 0.5,  
                        cov.model = "matern", messages = FALSE
      )
      sim_geor[i, ] <- c(fit_sim$sigmasq, fit_sim$phi)
    })
    # Save the simulated matrix as a backup and for further use
    write.csv(data.frame(sim_geor), "results/simulations_variofit.csv")
}

## Step 3: Calculate the 95% CI for each parameter

# 95% CI for sigma2
quantile(sim_geor[, "sigmasq"], probs = c(0.025, 0.975))

# 95% CI for phi 
quantile(sim_geor[, "phi"], probs = c(0.025, 0.975))

# 95% CI for kappa (1/phi)
quantile(1/sim_geor[, "phi"], probs = c(0.025, 0.975))

# 95% CI for the spatial range: rho = phi*sqrt(8)
quantile(sim_geor[, "phi"]*sqrt(8), probs = c(0.025, 0.975))

# The 10% Practical range 
solve_practical_range(param_val = est_phi, nu = 1, thresh = 0.1, engine = "geor")

# The 10% practical range for replicated data
prac_vec_geor <- sapply(sim_geor[,"phi"], function(p) {
  solve_practical_range(param_val = p, nu = 1, thresh = 0.1, engine = "geor")
})

# 95% CI estimates
quantile(prac_vec_geor, probs = c(0.025, 0.975))


#--- B) Simulation of LGCP covariance parameters ----

# Refit the classic LGCP if needed
if(!exists("fit_lgcp_spat")) {
  fit_lgcp_spat <- kppm(pp, cluster = "LGCP", model = "matern", nu = 1)
}

## Step 1: Simulate realisations from the fitted model

# settings
set.seed(234)
nsim <- 1000
sim_lgcp_list <- simulate(fit_lgcp_spat, nsim = nsim)

## Step 2: Re-fit the LGCP model to each simulated point pattern

if (file.exists("results/simulations_lgcp.csv")) {
  message("Loading pre-calculated replicates from outputs folder...")
  sim_lgcp <- read.csv("results/simulations_lgcp.csv") |> as.matrix() 
} else {
  message("No pre-calculated files found. Starting 1,000 simulations (this may take > 1 hour)...")
  sim_lgcp <- matrix(NA, nrow = nsim, ncol = 2) 
  colnames(sim_lgcp) <- c("sigma_sq", "alpha")
  
  # These simulations are highly time consuming !!!
   system.time(
    for(i in 1:nsim) { 
      fit_sim <- kppm(sim_lgcp_list[[i]], cluster = "LGCP", model = "matern", nu = 1)
      sim_lgcp[i, ] <- fit_sim$par 
    })
   # Save the simulated matrix as a backup or for further analyses
   write.csv(data.frame(sim_lgcp), "results/simulations_lgcp.csv")
}

## Step 3: Calculate the quantiles of derived parameters

# 95% CI for sigma2
quantile(sim_lgcp[, "sigma_sq"], probs = c(0.025, 0.975))

# 95% CI for alpha
quantile(sim_lgcp[, "alpha"], probs = c(0.025, 0.975))

# 95% CI for kappa
quantile(sqrt(8)/(sim_lgcp[, "alpha"] * 2), probs = c(0.025, 0.975))

# 95% CI for the spatial range: 2*alpha 
quantile(sim_lgcp[, "alpha"] * 2, probs = c(0.025, 0.975))

# 95% CI for the 10% practical range 
prac_vec_spat <- sapply(sim_lgcp[, "alpha"], function(p) {
  solve_practical_range(param_val = p,  nu = 1, thresh = 0.1, 
                        engine = "spatstat")
})

quantile(prac_vec_spat, probs = c(0.025, 0.975))

# The 10% practical range estimate
solve_practical_range(param_val = fit_lgcp_spat$par["alpha"], 
                 nu = 1, thresh = 0.10, engine = "spatstat"
                 )


#-------------------------------------------------------
#--- 6) Exploratory modelling for the abundance data
#-------------------------------------------------------

#--- A) Non-spatial models -----

# a) Poisson model 
glm_c <- glm(counts ~ bio1_wc30s +bio14_wc30s +srtm_slope +SLTPPT_d2 +CLYPPT_d6, 
             family = "poisson", offset = log(area), data = data_abund) 
summary(glm_c)

# Dispersion and Goodness-of-fit (GOF) test
glm_c$deviance/glm_c$df.residual  
1 - pchisq(glm_c$deviance, glm_c$df.residual) 

# Deviance explained
(glm_c$null.deviance - glm_c$deviance)/glm_c$null.deviance

# Simulation-based test (greater => overdispersion)
res_sims_pois <- simulateResiduals(glm_c)
plot(res_sims_pois)
testDispersion(res_sims_pois, alternative = "greater") 

# b) Negative binomial model 
glm_nb <- glm.nb(counts ~ bio1_wc30s +bio14_wc30s +srtm_slope +SLTPPT_d2 +CLYPPT_d6, 
                 link = "log", offset(log(area)), data = data_abund) 
summary(glm_nb)

# Dispersion & GOF test
glm_nb$deviance/glm_nb$df.residual
1 - pchisq(glm_nb$deviance, glm_nb$df.residual)  

# Deviance explained
(glm_nb$null.deviance - glm_nb$deviance)/glm_nb$null.deviance

# Simulation-based test
res_sims_nb <- simulateResiduals(glm_nb)
plot(res_sims_nb)
testDispersion(res_sims_nb, alternative = "greater")

# Check the model fit 
AIC(glm_c)
AIC(glm_nb)
BIC(glm_c)
BIC(glm_nb)

# The negative binomial family seems to correct the dispersion present in the data.  
# However, it shows a higher AIC and BIC values, likely due to additional model parameters. 
# Let's see what happen when we account for the unexplained spatial variation as 
# demonstrated by the variogram.

#--- B) Spatial GLMM with Matérn covariance -------

# a) Poisson model with Matern covariance
spamm_pois <- fitme(counts ~ bio1_wc30s +bio14_wc30s +srtm_slope +SLTPPT_d2 +CLYPPT_d6 + 
                            offset(log(area)) + Matern(1 | X + Y),
                    data = data_abund, 
                    family = poisson, 
                    fixed = list(nu = 1)
                    )
summary(spamm_pois)

# Dispersion - GOF - Deviance
deviance(spamm_pois)/df.residual(spamm_pois)
1-pchisq(deviance(spamm_pois), df.residual(spamm_pois))

spamm_pois0 <- update(spamm_pois, . ~ offset(log(area)))
(deviance(spamm_pois0) - deviance(spamm_pois))/deviance(spamm_pois0)

# Residuals diagnostics
set.seed(1234)
res_spamm <- simulateResiduals(spamm_pois, n = 1000) 
testResiduals(res_spamm)

# Figure A1
jpeg(file = "figures/fig_A1_spam_pois.jpeg", width = 800, height = 400)
par(mfrow = c(1,2))
plotQQunif(res_spamm)    
testDispersion(res_spamm)
dev.off()

# b) Negative binomial model with Matérn covariance
spamm_nb <- fitme(counts ~ bio1_wc30s +bio14_wc30s +srtm_slope +SLTPPT_d2 +CLYPPT_d6 + 
                    offset(log(area)) + Matern(1 | X + Y),
                  data = data_abund, 
                  family = spaMM::negbin, 
                  fixed = list(nu = 1),
                  )
summary(spamm_nb)

# Dispersion - GOF - deviance
deviance(spamm_nb)/df.residual(spamm_nb)
1-pchisq(deviance(spamm_nb), df.residual(spamm_nb))

spamm_nb0 <- update(spamm_nb, . ~ offset(log(area)))
(deviance(spamm_nb0) - deviance(spamm_nb))/deviance(spamm_nb0)

# Residuals diagnostics
set.seed(1234)
res_spamm_nb <- simulateResiduals(spamm_nb, n = 1000)
testResiduals(res_spamm_nb, plot=T)

# CONCLUSION: The two spatial models fit the observed abundance well. 
# However, we prefer the parsimonious model (i.e. the Poisson likelihood) for
# the Bayesian data integration framework
AIC(spamm_pois)
AIC(spamm_nb)

