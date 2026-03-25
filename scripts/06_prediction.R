# Copyright (c) 2026 SODE Akoeugnigan Idelphonse
# Licensed under the MIT License (see LICENSE file for details)
 

#-------------------------------------------------
#---- 1) Integrated Habitat Suitability Mapping 
#-------------------------------------------------

library(sf)
library(terra)
library(ggplot2)
library(patchwork)
library(PointedSDMs)
library(isdmtools)

source("scripts/08_utils.r")
fm_pix <- readRDS("data/pred_points.rds")

#-------------------------------
#--- (a) Integrated intensity 
#-------------------------------

# Baseline formula 
form_cov <- ~ bio1_wc30s +bio14_wc30s +srtm_slope +SLTPPT_d2 + CLYPPT_d6

# Prediction
jint <- predict(spat_shared, 
                data = fm_pix, 
                formula = update.formula(form_cov, ~. + shared_spatial +
                                           Presence_intercept),
                fun = "linear", n.samples = 1000, seed = 24)

jint <- format_predictions(jint)  
boxplot(jint$mean)

#--- log-intensity 
p <- generate_maps(jint , 
                    var_names = c("q0.025", "mean", "q0.975"), 
                    base_map = ben_utm,
                    color_gradient = col_grad, 
                    legend_title = bquote(log(lambda(x))),  
                    panel_labels = c("(a) 2.5% quantile", "(b) Mean", "(c) 97.5% quantile"),
                    xaxis_breaks = seq(0, 4, 1),
                    yaxis_breaks = seq(6, 13, 2)
                    )

# Figure 7 (log-intensity map)
(p <- p + ggplot2::geom_sf(data = zc_utm, fill = NA, color = "grey20", linewidth = 0.3))
ggsave("figures/fig7_shared_log-int.jpeg",  p, width = 8, height = 5, dpi = 300)

#--- Habitat suitability map
jt_prob <- suitability_index (jint, 
                              post_stat = c("q0.025", "mean", "q0.975"), 
                              output_format = "prob",
                              response_type = "joint.po",
                              has_offset = TRUE,
                              projection = proj)

p <- generate_maps(jt_prob, 
                    var_names = c("q0.025", "mean", "q0.975"), 
                    base_map = ben_utm,
                    color_gradient = col_grad, 
                    legend_title = , "suitability",  
                    panel_labels = c("(a) 2.5% quantile", "(b) Mean", "(c) 97.5% quantile"),
                    xaxis_breaks = seq(0, 4, 1),
                    yaxis_breaks = seq(6, 13, 2)
                    )

# Figure 8 (probability of presence)
p <- p + ggplot2::geom_sf(data = zc_utm, fill = NA, color = "grey20", linewidth = 0.3) 
print(p)
ggsave("figures/fig8_shared_suitability.jpeg",  p, width = 8, height = 5, dpi = 300)

#--- The shared spatial latent component 
jint_latent <- predict(spat_shared, data = fm_pix, formula = ~ shared_spatial, 
                    fun = "linear", n.samples = 1000, seed = 24)
jint_latent <- format_predictions(jint_latent)

# Visualisation 
p <- gg_map(jint_latent,
    vars_to_plot = c("mean", "sd", "IQR"),
    base_map = ben_utm,
    boundary_map = zc_utm,
    color_gradient = col_grad,
    x_axis_breaks = seq(0, 4, 1),
    y_axis_breaks = seq(6, 13, 2)
)

# Figure 9 (shared latent map)
p <- p + ggplot2::geom_sf(data = zc_utm, fill = NA, color = "grey20", linewidth = 0.3)
print(p)
ggsave("figures/fig9_shared_latent.jpeg",  p, width = 9, height = 7, dpi = 300)


#-------------------------------------------
#--- (b) Standalone count model (spat_pc) 
#-------------------------------------------

pred_pc <- predict(spat_pc, 
                  data = fm_pix, 
                  formula = update.formula(form_cov, ~. + shared_spatial + Count_intercept), 
                  fun = "linear", n.samples = 1000, seed = 24)

pred_pc <- format_predictions(pred_pc)
boxplot(pred_pc$mean)

resp_pc <- suitability_index(pred_pc, 
                            post_stat = c("q0.025", "mean", "q0.975"), 
                            output_format = "response", 
                            response_type = "count",
                            projection = proj
                            )

# Species relative abundance
p <- generate_maps (resp_pc, 
                    var_names = c("q0.025", "mean", "q0.975"), 
                    base_map = ben_utm,
                    color_gradient = col_grad, 
                    legend_title = "density", 
                    panel_labels = c("(a) 2.5% quantile", "(b) Mean", "(c) 97.5% quantile"),
                    xaxis_breaks = seq(0, 4, 1),
                    yaxis_breaks = seq(6, 13, 2)
                    )

# Figure B1 (density)
p <- p + ggplot2::geom_sf(data = zc_utm, fill = NA, color = "grey20", linewidth = 0.3)
print(p)
ggsave("figures/fig_B1_rel_abund.jpeg",  p, width = 7, height = 5, dpi = 300)


#-----------------------------------------
#--- (c) Standalone PO model (spat_pp) 
#-----------------------------------------

pred_po <- predict(spat_pp, 
                   data = fm_pix, 
                   formula = update.formula(form_cov, ~. + shared_spatial 
                                            + Presence_intercept), 
                   fun = "linear", n.samples = 1000, seed = 24)

pred_po <- format_predictions(pred_po)
boxplot(pred_po$mean)

prob_po  <- suitability_index(pred_po, 
                              post_stat = c("q0.025", "mean", "q0.975"), 
                              output_format = "prob",
                              response_type = "po", 
                              projection = proj)

# Probability of presence
p <- generate_maps(prob_po, 
                  var_names = c("q0.025", "mean", "q0.975"), 
                  base_map = ben_utm,
                  color_gradient = col_grad, 
                  legend_title = "suitability",  
                  panel_labels = c("(a) 2.5% quantile", "(b) Mean", "(c) 97.5% quantile"),
                  xaxis_breaks = seq(0, 4, 1),
                  yaxis_breaks = seq(6, 13, 2)
                  )

# Figure B2
p <- p + ggplot2::geom_sf(data = zc_utm, fill = NA, color = "grey20", linewidth = 0.3)
print(p)
ggsave("figures/fig_B2_po_suitability.jpeg",  p, width = 7, height = 5, dpi = 300)


#------------------------------------------------------
#--- 2) Prediction of response-specific contribution
#------------------------------------------------------

#--- a) Abundance data ----
jint_pc <- predict(spat_shared, 
                data = fm_pix, 
                covariates = names(covariates_pc), 
                spatial = TRUE,
                datasets = "Count",
                intercepts = TRUE,
                fun = "linear", 
                n.samples = 1000, 
                seed = 24) 
jint_ct <- format_predictions(jint_pc)

# On log scale 
var_names <- c("x", "y", "mean", "sd", "q0.025", "q0.5", "q0.975", "median")
jint_ct_r <- terra::rast(jint_ct[, var_names], type = "xyz", crs = proj) 
resp_val_jt <- extract(jint_ct_r, dataset$Count)

pred_ct_r <- terra::rast(pred_pc[, var_names], type = "xyz", crs = proj)
resp_val_ct <- extract(pred_ct_r, dataset$Count)

# Correlation (integrated counts, separate model counts)
cor.test(resp_val_jt$mean,  resp_val_ct$mean)

df_resp <- data.frame(resp_val_jt = resp_val_jt$mean, resp_val_ct = resp_val_ct$mean)
p_resp  <- ggplot(df_resp, aes(x = resp_val_jt, y = resp_val_ct)) +
  geom_point(alpha = 0.6, size = 0.6, color = "gray10") +
  geom_abline(intercept = 0, slope = 1, color = "blue", linewidth = 0.7, linetype = "dashed") + 
  annotate("text", x = -1, y = 3, label = "(b)", hjust = 0, fontface = "bold", size = 4) +
  labs(x = "log density (shared model)", y = "log density (count model)") +
  ylim(-1, 3) + xlim(-1, 3) + theme_bw() 
p_resp

#--- b) Presence-only ----
jint_po <- predict(spat_shared, 
                   data = fm_pix, 
                   covariates = names(covariates_pc), 
                   spatial = TRUE,
                   datasets = "Presence",
                   intercepts = TRUE,
                   fun = "linear", 
                   n.samples = 1000, 
                   seed = 24) 
jint_po <- format_predictions(jint_po)

# On log-intensity scale 
var_names <- c("x", "y", "mean", "sd", "q0.025", "q0.5", "q0.975", "median")
jint_po_r   <- terra::rast(jint_po[, var_names], type = "xyz", crs = proj) 
pred_jt_val <- extract(jint_po_r, dataset$Presence)

pred_po_r   <- terra::rast(pred_po[, var_names], type = "xyz", crs = proj) 
pred_po_val <- extract(pred_po_r, dataset$Presence)

# Correlation (joint prob, separate PO prob)
cor.test(pred_jt_val$mean,  pred_po_val$mean)

df_pred <- data.frame(pred_jt = pred_jt_val$mean, pred_po = pred_po_val$mean)
p_pred  <- ggplot(df_pred, aes(x = pred_jt, y = pred_po)) +
  geom_point(alpha = 0.6, size = 0.6, color = "gray10") +
  geom_abline(intercept = 0, slope = 1, color = "blue", linewidth = 0.7, linetype = "dashed") + 
  annotate("text", x = -7, y = 0, label = "(a)", hjust = 0, fontface = "bold", size = 4) +
  labs(x = "log intensity (shared model)", y = "log intensity (PO model)") +
  coord_cartesian(xlim = c(-7, 0), ylim = c(-7, 0)) +
  theme_bw()
p_pred

# Figure 10
(pj <- p_pred + p_resp)
ggsave("figures/fig10_relationship_response.jpeg", pj, width = 8, height = 5, dpi = 300)


#--------------------------------------------------
#--- 3) Posterior distribution of random effects 
#--------------------------------------------------

#--- a) Plot of random components ----
p <- ggplot()+
  gg(spat_shared$summary.random$bio14_wc30s)+
  theme_bw()

# Figure B3
ggsave("figures/fig_B3_random_bio14.jpeg", p, width = 5, height = 4, dpi = 300)

p <- ggplot()+
  gg(spat_shared$summary.random$shared_spatial)
ggsave("random_shared.jpeg", p, width = 5, height = 4, dpi = 300)

# Spatial range of bio14 on scaled level (corresponds to ID ~ 1)
bio14_rand <- spat_shared$summary.random$bio14_wc30s
bio14_est  <- bio14_rand[2, c("mean", "0.025quant", "0.975quant")]

#--- b) Posterior distribution of latent covariance function 
plot(spde.posterior(spat_shared, "shared_spatial", what = "range"))/
  plot(spde.posterior(spat_shared, "shared_spatial", what = "variance"))

plot(spde.posterior(spat_shared, "bio14_wc30s", what = "range"))/
  plot(spde.posterior(spat_shared, "bio14_wc30s", what = "variance"))

#--- c) Marginals of hyper parameters 
list_marginals <- list(
  "Range for shared spatial" = spat_shared$marginals.hyperpar$"Range for shared_spatial",
  "StDev for shared spatial" = spat_shared$marginals.hyperpar$"Stdev for shared_spatial",
  "Range for rainfall driest month" = spat_shared$marginals.hyperpar$"Range for bio14_wc30s",
  "StDev for rainfall driest month" = spat_shared$marginals.hyperpar$"Stdev for bio14_wc30s",
  "Beta for annual temperature" = spat_shared$marginals.fixed$"bio1_wc30s",
  "Beta for soil texture silt %" = spat_shared$marginals.fixed$"SLTPPT_d2",
  "Beta for SRTM slope" = spat_shared$marginals.fixed$"srtm_slope",
  "Presence Intercept" = spat_shared$marginals.fixed$"Presence_intercept",
  "Abundance Intercept" = spat_shared$marginals.fixed$"Count_intercept")

marginals <- data.frame(do.call(rbind, list_marginals))
marginals$parameter <- rep(names(list_marginals), times = sapply(list_marginals, nrow))

marginals$parameter <- as.factor(marginals$parameter)
marginals$parameter <- factor(marginals$parameter, levels = 
                                c("Abundance Intercept","Presence Intercept","Beta for annual temperature", 
                                  "Beta for SRTM slope", "Beta for soil texture silt %", "Range for rainfall driest month", 
                                  "StDev for rainfall driest month", "Range for shared spatial", "StDev for shared spatial")
)
levels(marginals$parameter)

# Figure 6 (parameters' marginals)
plot_hyper <- ggplot(marginals, aes(x = x, y = y)) + 
  geom_line() +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  labs(x = "", y = "Probability density") + theme_bw()
plot_hyper
ggsave("figures/fig_B4_shared_marginals_hyper.jpeg", plot_hyper, width = 7, height = 5, dpi = 200)


#------------------------------------------------------
#--- 4) Relationship of log-intensity vs covariates 
#------------------------------------------------------

var_names <- c("x", "y", "bio1_wc30s", "bio14_wc30s", "srtm_slope", "SLTPPT_d2", "CLYPPT_d6", 
             "mean", "sd", "q0.025", "q0.5", "q0.975", "median")

#jint_rast <- terra::rast(jint[, var_names], type = "xyz", crs = proj) 
#jint_df <- extract(jint_rast, xy_obs)

pred <- jint 
p1 <- ggplot(pred, aes(x = bio1_wc30s, y = mean)) +
  geom_point(alpha = 0.6, size = 0.1, color = "gray20") +
  #geom_smooth(method = "loess", formula = y ~ x, span = 0.6, col = "green", se = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x, span = 0.6, col = "blue", se = FALSE) +
  coord_cartesian(ylim = c(-8, max(pred$mean))) +
  labs(x = "Annual temperature", y = "Log intensity") +
  annotate("text", x = min(pred$bio1_wc30s), y = max(pred$mean), label = "(a)", 
           hjust = 0, fontface = "bold", size = 5) + 
  theme_gray()
p1

p3 <- ggplot(pred, aes(x = bio14_wc30s, y = mean)) +
  geom_point(alpha = 0.6, size = 0.1, color = "gray20") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "gp"), color = "blue", se = FALSE) + 
  labs(x = "Rainfall of the driest month", y = "Log intensity") +
  annotate("text", x = min(pred$bio14_wc30s), y = max(pred$mean), label = "(b)", 
           hjust = 0, fontface = "bold", size = 5) + 
  theme_gray()
p3

p4 <- ggplot(pred, aes(x = srtm_slope, y = mean)) +
  geom_point(alpha = 0.6, size = 0.1, color = "gray20") +
  geom_smooth(method = "lm", col = "blue", se = FALSE) +
  labs(x = "SRTM slope", y = "Log intensity") +
  annotate("text", x = min(pred$srtm_slope), y = max(pred$mean), label = "(c)", 
           hjust = 0, fontface = "bold", size = 5) +  
  theme_gray()
p4

p5 <- ggplot(pred, aes(x = SLTPPT_d2, y = mean)) +
  geom_point(alpha = 0.6, size = 0.1, color = "gray20") +
  geom_smooth(method = "lm", col = "blue", se = FALSE) +
  labs(x = "Soil texture silt fraction", y = "Log intensity") +
  annotate("text", x = min(pred$SLTPPT_d2), y = max(pred$mean), label = "(d)", 
           hjust = 0, fontface = "bold", size = 5) + 
  theme_gray()
p5

p6 <- ggplot(pred, aes(x = CLYPPT_d6, y = mean)) +
  geom_point(alpha = 0.6, size = 0.1, color = "gray20") +
  geom_smooth(method = "lm", col = "blue", se = FALSE) +
  labs(x = "Soil texture fraction clay", y = "Log intensity") +
  annotate("text", x = min(pred$CLYPPT_d6), y = max(pred$mean), label = "(e)", 
           hjust = 0, fontface = "bold", size = 5) + 
  theme_gray()
p6

# Visualisation
plot_cov <- p1 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)
print(plot_cov)
ggsave("figures/fig6_intensity_vs_covariate.jpeg",  plot_cov, width = 6, height = 6, dpi = 200)


