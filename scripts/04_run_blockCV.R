# Copyright (c) 2026 SODE Akoeugnigan Idelphonse
# Licensed under the MIT License (see LICENSE file for details)
 

#----------------------------------------------------
#--- Models evaluation via block cross-validation
#----------------------------------------------------

library(sf)
library(terra)
library(ggplot2)
library(inlabru)
library(PointedSDMs)
library(isdmtools)

library(dplyr)
library(data.table)
library(agricolae)

source("scripts/03_evaluation_pipeline.r")
source("scripts/08_utils.r")

#--- Visualise the CV folds (default seed = 23) 
datasets <- dataset
names(datasets) <- c("Presence", "Abundance")
folds <- create_folds(datasets, ben_utm, k = 5) 

p_cv <- plot(folds)
p_cv <- p_cv +
    scale_x_continuous(breaks = seq(0, 4, 1)) +
    scale_y_continuous(breaks = seq(6, 13, 2)) +
    labs(title = "")
print(p_cv)

# Figure 5
ggsave("figures/fig5_block_cv.jpeg",  p_cv, width = 8, height = 6, dpi = 300)
rm(datasets, p_cv)

#--- Notes !!! ------------------
# Prediction grids was created, and stored to disk for the prediction.
# saveRDS(fm_pixels(mesh, mask = ben_utm), "data/pred_points.RDS")
fm_pix <- readRDS("data/pred_points.rds")

# Cleaned rasters are also saved to disk to avoid serialization issues 
cov_path <- "data/covariates_pc.rds"

xy_obs <- rbind(st_coordinates(point_utm)[, c("X","Y")], 
                 st_coordinates(abund_utm)[abund_utm$counts > 0, c("X","Y")])

# Define the tasks configuration 
model_names <- c("shared", "shared_bias", "copy", "copy_bias", "individual", "shared_sc", 
                 "shared_bias_sc", "copy_sc", "copy_bias_sc", "individual_sc", "count", "po")

grid <- expand.grid(fold = 1:5, model = model_names)
grid$resp_type <- c(rep("joint.po", 50), rep("count", 5), rep("po", 5))

metrics <- c("auc", "tss", "recall", "precision", "f1", "specificity", "accuracy",  
             "npv", "fpr", "fnr", "mae", "rmse")

roc_composite <- c("auc", "tss", "accuracy") 

#--- a) Run a few jobs to test (optional)
system.time(
  res_test <- lapply(5:10, function(i){  
    evaluate_model_logged(i, grid, "results", dataset, mesh, proj, boundary = ben_utm, 
                          xy_excluded = xy_obs, cov_path = cov_path, bias_names = "Presence", 
                          metrics = metrics, roc_composite = roc_composite) 
  })
)

# View results 
res_df <- do.call(rbind, lapply(res_test, as.data.frame))
print(res_df)

#--- b) Run all jobs serially for ROC-based metrics
system.time(
  res_roc  <- lapply(seq_len(nrow(grid)), function(i){  
    evaluate_model_logged(i, grid, "results", dataset, mesh, proj, boundary = ben_utm, 
                   xy_excluded = xy_obs, cov_path = cov_path, bias_names = "Presence", 
                   metrics = metrics, roc_composite = roc_composite) 
     })
 )

#--- c) Resume missing jobs (if any job crashes or is not done yet) 
system.time(
   resume_tasks (grid, "results", dataset, mesh, proj, 
              boundary = ben_utm, xy_excluded = xy_obs, cov_path = cov_path, 
              bias_names = "Presence", metrics = metrics, roc_composite = roc_composite) 
 )

res_roc_df <- compile_results ("results")
write.csv(res_roc_df, "results/Evaluation_roc_final.csv", row.names = FALSE)

#--- d) Error-based metrics for count model alone ----
metrics_cont <- c("rmse", "mae")
system.time(
     res_cont  <- lapply(which(grid$resp_type == "count"), function(i) {
       evaluate_model_logged(i, grid, "result_cont", dataset, mesh, proj, boundary = ben_utm, 
                             xy_excluded = xy_obs, cov_path = cov_path, bias_names = "Presence", 
                             metrics = metrics_cont, roc_composite = NULL) 
     })
 )

res_count_df <- do.call(rbind, lapply(res_cont, as.data.frame))
write.csv(res_count_df, "result_cont/Evaluation_cont_final.csv", row.names = FALSE)

#--- f) Error-based metrics for the shared model ----
system.time(
  res_jt <- lapply(which(grid$model == "shared"), function(i) {
    evaluate_model_logged(i, grid, "result_cont_sh", dataset, mesh, proj, boundary = ben_utm, 
                          xy_excluded = xy_obs, cov_path = cov_path, bias_names = "Presence", 
                          metrics = metrics_cont, roc_composite = NULL) 
  })
)
res_jt_df <- do.call(rbind, lapply(res_jt, as.data.frame))
res_jt_df


#-------------------------------------------
#--- Analysis of results (metrics Table)
#-------------------------------------------

#--- 1) ROC-based metrics summary --------------
res_df <- compile_results ("results") |> 
  dplyr::slice(1:50)

# summary
summary_df <- res_df |>
  group_by(Model) |>
  summarise(
    across(
      ends_with("_Comp"),
      .fns = list(mean = ~ mean(.x, na.rm = TRUE), se = ~ std_err(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    TOT_SCORE_mean = mean(TOT_ROC_SCORE, na.rm = TRUE),
    TOT_SCORE_se = std_err(TOT_ROC_SCORE, na.rm = TRUE),
    DIC_mean = mean(DIC, na.rm = TRUE),
    DIC_se = std_err(DIC, na.rm = TRUE),
    .groups = "drop" 
  ) |>
  arrange(desc(AUC_Comp_mean)) |>
  dplyr::select( 
    Model,
    AUC_Comp_mean, AUC_Comp_se,
    TSS_Comp_mean, TSS_Comp_se,
    ACCURACY_Comp_mean, ACCURACY_Comp_se,
    TOT_SCORE_mean, TOT_SCORE_se,
    DIC_mean, DIC_se 
  )
summary_df

# The best model
which.max(summary_df$AUC_Comp_mean)
which.max(summary_df$TSS_Comp_mean)
which.max(summary_df$TOT_SCORE_mean)

#--- 2) Significance of the main composite metrics (AUC, TSS , etc.)

# a) Parametric and post-hoc for total composite score
mod_aov <- aov(TOT_ROC_SCORE ~ Model, data = res_df)
summary(mod_aov)
shapiro.test(resid(mod_aov))    # normality satisfied 
LSD.test(mod_aov, "Model", group = TRUE, p.adj="bonferroni")$groups

# b) Parametric and post-hoc tests for TSS
mod_aov <- aov(TSS_Comp ~ Model, data = res_df)
summary(mod_aov)
shapiro.test(resid(mod_aov))    # normality satisfied
LSD.test(mod_aov, "Model", group = TRUE, p.adj="bon")$groups

# c)  Parametric and post-hoc tests for AUC
mod_aov <- aov(AUC_Comp ~ Model, data = res_df)
summary(mod_aov)
shapiro.test(resid(mod_aov))    # normality satisfied
LSD.test(mod_aov, "Model", group = TRUE, p.adj="bon")$groups

# d)  Parametric test for Accuracy
mod_aov <- aov(ACCURACY_Comp ~ Model, data = res_df)
summary(mod_aov)
shapiro.test(resid(mod_aov))    
LSD.test(mod_aov, "Model", group = TRUE)$groups

# Nonparametric and post-hoc tests
kruskal.test(res_df$ACCURACY_Comp ~ res_df$Model)
kruskal(res_df$ACCURACY_Comp, res_df$Model, p.adj = "bon")$groups

# e)  Parametric and post-hoc tests for F1 score
mod_aov <- aov(F1_Comp ~ Model, data = res_df)
summary(mod_aov)
shapiro.test(resid(mod_aov))    
LSD.test(mod_aov, "Model", group = TRUE)$groups

# Nonparametric and post-hoc tests
kruskal.test(res_df$F1_Comp ~ res_df$Model)
kruskal(res_df$F1_Comp, res_df$Model, p.adj = "bon")$groups


#--- 3) Error-based metrics for count data ----
res_dfc <- read.csv("result_cont/Evaluation_cont_final.csv")
summary_dfc <- res_dfc |>
  group_by(Model) |>
  summarise(
    across(
      ends_with("_Comp"),
      .fns = list(mean = ~ mean(.x, na.rm = TRUE), se = ~ std_err(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    TOT_SCORE_mean = mean(TOT_ERROR_SCORE, na.rm = TRUE),
    TOT_SCORE_se = std_err(TOT_ERROR_SCORE, na.rm = TRUE),
    DIC_mean = mean(DIC, na.rm = TRUE),
    DIC_se = std_err(DIC, na.rm = TRUE),
    .groups = "drop" 
  ) |>
  arrange(desc(TOT_SCORE_mean)) |>
  select( 
    Model,
    MAE_Comp_mean, MAE_Comp_se,
    RMSE_Comp_mean, RMSE_Comp_se,
    TOT_SCORE_mean, TOT_SCORE_se,
    DIC_mean, DIC_se 
  )
summary_dfc


#--- 4) All metrics summary --------------------
summary_df2 <- res_df |>
  group_by(Model) |>
  summarise(
    across(
      ends_with("_Comp"),
      .fns = list(mean = ~ mean(.x, na.rm = TRUE), se = ~ std_err(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    TOT_SCORE_mean = mean(TOT_ROC_SCORE, na.rm = TRUE),
    TOT_SCORE_se = std_err(TOT_ROC_SCORE, na.rm = TRUE),
    DIC_mean = mean(DIC, na.rm = TRUE), 
    DIC_se = std_err(DIC, na.rm = TRUE),
    .groups = "drop" 
  ) |>
  arrange(desc(TOT_SCORE_mean))

# Transpose the table
t_summary_df2 <- transpose(summary_df2)
rownames(t_summary_df2) <- colnames(summary_df2)
colnames(t_summary_df2) <- rownames(summary_df2)

write.csv(t_summary_df2, "results/Average_composite_scores.csv")

