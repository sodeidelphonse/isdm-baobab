# Copyright (c) 2026 SODE Akoeugnigan Idelphonse
# Licensed under the MIT License (see LICENSE file for details)
 

#-----------------------
#--- Helper functions
#-----------------------

#--- 1) Function to compute the standard error
std_err <- function(x, na.rm = TRUE) {
  stopifnot(is.numeric(x))
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}

#--- 2) Map the model name to the available spatial model type 
set_points_spatial <- function(x = "") {
  x <- as.character(x)
  valid_points_spatial <- c("shared", "copy", "individual", "correlate", "count", "po")
  if (x == "") {
    return(NULL)
  }
  points_spatial_type <- unlist(strsplit(x, "_", fixed = TRUE))[1]
  if (!is.na(points_spatial_type) && points_spatial_type %in% valid_points_spatial) {
    return(points_spatial_type)
  } else {          
    stop(sprintf("Invalid 'model name' provided: '%s'. Expected a model name that starts with one of: %s, or an empty string for NULL.",
                 x, paste0("'", valid_points_spatial, "'", collapse = ", ")), call. = FALSE)
  }
}


#--- 3) Functions to generate maps for variables with different scales 
# Main wrapper
gg_map <- function(
    data_pred,
    vars_to_plot = c("mean", "q0.025", "q0.975"),
    base_map = NULL,
    boundary_map = NULL,
    color_gradient = NULL,
    x_axis_breaks = NULL,
    y_axis_breaks = NULL,
    label_x_coord = 580,  
    label_y_coord = 1400  
) {
  stopifnot(is.character(vars_to_plot))
  data_to_plot <- format_predictions(data_pred, base_map)
  
  # Calculate IQR directly on the prepared data, if requested.
  if ("IQR" %in% vars_to_plot) {
    if (!all(c("q0.975", "q0.025") %in% names(data_to_plot))) {
      stop("To plot IQR, columns 'q0.975' and 'q0.025' must be in the data.", call. = FALSE)
    }
    data_to_plot$IQR <- data_to_plot$q0.975 - data_to_plot$q0.025
  }
  
  # Ensure all variables to be plotted are present
  if (!all(vars_to_plot %in% names(data_to_plot))) {
    missing_vars <- setdiff(vars_to_plot, names(data_to_plot))
    stop(sprintf("Requested variables not found in data: '%s'.", paste(missing_vars, collapse = ", ")),
         call. = FALSE)
  }
  
  # Create a named vector to map variable names to user-friendly titles
  legend_titles <- c(
    "q0.025" = "q2.5%",
    "mean" = "Mean",
    "q0.975" = "q97.5%",
    "sd" = "StDev",
    "IQR" = "IQR"
  )
  
  # Generate a list of plots 
  plot_list <- lapply(seq_along(vars_to_plot), function(i) {
    var   <- vars_to_plot[i]
    title <- legend_titles[var]
    label <- paste0("(", letters[i], ")")
    plot_model_output(
      data_to_plot = data_to_plot,
      fill_variable = var,
      legend_title = title,
      panel_label = label,
      base_map = base_map,
      boundary_map = boundary_map,
      color_gradient = color_gradient,
      x_axis_breaks = x_axis_breaks,
      y_axis_breaks = y_axis_breaks,
      label_x_coord = label_x_coord,
      label_y_coord = label_y_coord
    )
  })
  patchwork::wrap_plots(plot_list)
}

# The plot_model_output() helper 
plot_model_output <- function(
    data_to_plot,
    fill_variable,
    legend_title,
    panel_label,
    base_map,
    boundary_map,
    color_gradient,
    x_axis_breaks,
    y_axis_breaks,
    label_x_coord,
    label_y_coord
) {
  # Dynamically get the variable name for the aesthetic
  fill_var_sym <- rlang::sym(fill_variable)
  is_sf_data <- inherits(data_to_plot, "sf")
  
  p <- ggplot(base_map) +
    labs(title = "", x = "Longitude", y = "Latitude") + # title = panel_label
    scale_x_continuous(breaks = x_axis_breaks) +
    scale_y_continuous(breaks = y_axis_breaks) +
    theme_bw() +
    annotate("text", x = label_x_coord, y = label_y_coord, label = panel_label, size = 3, fontface = "bold") +
    ggspatial::annotation_north_arrow(location = "tl", height = unit(0.8, "cm"), width = unit(0.3, "cm")) +
    ggspatial::annotation_scale(location = "br", bar_cols = c("grey60", "white"))
  
  if (is_sf_data) {
    p <- p + geom_sf(data = data_to_plot, aes(color = !!fill_var_sym))+
      scale_color_gradientn(colours = color_gradient, name = legend_title)
  } else {
    p <- p + geom_tile(data = data_to_plot, aes(x = x, y = y, fill = !!fill_var_sym))+
      scale_fill_gradientn(colours = color_gradient, name = legend_title)
  }
  p <- p + coord_sf(crs = st_crs(4326), expand = FALSE, label_axes = "SW") +
    geom_sf(data = boundary_map, fill = NA, color = "grey20")
  
  return(p)
}

