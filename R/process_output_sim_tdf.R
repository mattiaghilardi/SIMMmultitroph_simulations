#' Summarise performance of the TDF SD simulations
#'
#' @param sim_tdf_sd Output of [run_simulations_tdf_sd()] (one or more variants
#' combined with [dplyr::bind_rows()])
#'
#' @return A tibble with bias, absolute error, CI width and coverage of source
#' proportions, and median multiplicative residual error terms, for each
#' combination of TDF variation, error structure, sample size, TDF SD level and source
summarise_tdf_sd <- function(sim_tdf_sd) {
  
  sim_tdf_sd |>
    mutate(error = alpha_median - alpha_true,
           alpha_CI = alpha_upper - alpha_lower,
           covered = alpha_true >= alpha_lower & alpha_true <= alpha_upper) |>
    group_by(tdf_variation, error_structure, sample_size, tdf_sd_level, source) |>
    summarise(sd_C = mean(sd_C),
              sd_N = mean(sd_N),
              bias = mean(error),
              mae = mean(abs(error)),
              CI_width = mean(alpha_CI),
              coverage = mean(covered),
              eps_C = median(eps_C_median),
              eps_N = median(eps_N_median),
              n_rep = dplyr::n(),
              .groups = "drop")
}


#' Plot the effect of the TDF SD supplied to MixSIAR
#'
#' @param sim_tdf_sd Output of [run_simulations_tdf_sd()] for one TDF variation
#' @param alpha_points Point transparency
#' @param seed Seed for reproducibility of jittering
#' @param limits Optional data frame with columns `source` and `alpha_limit`:
#' proportions implied by the fixed source and baseline samples
#' ([get_sample_size_limits()] filtered to `simulation == "consumer"`, as the
#' consumer data are re-sampled and sources and baselines are fixed, as in the
#' consumer sample-size simulations). Drawn as dotted lines.
#'
#' Panels are faceted by error structure (rows) and TDF SD supplied to MixSIAR
#' (columns), with consumer sample size on the x axis. Sample sizes without a
#' given error structure (e.g. n = 1 with process x residual) are simply absent
#' from that row, so no empty panels are drawn.
#'
#' @return A `patchwork` object with (A) median source proportions,
#' (B) 95% CI width of source proportions, (C) coverage of the 95% CI and,
#' if any process x residual models were fitted, (D) multiplicative residual
#' error terms
plot_tdf_sd <- function(sim_tdf_sd, alpha_points = 0.2, seed = NULL, limits = NULL) {
  
  level_labels <- c(zero        = "SD = 0",
                    one_step    = "SD = &sigma;<sub>&lambda;</sub>",
                    scaled_sqrt = "SD = &radic;(TP<sub><i>c</i></sub>-1)&middot;&sigma;<sub>&lambda;</sub>",
                    scaled      = "SD = (TP<sub><i>c</i></sub>-1)&middot;&sigma;<sub>&lambda;</sub>",
                    scaled_TP   = "SD = &sigma;<sub>total</sub>")
  err_labels <- c(process_residual = "Proc. &times; resid.",
                  process = "Process only")
  
  dat <- sim_tdf_sd |>
    mutate(alpha_CI = alpha_upper - alpha_lower,
           source = gsub("source ", "", source),
           tdf_sd_level = factor(level_labels[as.character(tdf_sd_level)],
                                 levels = level_labels[levels(tdf_sd_level)]),
           error_structure = factor(err_labels[as.character(error_structure)],
                                    levels = err_labels[levels(error_structure)]),
           sample_size = factor(sample_size,
                                levels = sort(unique(sample_size)))) |>
    # drop unused levels (e.g. TDF SD levels not fitted)
    droplevels()
  
  facets <- facet_grid(error_structure ~ tdf_sd_level)
  dodge_width <- 0.6
  xlab <- "Consumer sample size"
  
  theme_tdf <- theme_bw() +
    theme(axis.title.y = ggtext::element_markdown(),
          strip.text = ggtext::element_markdown(),
          legend.position = "top")
  
  # Points for each replicate (posterior medians, dodged by source, not
  # jittered) + mean across replicates
  plot_response <- function(yvar, ylab) {
    dat |>
      ggplot(aes(x = sample_size, y = {{ yvar }}, color = source, fill = source,
                 group = source)) +
      geom_point(aes(fill = stage(source, after_scale = alpha(fill, alpha_points))),
                 position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75, seed = seed),
                 shape = 21, show.legend = FALSE) +
      ggdist::stat_pointinterval(.width = NA,
                                 point_interval = "mean_qi",
                                 position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75, seed = seed),
                                 shape = 21, size = 4, color = "black") +
      facets +
      labs(x = xlab,
           y = ylab,
           fill = "Source") +
      scale_fill_brewer(palette = "Dark2", aesthetics = c("color", "fill")) +
      guides(color = "none",
             fill = guide_legend(override.aes = list(linewidth = NA))) +
      theme_tdf
  }
  
  p1 <- plot_response(alpha_median, "Median <i>p<sub>k</sub></i>") +
    geom_hline(data = dat |> distinct(source, alpha_true),
               aes(yintercept = alpha_true, color = source),
               linetype = "dashed",
               show.legend = FALSE)
  
  # Expected limits given the fixed samples (dotted), as in the sample-size figures
  if (!is.null(limits)) {
    p1 <- p1 +
      geom_hline(data = limits |>
                   mutate(source = gsub("source ", "", source)),
                 aes(yintercept = alpha_limit, color = source),
                 linetype = "dotted",
                 show.legend = FALSE)
  }
  
  p2 <- plot_response(alpha_CI, "95% CI of <i>p<sub>k</sub></i>")
  
  # Coverage of 95% CI across replicates
  p3 <- dat |>
    group_by(error_structure, sample_size, tdf_sd_level, source) |>
    summarise(coverage = mean(alpha_true >= alpha_lower & alpha_true <= alpha_upper),
              .groups = "drop") |>
    ggplot(aes(x = sample_size, y = coverage, color = source, group = source)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line(position = position_dodge(width = dodge_width / 2),
              show.legend = FALSE) +
    geom_point(position = position_dodge(width = dodge_width / 2),
               size = 2, show.legend = FALSE) +
    facets +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = xlab,
         y = "Coverage of 95% CI") +
    scale_color_brewer(palette = "Dark2") +
    theme_tdf
  
  plots <- list(p1, p2, p3)
  
  # Multiplicative residual error terms (one value per replicate, not per source),
  # only when process x residual models were fitted
  eps <- dat |>
    filter(!is.na(eps_C_median)) |>
    distinct(sample_size, replicate, tdf_sd_level, eps_C_median, eps_N_median)
  
  if (nrow(eps) > 0) {
    p4 <- eps |>
      tidyr::pivot_longer(cols = c(eps_C_median, eps_N_median),
                          names_to = "isotope",
                          values_to = "epsilon") |>
      mutate(isotope = ifelse(isotope == "eps_C_median",
                              "&delta;<sup>13</sup>C",
                              "&delta;<sup>15</sup>N")) |>
      ggplot(aes(x = sample_size, y = epsilon)) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_point(alpha = alpha_points, shape = 21, fill = "grey50") +
      ggdist::stat_pointinterval(.width = NA,
                                 point_interval = "median_qi", 
                                 shape = 21, fill = "white", color = "black", 
                                 size = 4) +
      facet_grid(isotope ~ tdf_sd_level) +
      scale_y_log10() +
      labs(x = xlab,
           y = "Median &epsilon; (proc. &times; resid.)") +
      theme_tdf
    plots <- c(plots, list(p4))
  }
  
  # Legend in its own row above all panels; tags at the top-left of each panel
  if (length(plots) == 4) {
    layout <- patchwork::guide_area() /
      (plots[[1]] + plots[[2]]) /
      (plots[[3]] + plots[[4]])
    heights <- c(0.08, 1, 1)
  } else {
    layout <- patchwork::guide_area() /
      (plots[[1]] + plots[[2]] + plots[[3]])
    heights <- c(0.12, 1)
  }
  
  layout +
    patchwork::plot_layout(guides = "collect", heights = heights) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(legend.position = "top",
          plot.tag.position = c(0, 1),
          plot.tag = element_text(face = "bold", hjust = 0, vjust = 1))
}


#' Summarise performance of the TDF offset simulations (scenario B)
#'
#' @param sim_tdf_offset Output of [run_simulations_tdf_offset()]
#'
#' @return A tibble with bias, absolute error, CI width and coverage of source
#' proportions for each combination of TDF variation, error structure, sample
#' size, TDF SD level and source, plus coverage pooled across sources
summarise_tdf_offset <- function(sim_tdf_offset) {
  
  dat <- sim_tdf_offset |> 
    mutate(error = alpha_median - alpha_true,
           alpha_CI = alpha_upper - alpha_lower,
           covered = alpha_true >= alpha_lower & alpha_true <= alpha_upper)
  
  by_source <- dat |> 
    group_by(tdf_variation, error_structure, sample_size, tdf_sd_level, source) |> 
    summarise(bias = mean(error),
              mae = mean(abs(error)),
              CI_width = mean(alpha_CI),
              coverage = mean(covered),
              n_rep = dplyr::n(),
              .groups = "drop")
  
  pooled <- dat |> 
    group_by(tdf_variation, error_structure, sample_size, tdf_sd_level) |> 
    summarise(mae = mean(abs(error)),
              CI_width = mean(alpha_CI),
              coverage = mean(covered),
              n_intervals = dplyr::n(),
              .groups = "drop") |> 
    mutate(source = "pooled")
  
  bind_rows(by_source, pooled)
}


#' Plot the effect of a misspecified mean TDF (scenario B)
#'
#' @param sim_tdf_offset Output of [run_simulations_tdf_offset()]
#' @param sim_food_web The simulated food web used in [run_simulations_tdf_offset()],
#' to compute the expected TP error with [expected_tp_error()]
#' @param alpha_points Point transparency
#' @param seed Seed for reproducibility of jittering
#'
#' @return A `patchwork` object with (A) coverage of the 95% CI against sample
#' size, pooled across sources, (B) 95% CI width against sample size,
#' and (C) observed (points) and expected (background, see [expected_tp_error()])
#' error in estimated TP against the offsets in the mean C and N TDFs
plot_tdf_offset <- function(sim_tdf_offset, sim_food_web, alpha_points = 0.3, seed = NULL) {
  
  level_labels <- c(zero        = "SD = 0",
                    one_step    = "SD = &sigma;<sub>&lambda;</sub>",
                    scaled_sqrt = "SD = &radic;(TP<sub><i>c</i></sub>-1)&middot;&sigma;<sub>&lambda;</sub>",
                    scaled      = "SD = (TP<sub><i>c</i></sub>-1)&middot;&sigma;<sub>&lambda;</sub>",
                    scaled_TP   = "SD = &sigma;<sub>total</sub>")
  
  dat <- sim_tdf_offset |> 
    mutate(error = alpha_median - alpha_true,
           alpha_CI = alpha_upper - alpha_lower,
           covered = alpha_true >= alpha_lower & alpha_true <= alpha_upper,
           source = gsub("source ", "", source),
           tdf_sd_level = factor(level_labels[as.character(tdf_sd_level)],
                                 levels = level_labels[levels(tdf_sd_level)]),
           sample_size_f = factor(paste("n =", sample_size),
                                  levels = paste("n =", sort(unique(sample_size))))) |> 
    droplevels()
  
  theme_tdf <- theme_bw() +
    theme(axis.title = ggtext::element_markdown(),
          strip.text = ggtext::element_markdown(),
          legend.text = ggtext::element_markdown(),
          legend.title = ggtext::element_markdown(),
          legend.position = "top")
  
  dodge <- position_dodge(width = 0.4)
  
  # Coverage pooled across sources, with binomial 95% CI
  p1 <- dat |> 
    group_by(sample_size, tdf_sd_level) |> 
    summarise(k = sum(covered), n = dplyr::n(), .groups = "drop") |> 
    mutate(coverage = k / n,
           lower = purrr::map2_dbl(k, n, ~ stats::binom.test(.x, .y)$conf.int[1]),
           upper = purrr::map2_dbl(k, n, ~ stats::binom.test(.x, .y)$conf.int[2])) |> 
    ggplot(aes(x = factor(sample_size), y = coverage, 
               color = tdf_sd_level, shape = tdf_sd_level, group = tdf_sd_level)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line(position = dodge) +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = dodge) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Consumer sample size",
         y = "Coverage of 95% CI",
         color = "TDF SD supplied to the mixing model",
         shape = "TDF SD supplied to the mixing model") +
    # ordered levels: light to dark grey, distinct from source colours
    scale_color_grey(start = 0.65, end = 0.1) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
    theme_tdf
  
  # CI width
  p2 <- dat |> 
    ggplot(aes(x = factor(sample_size), y = alpha_CI, 
               color = tdf_sd_level, shape = tdf_sd_level)) +
    ggdist::stat_pointinterval(point_interval = "mean_qi", .width = 0.95,
                               position = dodge,
                               show.legend = FALSE) + # legend from panel A
    labs(x = "Consumer sample size",
         y = "95% CI of <i>p<sub>k</sub></i>") +
    scale_color_grey(start = 0.65, end = 0.1) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
    theme_tdf
  
  # TP error against both TDF offsets (one point per replicate and sample size).
  # Background: expected TP error (baseline mismatch + TDF offsets, no sampling
  # error), with the dashed line where it is zero.
  # Points: observed TP error. Fill: signed TP error (blue = underestimated,
  # red = overestimated); size: absolute observed TP error (no legend)
  tp_err <- dat |> 
    distinct(sample_size_f, replicate, offset_C, offset_N, TP_error)
  
  pad <- 0.1
  grid <- tidyr::expand_grid(
    offset_C = seq(min(tp_err$offset_C) - pad, max(tp_err$offset_C) + pad, length.out = 101),
    offset_N = seq(min(tp_err$offset_N) - pad, max(tp_err$offset_N) + pad, length.out = 101)
  ) |> 
    mutate(expected = expected_tp_error(sim_food_web, offset_C, offset_N))
  
  max_err <- max(abs(c(tp_err$TP_error, grid$expected)))
  
  p4 <- tp_err |> 
    ggplot(aes(x = offset_C, y = offset_N)) +
    geom_raster(data = grid, aes(fill = expected), alpha = 0.5, interpolate = TRUE) +
    geom_contour(data = grid, aes(z = expected), breaks = 0,
                 color = "grey20", linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_point(aes(fill = TP_error, size = abs(TP_error)),
               shape = 21, color = "grey20") +
    facet_wrap(~ sample_size_f, nrow = 1) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-max_err, max_err)) +
    scale_size_area(max_size = 5, guide = "none") +
    coord_cartesian(expand = FALSE, clip = "off") +
    labs(x = "Offset in mean &delta;<sup>13</sup>C TDF (&permil;)",
         y = "Offset in mean &delta;<sup>15</sup>N TDF (&permil;)",
         fill = "Error in TP<sub><i>c</i></sub>") +
    theme_tdf +
    theme(legend.position = "bottom")
  
  # TDF SD legend collected above A and B
  # TP error colour bar below C
  top <- (patchwork::guide_area() / (p1 + p2)) +
    patchwork::plot_layout(guides = "collect", heights = c(0.12, 1))
  
  patchwork::wrap_plots(top, p4, ncol = 1) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(face = "bold", hjust = 0, vjust = 1))
}