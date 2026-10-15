#' Model the effect of sample size on median proportions and intervals
#'
#' @param sim_sample_size Output of [run_simulations_consumer_sample_size()] or
#' [run_simulations_source_sample_size()]
#' @param response A string; the response variable
#' @param seed Seed for reproducibility
#'
#' @return A `brmsfit` object
fit_sample_size_effect <- function(sim_sample_size, response, seed = NA) {
  
  response <- rlang::arg_match(response, c("alpha_median", "alpha_CI"))
  
  # Fit distributional power function
  # Estimate separate parameters for each source
  fit <- brm(bf(paste(response, "~ a * sample_size^b"),
                a ~ 0 + source,
                b ~ 0 + source,
                sigma ~ 0 + source*sample_size,
                nl = TRUE),
             family = gaussian(),
             data = purrr::map(
               sim_sample_size[-1],
               ~ bind_rows(.x, .id = "replicate")) |> 
               bind_rows(.id = "sample_size") |> 
               mutate(sample_size = as.numeric(sample_size),
                      alpha_CI = alpha_upper - alpha_lower),
             prior = c(prior(beta(1, 1), nlpar = "a", lb = 0, ub = 1),
                       prior(normal(0, 1), nlpar = "b"),
                       prior(normal(0, 5), dpar = "sigma"),
                       prior(normal(0, 1), dpar = "sigma", coef = "sample_size"),
                       prior(normal(0, 1), dpar = "sigma", coef = "sourcesource2:sample_size"),
                       prior(normal(0, 1), dpar = "sigma", coef = "sourcesource3:sample_size")),
             cores = 4,
             backend = "cmdstan",
             seed = seed)
  
  plot(fit, ask = FALSE)
  
  return(fit)
}

#' Extract the slopes for the mean and for the variance of the power functions
#'
#' @param model Output of [fit_sample_size_effect()]
#' @param seed Seed for reproducibility
#'
#' @return A [brms::brmshypothesis] object
get_slopes <- function(model, seed = NULL) {
  
  hypothesis(model, 
             c("source 1_mean" = "b_sourcesource1 = 0", 
               "source 2_mean" = "b_sourcesource2 = 0",
               "source 3_mean" = "b_sourcesource3 = 0",
               "source 1_variance" = "sigma_sample_size = 0", 
               "source 2_variance" = "sigma_sample_size + sigma_sourcesource2:sample_size = 0", 
               "source 3_variance" = "sigma_sample_size + sigma_sourcesource3:sample_size = 0"),
             seed = seed)
  
}

#' Plot the effect of sample size
#'
#' @param sim_food_web Simulated food web with one baseline and one consumer species
#' @param sim_sample_size Output of [run_simulations_consumer_sample_size()] or
#' [run_simulations_source_sample_size()]
#' @param model_median_sample_size Output of [fit_sample_size_effect()]
#' @param slopes_median_sample_size Output of [get_slopes()]
#' @param model_CI_sample_size Output of [fit_sample_size_effect()]
#' @param slopes_CI_sample_size Output of [get_slopes()]
#' @param xlab The title of the x axis
#' @param .width The `.width` argument passed to [ggdist::point_interval()]
#' @param alpha_points Point transparency
#' @param limits Optional data frame with columns `source` and `alpha_limit`:
#' the proportions the estimates converge to given the fixed samples (output of
#' [get_sample_size_limits()] filtered to one simulation). Drawn as dotted lines.
#' If NULL (default), not drawn
#' @param seed Seed for reproducibility
#'
#' @return A `patchwork` object
plot_sample_size_effect <- function(sim_food_web, 
                                    sim_sample_size, 
                                    model_median_sample_size, 
                                    slopes_median_sample_size,
                                    model_CI_sample_size,
                                    slopes_CI_sample_size,
                                    xlab = "Sample size",
                                    .width = ppoints(20),
                                    alpha_points = 0.1,
                                    limits = NULL,
                                    seed = NULL) {
  
  plot_full_effect <- function(model, yvar, ylab) {
    
    sample_sizes <- names(sim_sample_size[-1]) |> as.numeric()
    
    preds <- expand.grid(source = paste("source", 1:3), 
                         sample_size = seq(min(sample_sizes), 
                                           max(sample_sizes), 
                                           length.out = 51)) |> 
      tidybayes::add_predicted_draws(model, seed = seed) |> 
      mutate(source = gsub("source ", "", source))
    
    purrr::map(
      sim_sample_size,
      ~ bind_rows(.x, .id = "replicate")) |> 
      bind_rows(.id = "sample_size") |> 
      mutate(sample_size = as.numeric(sample_size),
             alpha_CI = alpha_upper - alpha_lower,
             source = gsub("source ", "", source)) |> 
      group_by(sample_size, source) |> 
      ggplot(aes(x = sample_size,  
                 color = source, fill = source)) +
      ggdist::stat_lineribbon(data = preds, 
                              aes(y = .prediction,
                                  fill_ramp = after_stat(.width)), 
                              .width = .width,
                              alpha = 0.3) +
      ggdist::scale_fill_ramp_continuous(range = c(1, 0), 
                                         guide = ggdist::guide_rampbar(title = "CI",
                                                                       theme = theme(legend.title = element_text(vjust = 0.8)))) +
      geom_point(aes(y = {{ yvar }}, 
                     fill = stage(source, after_scale = alpha(fill, alpha_points))),
                 position = position_jitterdodge(jitter.width = 0, dodge.width = 1.25, seed = seed),
                 shape = 21) +
      ggdist::stat_pointinterval(aes(y = {{ yvar }}), 
                                 .width = NA,
                                 point_interval = "mean_qi",
                                 position = position_jitterdodge(jitter.width = 0, dodge.width = 1.25, seed = seed),
                                 shape = 21, color = "black", show.legend = FALSE) +
      labs(x = xlab,
           y = ylab,
           color = "Source",
           fill = "Source") +
      theme_bw() +
      scale_fill_brewer(palette = "Dark2", aesthetics = c("color", "fill"))
  }
  
  plot_slopes <- function(slopes) {
    slopes$samples |> 
      rlang::set_names(slopes$hypothesis$Hypothesis) |> 
      tidyr::pivot_longer(cols = everything(),
                          names_to = "slope") |> 
      tidyr::separate_wider_delim(slope, 
                                  delim = "_", 
                                  names = c("source", "estimate")) |> 
      mutate(source = gsub("source ", "", source),
             estimate = ifelse(estimate == "mean", "&mu;", "&sigma;")) |> 
      group_by(source, estimate) |> 
      ggplot(aes(y = source, fill = source)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      ggdist::stat_slabinterval(aes(x = value), 
                                .width = c(0.5, 0.95), 
                                normalize = "panels", 
                                slab_alpha = 0.5) +
      facet_grid(~estimate, scales = "free") +
      labs(x = "Sample size effect",
           y = "Source",
           fill = "Source") +
      theme_bw() +
      scale_fill_brewer(palette = "Dark2", guide = "none") +
      scale_x_continuous(n.breaks = 3) +
      theme(strip.text = ggtext::element_markdown(face = "italic"),
            panel.spacing.x = unit(5, "mm"))
  }
  
  p1 <- plot_full_effect(model_median_sample_size, 
                         yvar = alpha_median, 
                         ylab = "Median <i>p<sub>k</sub></i>") +
    geom_hline(data = sim_food_web$alpha_species |> 
                 mutate(source = gsub("source ", "", source)),
               aes(yintercept = alpha, color = source),
               linetype = "dashed", 
               show.legend = FALSE) +
    theme(axis.title.y = ggtext::element_markdown())
  
  # Expected limits given the fixed samples (dotted)
  if (!is.null(limits)) {
    p1 <- p1 +
      geom_hline(data = limits |>
                   mutate(source = gsub("source ", "", source)),
                 aes(yintercept = alpha_limit, color = source),
                 linetype = "dotted",
                 show.legend = FALSE)
  }
  p2 <- plot_slopes(slopes_median_sample_size)
  p3 <- plot_full_effect(model_CI_sample_size,
                         yvar = alpha_CI,
                         ylab = "95% CI of <i>p<sub>k</sub></i>") +
    theme(axis.title.y = ggtext::element_markdown())
  p4 <- plot_slopes(slopes_CI_sample_size)
  
  (patchwork::guide_area() /
      (p1 + p3) / (p2 + p4)) +
    patchwork::plot_layout(guides = "collect",
                heights = c(1/10, 1, 3/4)) +
    patchwork::plot_annotation(tag_levels = "A") & 
    theme(legend.position = "top")
}

#' Proportions the estimates converge to in the sample-size simulations
#'
#' In each sample-size simulation only one data set is re-sampled while the
#' other two are kept fixed, so as the sample size increases the estimated
#' proportions converge to those implied by the fixed samples rather than to
#' the true proportions:
#' - consumer sample size ([run_simulations_consumer_sample_size()]): consumers
#'   are re-sampled; sources and baselines are fixed. The limit combines the
#'   expected consumer mixture (population source means), the TP error from the
#'   observed baselines ([estimate_tp_error()]), and the observed source means.
#' - source sample size ([run_simulations_source_sample_size()]): sources are
#'   re-sampled; consumers and baselines are fixed. The limit combines the
#'   observed consumer and baseline means, from which TP is solved once, and
#'   the population source means.
#'
#' Limits are obtained by solving the mixing equations for the mean isotope
#' values, with TP obtained from the two-baseline equations. With three sources
#' and two isotopes the solution is unique.
#'
#' @param sim_food_web Simulated food web with two baselines and one consumer species
#'
#' @return A tibble with columns `simulation` ("consumer" or "source"),
#' `source`, `alpha_true`, `alpha_limit`, plus the TP error or TP estimate used
get_sample_size_limits <- function(sim_food_web) {
  
  src_pop <- sim_food_web$sources_summary
  TDF_C <- mean(sim_food_web$TDF$deltaC)
  TDF_N <- mean(sim_food_web$TDF$deltaN)
  TDF <- c(TDF_C, TDF_N)
  p_true <- sim_food_web$alpha_species$alpha
  source_names <- as.character(src_pop$name)
  
  # Observed means of the fixed samples (ordered as factor levels)
  obs_means <- function(type_) {
    sim_food_web$isotopes |> 
      dplyr::filter(type == type_) |> 
      dplyr::group_by(name) |> 
      dplyr::summarise(d13C_mean = mean(d13C), 
                       d15N_mean = mean(d15N), 
                       .groups = "drop") |> 
      dplyr::arrange(name)
  }
  src_obs <- obs_means("source")
  base_obs <- obs_means("baseline")
  cons_obs <- obs_means("higher consumer")
  
  S_pop <- rbind(src_pop$d13C_mean, src_pop$d15N_mean, 1)
  S_obs <- rbind(src_obs$d13C_mean, src_obs$d15N_mean, 1)
  
  # Consumer sample-size simulation
  tp_err <- estimate_tp_error(src_pop, base_obs, TP_b = 2, TDF_N, TDF_C, p_true)
  mix_true <- c(sum(p_true * src_pop$d13C_mean), sum(p_true * src_pop$d15N_mean))
  lim_consumer <- solve(S_obs, c(mix_true - TDF * tp_err, 1))
  
  # Source sample-size simulation: TP estimated once from fixed consumers and baselines
  # Two-baseline equations with observed means (baselines at TP = 2):
  # d13C_c - d13C_b2 = TDF_C (TP - 2) + a (d13C_b1 - d13C_b2), same for d15N
  x <- cons_obs$d13C_mean - base_obs$d13C_mean[2]
  y <- cons_obs$d15N_mean - base_obs$d15N_mean[2]
  dC_b <- base_obs$d13C_mean[1] - base_obs$d13C_mean[2]
  dN_b <- base_obs$d15N_mean[1] - base_obs$d15N_mean[2]
  r <- TDF_C / TDF_N
  a <- min(max((x - r * y) / (dC_b - r * dN_b), 0), 1)
  TP <- 2 + (y - a * dN_b) / TDF_N
  cons_rescaled <- c(cons_obs$d13C_mean, cons_obs$d15N_mean) - TDF * (TP - 1)
  lim_source <- solve(S_pop, c(cons_rescaled, 1))
  
  dplyr::bind_rows(
    tibble::tibble(simulation = "consumer",
                   source = source_names,
                   alpha_true = p_true,
                   alpha_limit = lim_consumer,
                   TP_error = tp_err),
    tibble::tibble(simulation = "source",
                   source = source_names,
                   alpha_true = p_true,
                   alpha_limit = lim_source,
                   TP_estimate = TP)
  )
}