#' Compute the absolute error in *p~k~* and the proportion of true *p~k~* within the 95% CI
#'
#' @param sim_baselines Output of [run_simulations_baselines()]
#' @param nrep Number of simulations run
#'
#' @return A data frame
get_performance_pk <- function(sim_baselines, nrep = 100) {
  
  purrr::map(
  1:2,
  function(i)
    purrr::map(1:nrep,
               function(j)
                 sim_baselines[[i]][[j]][["estimated_pk"]] |>
                 mutate(n_baselines = i, replicate = j))) |>
  bind_rows() |>
  mutate(n_baselines = as.factor(n_baselines),
         within_95CI = ifelse(alpha_observed >= alpha_lower & alpha_observed <= alpha_upper, 1, 0),
         abs_error = abs(alpha_median - alpha_observed))
}

# performance_pk |> 
#   group_by(n_baselines) |> 
#   summarise(prop_within_95CI = sum(within_95CI) / n())

#' Compute the absolute error in TP and the proportion of true TP within the 95% CI
#'
#' @param sim_baselines Output of [run_simulations_baselines()]
#' @param nrep Number of simulations run
#'
#' @return A data frame
get_performance_TP <- function(sim_baselines, nrep = 100) {
  
  purrr::map(
  1:2,
  function(i)
    purrr::map(1:nrep,
               function(j)
                 sim_baselines[[i]][[j]][["estimated_TP"]] |>
                 mutate(n_baselines = i, replicate = j))) |>
  bind_rows() |>
  mutate(n_baselines = as.factor(n_baselines),
         within_95CI = ifelse(TP_observed >= lower & TP_observed <= upper, 1, 0),
         abs_error = abs(mode - TP_observed))
}

# performance_TP |> 
#   group_by(n_baselines) |> 
#   summarise(prop_within_95CI = sum(within_95CI) / n())

#' Model the maximum absolute error in *p~k~*
#'
#' @param performance_pk Output of [get_performance_pk()]
#' @param performance_TP Output of [get_performance_TP()]
#' @param seed Seed for reproducibility
#'
#' @return A `brmsfit` object
fit_baseline_effect <- function(performance_pk, performance_TP, seed = NA) {
  
  data <- performance_pk |> 
    group_by(n_baselines, replicate) |> 
    summarise(max_abs_error = max(abs_error)) |>
    ungroup() |> 
    left_join(performance_TP |> 
                select(n_baselines, replicate, abs_error_TP = abs_error)) |> 
    mutate(replicate = as.factor(replicate))
  
  fit <- brm(bf(max_abs_error ~ 1 + n_baselines*log(abs_error_TP) + (1 | replicate),
                phi ~ 1 + n_baselines + (1 | replicate)),
             family = Beta(),
             data = data,
             prior = prior(normal(0, 1), class = b),
             cores = 4,
             backend = "cmdstan",
             control = list(adapt_delta = 0.999),
             seed = seed)
  
  plot(fit, ask = FALSE)
  
  return(fit)
}

#' Extract contrasts in mean and precision
#'
#' @param model_baseline Output of [fit_baseline_effect()]
#' @param seed Seed for reproducibility
#'
#' @return A [brms::brmshypothesis] object
get_baseline_contrasts <- function(model_baseline, seed = NULL) {
  
  brms::hypothesis(model_baseline,
                   c("&mu;" = "plogis(Intercept) > plogis(Intercept + n_baselines2)",
                     "&varphi;" = "exp(phi_Intercept) > exp(phi_Intercept + phi_n_baselines2)"),
                   seed = seed)
}

#' Plot the effect of the number of baselines on maximum absolute error in *p~k~* and contrasts
#'
#' @param performance_pk Output of [get_performance_pk()]
#' @param performance_TP Output of [get_performance_TP()]
#' @param model_baseline Output of [fit_baseline_effect()]
#' @param baseline_contrasts Output of [get_baseline_contrasts()]
#' @param seed Seed for reproducibility
#'
#' @return A `patchwork` object
plot_baseline_effect <- function(performance_pk, 
                                 performance_TP,
                                 model_baseline, 
                                 baseline_contrasts, 
                                 seed = NULL) {
  
  # Get predictions at the intercept (i.e. log(abs_error_TP) = 0)
  preds <- data.frame(n_baselines = 1:2, 
                      abs_error_TP = exp(0)) |>
    tidybayes::add_predicted_draws(model_baseline,
                                   re_formula = NA,
                                   seed = seed)
  
  # Plot effect
  p1 <- performance_pk |> 
    group_by(n_baselines, replicate) |> 
    summarise(max_abs_error = max(abs_error)) |>
    ungroup() |> 
    ggplot(aes(x = n_baselines, y = max_abs_error)) +
    geom_point(position = position_jitter(width = 0.3, seed = seed), 
               alpha = 0.1) +
    ggdist::stat_pointinterval(data = preds,
                               aes(y = .prediction),
                               .width = c(0.5, 0.95)) +
    ylim(0, 1) +
    labs(x = "Number of baselines",
         y = "Max absolute error in <i>p<sub>k</sub></i>") +
    theme_bw() +
    theme(axis.title.y = ggtext::element_markdown())
  
  # Plot contrasts
  p2 <- baseline_contrasts$samples |>
    rlang::set_names(baseline_contrasts$hypothesis$Hypothesis) |>
    tidyr::pivot_longer(cols = everything(),
                        names_to = "hyp") |>
    ggplot(aes(x = value)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggdist::stat_slabinterval(.width = c(0.5, 0.95),
                              normalize = "xy",
                              slab_alpha = 0.5) +
    facet_grid(~hyp, scales = "free") +
    labs(x = "Contrast: 1 baseline - 2 baselines",
         y = "Normalised density") +
    theme_bw() +
    theme(strip.text.x = ggtext::element_markdown(face = "italic"))
  
  p1 / p2 +
    patchwork::plot_annotation(tag_levels = "A")
}

#' #' Compute TP-corrected difference in \u03b4^15^N between consumer and baseline(s)
#' #'
#' #' @param sim_results Output of [run_simulations_baselines()] for a single simulation
#' #'
#' #' @return A value
#' get_d15N_diff <- function(sim_results) {
#'   
#'   # isotope summary at source TP (TP = 1)
#'   summary_tp1 <- sim_results$data$alpha_baselines |> 
#'     bind_rows(sim_results$data$alpha_species |> 
#'                 select(name = species, source, alpha)) |> 
#'     left_join(sim_results$data$sources_summary, by = c("source" = "name")) |> 
#'     mutate(d13C = alpha * d13C_mean, d15N = alpha * d15N_mean) |> 
#'     group_by(name) |> 
#'     summarise(d13C = sum(d13C),
#'               d15N = sum(d15N)) |> 
#'     ungroup()
#'   
#'   # get baselines proportions
#'   # for two baselines
#'   if (nrow(summary_tp1) == 3) {
#'     # d13C averages at actual TP 
#'     d13C <- sim_results$data$isotopes |> 
#'       filter(type != "source") |> 
#'       group_by(name) |> 
#'       summarise(d13C = mean(d13C)) |> 
#'       ungroup() |> 
#'       pull(d13C)
#'     
#'     # proportion of baseline 1 with Carbon TDF = 0
#'     # alpha <- (d13C[3] - d13C[2]) / (d13C[1] - d13C[2])
#'     
#'     # proportion of baseline 1 including Carbon TDF
#'     mean_deltaC <- mean(sim_results$data$TDF$deltaC)
#'     TPc <- sim_results$estimated_TP$mode
#'     alpha <- (d13C[3] - d13C[2] - mean_deltaC * (TPc - 2)) / (d13C[1] - d13C[2])
#'     # correct cases > 1 (i.e. d13C consumer < d13C baseline 1) and < 1 (i.e. d13C consumer > d13C baseline 2)
#'     alpha <- ifelse(alpha > 1, 1, ifelse(alpha < 0, 0, alpha))
#'     # add proportion of baseline 2
#'     alpha <- c(alpha, 1 - alpha)
#'   } else {
#'     # for one baseline is 1
#'     alpha <- 1
#'   }
#'   
#'   # multiply baselines d15N at source TP by their proportions and sum
#'   d15N_b <- sum(summary_tp1 |> filter(name != "species 1") |> pull(d15N) * alpha)
#'   # get consumer d15N at source TP
#'   d15N_c <- summary_tp1 |> filter(name == "species 1") |> pull(d15N)
#'   # compute difference
#'   dN <- d15N_c - d15N_b
#'   # divide by D15N
#'   mean_deltaN <- mean(sim_results$data$TDF$deltaN)
#'   dN <- dN / mean_deltaN
#'   
#'   return(dN)
#'   
#' }

#' Estimate TP error
#'
#' @param sim_results Output of [run_simulations_baselines()] for a single simulation
#'
#' @return A value
get_tp_error <- function(sim_results) {
  
  # sources summary
  sources <- sim_results$data$sources_summary
  # baselines summary
  baselines <- sim_results$data$isotopes |> 
    filter(type == "baseline") |> 
    group_by(name) |> 
    summarise(d13C_mean = mean(d13C),
              d15N_mean = mean(d15N))
  # mean TDFs
  TDF_C <- mean(sim_results$data$TDF$deltaC)
  TDF_N <- mean(sim_results$data$TDF$deltaN)
  # source relative contributions to consumer
  consumer_p <- sim_results$data$alpha_species$alpha
  
  estimate_tp_error(sources, baselines, TP_b = 2, TDF_N, TDF_C, consumer_p)
  
}

#' Plot the effect of absolute error in TP
#'
#' @param sim_baselines Output of [run_simulations_baselines()]
#' @param performance_pk Output of [get_performance_pk()]
#' @param performance_TP Output of [get_performance_TP()]
#' @param model_baseline Output of [fit_baseline_effect()]
#' @param .width The `.width` argument passed to [ggdist::point_interval()]
#' @param alpha_points Point transparency
#' @param seed Seed for reproducibility
#'
#' @return A `patchwork` object
plot_TP_error_effect <- function(sim_baselines,
                                 performance_pk, 
                                 performance_TP, 
                                 model_baseline, 
                                 .width = ppoints(20), 
                                 alpha_points = 0.2, 
                                 seed = NULL) {
  
  # Get predictions
  preds <- expand.grid(n_baselines = 1:2, 
                       abs_error_TP = seq(min(performance_TP$abs_error),
                                          max(performance_TP$abs_error),
                                          length.out = 51)) |>
    tidybayes::add_epred_draws(model_baseline,
                               re_formula = NA,
                               seed = seed) |> 
    mutate(n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("1 baseline", "2 baselines")))
  
  # Plot effect
  p1 <- performance_pk |> 
    group_by(n_baselines, replicate) |> 
    summarise(max_abs_error = max(abs_error)) |>
    ungroup() |> 
    left_join(performance_TP |> 
                select(n_baselines, replicate, abs_error_TP = abs_error)) |> 
    mutate(replicate = as.factor(replicate),
           n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("1 baseline", "2 baselines"))) |> 
    ggplot(aes(x = abs_error_TP, y = max_abs_error)) + 
    ggdist::stat_lineribbon(data = preds,
                            aes(y = .epred, 
                                fill_ramp = after_stat(.width)),
                            .width = .width,
                            alpha = 0.3,
                            fill = "#2171b5") +
    ggdist::scale_fill_ramp_continuous(range = c(1, 0), 
                                       guide = ggdist::guide_rampbar(title = "CI",
                                                                     to = "#2171b5")) +
    geom_point(alpha = alpha_points) + 
    facet_grid(~ n_baselines, scales = "free") +
    labs(x = "Absolute error in TP",
         y = "Max absolute error in <i>p<sub>k</sub></i>") +
    theme_bw() +
    theme(axis.title.y = ggtext::element_markdown())
  
  # Get TP error
  tp_error <- purrr::map(
    1:2,
    function(i) {
      purrr::map(
        1:100,
        function(j) {
          data.frame(n_baselines = as.factor(i), 
                     replicate = j,
                     estimated_tp_error = get_tp_error(sim_baselines[[i]][[j]]))
        }
      )
    }
  ) |>
    bind_rows()
  
  # Plot relationship between error in TP and difference in d15N
  p2 <- performance_TP |> 
    mutate(tp_error = mode - TP_observed) |> 
    select(n_baselines, replicate, tp_error) |> 
    left_join(tp_error) |> 
    mutate(n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("1 baseline", "2 baselines"))) |> 
    ggplot(aes(x = estimated_tp_error, y = tp_error)) + 
    geom_abline() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = alpha_points) + 
    facet_grid(~ n_baselines, scales = "free") +
    labs(x = "TP-corrected &Delta;&delta;<sup>15</sup>N between consumer and baseline(s)",
         y = "Error in TP") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown())
  
  # Combine plots
  p1 / p2 +
    patchwork::plot_annotation(tag_levels = "A")
}
