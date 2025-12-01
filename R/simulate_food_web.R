#' Simulate consumer isotope data
#' 
#' Simulate isotope data for a consumer based on known source proportions and TP
#'
#' @param n Number of observations to simulate
#' @param sources A data frame of source isotope data with three columns:
#' - `name` : source names
#' - `d13C_mean`: mean carbon stable isotope ratio
#' - `d13C_sd`: standard deviation of carbon stable isotope ratio
#' - `d15N_mean`: mean nitrogen stable isotope ratio
#' - `d15N_sd`: standard deviation of nitrogen stable isotope ratio
#' @param alpha A numeric vector of known source proportions
#' @param C Total consumption rate (i.e. biomass that the consumer 
#' takes from all sources in the tissue turnover period). Default to 10
#' @param TP_mean Mean trophic position
#' @param TP_sd Standard deviation for trophic position. Default to 0.05
#' @param deltaC A numeric vector. Trophic discrimination factor for C
#' @param deltaN A numeric vector. Trophic discrimination factor for N
#' @param seed The seed passed to [base::set.seed()] to make simulations reproducible.
#' If `NA` (the default), `set.seed` will not be called
#'
#' @return A data frame with `n` rows and two columns, "d13C" and "d15N".
#'
#' @examples 
#' sources <- data.frame(name = paste("source", 1:3),
#'                       d13C_mean = rnorm(3, 0, 5),
#'                       d13C_sd = runif(3, 0.5, 1.5),
#'                       d15N_mean = rnorm(3, 0, 5),
#'                       d15N_sd = runif(3, 0.5, 1.5))
#' alpha <- brms::rdirichlet(1, c(1, 1, 1))
#' # TDFs from Post 2002
#' deltaC <- rnorm(59, 0.39, 1.69)
#' deltaN <- rnorm(107, 3.4, 0.96)
#' simulate_consumer_data(10, sources, alpha, C = 10, TP_mean = 3.5, TP_sd = 0.05, deltaC, deltaN)
#' 
simulate_consumer_data <- function(n, sources, alpha, C = 10, TP_mean, TP_sd = 0.05, deltaC, deltaN, seed = NA) {
  
  if (!is.na(seed)) set.seed(seed)
  
  # Matrix of sample sizes, Z, taken by each consumer i from each source k
  Zik <- rmultinom(n, C, alpha)
  
  # Function to get consumer values for each tracer j
  # j: a string, "d13C" or "d15N"
  get_tracers_values <- function(j) {
    # Consumer i draws Zik samples from source k, resulting in sample means Yijk for each tracer j
    j_mean <- sources |> pull(paste(j, "mean", sep = "_"))
    j_sd <- sources |> pull(paste(j, "sd", sep = "_"))
    Yijk <- purrr::map(
      .x = 1:dim(Zik)[2], 
      function(i) 
        purrr::map(1:dim(Zik)[1],
                   ~ sum(rnorm(Zik[.x, i], j_mean[.x], j_sd[.x])) / Zik[.x, i]
        ) |> 
        purrr::as_vector()
    )
    # Each consumer i's tracer values, Xij, are means of all source samples
    Xij <- purrr::map_dbl(1:length(Yijk), ~ sum(Yijk[[.x]] * Zik[, .x], na.rm = TRUE) / C)
    
    # Add residual error
    res_err <- rnorm(n, 0, 0.1)
    
    return(Xij + res_err)
  }
  
  # Get values and add total fractionation (TDF * number of trophic steps)
  # Sample TP around its mean
  TP_reps <- rnorm(n, TP_mean, TP_sd)
  
  d13C <- get_tracers_values("d13C") + mean(deltaC) * (TP_reps - 1)
  d15N <- get_tracers_values("d15N") + mean(deltaN) * (TP_reps - 1)
  
  data.frame(d13C = d13C, d15N = d15N)
}

#' Simulate a food web
#' 
#' Simulate isotope data for a food web including
#' basal sources, baselines (i.e. primary consumers), 
#' and consumers.
#'
#' @param n_sources Number of sources to simulate. Default to 3
#' @param source_C_range The range in source \u03b4^13^C. A numeric vector of length 2 (minimum and maximum).
#' Ignored if `source_C_mean` is provided
#' @param source_N_range The range in source \u03b4^15^N. A numeric vector of length 2 (minimum and maximum).
#' Ignored if `source_N_mean` is provided
#' @param source_C_mean Numeric vector of length `n_sources` with the mean \u03b4^13^C for each source. 
#' If NULL (default), mean \u03b4^13^C values are sampled from equal (n_sources * 2) intervals along the given `source_C_range`
#' @param source_N_mean Numeric vector of length `n_sources` with the mean \u03b4^15^N for each source.
#' If NULL (default), mean \u03b4^15^N values are sampled from equal (n_sources * 2) intervals along the given `source_N_range`
#' @param source_C_sd Numeric vector of length `n_sources` with the \u03b4^13^C standard deviation for each source
#' @param source_N_sd Numeric vector of length `n_sources` with the \u03b4^15^N standard deviation for each source
#' @param n_obs_sources Numeric vector of length `n_sources` with the number of observations 
#' for each source. Can be a single value if all source have the same number of observations. 
#' Default to 10
#' @param n_baselines Number of baselines (i.e. primary consumers, TP = 2) to simulate. 
#' Must be 1 (default) or 2
#' @param n_obs_baselines Numeric vector of length `n_baselines` with number of observations 
#' for each baseline. Can be a single value if all baselines have the same number of observations. 
#' Default to 10
#' @param n_families Number of consumer families to simulate. Default to 3
#' @param n_species Number of consumer species to simulate. Default to 10.
#' These are randomly split among the families
#' @param n_obs_consumers Numeric vector of length `n_species` with number of observations 
#' for each consumer species. Can be a single value if all consumers have the same number of 
#' observations. Default to 10
#' @param sd_families Variation in diet across families. Default to 1.5
#' @param sd_species Variation in diet across species. Default to 0.5
#' @param TP_sd Standard deviation of trophic position for baselines and consumer species. Default to 0.05
#' @param n_TDF_C Number of observations for carbon's trophic discrimination factor
#' @param mean_TDF_C Mean for carbon's trophic discrimination factor
#' @param sd_TDF_C Standard deviation for carbon's trophic discrimination factor
#' @param n_TDF_N Number of observations for nitrogen's trophic discrimination factor
#' @param mean_TDF_N Mean for nitrogen's trophic discrimination factor
#' @param sd_TDF_N Standard deviation for nitrogen's trophic discrimination factor
#' @param C_baselines,C_consumers Total consumption rate of baselines and consumers
#' (i.e. biomass taken from all sources in the tissue turnover period). Default to 20
#' @param seed The seed passed to [base::set.seed()] to make simulations reproducible.
#' If `NA` (the default), `set.seed` will not be called
#'
#' @return A list with the following elements:
#' - `isotopes`: data frame with isotope data for all observations
#'- `sources_summary`: data frame with isotopic means and standard deviations for each source
#' - `TP_families`: data frame with number of species, mean and
#' standard deviation of trophic position for each family
#' - `TP_species`: data frame with trophic position for each species
#' - `alpha_baselines`: data frame with source proportions for each baseline
#' - `alpha_global`: data frame with global source proportions
#' - `alpha_families`: data frame with source proportions for each family
#' - `alpha_species`: data frame with source proportions for each species
#' - `TDF`: list of length two including the simulated values of 
#' C and N trophic discrimination factors
#' - `TDF_summary`: data frame with mean and standard deviation of trophic 
#' discrimination factors for each source
#' - `isospace_plot`: ggplot object displaying the food web in isotopic space
#' 
#' @examples simulate_food_web()
#' 
simulate_food_web <- function(n_sources = 3,
                              source_C_range = c(-10, 10),
                              source_N_range = c(-10, 10),
                              source_C_mean = NULL,
                              source_N_mean = NULL,
                              source_C_sd = runif(n_sources, 0.5, 1.5),
                              source_N_sd = runif(n_sources, 0.5, 1.5),
                              n_obs_sources = 10,
                              n_baselines = 1,
                              n_obs_baselines = 10,
                              n_families = 3,
                              n_species = 10,
                              n_obs_consumers = 10,
                              sd_families = 1.5,
                              sd_species = 0.5,
                              TP_sd = 0.05,
                              n_TDF_C = 56, 
                              mean_TDF_C = 0.39, 
                              sd_TDF_C = 1.3,
                              n_TDF_N = 107,
                              mean_TDF_N = 3.4, 
                              sd_TDF_N = 0.98,
                              C_baselines = 20,
                              C_consumers = 20,
                              seed = NA) {
  
  # TO DO: implement checks for all arguments
  if (!is.null(source_C_mean) & 
      (!is.numeric(source_C_mean) | 
       length(source_C_mean) != n_sources)) stop("`source_C_mean` must be NULL or a numeric vector of length `n_sources`")
  if (!is.null(source_N_mean) & 
      (!is.numeric(source_N_mean) | 
       length(source_N_mean) != n_sources)) stop("`source_N_mean` must be NULL or a numeric vector of length `n_sources`")
  if (!is.numeric(source_C_sd) | 
       length(source_N_sd) != n_sources) stop("`source_C_sd` must be a numeric vector of length `n_sources`")
  if (!is.numeric(source_N_sd) | 
      length(source_C_sd) != n_sources) stop("`source_N_sd` must be a numeric vector of length `n_sources`")
  if (!is.integer(n_baselines) & !n_baselines %in% 1:2) stop("`n_baselines` must be an integer, either 1 or 2")
  if (n_species < n_families) stop("`n_species` must be higher then or equal to `n_families`")
  if (!length(n_obs_sources %in% c(1, n_sources))) stop("`n_obs_sources` must be a numeric vector of length `n_sources` or 1, 
                                                        if all sources have the same number of observations")
  if (!length(n_obs_baselines %in% c(1, n_baselines))) stop("`n_obs_baselines` must be a numeric vector of length `n_baselines`. 
                                                            It can be a single value when `n_baselines = 2` and the baselines 
                                                            have the same number of observations")
  if (!length(n_obs_consumers %in% c(1, n_species))) stop("`n_obs_consumers` must be a numeric vector of length `n_species`or 1,
                                                          if all species have the same number of observations")
  
  if (!is.na(seed)) set.seed(seed)

  # TDFs
  deltaC <- rnorm(n_TDF_C, mean_TDF_C, sd_TDF_C)
  deltaN <- rnorm(n_TDF_N, mean_TDF_N, sd_TDF_N)
  
  # Sources
  source_names <- paste("source", 1:n_sources, sep = " ")
  if (is.null(source_C_mean)) source_C_mean <- sample(seq(source_C_range[1], source_C_range[2], length = n_sources * 2), 
                                                      size = n_sources, 
                                                      replace = FALSE)
  if (is.null(source_N_mean)) source_N_mean <- sample(seq(source_N_range[1], source_N_range[2], length = n_sources * 2), 
                                                      size = n_sources,
                                                      replace = FALSE)
  # While loop to avoid collinearity
  while (abs(cor(source_C_mean, source_N_mean)) == 1) {
    source_N_mean <- sample(seq(source_N_range[1], source_N_range[2], length = n_sources * 2), 
                            size = n_sources,
                            replace = FALSE)
  }
  
  if (length(n_obs_sources) == 1) n_obs_sources <- rep(n_obs_sources, n_sources)
  
  sources <- purrr::map(
    1:n_sources,
    function(x) {
      data.frame(type = "source",
                 name = factor(source_names[x]),
                 d13C = rnorm(n_obs_sources[x], source_C_mean[x], source_C_sd[x]),
                 d15N = rnorm(n_obs_sources[x], source_N_mean[x], source_N_sd[x]))
    }
  ) |> 
    bind_rows()
  
  sources_summary <- data.frame(name = factor(source_names),
                                d13C_mean = source_C_mean,
                                d13C_sd = source_C_sd,
                                d15N_mean = source_N_mean,
                                d15N_sd = source_N_sd)
  
  # Baselines - Primary consumers (TP = 2)
  
  if (!is.na(seed)) set.seed(seed)
  
  if (n_baselines == 1) {
    # Assume a relatively even source contribution (equal relatively high shape parameter for all sources)
    alpha_baselines <- brms::rdirichlet(1, rep(10, n_sources)) |> 
      as.data.frame() |> 
      rlang::set_names(source_names) |> 
      mutate(name = factor("primary consumer 1"),
             n_obs_baselines = n_obs_baselines) |> 
      tidyr::pivot_longer(cols = starts_with("source"), 
                          names_to = "source", 
                          values_to = "alpha")
  } else {
    # Two distinct baselines covering the range in C and N
    source_minC <- which.min(source_C_mean)
    source_maxC <- which.max(source_C_mean)
    source_minN <- which.min(source_N_mean)
    source_maxN <- which.max(source_N_mean)
    alpha_params_min <- alpha_params_max <- rep(1, n_sources)
    if (source_minC == source_minN | source_minC == source_maxN) {
      alpha_params_min[source_minC] <- 20
    } else if (source_maxC == source_minN) {
      alpha_params_min[c(source_minC, source_maxN)] <- 10
    } else if (source_maxC == source_maxN) {
      alpha_params_min[c(source_minC, source_minN)] <- 10
    } else {
      alpha_params_min[c(source_minC, source_minN)] <- 10
    }
    
    if (source_maxC == source_minN | source_maxC == source_maxN) {
      alpha_params_max[source_maxC] <- 20
    } else if (source_minC == source_minN) {
      alpha_params_max[c(source_maxC, source_maxN)] <- 10
    } else if (source_minC == source_maxN) {
      alpha_params_max[c(source_maxC, source_minN)] <- 10
    } else {
      alpha_params_max[c(source_maxC, source_maxN)] <- 10
    }
    
    if (length(n_obs_baselines) == 1) n_obs_baselines <- rep(n_obs_baselines, 2)
    
    alpha_baselines <- purrr::map(
      list(alpha_params_min, 
           alpha_params_max), 
      ~ brms::rdirichlet(1, .x) |> 
        as.data.frame() |> 
        rlang::set_names(source_names)) |>
      bind_rows() |> 
      rescale_dirichlet(d = n_sources) |> 
      mutate(name = factor(paste("primary consumer", 1:n_baselines)),
             n_obs_baselines = n_obs_baselines) |> 
      tidyr::pivot_longer(cols = starts_with("source"), 
                          names_to = "source", 
                          values_to = "alpha")
  }
  
  # Isotope data for baselines
  baselines <- alpha_baselines |> 
    split(~name) |> 
    purrr::map(
      ~ simulate_consumer_data(n = unique(.x$n_obs_baselines), 
                               sources = sources_summary,
                               alpha = .x$alpha,
                               C = C_baselines,
                               TP_mean = 2, 
                               TP_sd = TP_sd,
                               deltaC = deltaC, 
                               deltaN = deltaN, 
                               seed = NA) # seed already set upfront
    ) |> 
    bind_rows(.id = "name") |> 
    mutate(type = "baseline",
           name = factor(name, levels = levels(alpha_baselines$name))) |> 
    select(type, name, d13C, d15N)
  
  # Higher consumers
  
  if (!is.na(seed)) set.seed(seed)
  
  # Global source proportions
  alpha_global <- brms::rdirichlet(1, rep(1, n_sources)) |> 
    as.data.frame() |> 
    rlang::set_names(source_names)
  
  alpha_global <- rescale_dirichlet(alpha_global, d = n_sources)
  
  # Source proportions at family level
  if (n_families == 1) {
    # if only 1 family use global proportions
    alpha_families <- tibble(family = factor("family 1"),
                             n_species = n_species,
                             source = source_names,
                             alpha = as.numeric(alpha_global))
  } else {
    B_families <- rnorm(n_families, 0, sd_families)
    ilr_global <- compositions::ilr(alpha_global)
    alpha_families <- sapply(B_families, function(i) ilr_global + i) |> 
      t() |> 
      compositions::ilrInv() |> 
      as.data.frame() |> 
      rescale_dirichlet(d = n_sources) |> 
      rlang::set_names(source_names) |> 
      mutate(family = factor(paste("family", 1:n_families),
                             levels = paste("family", 1:n_families)),
             # randomly split the total number of species across families
             n_species = rmultinom(n = 1, 
                                   size = n_species, 
                                   prob = rep.int(1 / 10, n_families)) |> 
               as.vector())
    
    # While loop to avoid families with 0 species
    while(any(alpha_families$n_species == 0)) {
      alpha_families$n_species = rmultinom(n = 1, 
                                           size = n_species, 
                                           prob = rep.int(1 / 10, n_families)) |> 
        as.vector()
    }
    
    alpha_families <- alpha_families |> 
      tidyr::pivot_longer(cols = starts_with("source"), 
                          names_to = "source", 
                          values_to = "alpha")
  }
  
  # Source proportions at species level
  if (n_species == 1) {
    # if only 1 species use global proportions
    alpha_species <- tibble(family = factor("family 1"),
                            species = factor("species 1"),
                            n_obs_consumers = n_obs_consumers,
                            source = source_names,
                            alpha = as.numeric(alpha_global))
  } else {
    family_n_species <- alpha_families |> 
      select(family, n_species) |> 
      distinct()
    B_species <- data.frame(family = rep(family_n_species$family, 
                                         family_n_species$n_species),
                            B = rnorm(n_species, 0, sd_species))
    if (length(n_obs_consumers) == 1) n_obs_consumers <- rep(n_obs_consumers, n_species)
    alpha_species <- purrr::map(
      1:n_families, 
      function(i) {
        ilr_family_i <- alpha_families |> 
          filter(family == paste("family", i)) |> 
          pull(alpha) |> 
          compositions::ilr()
        B_species_family_i <- B_species |> 
          filter(family == paste("family", i)) |> 
          pull(B)
        sapply(B_species_family_i, function(i) compositions::ilrInv(ilr_family_i + i)) |> 
          t() |> 
          as.data.frame() |> 
          rescale_dirichlet(d = n_sources) |> 
          rlang::set_names(source_names) |> 
          mutate(family = paste("family", i))
      }) |> 
      bind_rows() |> 
      mutate(family = factor(family, 
                             levels = levels(alpha_families$family)),
             species = factor(paste("species", 1:n_species),
                              levels = paste("species", 1:n_species)),
             n_obs_consumers = n_obs_consumers) |> 
      tidyr::pivot_longer(cols = starts_with("source"), 
                          names_to = "source", 
                          values_to = "alpha")
  }
  
  # Trophic positions
  # Average TP of families
  TP_families <- alpha_families |> 
    select(family, n_species) |> 
    distinct() |> 
    # constrain lower and upper bound as family averages are unlikely to be at the extremes
    mutate(TP_mean = runif(n = n_families, 2.5, 4.5),
           TP_sd = runif(n = n_families, 0.1, 0.3))
  
  # Randomly sample around families averages
  TP_species <- alpha_species |> 
    select(family, species) |> 
    distinct() |> 
    bind_cols(TP = TP_families |>
                split(~family) |> 
                # use truncated normal distribution to avoid species with TP < 2 or > 5
                purrr::map(~ EnvStats::rnormTrunc(.x$n_species, .x$TP_mean, .x$TP_sd, min = 2, max = 5)) |> 
                purrr::as_vector())
  
  # Generate data for higher consumers
  higher_consumers <- alpha_species |> 
    left_join(TP_species) |> 
    split(~species) |> 
    purrr::map(
      ~ simulate_consumer_data(unique(.x$n_obs_consumers), 
                               sources = sources_summary,
                               C = C_consumers,
                               alpha = .x$alpha,
                               TP_mean = .x$TP[1], 
                               TP_sd = TP_sd,
                               deltaC = deltaC, 
                               deltaN = deltaN, 
                               seed = NA) # seed already set upfront
    ) |> 
    bind_rows(.id = "species") |> 
    mutate(species = factor(species,
                            levels = levels(alpha_species$species)))
  
  # Add family column to isotopes data
  higher_consumers <- alpha_species |> 
    select(family, species) |> 
    distinct()  |> 
    mutate(type = "higher consumer", 
           name = paste(family, species, sep = "_"),
           name = factor(name, 
                         levels = unique(name))) |> 
    left_join(higher_consumers, by = "species")
  
  # Plot
  data <- bind_rows(sources,
                    baselines, 
                    higher_consumers)
  
  p <- data |> 
    group_by(type, name) |> 
    summarise(across(where(is.double), list(mean = mean, sd = sd))) |> 
    ggplot(aes(x = d13C_mean, y = d15N_mean, shape = type, color = name)) +
    geom_point(data = data, 
               aes(x = d13C, y = d15N),
               alpha = 0.5) +
    geom_linerange(aes(xmin = d13C_mean - d13C_sd,
                       xmax = d13C_mean + d13C_sd)) +
    geom_pointrange(aes(ymin = d15N_mean - d15N_sd,
                        ymax = d15N_mean + d15N_sd)) +
    labs(x = "&delta;<sup>13</sup>C (&permil;)",
         y = "&delta;<sup>15</sup>N (&permil;)") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown())
  
  return(list(isotopes = data,
              sources_summary = sources_summary,
              TP_families = TP_families,
              TP_species = TP_species,
              alpha_baselines = alpha_baselines,
              alpha_global = alpha_global,
              alpha_families = alpha_families,
              alpha_species = alpha_species,
              TDF = list(deltaC = deltaC,
                         deltaN = deltaN),
              TDF_summary = data.frame(source = source_names,
                                       Meand13C = mean(deltaC),
                                       SDd13C = sd(deltaC),
                                       Meand15N = mean(deltaN),
                                       SDd15N = sd(deltaN)),
              isospace_plot = p))
}


#' Rescale Dirichlet proportions to avoid values < 0.05
#'
#' @param data A data frame with one column for each category and one row for each mixture. 
#' Each row must sum up to 1
#' @param d Number of categories
#'
#' @return A data frame
rescale_dirichlet <- function(data, d) {
  scale <- (1/d) / 0.05
  ((data * (scale - 1)) + 1/d) / scale
}
