library(targets)
library(ggplot2)
library(ggtext)
library(patchwork)
library(purrr)
library(brms)
library(withr)

# create figures folder
if (!dir.exists("figures")) dir.create("figures")

# Maximum absolute error in estimated source proportions
tar_read(plot_baseline) &
  theme_bw(base_size = 7) +
  theme(strip.text.x = ggtext::element_markdown(face = "italic"),
        axis.title.y = ggtext::element_markdown())
withr::with_seed(25, # seed to reproduce jittered positions
                 ggsave("figures/sim_baseline.png", width = 8, height = 10, units = "cm"))

# Fig 2
tar_read(plot_pk_vs_tp)
ggsave("figures/pk_vs_tp_error.png", width = 18, height = 18, units = "cm")

# Fig 3
tar_read(plot_sample_size_consumer)
ggsave("figures/sim_sample_size_consumer.png", width = 18, height = 18, units = "cm")

# Fig 4
tar_read(plot_sample_size_source)
ggsave("figures/sim_sample_size_source.png", width = 18, height = 18, units = "cm")
# Warnings in fig 3 and 4 because we show the medians but not the intervals

# Example food webs with one or two baselines
purrr::map(
  1:2, 
  ~ tar_read(sim_baselines)[[.x]][[1]][["data"]]$isospace_plot +
    labs(shape = "Sample type",
         color = paste0("Sample name (",
                        .x, 
                        ifelse(.x == 1, " baseline", " baselines"),
                        ")")) +
    ggtitle(paste(.x, ifelse(.x == 1, "baseline", "baselines"))) + 
    theme(plot.title = element_text(hjust = 0.5))
) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.position = "bottom", 
        legend.direction = "vertical")
ggsave("figures/example_food_web_baselines.png", width = 18, height = 16, units = "cm")

# All 100 simulated source polygons
purrr::map(
  1:100,
  function(i) 
    ggplot(tar_read(sim_baselines)[[1]][[i]]$data$sources_summary, 
           aes(x = d13C_mean, y = d15N_mean,
               xmin = d13C_mean - d13C_sd, xmax = d13C_mean + d13C_sd,
               ymin = d15N_mean - d15N_sd, ymax = d15N_mean + d15N_sd)) + 
    geom_polygon(color = "black", fill = NA) +
    geom_linerange(orientation = "x") + 
    geom_linerange(orientation = "y") +
    geom_point() +
    xlim(-12, 12) + 
    ylim(-12, 12) + 
    theme_bw(base_size = 8) +
    labs(x = "&delta;<sup>13</sup>C (&permil;)",
         y = "&delta;<sup>15</sup>N (&permil;)",
         title = i) +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          plot.title = element_text(hjust = 0.5, size = 6),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))) |> 
  patchwork::wrap_plots() + 
  patchwork::plot_layout(axes = "collect")
ggsave("figures/100_polygons.png", width = 18, height = 18, units = "cm")

# Simulated food web
tar_read(sim_food_web)$isospace_plot
ggsave("figures/simulated_food_web.png", width = 16, height = 10, units = "cm")

# PP checks brms models
# Baselines
"model_baseline"|>
  tar_read_raw() |>
  brms::pp_check(ndraws = 100) +
  xlab("Max absolute error")
ggsave("figures/pp_check_model_baselines.png", width = 12, height = 10, units = "cm")
# Sample size
purrr::map2(
  .x = c("model_median_sample_size_consumer",
         "model_CI_sample_size_consumer",
         "model_median_sample_size_source",
         "model_CI_sample_size_source"),
  .y = c("Median <i>p<sub>k</sub></i>",
         "95% CI of <i>p<sub>k</sub></i>",
         "Median <i>p<sub>k</sub></i>",
         "95% CI of <i>p<sub>k</sub></i>"),
  ~ .x |>
    tar_read_raw() |>
    brms::pp_check(ndraws = 100) +
    xlab(.y) +
    theme(axis.title.x = ggtext::element_markdown())
) |>
  patchwork::wrap_plots(ncol = 2) +
  patchwork::plot_layout(guides = "collect", axis_titles = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
ggsave("figures/pp_check_models_sample_size.png", width = 18, height = 14, units = "cm")
