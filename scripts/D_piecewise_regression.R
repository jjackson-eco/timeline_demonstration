################################################################################
##                                                                            ##
##           Multivariate signals of population collapse in a                 ##
##                high-throughput ecological experiment                       ##
##                                                                            ##
##            Threshold regression with brms - November 2024                  ##
##                                                                            ##
##                                                                            ##
################################################################################
rm(list = ls())

library(tidyverse)
library(cmdstanr)
library(brms)
library(patchwork)
library(ggridges)
library(tidybayes)

#_______________________________________________________________________________
#### 1. Load and wrangle ####

# load the detrended data
load("data/multivariate_caudatum_no_cycle_for_analysis.RData")

ttc_caudatum_long <- caudatum_no_cycle_for_analysis %>%
  ungroup() %>% 
  # only data from the start point
  filter(time_point >= start_stress_time_point) %>% 
  # this control crashed when it shouldn't have
  filter(!(treatment == "control" & replicate == 7)) %>% 
  # these two replicates didn't work as expected
  filter(!(treatment == "predator" & replicate == 1)) %>% 
  filter(!(treatment == "predator" & replicate == 2)) %>% 
  dplyr::select(id:treatment, 
                behaviour = mean_speed_rep_nocycle,
                morphology = mean_length_rep_nocycle,
                abundance = max_abundance_rep_nocycle) %>% 
  mutate(across(behaviour:abundance, ~ as.numeric(scale(.x)))) %>% 
  pivot_longer(cols = behaviour:abundance, 
               names_to = "timeline_component", 
               values_to = "response") %>% 
  mutate(id_full = paste0(timeline_component, gsub("caudatum", "", id)),
         id_trt = paste0(timeline_component, "_", treatment),
         across(c(id_full, id_trt, replicate, treatment), as.factor)) %>% 
  dplyr::select(id_trt, id_full, time_point, time_date, 
                treatment, timeline_component, response)

# save(ttc_caudatum_long, file = "data/ttc_caudatum_data_long.RData")

## Data for each treatment - Adding in start and end points per ID
tt_control <- ttc_caudatum_long %>% 
  filter(treatment == "control") %>% 
  mutate(start_time_point = min(time_point + 0.1),
         end_time_point = max(time_point - 0.1), .by = id_full)

tt_pollution <- ttc_caudatum_long %>% 
  filter(treatment == "pollution") %>% 
  mutate(start_time_point = min(time_point + 0.1),
         end_time_point = max(time_point - 0.1), .by = id_full)

tt_predator <- ttc_caudatum_long %>% 
  filter(treatment == "predator") %>% 
  mutate(start_time_point = min(time_point + 0.1),
         end_time_point = max(time_point - 0.1), .by = id_full)

#_______________________________________________________________________________
#### 2. Set up piece wise regression ####

timeline_piecewise_formula <- bf(
  response ~ b0 + b1 * (time_point - omega) * inv_logit(omega - time_point) + 
    b2 * (time_point - omega) * inv_logit(time_point - omega),
  alpha ~ -1 + timeline_component + (1|id_full),
  b0 + b1 + b2 ~ 1 + timeline_component + (1|id_full),
  nlf(omega ~ inv_logit(alpha)*end_time_point + start_time_point),
  nl = TRUE
)

timeline_piecewise_priors <- 
  prior(normal(0, 0.2), class = "b", nlpar = "b0") +
  prior(normal(0, 0.2), class = "b", nlpar = "b1") +
  prior(normal(0, 0.2), class = "b", nlpar = "b2") +
  prior(normal(0, 0.2), class = "b", nlpar = "alpha") +
  prior(exponential(5), class = "sigma")

#_______________________________________________________________________________
#### 3. Run the regressions ####

set.seed(420)
timeline_piecewise_control  <- 
  brm(timeline_piecewise_formula, 
      data = tt_control, prior = timeline_piecewise_priors, 
      chains = 3, cores = 3, iter = 3000, warmup = 1500,
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      backend = "cmdstanr")

set.seed(420)
timeline_piecewise_pollution  <- 
  brm(timeline_piecewise_formula, 
      data = tt_pollution, prior = timeline_piecewise_priors, 
      chains = 3, cores = 3, iter = 3000, warmup = 1500,
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      backend = "cmdstanr")

set.seed(420)
timeline_piecewise_predator  <- 
  brm(timeline_piecewise_formula, 
      data = tt_predator, prior = timeline_piecewise_priors, 
      chains = 3, cores = 3, iter = 3000, warmup = 1500,
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      backend = "cmdstanr")

# save(timeline_piecewise_control, timeline_piecewise_pollution, timeline_piecewise_predator,
#      file = "output/piecewise_regression_models_dec24.RData")

#_______________________________________________________________________________
#### 4. Turn point plots ####

load("output/piecewise_regression_models_dec24.RData", verbose = T)

## getting a inverse logit function
inv_logit <- plogis

## Start and end points for equation 
start_ends <- ttc_caudatum_long %>% 
  group_by(id_full) %>% 
  summarise(start_time_point = min(time_point + 0.1),
         end_time_point = max(time_point - 0.1))

start_ends_treatments <- start_ends %>% 
  separate_wider_delim(cols = id_full, delim =  "_", 
                       names = c("timeline_component",
                                 "treatment", "replicate"), 
                       too_few  = "align_start",
                       cols_remove = FALSE) %>% 
  summarise(across(start_time_point:end_time_point, mean),.by = treatment)

## Extracting turn points per replicate - Control
alpha_preds_replicate_control <- 
  as_draws_df(timeline_piecewise_control) %>% 
  dplyr::select(18:44) %>% 
  pivot_longer(cols = 1:27, names_to = "id_full", values_to = "alpha_post") %>% 
  mutate(id_full = gsub("r_id_full__alpha[[]", "", id_full),
         id_full = gsub(",Intercept[]]", "", id_full)) %>% 
  separate_wider_delim(cols = id_full, delim =  "_", 
                       names = c("timeline_component",
                                 "treatment", "replicate"), 
                       too_few  = "align_start",
                       cols_remove = FALSE) %>% 
  left_join(x = ., y = start_ends, by = "id_full") %>% 
  mutate(turn_point = inv_logit_scaled(alpha_post)*end_time_point + start_time_point,
         timeline_component = str_to_title(timeline_component))

alpha_preds_replicate_control$timeline_component <- factor(alpha_preds_replicate_control$timeline_component, 
                                                             levels = c("Behaviour", "Morphology", "Abundance"))
## Extracting turn points per replicate - Pollution
alpha_preds_replicate_pollution <- as_draws_df(timeline_piecewise_pollution) %>% 
  dplyr::select(18:47) %>% 
  pivot_longer(cols = 1:30, names_to = "id_full", values_to = "alpha_post") %>% 
  mutate(id_full = gsub("r_id_full__alpha[[]", "", id_full),
         id_full = gsub(",Intercept[]]", "", id_full)) %>% 
  separate_wider_delim(cols = id_full, delim =  "_", 
                       names = c("timeline_component",
                                 "treatment", "replicate"), 
                       too_few  = "align_start",
                       cols_remove = FALSE) %>% 
  left_join(x = ., y = start_ends, by = "id_full") %>% 
  mutate(turn_point = inv_logit_scaled(alpha_post)*end_time_point + start_time_point,
         timeline_component = str_to_title(timeline_component))

alpha_preds_replicate_pollution$timeline_component <- factor(alpha_preds_replicate_pollution$timeline_component, 
                                                           levels = c("Behaviour", "Morphology", "Abundance"))

## Extracting turn points per replicate - Predator
alpha_preds_replicate_predator <- as_draws_df(timeline_piecewise_predator) %>% 
  dplyr::select(18:41) %>% 
  pivot_longer(cols = 1:24, names_to = "id_full", values_to = "alpha_post") %>% 
  mutate(id_full = gsub("r_id_full__alpha[[]", "", id_full),
         id_full = gsub(",Intercept[]]", "", id_full)) %>% 
  separate_wider_delim(cols = id_full, delim =  "_", 
                       names = c("timeline_component",
                                 "treatment", "replicate"), 
                       too_few  = "align_start",
                       cols_remove = FALSE) %>% 
  left_join(x = ., y = start_ends, by = "id_full") %>% 
  mutate(turn_point = inv_logit_scaled(alpha_post)*end_time_point + start_time_point,
         timeline_component = str_to_title(timeline_component))

alpha_preds_replicate_predator$timeline_component <- factor(alpha_preds_replicate_predator$timeline_component, 
                                         levels = c("Behaviour", "Morphology", "Abundance"))

## Summaries per timeline component
tc_summary_control <- alpha_preds_replicate_control %>% 
  group_by(timeline_component) %>% 
  summarise(md_plt = quantile(turn_point, prob = 0.5),
            md = quantile(turn_point, prob = 0.5) - 329.9333,
            upr = quantile(turn_point, prob = 0.75) - 329.9333,
            lwr = quantile(turn_point, prob = 0.25) - 329.9333)

tc_summary_pollution <- alpha_preds_replicate_pollution %>% 
  group_by(timeline_component) %>% 
  summarise(md_plt = quantile(turn_point, prob = 0.5),
            md = quantile(turn_point, prob = 0.5) - 329.9333,
            upr = quantile(turn_point, prob = 0.75) - 329.9333,
            lwr = quantile(turn_point, prob = 0.25) - 329.9333) %>% 
  mutate(gens = md/22.9798)
  
tc_summary_predator <- alpha_preds_replicate_predator %>% 
  group_by(timeline_component) %>% 
  summarise(md_plt = quantile(turn_point, prob = 0.5),
            md = quantile(turn_point, prob = 0.5) - 329.9333,
            upr = quantile(turn_point, prob = 0.75) - 329.9333,
            lwr = quantile(turn_point, prob = 0.25) - 329.9333)

#_______________________________________________________________________________
#### 5. Turnpoint summary plots ####

pollution_turnpoint_halfeye <- ggplot() +
  stat_halfeye(data = alpha_preds_replicate_pollution, aes(x = turn_point), 
               fill = "#d95f02") +
  geom_vline(data = tc_summary_pollution, aes(group = timeline_component,
                                              xintercept = md_plt)) +
  facet_wrap(~ timeline_component, ncol = 1) +
  scale_x_continuous(limits = c(350,600)) +
  labs(y = "Posterior density", x = "Regression turning point (hours)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())


control_turnpoint_halfeye <- ggplot() +
  stat_halfeye(data = alpha_preds_replicate_control, aes(x = turn_point), 
               fill = "#1b9e77") +
  geom_vline(data = tc_summary_control, aes(group = timeline_component,
                                              xintercept = md_plt)) +
  facet_wrap(~ timeline_component, ncol = 1) +
  scale_x_continuous(limits = c(350,600)) +
  labs(y = "Posterior density", x = "Regression turning point (hours)", tag = "a)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())


predator_turnpoint_halfeye <- ggplot() +
  stat_halfeye(data = alpha_preds_replicate_predator, aes(x = turn_point), 
               fill = "#6654BF") +
  geom_vline(data = tc_summary_predator, aes(group = timeline_component,
                                              xintercept = md_plt)) +
  facet_wrap(~ timeline_component, ncol = 1) +
  scale_x_continuous(limits = c(350,600)) +
  labs(y = "Posterior density", x = "Regression turning point (hours)", tag = "c)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave(pollution_turnpoint_halfeye, file = "output/manuscript/pollution_turnpoint_posterior.jpeg",
       width = 15, height = 24, units = "cm", dpi = 600)

ggsave(control_turnpoint_halfeye | pollution_turnpoint_halfeye + labs(tag = "b)") |
         predator_turnpoint_halfeye, 
       file = "output/manuscript/turnpoint_posterior_all.jpeg",
       width = 27, height = 24, units = "cm", dpi = 600)

