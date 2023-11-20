################################################################################
##                                                                            ##
##           Multivariate signals of population collapse in a                 ##
##                high-throughput ecological experiment                       ##
##                                                                            ##
##          Cerini, F., Jackson, J., Childs, D.Z., & Clements, C.F.           ##
##                                                                            ##
##                  Analysis pipeline - November 2023                         ##
##                                                                            ##
##                                                                            ##
################################################################################
rm(list = ls())

library(tidyverse)
library(EWSmethods)
library(mgcv)
library(MASS)
library(patchwork)

## set the time point of the stress treatment in the time series
start_stress_time_point <- 329.9333

#_______________________________________________________________________________
#### 1. Load and wrangle ####

# load the detrended data
load("Data//multivariate_caudatum_no_cycle_for_analysis.RData")

# wrangle
ttc_caudatum <- caudatum_no_cycle_for_analysis %>%
  # only data from the start point
  filter(time_point >= start_stress_time_point) %>% 
  # this control crashed when it shouldn't have
  filter(!(treatment == "control" & replicate == 7)) %>% 
  # these two replicates didn't work as expected
  filter(!(treatment == "predator" & replicate == 1)) %>% 
  filter(!(treatment == "predator" & replicate == 2)) %>% 
  mutate(replicate = as.factor(replicate)) 

# raw plot - check
ggplot(ttc_caudatum, aes(x = time_point, 
                         y = mean_speed_rep_nocycle, 
                         colour = treatment, 
                         group = interaction(treatment, replicate))) +
  geom_line(alpha = 0.6,  show.legend = F) +
  facet_wrap(~ treatment) +
  theme_bw()

# min and max timepoint
mintp <- min(ttc_caudatum$time_point)
maxtp <- max(ttc_caudatum$time_point)

#_______________________________________________________________________________
#### 2. Autocorrelation ####

### iterating through ids and pulling out autocorrelation
autocorr_id <- expand.grid(id = unique(ttc_caudatum$id),
                           lag = 1:10, 
                           cor_ab = NA, cor_ms = NA, cor_bs = NA,
                           ci = NA,
                           stringsAsFactors = F)

for(i in 1:nrow(autocorr_id)){
  
  #  right data
  cid = autocorr_id[i,"id"]
  cdat = filter(ttc_caudatum, id == cid)
  
  # time series for each timeline
  n_id = nrow(cdat)
  ts_ab = ts(cdat$max_abundance_rep_nocycle)
  ts_ms = ts(cdat$mean_speed_rep_nocycle)
  ts_bs = ts(cdat$mean_length_rep_nocycle)
  
  # partial autocorrelation for each time-series
  acf_ab = pacf(ts_ab, lag.max = 10, type = "correlation", plot = F)
  acf_ms = pacf(ts_ms, lag.max = 10, type = "correlation", plot = F)
  acf_bs = pacf(ts_bs, lag.max = 10, type = "correlation", plot = F)
  
  # fill in with auto-correllation results
  autocorr_id[which(autocorr_id$id == cid), "cor_ab"] <- as.numeric(acf_ab$acf)
  autocorr_id[which(autocorr_id$id == cid), "cor_ms"] <- as.numeric(acf_ms$acf)
  autocorr_id[which(autocorr_id$id == cid), "cor_bs"] <- as.numeric(acf_bs$acf)
  autocorr_id[which(autocorr_id$id == cid), "ci"] <- qnorm((1 + 0.95)/2)/sqrt(n_id)
  
  cat("\r", "Your job is", round((i/nrow(autocorr_id)), 3)*100, "% Percent complete       ")
}

## Visualize
autocorr_summary <- autocorr_id %>% 
  mutate(mean_ci = mean(ci)) %>% 
  group_by(lag) %>% 
  summarise(mean_ci = mean_ci[1],
            mn_ab = mean(cor_ab), mn_ms = mean(cor_ms), mn_bs = mean(cor_bs),
            se_ab = sd(cor_ab)/sqrt(n()), 
            se_ms = sd(cor_ms)/sqrt(n()), 
            se_bs = sd(cor_bs)/sqrt(n())) 

#abundance time series
ab_auto_plot <- ggplot(autocorr_summary, aes(x = lag, y = mn_ab)) +
  geom_hline(yintercept = autocorr_summary$mean_ci[1], linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -autocorr_summary$mean_ci[1], linetype = "dashed") +
  geom_errorbar(width = 0, aes(ymax = mn_ab + se_ab, ymin = mn_ab - se_ab)) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 0:15) +
  labs(x = "Lag", y = "Mean partial autocorrelation", title = "Abundance") +
  theme_classic(base_size = 12)
ab_auto_plot

# speed time series
ms_auto_plot <- ggplot(autocorr_summary, aes(x = lag, y = mn_ms)) +
  geom_hline(yintercept = autocorr_summary$mean_ci[1], linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -autocorr_summary$mean_ci[1], linetype = "dashed") +
  geom_errorbar(width = 0, aes(ymax = mn_ms + se_ms, ymin = mn_ms - se_ms)) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 0:15) +
  labs(x = "Lag", y = "Mean partial autocorrelation", title = "Movement speed") +
  theme_classic(base_size = 12)
ms_auto_plot

# body size time series
bs_auto_plot <- ggplot(autocorr_summary, aes(x = lag, y = mn_bs)) +
  geom_hline(yintercept = autocorr_summary$mean_ci[1], linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -autocorr_summary$mean_ci[1], linetype = "dashed") +
  geom_errorbar(width = 0, aes(ymax = mn_bs + se_bs, ymin = mn_bs - se_bs)) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 0:15) +
  labs(x = "Lag", y = "Mean partial autocorrelation", title = "Body size") +
  theme_classic(base_size = 12)
bs_auto_plot

# ggsave(ab_auto_plot / ms_auto_plot / bs_auto_plot, filename = "temporal_autocorrelation_averages.jpeg",
#        width = 14, height = 30, units = "cm", dpi = 800)

#_______________________________________________________________________________
#### 3. Model selection for k ####

## Simple model with just a treatment by time point effect and random id effect (as we have all the treatments). 
iterdat_k <- expand_grid(k = 4:30, AIC_ab = NA, BIC_ab = NA, 
                         AIC_ms  = NA, BIC_ms = NA,
                         AIC_bs = NA, BIC_bs = NA)

# turn the id and treatment variables to factors
ttc_caudatum$id  <- as.factor(ttc_caudatum$id)
ttc_caudatum$treatment  <- as.factor(ttc_caudatum$treatment)

for(i in 1:nrow(iterdat_k)){
  
  # current k
  ck = iterdat_k[i,]$k
  
  # models
  abmod = gam(max_abundance_rep_nocycle ~ s(time_point, by = treatment, k = ck, bs = "tp") + s(id, bs = "re"),
              data = ttc_caudatum, family = "gaussian", method = "REML")
  
  msmod = gam(mean_speed_rep_nocycle ~ s(time_point, by = treatment, k = ck, bs = "tp") + s(id, bs = "re"),
              data = ttc_caudatum, family = "gaussian", method = "REML")
  
  bsmod = gam(mean_length_rep_nocycle ~ s(time_point, by = treatment, k = ck, bs = "tp") + s(id, bs = "re"),
              data = ttc_caudatum, family = "gaussian", method = "REML")
  
  
  # filling in the data
  iterdat_k[i,"AIC_ab"] <- AIC(abmod)
  iterdat_k[i,"BIC_ab"] <- BIC(abmod)
  iterdat_k[i,"AIC_ms"] <- AIC(msmod)
  iterdat_k[i,"BIC_ms"] <- BIC(msmod)
  iterdat_k[i,"AIC_bs"] <- AIC(bsmod)
  iterdat_k[i,"BIC_bs"] <- BIC(bsmod)
  
  cat("\r", "Your job is", round((i/nrow(iterdat_k)), 2)*100, "% Percent complete       ")
}

## Visualize change in AIC based on k values
k_mod_plot <- iterdat_k %>% 
  pivot_longer(-k) %>% 
  mutate(resp = gsub(".IC[_]", "", name),
         quant = gsub("[_].*", "", name),
         resp = case_when(
           resp == "ab" ~ "Abundance",
           resp == "ms" ~ "Movement speed",
           resp == "bs" ~ "Body size")) %>% 
  filter(quant == "AIC") %>% 
  ggplot(aes(x = k, y = value)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(4, 40, by = 1)) +
  labs(x = "Basis dimension (k)", y = "Information Criterion value") +
  facet_wrap(quant ~ resp, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))
k_mod_plot

## We can use for Abundance k = 17 
## for Body size we go k = 11 
## for speed k = 12

#_______________________________________________________________________________
#### 4. Models and predictions for the timeline ####

# 1. create all the variables to iterate over
iteration_vars <- expand_grid(response = c("max_abundance_rep_nocycle", 
                                           "mean_length_rep_nocycle",
                                           "mean_speed_rep_nocycle"),
                              treatment = c("control", "pollution", "predator")) %>% 
  mutate(timeline_component = 
           case_when(response == "max_abundance_rep_nocycle" ~ "abundance",
                     response == "mean_length_rep_nocycle" ~ "morphology",
                     response == "mean_speed_rep_nocycle" ~ "behaviour"))

## list to store posterior simulations
iteration_list <- vector(mode = "list")

# 2. iterate through and do the model + predictions for each
prediction_data_all <- lapply(X= 1:nrow(iteration_vars), function(x){
  
  # current variables
  cresp = iteration_vars$response[x]
  ctrt = iteration_vars$treatment[x]
  ctml = iteration_vars$timeline_component[x]
  
  # correct data
  cdat = filter(ttc_caudatum, treatment == ctrt)
  
  # basis dimension
  cK = case_when(ctml == "abundance" ~ 17,
                 ctml == "morphology" ~ 11,
                 ctml == "behaviour" ~ 12)
  
  # create the formula
  cform = formula(paste0(cresp, "~ s(time_point, bs = 'tp', k = ", cK,
                         ") + s(time_point, by = replicate, m = 1) + s(replicate, bs = 're')"))
  
  # run the model
  cmod = gamm(cform, correlation = corARMA(form = ~1 | time_point, p = 1),
              data = cdat, family = "gaussian", method = "REML")
  
  # excluding all but time-point effect for some predictions
  terms_exclude <- names(cmod$gam$sp)[2:11]
  
  # create newdata
  cnewdata <- expand_grid(timeline_component = ctml, treatment = ctrt,
                          time_point  = seq(from = mintp, to  = maxtp, length.out = 500), 
                          replicate  = unique(cdat$replicate)) %>%
    filter(time_point <= max(cdat$time_point))
  
  # predict and add to data
  cfit = predict(cmod$gam, newdata = cnewdata, se.fit = T, type = "response", exclude = terms_exclude)
  
  # posterior simulation
  # Xp matrix
  cXp = predict(cmod$gam, cnewdata, type="lpmatrix", exclude = terms_exclude) 
  
  # coefficients
  cbeta = coef(cmod$gam)
  cVb = vcov(cmod$gam)
  
  # simulate 1000 posterior coefficients
  set.seed(420)
  cmrand = mvrnorm(1000, mu = cbeta, Sigma = cVb) 
  cpred_sim = matrix(nrow = 1000, ncol = nrow(cnewdata))
  
  # from prediction matrix, multiply to get predictions under uncertainty
  for (i in seq_len(1000)) {
    cpred_post  =  c(cXp %*% cmrand[i,])
    cpred_sim[i,] = cpred_post
  }
  
  # posterior simulation output
  post_sim_all <- cpred_sim %>% as_tibble() %>% 
    slice(1:100) %>% mutate(sim = 1:n()) %>% 
    pivot_longer(-sim, names_to = "pred_row", values_to = "posterior_simulation") %>% 
    bind_cols(., slice(cnewdata, rep(1:n(), 100)))
  
  # create predictions
  cpred = cnewdata %>% 
    mutate(fit = as.numeric(cfit$fit), 
           upr = apply(cpred_sim, 2, quantile, prob = 0.975),
           lwr = apply(cpred_sim, 2, quantile, prob = 0.025)) %>% 
    left_join(post_sim_all, by = c("timeline_component", "treatment", 
                                   "time_point", "replicate"))
  
  print(paste0(ctml, " - ", ctrt, " - modelling completed [", x, " out of ", nrow(iteration_vars), "]"))
  return(cpred)}) %>% 
  bind_rows()

# 4. Aggregating predictions for each time point
prediction_data <- prediction_data_all %>% 
  group_by(timeline_component, treatment) %>% 
  mutate(posterior_simulation = as.numeric(scale(posterior_simulation, scale = F))) %>% 
  ungroup() %>% 
  group_by(timeline_component, treatment, time_point, replicate) %>% 
  summarise(fit = mean(posterior_simulation), 
            lwr = quantile(posterior_simulation, 0.025),
            upr = quantile(posterior_simulation, 0.975)) %>% 
  ungroup() %>% 
  group_by(timeline_component, treatment, time_point) %>% 
  summarise(fit  = mean(fit), upr = max(upr), lwr = min(lwr)) %>%
  ungroup() %>% 
  group_by(timeline_component, treatment) %>% 
  mutate(fit1 = fit[1],
         across(fit:lwr, ~ .x + (1 - fit1)))

#_______________________________________________________________________________
#### 5. Assessing CI splits for each timeline component ####

split_data <- expand_grid(treatment = c("pollution", "predator"),
                          timeline_component = c("behaviour", "morphology", "abundance"),
                          split_point = NA, split_point_lwr = NA, split_point_upr = NA,
                          split_point_diff = NA)

for(i in 1:nrow(split_data)){
  
  # current vars
  ctlc = split_data[i,]$timeline_component
  ctrt = split_data[i,]$treatment
  
  # current data
  cntrl_data = filter(prediction_data, timeline_component == ctlc & 
                        treatment == "control")
  test_data = filter(prediction_data, timeline_component == ctlc &
                       treatment == ctrt)
  
  # join add difference and take split point
  comb_data = left_join(cntrl_data, test_data, by = "time_point") %>% 
    drop_na() %>% 
    mutate(abs_diff = abs(lwr.x - upr.y)) 
  
  # split_point not excluding replicate effects
  min_abs_diff = min(comb_data$abs_diff)
  split_point_data = filter(comb_data, abs_diff == min_abs_diff)
    
  # timepoint in original data
  split_data[i,]$split_point <- split_point_data$time_point
  
  # current data 2
  cntrl_data2 = filter(prediction_data_all, timeline_component == ctlc & 
                        treatment == "control")
  test_data2 = filter(prediction_data_all, timeline_component == ctlc &
                       treatment == ctrt)
  
  # join two big difference and calculate differences in posterior
  comb_data2 = left_join(cntrl_data2, test_data2, 
                         by = c("time_point","replicate","sim")) %>% 
    drop_na() %>% 
    mutate(diff = posterior_simulation.x - posterior_simulation.y) 
  
  # credible interval of posterior mean differences per timepoint crosses 0?
  split_data[i,]$split_point_diff <- 
    comb_data2 %>% 
    group_by(time_point) %>% 
    summarise(lwr = quantile(diff, 0.025),
              upr = quantile(diff, 0.975)) %>% 
    ungroup() %>% 
    mutate(crosses_0 = if_else(lwr > 0 | upr < 0, 1, 0)) %>% 
    filter(crosses_0 == 1) %>% 
    pull(time_point) %>% 
    min()
  
  ## Assessing the upr and lwr bounds of the splits based on posterior resampling (100 samples)
  split_points_resample <- vector(length = 100)
  set.seed(420)
  for(x in 1:100){
    # join add difference and take split point
    comb_data_rs = left_join(cntrl_data, test_data, by = "time_point") %>% 
      drop_na() %>% 
      # sample 10% of the data
      slice_sample(prop = 0.1) %>% 
      mutate(abs_diff = abs(lwr.x - upr.y)) 
    
    # split_point not excluding replicate effects
    min_abs_diff_rs = min(comb_data_rs$abs_diff)
    split_point_data_rs = filter(comb_data_rs, abs_diff == min_abs_diff_rs)
    
    # timepoint in original data
    split_points_resample[x] <- split_point_data_rs$time_point
  }
  
  split_data[i,]$split_point_lwr <- quantile(split_points_resample, prob = 0.025)
  split_data[i,]$split_point_upr <- quantile(split_points_resample, prob = 0.975)
  
}

## split data in generations after establishment to results
split_data %>% 
  mutate(across(split_point:split_point_upr, ~ .x/23)) %>% 
  write_csv(x = ., file = "Results/manuscript/split_results.csv")

#_______________________________________________________________________________
#### 6. Plotting - Main Pollution plot ####

rawdat_pol <- filter(ttc_caudatum, treatment != "predator") %>% 
  group_by(id) %>% 
  mutate(across(max_abundance_rep_nocycle:mean_length_rep_nocycle, 
                ~ scale(.x, scale= F)  + (1 - scale(.x, scale = F)[1])))
pred_pol <- filter(prediction_data, treatment != "predator")
split_pol <- filter(split_data, treatment == "pollution")

## 6a. Behaviour ##
ms_plt_pol <- ggplot(rawdat_pol, aes(x = time_point, y = mean_speed_rep_nocycle,
                                     colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(pred_pol, 
                            timeline_component == "behaviour"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  annotate("rect", alpha = 0.3, ymin = 0.92, ymax = 1.025,
           xmin = split_data$split_point_lwr[1], 
           xmax = split_data$split_point_upr[1]) +
  geom_vline(xintercept = split_data$split_point[1]) +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", "control" = "#1b9e77"),
                     labels = c("Control", "Pollution"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(330, 600, by = 30),
                     limits = c(mintp-1, 455)) +
  coord_cartesian(ylim = c(0.925,1.02)) +
  labs(x= NULL, y = "Scaled Movement Speed", tag = "a)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

## 6b. Morphology ##
bs_plt_pol <- ggplot(rawdat_pol, aes(x = time_point, y = mean_length_rep_nocycle,
                                     colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(pred_pol, 
                            timeline_component == "morphology"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  annotate("rect", alpha = 0.3, ymin = -80, ymax = 25,
           xmin = split_data$split_point_lwr[2], 
           xmax = split_data$split_point_upr[2]) +
  geom_vline(xintercept = split_data$split_point[2]) +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", "control" = "#1b9e77"),
                     labels = c("Control", "Pollution"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(330, 600, by = 30),
                     limits = c(mintp - 1, 455)) +
  coord_cartesian(ylim = c(-62,19)) +
  labs(x= NULL, 
       y = "Scaled Body Length", tag = "b)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

## 6c. Abundance ##
ab_plt_pol <- ggplot(rawdat_pol, aes(x = time_point, y = max_abundance_rep_nocycle,
                                     colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(pred_pol, 
                            timeline_component == "abundance"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  annotate("rect", alpha = 0.3, ymin = -200, ymax = 50,
           xmin = split_data$split_point_lwr[3], 
           xmax = split_data$split_point_upr[3]) +
  geom_vline(xintercept = split_data$split_point[3]) +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", "control" = "#1b9e77"),
                     labels = c("Control", "Pollution"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(330, 600, by = 30),
                     limits = c(mintp - 1, 455)) +
  coord_cartesian(ylim = c(-150,40)) +
  labs(x= "Experimentation time (hours)", 
       y = "Scaled abundance", tag = "c)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave(ms_plt_pol / bs_plt_pol / ab_plt_pol + 
         plot_layout(guides = "collect") & theme(legend.position = 'top'), 
       filename = "Results/predictions_pollution.jpeg",
       width = 15, height = 23, units ="cm", dpi = 600)

#_______________________________________________________________________________
#### 7. Plotting - Predator plot ####

rawdat_pred <- filter(ttc_caudatum, treatment != "pollution") %>% 
  group_by(id) %>% 
  mutate(across(max_abundance_rep_nocycle:mean_length_rep_nocycle, 
                ~ scale(.x, scale= F)  + (1 - scale(.x, scale = F)[1])))
pred_pred <- filter(prediction_data, treatment != "pollution")
split_pred <- filter(split_data, treatment == "predator")

## 6a. Behaviour ##
ms_plt_pred <- ggplot(rawdat_pred, aes(x = time_point, y = mean_speed_rep_nocycle,
                                       colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(pred_pred, 
                            timeline_component == "behaviour"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  scale_color_manual(name = "Treatment",
                     values = c("predator" ="#6654BF", "control" = "#1b9e77"),
                     labels = c("Control", "Predator"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(330, 600, by = 30),
                     limits = c(mintp-1, 455)) +
  coord_cartesian(ylim = c(0.98,1.03)) +
  labs(x= NULL, y = "Scaled Movement Speed", tag = "a)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

## 6b. Morphology ##
bs_plt_pred <- ggplot(rawdat_pred, aes(x = time_point, y = mean_length_rep_nocycle,
                                       colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(pred_pred, 
                            timeline_component == "morphology"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  scale_color_manual(name = "Treatment",
                     values = c("predator" ="#6654BF", "control" = "#1b9e77"),
                     labels = c("Control", "Predator"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(330, 600, by = 30),
                     limits = c(mintp - 1, 455)) +
  coord_cartesian(ylim = c(-20,10)) +
  labs(x= NULL, 
       y = "Scaled Body Length", tag = "b)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

## 6c. Abundance ##
ab_plt_pred <- ggplot(rawdat_pred, aes(x = time_point, y = max_abundance_rep_nocycle,
                                       colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(pred_pred, 
                            timeline_component == "abundance"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  annotate("rect", alpha = 0.3, ymin = -100, ymax = 40,
           xmin = split_data$split_point_lwr[6], 
           xmax = split_data$split_point_upr[6]) +
  geom_vline(xintercept = split_data$split_point[6]) +
  scale_color_manual(name = "Treatment",
                     values = c("predator" ="#6654BF", "control" = "#1b9e77"),
                     labels = c("Control", "Predator"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(330, 600, by = 30),
                     limits = c(mintp - 1, 455)) +
  coord_cartesian(ylim = c(-90,30)) +
  labs(x= "Experimentation time (hours)", 
       y = "Scaled abundance", tag = "c)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave(ms_plt_pred / bs_plt_pred / ab_plt_pred + 
         plot_layout(guides = "collect") & theme(legend.position = 'top'), 
       filename = "Results/predictions_predator.jpeg",
       width = 15, height = 23, units ="cm", dpi = 600)

#_______________________________________________________________________________
#### 8. EWS signals for timeline components ####

## the maximum possible time-point of experiment (end of pollution)
end_pol <- ttc_caudatum %>% 
  filter(treatment  == "pollution") %>% 
  pull(time_point) %>% 
  max()

ttc_long <- ttc_caudatum %>% ungroup() %>% 
  dplyr::select(id:treatment, max_abundance_rep_nocycle, 
                mean_length_rep_nocycle,
                mean_speed_rep_nocycle) %>% 
  pivot_longer(max_abundance_rep_nocycle:mean_speed_rep_nocycle, 
               names_to = "response", values_to = "value") %>% 
  mutate(timeline_component = 
           case_when(response == "max_abundance_rep_nocycle" ~ "abundance",
                     response == "mean_length_rep_nocycle" ~ "morphology",
                     response == "mean_speed_rep_nocycle" ~ "behaviour"),
         id = paste0(id, "_", timeline_component)) %>% 
  filter(time_point <= end_pol) %>% 
  group_by(id) %>% 
  mutate(time_before_end = seq((n()-1)*3, 0, by = -3)) %>% 
  ungroup()

## iterate over IDs and extract EWS metrics
ids <- unique(ttc_long$id)
ews_ttc_rawdat <- lapply(ids, function(x){
  
  # correct data
  cdat = filter(ttc_long, id == x) %>% 
    dplyr::select(time_point, value)
  
  # pull out ews metrics
  cews = uniEWS(data = cdat,
                metrics =  c("cv", "SD", "acf"),
                method = "rolling", winsize = 50) 
  
  return(tibble(id = x,
                cews$EWS$raw))}) %>% 
  bind_rows() %>% 
  rename(time_point = time,
         count_used = count.used)

# bind with original data
ews_ttc <- ttc_long %>% 
  dplyr::select(id, replicate, treatment, timeline_component,
               time_date, time_point, time_before_end) %>% 
  left_join(x =., y = ews_ttc_rawdat, by = c("id", "time_point")) %>% 
  na.omit()

#_______________________________________________________________________________
#### 9. EWS models, predictions and splits ####

# 1. create all the variables to iterate over
iteration_ews <- expand_grid(response = c("cv", "SD", "acf"), 
                             treatment = c("control", "pollution", "predator"),
                             timeline_component = c("behaviour", "morphology", "abundance"))

## list to store posterior simulations
iteration_list <- vector(mode = "list")

## time points for this analysis
minews <- min(ews_ttc$time_before_end)
maxews <- max(ews_ttc$time_before_end)

# 2. iterate through and do the model + predictions for each
prediction_ews_all <- lapply(X= 1:nrow(iteration_ews), function(x){
  
  # current variables
  cresp = iteration_ews$response[x]
  ctrt = iteration_ews$treatment[x]
  ctml = iteration_ews$timeline_component[x]
  
  # correct data
  cdat = filter(ews_ttc, treatment == ctrt & timeline_component == ctml)
  
  # basis dimension  - just going with 10 here to be consistent
  cK = 10
  
  # create the formula
  cform = formula(paste0(cresp, "~ s(time_before_end, bs = 'tp', k = ", cK,
                         ") + s(time_before_end, by = replicate, m = 1) + s(replicate, bs = 're')"))
  
  # run the model
  cmod = gamm(cform, correlation = corARMA(form = ~1 | time_before_end, p = 1),
              data = cdat, family = "gaussian", method = "REML")
  
  # excluding all but time-point effect for some predictions
  terms_exclude <- names(cmod$gam$sp)[2:11]
  
  # create newdata
  cnewdata <- expand_grid(timeline_component = ctml, treatment = ctrt,
                          time_before_end  = seq(from = minews, to  = maxews, length.out = 500), 
                          replicate  = unique(cdat$replicate)) %>%
    filter(time_before_end <= max(cdat$time_before_end)) %>% 
    mutate(response = cresp)
  
  # predict and add to data
  cfit = predict(cmod$gam, newdata = cnewdata, se.fit = T, 
                 type = "response", exclude = terms_exclude)
  
  # posterior simulation
  # Xp matrix
  cXp = predict(cmod$gam, cnewdata, type="lpmatrix", exclude = terms_exclude) 
  
  # coefficients
  cbeta = coef(cmod$gam)
  cVb = vcov(cmod$gam)
  
  # simulate 1000 posterior coefficients
  set.seed(420)
  cmrand = mvrnorm(1000, mu = cbeta, Sigma = cVb) 
  cpred_sim = matrix(nrow = 1000, ncol = nrow(cnewdata))
  
  # from prediction matrix, multiply to get predictions under uncertainty
  for (i in seq_len(1000)) {
    cpred_post  =  c(cXp %*% cmrand[i,])
    cpred_sim[i,] = cpred_post
  }
  
  # posterior simulation output
  post_sim_all <- cpred_sim %>% as_tibble() %>% 
    slice(1:100) %>% mutate(sim = 1:n()) %>% 
    pivot_longer(-sim, names_to = "pred_row", values_to = "posterior_simulation") %>% 
    bind_cols(., slice(cnewdata, rep(1:n(), 100))) %>% 
    dplyr::select(-response)
  
  # create predictions
  cpred = cnewdata %>% 
    mutate(fit = as.numeric(cfit$fit), 
           upr = apply(cpred_sim, 2, quantile, prob = 0.975),
           lwr = apply(cpred_sim, 2, quantile, prob = 0.025)) %>% 
    left_join(post_sim_all, by = c("timeline_component", "treatment", 
                                   "time_before_end", "replicate"))
  
  print(paste0(cresp, " - ", ctml, " - ", ctrt, " - modelling completed [", x, " out of ", nrow(iteration_ews), "]"))
  return(cpred)}) %>% 
  bind_rows() 

# 3. Aggregating predictions for each time point
prediction_ews <- prediction_ews_all %>% 
  group_by(response, timeline_component, treatment, time_before_end, replicate) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(response, timeline_component, treatment, time_before_end) %>% 
  summarise(time_before_end = time_before_end[1],
            fit  = mean(fit), upr = max(upr), lwr = min(lwr)) %>%
  ungroup() 

# 4. Confidence interval splits
split_ews <- expand_grid(treatment = c("pollution", "predator"),
                         response = c("acf", "cv", "SD"),
                          timeline_component = c("behaviour", "morphology", "abundance"),
                          split_point = NA, split_point_lwr = NA, split_point_upr = NA)

for(i in 1:nrow(split_ews)){
  
  # current vars
  ctrt = split_ews[i,]$treatment
  cresp = split_ews[i,]$response
  ctlc = split_ews[i,]$timeline_component
  
  # current data
  cntrl_data = filter(prediction_ews, response == cresp & 
                        timeline_component == ctlc & 
                        treatment == "control")
  test_data = filter(prediction_ews, response == cresp & 
                       timeline_component == ctlc &
                       treatment == ctrt)
  
  # join add difference and take split point
  comb_data = left_join(cntrl_data, test_data, by = "time_before_end") %>% 
    drop_na() %>% 
    mutate(abs_diff = abs(upr.x - lwr.y)) 
  
  # split_point not excluding replicate effects
  min_abs_diff = min(comb_data$abs_diff)
  split_point_data = filter(comb_data, abs_diff == min_abs_diff)
  
  # timepoint in original data
  split_ews[i,]$split_point <- split_point_data$time_before_end
  
  ## Assessing the upr and lwr bounds of the splits based on posterior resampling (100 samples)
  split_points_resample <- vector(length = 100)
  set.seed(420)
  for(x in 1:100){
    # join add difference and take split point
    comb_data_rs = left_join(cntrl_data, test_data, by = "time_before_end") %>% 
      drop_na() %>% 
      # sample 10% of the data
      slice_sample(prop = 0.1) %>% 
      mutate(abs_diff = abs(lwr.x - upr.y)) 
    
    # split_point not excluding replicate effects
    min_abs_diff_rs = min(comb_data_rs$abs_diff)
    split_point_data_rs = filter(comb_data_rs, abs_diff == min_abs_diff_rs)
    
    # timepoint in original data
    split_points_resample[x] <- split_point_data_rs$time_before_end
  }
  
  split_ews[i,]$split_point_lwr <- quantile(split_points_resample, prob = 0.025)
  split_ews[i,]$split_point_upr <- quantile(split_points_resample, prob = 0.975)
  
}

#_______________________________________________________________________________
#### 10. Plotting EWS ####

## Getting all data together
ews_rawdat <- ews_ttc %>% 
  filter(treatment != "predator") %>% 
  dplyr::select(timeline_component, treatment, replicate, 
                time_date, time_point, time_before_end, 
                cv:acf) %>% 
  pivot_longer(cols  = cv:acf,names_to = "response") 

ews_preddat <- prediction_ews %>%
  filter(treatment != "predator")

## 6a. Behaviour ##
ms_plt_ews_pol <- filter(ews_rawdat, timeline_component == "behaviour") %>% 
  ggplot(aes(x = -time_before_end, y = value, 
             group = interaction(treatment, response), 
             colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(ews_preddat, 
                            timeline_component == "behaviour"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  facet_wrap(~ response, ncol = 3, scales = "free_y") +
  geom_vline(data = filter(split_ews, treatment == "pollution" & 
                             timeline_component == "behaviour"),
             aes(xintercept = -split_point, 
                 group = response, colour = NULL)) +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", "control" = "#1b9e77"),
                     labels = c("Control", "Pollution"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(-80,0, by = 10),
                     labels = seq(80,0, by = -10)) +
  labs(x= NULL, y = "Movement speed early warning signal", tag = "a)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        strip.background = element_blank())

## 6b. Morphology ##
bs_plt_ews_pol <- filter(ews_rawdat, timeline_component == "morphology") %>% 
  ggplot(aes(x = -time_before_end, y = value, 
             group = interaction(treatment, response), 
             colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(ews_preddat, 
                            timeline_component == "morphology"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  facet_wrap(~ response, ncol = 3, scales = "free_y") +
  geom_vline(data = filter(split_ews, treatment == "pollution" & 
                             timeline_component == "morphology"),
             aes(xintercept = -split_point, 
                 group = response, colour = NULL)) +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", "control" = "#1b9e77"),
                     labels = c("Control", "Pollution"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(-80,0, by = 10),
                     labels = seq(80,0, by = -10)) +
  labs(x= NULL, y = "Mean length early warning signal", tag = "b)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        strip.background = element_blank())

## 6c. Abundance ##
ab_plt_ews_pol <- filter(ews_rawdat, timeline_component == "abundance") %>% 
  ggplot(aes(x = -time_before_end, y = value, 
             group = interaction(treatment, response), 
             colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(ews_preddat, 
                            timeline_component == "abundance"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  facet_wrap(~ response, ncol = 3, scales = "free_y") +
  # geom_vline(data = filter(split_ews, treatment == "pollution" & 
  #                            timeline_component == "abundance"),
  #            aes(xintercept = -split_point, 
  #                group = response, colour = NULL)) +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", "control" = "#1b9e77"),
                     labels = c("Control", "Pollution"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(-80,0, by = 10),
                     labels = seq(80,0, by = -10)) +
  labs(x= "Time before collapse (hours)", y = "Abundance early warning signal", tag = "c)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        strip.background = element_blank())

ggsave(ms_plt_ews_pol / bs_plt_ews_pol / ab_plt_ews_pol + 
         plot_layout(guides = "collect") & 
         theme(legend.position = 'top'), 
       filename = "Results/predictions_ews_pollution_all.jpeg",
       width = 23, height = 28, units ="cm", dpi = 600)

ggsave(ab_plt_ews_pol & 
         theme(legend.position = 'top'), 
       filename = "Results/predictions_ews_pollution.jpeg",
       width = 23, height = 10, units ="cm", dpi = 600)

#_______________________________________________________________________________
#### 11. Plotting EWS - predator ####

## Getting all data together
ews_rawdat <- ews_ttc %>% 
  filter(treatment != "pollution") %>% 
  dplyr::select(timeline_component, treatment, replicate, 
                time_date, time_point, time_before_end, 
                cv:acf) %>% 
  pivot_longer(cols  = cv:acf,names_to = "response") 

ews_preddat <- prediction_ews %>%
  filter(treatment != "pollution")

## 6a. Behaviour ##
ms_plt_ews_pred <- filter(ews_rawdat, timeline_component == "behaviour") %>% 
  ggplot(aes(x = -time_before_end, y = value, 
             group = interaction(treatment, response), 
             colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(ews_preddat, 
                            timeline_component == "behaviour"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  facet_wrap(~ response, ncol = 3, scales = "free_y") +
  geom_vline(data = filter(split_ews, treatment == "predator" & 
                             timeline_component == "behaviour"),
             aes(xintercept = -split_point, 
                 group = response, colour = NULL)) +
  scale_color_manual(name = "Treatment",
                     values = c("predator" ="#6654BF", "control" = "#1b9e77"),
                     labels = c("Control", "Predator"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(-80,0, by = 10),
                     labels = seq(80,0, by = -10)) +
  labs(x= NULL, y = "Movement speed early warning signal", tag = "a)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        strip.background = element_blank())

## 6b. Morphology ##
bs_plt_ews_pred <- filter(ews_rawdat, timeline_component == "morphology") %>% 
  ggplot(aes(x = -time_before_end, y = value, 
             group = interaction(treatment, response), 
             colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(ews_preddat, 
                            timeline_component == "morphology"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  facet_wrap(~ response, ncol = 3, scales = "free_y") +
  geom_vline(data = filter(split_ews, treatment == "predator" & 
                             timeline_component == "morphology"),
             aes(xintercept = -split_point, 
                 group = response, colour = NULL)) +
  scale_color_manual(name = "Treatment",
                     values = c("predator" ="#6654BF", "control" = "#1b9e77"),
                     labels = c("Control", "Predator"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(-80,0, by = 10),
                     labels = seq(80,0, by = -10)) +
  labs(x= NULL, y = "Mean length early warning signal", tag = "b)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        strip.background = element_blank())

## 6c. Abundance ##
ab_plt_ews_pred <- filter(ews_rawdat, timeline_component == "abundance") %>% 
  ggplot(aes(x = -time_before_end, y = value, 
             group = interaction(treatment, response), 
             colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = filter(ews_preddat, 
                            timeline_component == "abundance"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  facet_wrap(~ response, ncol = 3, scales = "free_y") +
  # geom_vline(data = filter(split_ews, treatment == "predator" & 
  #                            timeline_component == "abundance"),
  #            aes(xintercept = -split_point, 
  #                group = response, colour = NULL)) +
  scale_color_manual(name = "Treatment",
                     values = c("predator" ="#6654BF", "control" = "#1b9e77"),
                     labels = c("Control", "Predator"),
                     aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(-80,0, by = 10),
                     labels = seq(80,0, by = -10)) +
  labs(x= "Time before collapse (hours)", y = "Abundance early warning signal", tag = "c)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        strip.background = element_blank())

ggsave(ms_plt_ews_pred / bs_plt_ews_pred / ab_plt_ews_pred + plot_layout(guides = "collect") & 
         theme(legend.position = 'top'), 
       filename = "Results/predictions_ews_predator_all.jpeg",
       width = 23, height = 28, units ="cm", dpi = 600)

ggsave(ab_plt_ews_pred & theme(legend.position = 'top'), 
       filename = "Results/predictions_ews_predator.jpeg",
       width = 23, height = 10, units ="cm", dpi = 600)
#_______________________________________________________________________________
#### 12. Manuscript figures ####

stress_start_data <- tibble(treatment = c("control", "pollution", "predator"),
                            start = c(NA, 329.9333, 329.9333))

## 12a. Raw data plot 
# wrangle
ttc_raw <- ttc_caudatum %>%
  mutate(replicate = as.factor(replicate))

ms_plt_raw <- ggplot(ttc_raw, aes(x = time_point, 
                                       y = mean_speed_rep_nocycle,
                                       colour = treatment,
                                       group = id)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  stat_summary(geom = "line", aes(group = treatment), 
               fun = "mean", colour = "black", linewidth = 2) +
  geom_vline(data = stress_start_data, aes(xintercept = start),  
             linetype = "dashed") +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", 
                                "control" = "#1b9e77",
                                "predator" = "#6654BF"),
                     labels = c("Control", "Pollution", "Predator"),
                     aesthetics = c("colour","fill")) +
  facet_grid(~ treatment) +
  labs(x= NULL, y = "Mean Speed\n(mm/s - detrended)", tag = "a)") +
  theme_bw(base_size = 16) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_blank())

bs_plt_raw <- ggplot(ttc_raw, aes(x = time_point, 
                                  y = mean_length_rep_nocycle,
                                  colour = treatment,
                                  group = id)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  stat_summary(geom = "line", aes(group = treatment), 
               fun = "mean", colour = "black", linewidth = 2) +
  geom_vline(data = stress_start_data, aes(xintercept = start),  
             linetype = "dashed") +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", 
                                "control" = "#1b9e77",
                                "predator" = "#6654BF"),
                     labels = c("Control", "Pollution", "Predator"),
                     aesthetics = c("colour","fill")) +
  facet_grid(~ treatment) +
  labs(x= NULL, 
       y = "Mean Body Length\n(\u03BCm - detrended)", tag = "b)") +
  theme_bw(base_size = 16) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_blank())

ab_plt_raw <- ggplot(ttc_raw, aes(x = time_point, 
                                  y = max_abundance_rep_nocycle,
                                  colour = treatment,
                                  group = id)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  stat_summary(geom = "line", aes(group = treatment), 
               fun = "mean", colour = "black", linewidth = 2) +
  geom_vline(data = stress_start_data, aes(xintercept = start),  
             linetype = "dashed") +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", 
                                "control" = "#1b9e77",
                                "predator" = "#6654BF"),
                     labels = c("Control", "Pollution", "Predator"),
                     aesthetics = c("colour","fill")) +
  facet_grid(~ treatment) +
  labs(x= "Experimentation time (hours)", y = "Abundance (detrended)", tag = "c)") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank())

ggsave(ms_plt_raw / bs_plt_raw / ab_plt_raw + plot_layout(guides = "collect") & 
         theme(legend.position = 'top'), 
       filename = "Results/manuscript/figure_1.jpeg",
       width = 26, height = 28, units ="cm", dpi = 600)

## 12b. Pollution timeline plot full
ews_plotdat <- ews_ttc %>% 
  filter(treatment != "predator" &
           timeline_component == "abundance") %>% 
  dplyr::select(timeline_component, treatment, replicate, 
                time_date, time_point, time_before_end, cv)

ews_preddat_plot <- prediction_ews %>%
  filter(treatment != "predator" & response == "cv" & 
           timeline_component == "abundance")

ab_plt_ews_pol_cv <- ggplot(ews_plotdat, aes(x = -time_before_end, y = cv, 
                                          group = treatment, colour = treatment)) +
  stat_summary(fun.data = "mean_se") +
  geom_smooth(stat = "identity", 
              data = ews_preddat_plot,
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  geom_vline(data = filter(split_ews, treatment == "pollution" & 
                             timeline_component == "abundance" &
                             response == "cv"),
             aes(xintercept = -split_point, 
                 group = response, colour = NULL)) +
  scale_color_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", "control" = "#1b9e77"),
                     labels = c("Control", "Pollution"),
                     aesthetics = c("colour","fill"), guide = "none") +
  scale_x_continuous(breaks = seq(-80,0, by = 10),
                     labels = seq(80,0, by = -10)) +
  labs(x= "Time before collapse (hours)", y = "Abundance CV", 
       tag = "d)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        strip.background = element_blank())

ggsave(ms_plt_pol / bs_plt_pol / ab_plt_pol / ab_plt_ews_pol_cv + 
         plot_layout(guides = "collect") & theme(legend.position = 'top'), 
       filename = "Results/timeline_figure.jpeg",
       width = 15, height = 30, units ="cm", dpi = 600)

ggsave(ms_plt_pol / bs_plt_pol / ab_plt_pol / ab_plt_ews_pol_cv + 
         plot_layout(guides = "collect") & theme(legend.position = 'top'), 
       filename = "Results/figure2.pdf",
       width = 15, height = 30, units ="cm")

## 12c. Predator plot
ggsave(ms_plt_pred / bs_plt_pred / ab_plt_pred + 
         plot_layout(guides = "collect") & theme(legend.position = 'top'), 
       filename = "Results/predator_figure.jpeg",
       width = 15, height = 23, units ="cm", dpi = 600)

## 12d. EWS SI plot
ggsave((ab_plt_ews_pol + labs(tag = "a)")) / (ab_plt_ews_pred + labs(tag = "b)")),
  filename = "Results/EWS_figure_SI.jpeg",
  width = 25, height = 23, units ="cm", dpi = 600)

#_______________________________________________________________________________
#### 13. Looking for inflection points in original timeline ####

inflection <- prediction_data_all %>% 
  group_by(timeline_component, treatment, time_point) %>% 
  summarise(fit  = mean(fit), upr = mean(upr), lwr = mean(lwr)) %>%
  ungroup(timeline_component, treatment) %>% 
  group_by() %>% 
  mutate(inflection = c(FALSE, FALSE, diff(diff(fit) > 0) != 0)) %>% 
  ungroup() %>% 
  filter(treatment != "predator", inflection == TRUE)

max_inflection <- inflection %>% 
  filter(treatment != "control") %>% 
  group_by(timeline_component) %>% 
  filter(time_point == max(time_point))

ms_plt <- ggplot(rawdat_pol, aes(x = time_point, y = mean_speed_rep_nocycle, colour = treatment)) +
  geom_smooth(stat = "identity", 
              data = filter(prediction_data,
                            timeline_component == "behaviour"),
              aes(y = fit, ymax = upr, ymin = lwr, 
                  fill = treatment), alpha = 0.3) +
  geom_vline(xintercept = start_stress_time_point, linetype = "dashed") +
  geom_vline(data = filter(max_inflection, timeline_component == "behaviour"),
             aes(xintercept = time_point)) +
  geom_point(data = filter(inflection, timeline_component == "behaviour" & 
                             treatment != "control"), 
             aes(y = fit), colour = "black", size = 4) +
  scale_color_manual(name = "Treatment", values = c("#1b9e77", "#d95f02"),
                     labels = c("Control", "Pollution"), aesthetics = c("colour","fill")) +
  scale_x_continuous(breaks = seq(330, 600, by = 30),
                     limits = c(mintp - 1, 455)) +
  labs(x= NULL, y = "Mean Speed (mm/s)", tag = "a)") +
  theme_bw(base_size = 13) 

ggsave(ms_plt / bs_plt / ab_plt, filename = "Results/predictions_inflection.jpeg",
       width = 15, height = 23, units ="cm", dpi = 600)


