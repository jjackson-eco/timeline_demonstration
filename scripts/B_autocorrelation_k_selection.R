################################################################################
##                                                                            ##
##           Multivariate signals of population collapse in a                 ##
##                high-throughput ecological experiment                       ##
##                                                                            ##
##    Cerini, F., Jackson, J., O'Brien, D., Childs, D.Z., & Clements, C.F.    ##
##                                                                            ##
##      Exploring autocorrelation and k model selection - February 2025       ##
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
load("data/multivariate_caudatum_no_cycle_for_analysis.RData")

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
#### 2. Histograms ####

speed_hist <- ggplot(ttc_caudatum, aes(x = mean_speed_rep_nocycle, fill = treatment)) +
  geom_histogram() +
  scale_fill_manual(name = "Treatment",
                     values = c("pollution" ="#d95f02", 
                                "control" = "#1b9e77",
                                "predator" = "#6654BF"),
                     labels = c("Control", "Pollution", "Predator")) +
  labs(y = "Frequency", x = "Mean Speed (mm/s)") +
  theme_bw(base_size = 12)

length_hist <- ggplot(ttc_caudatum, aes(x = mean_length_rep_nocycle, fill = treatment)) +
  geom_histogram(bins = 20) +
  scale_fill_manual(name = "Treatment",
                    values = c("pollution" ="#d95f02", 
                               "control" = "#1b9e77",
                               "predator" = "#6654BF"),
                    labels = c("Control", "Pollution", "Predator")) +
  labs(y = "Frequency", x = "Mean Body Length (\u03BCm)") +
  theme_bw(base_size = 12)

abundance_hist <- ggplot(ttc_caudatum, aes(x = max_abundance_rep_nocycle, fill = treatment)) +
  geom_histogram(bins = 20) +
  scale_fill_manual(name = "Treatment",
                    values = c("pollution" ="#d95f02", 
                               "control" = "#1b9e77",
                               "predator" = "#6654BF"),
                    labels = c("Control", "Pollution", "Predator")) +
  labs(y = "Frequency", x = "Abundance") +
  theme_bw(base_size = 12)


ggsave(speed_hist / length_hist / abundance_hist, filename = "output/supplementary/figS4_histograms.jpeg",
       width = 13, height = 15, units = "cm", dpi = 300)

#_______________________________________________________________________________
#### 3. Autocorrelation ####

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

ggsave(ab_auto_plot / ms_auto_plot / bs_auto_plot, filename = "output/supplementary/figS5_temporal_autocorrelation_averages.jpeg",
       width = 14, height = 30, units = "cm", dpi = 800)

#_______________________________________________________________________________
#### 4. Model selection for k ####

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
  labs(x = "Basis dimension (k)", y = "Akaike Information Criterion") +
  facet_wrap( ~ resp, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave(k_mod_plot, filename = "output/supplementary/figS6_k_selection_averages.jpeg",
       width = 14, height = 12, units = "cm", dpi = 800)

## We can use for Abundance k = 17 
## for Body size we go k = 11 
## for speed k = 12

