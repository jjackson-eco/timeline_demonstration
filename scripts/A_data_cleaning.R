################################################################################
##                                                                            ##
##           Multivariate signals of population collapse in a                 ##
##                high-throughput ecological experiment                       ##
##                                                                            ##
##          Cerini, F., Jackson, J., Childs, D.Z., & Clements, C.F.           ##
##                                                                            ##
##                    Data cleaning - November 2023                           ##
##                                                                            ##
##                                                                            ##
################################################################################

## This script is used to clean the raw data for Paramecium caudatum and:
##
##  - Formulate the collapse column metric for the timeseries
##  - Do seasonal decomposition by loess for daily fluctuations
##  - Supplementary figures

rm(list = ls())
options(width = 100)

# packages
library(tidyverse)
library(lubridate) 
library(patchwork)
library(ggpubr)

#_______________________________________________________________________________
#### 1. Load and wrangle ####

## Load the data
load("Data//data_complete_rep_caudatum.RData")

## Temporal events
t_init <- as.POSIXct("2022-10-08 06:00:00", tz = "UTC") #reference point: we begin analysing after 5 days of growth and acclimation
treatment_start <- as.POSIXct("2022-10-17 12:00:00", tz = "UTC") #this is the time point when the stressors began
start_stress_time_point<-329.9333 #time point of the onset of stressors in the time series (in hours scale)

#### 1a. Adding collapsed collumn to the data

caudatum_collapsed<-data_complete_rep_caudatum %>% 
  
  ## 2a) Adding in an initial collapsed category for abundance below 10% of reference point
  group_by(id) %>% 
  group_modify(~{
    
    # only max abundance from 8th October as an established reference
    ab_reference = filter(., time_date >= t_init & 
                            time_date < t_init + hours(24))
    max_ab = max(ab_reference$max_abundance_rep)
    
    mutate(., collapsed = if_else(max_abundance_rep/max_ab <= 0.1, "Yes", "No"))}) %>% 
  ungroup() %>% 
  
  ## 2b) Binary collapse column with only data up to collapse - only after treatment to account for odd daily fluctuations
  group_by(id) %>% 
  group_modify(~ {
    # where do collapses happen (after treatment was applied)
    collapse_pos = which(.$time_date > treatment_start & .$collapsed == "Yes")
    # if no collapses, then just have a zero vector the length of the timeseries
    if(length(collapse_pos) == 0){
      col_bin = rep(0, length = nrow(.))
      mutate(., collapse_binary = col_bin)}
    # otherwise, binary collapse variable and only data up to the first post-treatment collapse
    else{
      col_bin = c(rep(0, length = collapse_pos[1] - 1), 1)
      dplyr::slice(., 1:collapse_pos[1]) %>% 
        mutate(collapse_binary = col_bin)}}) %>% 
  ungroup()

#_______________________________________________________________________________
#### 2. STL decomposition ####

caudatum_stl <- caudatum_collapsed %>% 
  
  ## 2b. setting movement speed/body size to 0 if the population has already collapsed
  mutate(mean_speed_rep = if_else(is.nan(mean_speed_rep ) == TRUE, 0, mean_speed_rep),
         mean_length_rep = if_else(is.nan(mean_length_rep) == TRUE, 0, mean_length_rep),
         max_abundance_rep  = if_else(is.na(max_abundance_rep) == TRUE, 0, max_abundance_rep)) %>% 
  mutate(mean_speed_rep = if_else(collapse_binary == 1, 0, mean_speed_rep),
         mean_length_rep = if_else(collapse_binary == 1, 0, mean_length_rep)) %>% 
  filter(time_date >= t_init & 
           time_date < as.POSIXct("2022-10-27 06:00:00") & # reasonable end to the experiment
           collapse_binary == 0) %>% 
  group_by(id) %>% 
  group_modify(~ {
    
    ## 2d. converting data to time-series with eight -observations per day (sampling frequency was every 3 horus)
    ms_ts = ts(.$mean_speed_rep, frequency = 8)
    bs_ts = ts(.$mean_length_rep, frequency = 8)
    ab_ts = ts(.$max_abundance_rep, frequency = 8)
    
    ## 2e. seasonal loess decomposition for each. daily cycle acting as 'season'
    ms_stl = stl(ms_ts, s.window = 9, t.window = 20)
    bs_stl = stl(bs_ts, s.window = 9, t.window = 20)
    ab_stl = stl(ab_ts, s.window = 9, t.window = 20)
    
    ## 2f. add new columns for each of the three additive components (which are columns of the stl object created) of the decomposition
    mutate(., 
           mean_speed_rep_season = ms_stl$time.series[,1],
           mean_speed_rep_trend =  ms_stl$time.series[,2],
           mean_speed_rep_anomaly =  ms_stl$time.series[,3],
           mean_length_rep_season =  bs_stl$time.series[,1],
           mean_length_rep_trend =  bs_stl$time.series[,2],
           mean_length_rep_anomaly =  bs_stl$time.series[,3],
           max_abundance_rep_season =  ab_stl$time.series[,1],
           max_abundance_rep_trend =  ab_stl$time.series[,2],
           max_abundance_rep_anomaly =  ab_stl$time.series[,3])}) %>% 
  
  ## 2g. Add the trend and anomaly component together to remove only the daily cycle component
  mutate(mean_speed_rep_nocycle = mean_speed_rep_trend + mean_speed_rep_anomaly,
         mean_length_rep_nocycle = mean_length_rep_trend + mean_length_rep_anomaly,
         max_abundance_rep_nocycle = max_abundance_rep_trend + max_abundance_rep_anomaly)


#Save dataset with just de-trended data on speed, length and abundance
caudatum_no_cycle_for_analysis<-caudatum_stl %>% select(replicate,time_date, time_point, treatment,mean_speed_rep,
                                                        mean_length_rep, max_abundance_rep, max_abundance_rep_nocycle,
                                                        mean_speed_rep_nocycle, mean_length_rep_nocycle)


save(caudatum_no_cycle_for_analysis, file = "Data/multivariate_caudatum_no_cycle_for_analysis.RData")

#_______________________________________________________________________________
#### 3. Compare plots raw vs de-trended data ####

## 3.1. Plot raw time series CONTROL treatment####

raw_ab_plot_control<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "control",
                                                                         time_date >= t_init), 
                            aes(x = time_point, y = max_abundance_rep, group = replicate))+
  geom_line(col = "#1b9e77", alpha = .7, linewidth = .5)+
  xlab("Time point (hours)")+
  ylab("Abundance")+
  guides(color = "none")+ 
  theme_bw()

raw_moprh_plot_control<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "control",
                                                                            time_date >= t_init), 
                               aes(x = time_point, y = mean_length_rep, group = replicate))+
  geom_line(col = "#1b9e77", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(expression(paste("Mean body length (", mu,"m)")))+
  guides(color = "none")+ 
  theme_bw()

raw_speed_plot_control<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "control",
                                                                            time_date >= t_init), 
                               aes(x = time_point, y = mean_speed_rep , group = replicate))+
  geom_line(col = "#1b9e77", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(("Mean swimming speed (mm/s)"))+
  guides(color = "none")+ 
  theme_bw()

raw_stack_plots_control<-ggarrange(raw_speed_plot_control,raw_moprh_plot_control,raw_ab_plot_control,
                                   nrow = 3, ncol = 1)


## 3.2. Plot detrended time series CONTROL treatment####
det_ab_plot_control<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "control"), 
                            aes(x = time_point, y = max_abundance_rep_nocycle, group = replicate))+
  geom_line(col = "#1b9e77", alpha = .7, linewidth = .5)+
  xlab("Time point (hours)")+
  ylab("Detrended abundance")+
  guides(color = "none")+ 
  theme_bw()

det_moprh_plot_control<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "control"), 
                               aes(x = time_point, y = mean_length_rep_nocycle, group = replicate))+
  geom_line(col = "#1b9e77", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(expression(paste("Detrended mean body length (", mu,"m)")))+
  guides(color = "none")+ 
  theme_bw()

det_speed_plot_control<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "control"), 
                               aes(x = time_point, y = mean_speed_rep_nocycle , group = replicate))+
  geom_line(col = "#1b9e77", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(("Detrended mean swimming speed (mm/s)"))+
  guides(color = "none")+ 
  theme_bw()

det_stack_plots_control<-ggarrange(det_speed_plot_control,det_moprh_plot_control,det_ab_plot_control,
                                   nrow = 3, ncol = 1)

# 3.3. Print Figure S1####

label_text <- c("A","B")

fig.S1<-ggarrange(raw_stack_plots_control, det_stack_plots_control, ncol = 2, labels = c("A","B"))
#save
ggsave(fig.S1, filename = "Results/Fig.S1_raw_vs_detrended_control.pdf",
       width = 45, height = 30, units = "cm", dpi = 800)




## 3.4. Plot raw time series pollution treatment ####
raw_ab_plot_pollution<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "pollution",
                                                                         time_date >= t_init), 
                            aes(x = time_point, y = max_abundance_rep, group = replicate))+
  geom_line(col = "#d95f02", alpha = .7, linewidth = .5)+
  xlab("Time point (hours)")+
  ylab("Abundance")+
  scale_x_continuous(limits = c(100,500))+
  guides(color = "none")+ 
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  theme_bw()

raw_moprh_plot_pollution<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "pollution",
                                                                            time_date >= t_init,
                                                                            time_point <= 450
                                                                            ),
                               aes(x = time_point, y = mean_length_rep, group = replicate))+
  geom_line(col = "#d95f02", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(expression(paste("Mean body length (", mu,"m)")))+
  scale_x_continuous(limits = c(100,500))+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  guides(color = "none")+ 
  theme_bw()

raw_speed_plot_pollution<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "pollution",
                                                                            time_date >= t_init,
                                                                            time_point <= 450
                                                                            ), 
                               aes(x = time_point, y = mean_speed_rep , group = replicate))+
  geom_line(col = "#d95f02", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(("Mean swimming speed (mm/s)"))+
  # scale_x_continuous(limits = c(100,500))+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  guides(color = "none")+ 
  theme_bw()

raw_stack_plots_pollution<-ggarrange(raw_speed_plot_pollution,raw_moprh_plot_pollution,raw_ab_plot_pollution,
                                   nrow = 3, ncol = 1)

## 3.5. Plot detrended time series pollution treatment####
det_ab_plot_pollution<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "pollution"), 
                            aes(x = time_point, y = max_abundance_rep_nocycle, group = replicate))+
  geom_line(col = "#d95f02", alpha = .7, linewidth = .5)+
  xlab("Time point (hours)")+
  scale_x_continuous(limits = c(100,500))+
  ylab("Detrended abundance")+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  guides(color = "none")+ 
  theme_bw()

det_moprh_plot_pollution<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "pollution",
                                                                                  time_point <= 450), 
                               aes(x = time_point, y = mean_length_rep_nocycle, group = replicate))+
  geom_line(col = "#d95f02", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(expression(paste("Detrended mean body length (", mu,"m)")))+
  scale_x_continuous(limits = c(100,500))+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  guides(color = "none")+ 
  theme_bw()

det_speed_plot_pollution<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "pollution",
                                                                                  time_point <= 450), 
                               aes(x = time_point, y = mean_speed_rep_nocycle , group = replicate))+
  geom_line(col = "#d95f02", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  scale_x_continuous(limits = c(100,500))+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  ylab(("Detrended mean swimming speed (mm/s)"))+
  guides(color = "none")+ 
  theme_bw()

det_stack_plots_pollution<-ggarrange(det_speed_plot_pollution,det_moprh_plot_pollution,det_ab_plot_pollution,
                                   nrow = 3, ncol = 1)

# 3.7. Print Figure S2 ####

label_text <- c("A","B")

fig.S2<-ggarrange(raw_stack_plots_pollution, det_stack_plots_pollution, ncol = 2, labels = c("A","B"))
#save
ggsave(fig.S2, filename = "Results/Fig.S2_raw_vs_detrended_pollution.pdf",
       width = 45, height = 30, units = "cm", dpi = 800)



## 3.8. Plot raw time series PREDATOR treatment####

raw_ab_plot_predator<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "predator",
                                                                         time_date >= t_init), 
                            aes(x = time_point, y = max_abundance_rep, group = replicate))+
  geom_line(col = "purple", alpha = .7, linewidth = .5)+
  xlab("Time point (hours)")+
  ylab("Abundance")+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  scale_x_continuous(limits = c(100,500))+
  guides(color = "none")+ 
  theme_bw()

raw_moprh_plot_predator<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "predator",
                                                                            time_date >= t_init), 
                               aes(x = time_point, y = mean_length_rep, group = replicate))+
  geom_line(col = "purple", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(expression(paste("Mean body length (", mu,"m)")))+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  scale_x_continuous(limits = c(100,500))+
  guides(color = "none")+ 
  theme_bw()

raw_speed_plot_predator<-ggplot(data = data_complete_rep_caudatum %>% filter(treatment == "predator",
                                                                            time_date >= t_init), 
                               aes(x = time_point, y = mean_speed_rep , group = replicate))+
  geom_line(col = "purple", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(("Mean swimming speed (mm/s)"))+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  scale_x_continuous(limits = c(100,500))+
  guides(color = "none")+ 
  theme_bw()

raw_stack_plots_predator<-ggarrange(raw_speed_plot_predator,raw_moprh_plot_predator,raw_ab_plot_predator,
                                   nrow = 3, ncol = 1)


## 3.9. Plot detrended time series PREDATOR treatment####
det_ab_plot_predator<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "predator"), 
                            aes(x = time_point, y = max_abundance_rep_nocycle, group = replicate))+
  geom_line(col = "purple", alpha = .7, linewidth = .5)+
  xlab("Time point (hours)")+
  ylab("Detrended abundance")+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  scale_x_continuous(limits = c(100,500))+
  guides(color = "none")+ 
  theme_bw()

det_moprh_plot_predator<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "predator"), 
                               aes(x = time_point, y = mean_length_rep_nocycle, group = replicate))+
  geom_line(col = "purple", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(expression(paste("Detrended mean body length (", mu,"m)")))+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  scale_x_continuous(limits = c(100,500))+
  guides(color = "none")+ 
  theme_bw()

det_speed_plot_predator<-ggplot(data = caudatum_no_cycle_for_analysis %>% filter(treatment == "predator"), 
                               aes(x = time_point, y = mean_speed_rep_nocycle , group = replicate))+
  geom_line(col = "purple", alpha = .7, linewidth = .5)+
  xlab(NULL)+
  ylab(("Detrended mean swimming speed (mm/s)"))+
  geom_vline(xintercept=start_stress_time_point, linetype="dashed", color = "red", linewidth = 0.5)+
  scale_x_continuous(limits = c(100,500))+
  guides(color = "none")+ 
  theme_bw()

det_stack_plots_predator<-ggarrange(det_speed_plot_predator,det_moprh_plot_predator,det_ab_plot_predator,
                                   nrow = 3, ncol = 1)

# 3.10. Print Figure S3####

label_text <- c("A","B")

fig.S3<-ggarrange(raw_stack_plots_predator, det_stack_plots_predator, ncol = 2, labels = c("A","B"))
#save
ggsave(fig.S3, filename = "Results/Fig.S3_raw_vs_detrended_predator.pdf",
       width = 45, height = 30, units = "cm", dpi = 800)
