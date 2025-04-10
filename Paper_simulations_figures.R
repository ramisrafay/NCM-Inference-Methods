# load functions for simulation
source("Function_definitions.R")  ## MAKE SURE THIS IS IN YOUR FILEPATH 


# Make source community abundances  
IWW <- rsad(sad = "ls", frac = 1, coef = list(alpha = 250, N = 1E10), ssize = 1)
saveRDS(IWW, "source_sim_1E10.RDS")

# Once initially created, its faster to load it, uncomment below if needed
# IWW <- readRDS("source_sim_1E10.RDS")
# IWW <- sort(IWW, decreasing = T)

# Scenario 1: Source community abundances are known ####
# Parameters
NTms <- c(5E2, 5E3, 5E4)
n_local <- logspace(log10(2), log10(100), n = 25) #seq(2, 100)
n_source <- NA
sample_sizes <- logspace(log10(1E2), log10(1E5), n = 25)

nsims <- length(NTms) * length(sample_sizes) * length(n_source) * length(n_local)

params_scenario_1 <- tidyr::crossing(NTms, sample_sizes, n_source, n_target)
colnames(params_scenario_1) <- c("sim_ntm", "samplesize", "n_source", "n_local")

# Occupancy inference
occ_results_scenario_1 <- sample_and_fit_NCMs(params = params_scenario_1, true_source = IWW, s = 1, 
                                              method = "occupancy", sample_source = F)

# Variance inference
var_results_scenario_1 <- sample_and_fit_NCMs(params = params_scenario_1, true_source = IWW, s = 1, 
                                              method = "variance", sample_source = F)

# Loglikelihood inference
LL_results_scenario_1 <- sample_and_fit_NCMs(params = params_scenario_1, true_source = IWW, s = 1,
                                             method = "loglikelihood", sample_source = F)


# combine and visualize 
scenario1_all <- bind_rows(occ_results_scenario_1, var_results_scenario_1, LL_results_scenario_1)
scenario1_all$method <- rep(c("occupancy", "variance", "loglikelihood"), each = dim(occ_results_scenario_1)[1])
scenario1_all$method <- factor(scenario1_all$method, levels = c("occupancy", "variance", "loglikelihood"),
                               labels = c("occupancy", "variance", "DM-LL"))

Fig_2 <- ggplot(scenario1_all[scenario1_all$Input_NTm == 500 & scenario1_all$target_samples %in% c(3, 10, 100), ],
                aes(x = sample_size, y = f_NTm, color = method)) +
  geom_point(aes(shape = as.factor(target_samples)), size = 1.5, alpha = 0.9) +
  scale_x_continuous(trans = "log10",
                     breaks = breaks_log(base = 10),
                     labels = label_log(base = 10)) + 
  annotation_logticks(base = 10, outside = T, scaled = T, sides = "bl",
                      short = unit(0.05, "cm"), mid = unit(0.05, "cm"), 
                      long = unit(0.1, "cm")) +
  scale_shape_manual(values = c("3" = 21, "10" = 15, "100" = 25)) +
  geom_hline(yintercept = 500, linetype = "dashed", color = "#000000", size = 0.75) +
  scale_color_manual(values = c("occupancy" = "#C00000", "variance" = "#156082", "DM-LL" = "#CE4ABE"))+
  scale_y_continuous(limits = c(0, 500/4*9), breaks = seq(0, 500/4*9, 500/4)) +
  facet_wrap(~method) +
  labs(x = "read depth", 
       y = expression(inferred~N[T]*'m'),
       shape = "samples taken") + 
  theme_bw() + guides(color = "none")

# Scenario 2: Source community abundances needs to be sampled ####
# Parameters
NTms <- c(5E2, 5E3, 5E4)
n_target <- 10 
n_source <- logspace(log10(2), log10(100), n = 25)
sample_sizes <- logspace(log10(1E2), log10(1E5), n = 25)

nsims <- length(NTms) * length(sample_sizes) * length(n_source) * length(n_target)

params_scenario_2 <- tidyr::crossing(NTms, sample_sizes, n_source, n_target)
colnames(params_scenario_2) <- c("sim_ntm", "samplesize", "n_source", "n_target")

# Occupancy inference
occ_results_scenario_2 <- sample_and_fit_NCMs(params = params_scenario_2, true_source = IWW, s = 1,
                                              method = "occupancy", sample_source = T)
# Variance inference
var_results_scenario_2 <- sample_and_fit_NCMs(params = params_scenario_2, true_source = IWW, s = 1,
                                              method = "variance", sample_source = T)

# Loglikelihood inference
LL_results_scenario_2 <- sample_and_fit_NCMs(params = params_scenario_2, true_source = IWW, s = 1,
                                             method = "loglikelihood", sample_source = T)

# Combine and visualize 
scenario2_all <- bind_rows(occ_results_scenario_2, var_results_scenario_2, LL_results_scenario_2)
scenario2_all$method <- rep(c("occupancy", "variance", "loglikelihood"), each = dim(occ_results_scenario_1)[1])
scenario2_all$method <- factor(scenario1_all$method, levels = c("occupancy", "variance", "loglikelihood"),
                               labels = c("occupancy", "variance", "DM-LL"))

Fig_3 <- ggplot(scenario2_all[scenario2_all$Input_NTm == 500 & scenario2_all$source_samples %in% c(3, 10, 100), ],
       aes(x = sample_size, y = f_NTm, color = method)) +
  geom_point(aes(shape = as.factor(source_samples)), size = 1.5, alpha = 0.9) +
  scale_x_continuous(trans = "log10",
                     breaks = breaks_log(base = 10),
                     labels = label_log(base = 10)) + 
  annotation_logticks(base = 10, outside = T, scaled = T, sides = "b",
                      short = unit(0.05, "cm"), mid = unit(0.05, "cm"), 
                      long = unit(0.1, "cm")) +
  scale_shape_manual(values = c("3" = 21, "10" = 15, "100" = 25)) +
  geom_hline(yintercept = 500, linetype = "dashed", color = "#000000", size = 0.75) +
  scale_color_manual(values = c("occupancy" = "#C00000", "variance" = "#156082", "DM-LL" = "#CE4ABE"))+
  scale_y_continuous(limits = c(0, 500*2), breaks = seq(0, 500*2, 500/4)) + 
  facet_wrap(~method) +
  labs(x = "read depth", 
       y = expression(inferred~N[T]*'m'),
       shape = "samples taken") + 
  theme_bw() + guides(color = "none")

# Scenario 3: Deviations from neutrality ####
# Parameters
NTms <- c(5E2, 5E3, 5E4)
n_source <- c(10, 30, 100)
sample_sizes <- logspace(log10(1E2), log10(1E5), n = 10)

params_scenario_3 <- tidyr::crossing(NTms, sample_sizes, n_source)
params_scenario_3$n_target <- params_scenario_3$n_source
colnames(params_scenario_3) <- c("sim_ntm", "samplesize", "n_source", "n_local")

s3_results_1 <- sample_and_fit_non_neutral_NCMs(params = params_scenario_3, true_source = IWW, s = 1, sample_source = T, 
                                      alpha_deviation = 0.00005, remove_violations = F, reps = 10)

s3_results_2 <- sample_and_fit_non_neutral_NCMs(params = params_scenario_3, true_source = IWW, s = 1, sample_source = T, 
                                      alpha_deviation = 0.0005, remove_violations = F, reps = 10)

s3_results_3 <- sample_and_fit_non_neutral_NCMs(params = params_scenario_3, true_source = IWW, s = 1, sample_source = T, 
                                      alpha_deviation = 0.005, remove_violations = F, reps = 10)

s3_results_4 <- sample_and_fit_non_neutral_NCMs(params = params_scenario_3, true_source = IWW, s = 1, sample_source = T, 
                                      alpha_deviation = 0, remove_violations = F, reps = 10)

s3_results <- rbind(s3_results_1, s3_results_2, 
                    s3_results_3, s3_results_4)


Fig_4 <- vis_scenario_3(s3_results, 10, 500, T) # see function definitions file for 'vis_scenario_3' 

# # Zheng et al dataset inference ####
# ps_zheng <- readRDS("ps_zheng2019.RDS") # this is a phyloseq object 
zheng_source <- subset_samples(ps_zheng, sample_type == "IWW")
zheng_local <- subset_samples(ps_zheng, sample_type == "AS")
Ns <- round(logspace(log10(64000), log10(500), 10))


zheng_results <- depth_profile_ncm(ps_local = zheng_local, ps_source = zheng_source, 
                                       samplesizes = round(logspace(log10(90000), log10(500), 10)), 
                                       reps = 10)

# Visualize 
Fig5 <- ggplot(zheng_results, aes(x = sample_size, y = fit)) + 
  geom_line(aes(color = method), size = 0.75, alpha = 0.7) +
  geom_point(aes(color = method), size = 2, shape = 21) + 
  geom_errorbar(aes(ymax = upper, ymin = lower, color = method), width = 0) +
  scale_color_manual(values = c("occupancy" = "#C00000", "variance" = "#156082", "DM-LL" = "#CE4ABE")) +
  scale_fill_manual(values = c("occupancy" = "#C0000050", "variance" = "#15608250", "DM-LL" = "#CE4ABE50")) +
  scale_x_continuous(trans = "log10",
                     breaks = breaks_log(base = 10),
                     labels = label_log(base = 10),
                     limits = c(5E2, 1E5)) + 
  scale_y_continuous(trans = "log10",
                     breaks = breaks_log(base = 10, n = 6),
                     labels = label_log(base = 10),
                     limits = c(1E1, 1E20)) +
  coord_cartesian(ylim = c(1e1, 1e6)) + 
  annotation_logticks(base = 10, outside = T, scaled = T, sides = "bl",
                      short = unit(0.05, "cm"), mid = unit(0.05, "cm"), 
                      long = unit(0.1, "cm")) +
  labs(y = expression(inferred~N[T]*'m'),
       x = "read depth") +
  guides(fill = "none",
         color = "none") +
  theme(axis.title.x = element_blank()) + 
  theme_bw()



