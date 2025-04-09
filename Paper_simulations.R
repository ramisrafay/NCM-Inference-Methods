# load functions for simulation ## ADD PATH TO FILES IF NECESSARY
source("Function_definitions.R") 


# Make source community abundances  
IWW <- rsad(sad = "ls", frac = 1, coef = list(alpha = 250, N = 1E10), ssize = 1)
saveRDS(IWW, "source_sim_1E10.RDS")

# Once initially created, its faster to load it, uncomment below if needed
# IWW <- readRDS("source_sim_1E10.RDS")
# IWW <- sort(IWW, decreasing = T)

# Scenario 1: Source community abundances are known ####
# Parameters
NTms <- c(5E2, 5E3, 5E4)
n_target <- logspace(log10(2), log10(100), n = 25) #seq(2, 100)
n_source <- NA
sample_sizes <- logspace(log10(1E2), log10(1E5), n = 25)

nsims <- length(NTms) * length(sample_sizes) * length(n_source) * length(n_target)

params_scenario_1 <- tidyr::crossing(NTms, sample_sizes, n_source, n_target)
colnames(params_scenario_1) <- c("sim_ntm", "samplesize", "n_source", "n_target")

# Occupancy inference
occ_results_scenario_1 <- sample_and_fit_NCMs(params = params_scenario_1, true_source = IWW, s = 1, 
                                              method = "occupancy", sample_source = F)
# occ_results_scenario_1$x_coord <- params_scenario_1$n_target
# occ_results_scenario_1$y_coord <- params_scenario_1$samplesize
saveRDS(occ_results_scenario_1, "occ_results_scenario_1.RDS")

# Variance inference
var_results_scenario_1 <- sample_and_fit_NCMs(params = params_scenario_1, true_source = IWW, s = 1, 
                                              method = "variance", sample_source = F)
# var_results_scenario_1$x_coord <- params_scenario_1$n_target
# var_results_scenario_1$y_coord <- params_scenario_1$samplesize
saveRDS(var_results_scenario_1, "var_results_scenario_1.RDS")

# Loglikelihood inference
LL_results_scenario_1 <- sample_and_fit_NCMs(params = params_scenario_1, true_source = IWW, s = 1,
                                             method = "loglikelihood", sample_source = F)
# LL_results_scenario_1$x_coord <- params_scenario_1$n_target
# LL_results_scenario_1$y_coord <- params_scenario_1$samplesize
saveRDS(LL_results_scenario_1, "LL_results_scenario_1.RDS")



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
# occ_results_scenario_2$x_coord <- params_scenario_2$n_source
# occ_results_scenario_2$y_coord <- params_scenario_2$samplesize
saveRDS(occ_results_scenario_2, "occ_results_scenario_2.RDS")

# Variance inference
var_results_scenario_2 <- sample_and_fit_NCMs(params = params_scenario_2, true_source = IWW, s = 1,
                                              method = "variance", sample_source = T)
# var_results_scenario_2$x_coord <- params_scenario_2$n_source
# var_results_scenario_2$y_coord <- params_scenario_2$samplesize
saveRDS(var_results_scenario_2, "var_results_scenario_2.RDS")

# Loglikelihood inference
LL_results_scenario_2 <- sample_and_fit_NCMs(params = params_scenario_2, true_source = IWW, s = 1,
                                             method = "loglikelihood", sample_source = T)
# LL_results_scenario_2$x_coord <- params_scenario_2$n_source
# LL_results_scenario_2$y_coord <- params_scenario_2$samplesize
saveRDS(LL_results_scenario_2, "LL_results_scenario_2.RDS")



# Scenario 3: Deviations from neutrality ####
# Parameters
NTms <- c(5E2, 5E3, 5E4)
n_source <- c(10, 30, 100)
sample_sizes <- logspace(log10(1E2), log10(1E5), n = 10)

params_scenario_3 <- tidyr::crossing(NTms, sample_sizes, n_source)
params_scenario_3$n_target <- params_scenario_3$n_source
colnames(params_scenario_3) <- c("sim_ntm", "samplesize", "n_source", "n_target")

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

saveRDS(s3_results, "All_results_scenario_3.RDS")



# # Zheng et al dataset inference ####
# ps_zheng <- readRDS("ps_zheng2019.RDS")
# zheng_source <- subset_samples(ps_zheng, sample_type == "IWW")
# zheng_target <- subset_samples(ps_zheng, sample_type == "AS")
# Ns <- round(logspace(log10(64000), log10(500), 10))
# 
# 
# zheng_results_more <- depth_profile_ncm(ps_target = zheng_target, ps_source = zheng_source, 
#                                         samplesizes = round(logspace(log10(90000), log10(500), 10)), 
#                                         reps = 10)



