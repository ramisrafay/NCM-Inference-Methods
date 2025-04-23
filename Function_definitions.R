# This file contains all of the functions used in the fitting of the neutral community model (NCM) to abundance data

# required packages 
library(sads)
library(tidyr) # for tidyr::crossing 
library(phyloseq)
library(vegan)
library(gtools)
library(dplyr)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(svglite)
library(patchwork)
library(pracma)
library(scales)

# Functions for inference ####
# Preparing abundance data to be analyzed via NCM # 
NCM_sample_prep <- function(samplesize, ps_local, ps_source, correction = F, 
                            verbose = T, meta = F, cutoff = 1){
  # Inputs # 
  # samplesize: sample depth (# of reads in each sample) used to normalize the abundance tables
  # ps_local: subsetted phyloseq object with abundance table for the local/island community normalized to 'samplesize' 
  # ps_source: subsetted phyloseq object with abundance table for the source/mainland community normalized to 'samplesize' 
  # correction: averaging the abundances for very rare taxa can result in average relative abundances below  
  #             detection limit (1/sample size), this may add noise to NCM inference. Setting to TRUE 
  #             removes those noisy OTUs/ASVs from the analysis. 
  # verbose: if TRUE, prints summary of sample prep methods to console
  # meta: Setting to TRUE, the metacommunity (local community) is considered to be the source community of the system
  
  # Outputs # 
  # a list consisting of the following elements: 
  # $Neutral_df: a table containing the means and variances of relative abundances, and their sample occupancies
  #              for the shared species between input source and local communities
  # $summary: summary of the sample prep process 
  # $source_only: OTUs/ASVs only seen in the source community
  # $local_only: OTUs/ASVs only seen in the local communty
  # $df: data for all species
  # $df_c: data for all species with relative abundance above detection limit d (d = 1/sample size)
  
  # convert abundance to relative abundance
  ra_local <- transform_sample_counts(ps_local, function(x) x/sum(x))
  ra_source <- transform_sample_counts(ps_source, function(x) x/sum(x))
  
  # get abundance tables from phyloseq object
  localtable <- as.matrix(otu_table(ra_local))
  sourcetable <- as.matrix(otu_table(ra_source))
  
  if(meta){
    sourcetable <- localtable
  }
  
  # determine average relative abundances, variances in RA, sample occurrence frequencies for each taxa in each sample
  allOTUs      <- colnames(localtable) 
  
  t_freq <- apply((1 * (localtable > 0)), 2, mean) # occupancy/occurrence frequency in local community
  t_abun <- apply(localtable, 2, mean) # observed mean relative abundance in local community
  t_var  <- apply(localtable, 2, var) # observed variance in local community relative abundances
  
  
  s_freq <- apply((1 * (sourcetable > 0)), 2, mean) # # occupancy/occurrence frequency in source community
  s_abun <- apply(sourcetable, 2, mean) # observed mean relative abundance in source community
  s_var  <- apply(sourcetable, 2, var) # observed variance in source community relative abundances
  
  col_ind <- seq(1:ncol(localtable))
  
  # Make a dataframe containing all the data required for fitting the NCM 
  df <- data.frame(column_index = col_ind,
                   otu_name = allOTUs,
                   local_frequency = t_freq,
                   local_abundance = t_abun,
                   local_variance = t_var,
                   source_frequency = s_freq,
                   source_abundance = s_abun,
                   source_variance = s_var)
  
  shared <- df[which(df$source_frequency > 0 & df$local_frequency > 0), ]
  s_only <- df[which(df$source_frequency > 0 & df$local_frequency == 0), ]
  t_only <- df[which(df$source_frequency == 0 & df$local_frequency > 0), ]
  
  
  if(correction){
    shared_low <- shared[shared$source_abundance < cutoff/samplesize, ]
    shared <- shared[shared$source_abundance >= cutoff/samplesize, ]
  }
  
  if(verbose){
    if(correction){
      print("taxa below detection limit (1/sample size) are not considered for analysis")
      print(paste("# of shared taxa:", dim(shared)[1] + dim(shared_low)[1]))
    } else {
      print("no abundance correction is being used")
      print(paste("# of shared taxa:", dim(shared)[1]))}
    
    print(paste("number of shared taxa above detection limit:", dim(shared)[1]))
    print(paste("Abundance in shared community (%)", round(sum(shared$local_abundance) * 100, 2)))
    print(paste("number of source-only taxa:", dim(s_only)[1]))
    print(paste("number of local-only taxa:", dim(t_only)[1]))
  }
  
  # Summary for shared taxa 
  summ <- setNames(c(ifelse(correction, dim(shared)[1] + dim(shared_low)[1], dim(shared)[1]), 
                     dim(shared)[1], 
                     round(sum(shared$local_abundance) * 100, 2), 
                     dim(s_only)[1], 
                     dim(t_only)[1]),
                   c("all shared taxa", 
                     "shared taxa corrected", 
                     "abundance shared", 
                     "source only taxa", 
                     "local only taxa"))
  
  out <- list(shared, summ, s_only, t_only, df)
  
  names(out) <- c("neutral_df", "summary", "source_only", "local_only", "all_data")
  
  return(out)
}

# neutral dataframe for simulation, does not use phyloseq  
NCM_sample_prep_sim <- function(samplesize, local, source, correction = F, 
                                verbose = T, meta = F, cutoff = 1){
  # Inputs # 
  # samplesize: read depth (# of reads in each sample) used to normalize the abundance tables
  # local: abundance table for the local/island community normalized to 'samplesize'
  # source: abundance table for the source/mainland community normalized to 'samplesize'
  # correction: averaging the abundances for very rare taxa can result in average relative abundances below  
  #             detection limit (1/sample size), this may add noise to NCM inference. Setting to TRUE 
  #             removes those noisy OTUs/ASVs from the analysis. 
  # verbose: if TRUE, prints summary of sample prep to console, steps like rarefaction are also verbose
  # meta: the metacommunity is considered to be the source community of the system
  
  # Outputs # 
  # a list consisting of the following elements: 
  # $Neutral_df: a table containing the means and variances of relative abundances, and their occurrence
  #                   frequences for the shared species between input source and local communities
  # $summary: summary of the sample prep process 
  # $source_only: OTUs/ASVs only seen in the source community
  # $local_only: OTUs/ASVs only seen in the local communty
  # $df: data for all species
  # $df_c: data for all species with relative abundance above detection limit d (d = 1/sample size)
  
  # convert data to relative abundance
  ra_local <- t(apply(local, 1, function(x) x/sum(x))) # 
  ra_source <- t(apply(source, 1, function(x) x/sum(x)))
  
  # get abundance tables from phyloseq object
  localtable <- as.matrix(ra_local)
  sourcetable <- as.matrix(ra_source)
  
  if(meta){
    sourcetable <- localtable 
  }
  
  # determine average relative abundances, variances in RA, sample occurrence frequencies for each taxa in each sample
  allOTUs      <- colnames(localtable) 
  
  t_freq <- apply((1 * (localtable > 0)), 2, mean) # occupancy/occurrence frequency in local community
  t_abun <- apply(localtable, 2, mean) # mean relative abundance in local
  t_var  <- apply(localtable, 2, var) # observed variance in local community relative abundances
  
  
  s_freq <- apply((1 * (sourcetable > 0)), 2, mean) # occupancy/occurrence frequency in source community
  s_abun <- apply(sourcetable, 2, mean) # mean relative abundance in source
  s_var  <- apply(sourcetable, 2, var) # observed variance in source community relative abundances
  
  col_ind <- seq(1:ncol(localtable))
  
  # Make a dataframe out of the whole thing 
  df <- data.frame(column_index = col_ind,
                   otu_name = allOTUs,
                   local_frequency = t_freq,
                   local_abundance = t_abun,
                   local_variance = t_var,
                   source_frequency = s_freq,
                   source_abundance = s_abun,
                   source_variance = s_var,
                   alpha = t_abun - s_abun)
  
  shared <- df[which(df$source_frequency > 0 & df$local_frequency > 0), ]
  s_only <- df[which(df$source_frequency > 0 & df$local_frequency == 0), ]
  t_only <- df[which(df$source_frequency == 0 & df$local_frequency > 0), ]
  
  
  if(correction){
    shared_low <- shared[shared$source_abundance < cutoff/samplesize, ]
    shared <- shared[shared$source_abundance >= cutoff/samplesize, ]
  }
  
  
  if(verbose){
    if(correction){
      print("taxa below detection limit (1/sample size) are not considered for analysis")
      print(paste("# of shared taxa:", dim(shared)[1] + dim(shared_low)[1]))
    } else {
      print("no abundance correction is being used")
      print(paste("# of shared taxa:", dim(shared)[1]))}
    
    print(paste("number of shared taxa above detection limit:", dim(shared)[1]))
    print(paste("Abundance in shared community (%)", round(sum(shared$local_abundance) * 100, 2)))
    print(paste("number of source-only taxa:", dim(s_only)[1]))
    print(paste("number of local-only taxa:", dim(t_only)[1]))
  }
  
  summ <- setNames(c(ifelse(correction, dim(shared)[1] + dim(shared_low)[1], dim(shared)[1]), 
                     dim(shared)[1], 
                     round(sum(shared$local_abundance) * 100, 2), 
                     dim(s_only)[1], 
                     dim(t_only)[1]),
                   c("all shared taxa", 
                     "shared taxa corrected", 
                     "abundance shared", 
                     "source only taxa", 
                     "local only taxa"))
  
  out <- list(shared, summ, s_only, t_only, df)
  
  names(out) <- c("neutral_df", "summary", "source_only", "local_only", "all_data")
  
  return(out)
}

# Functions used to fit NCM # 
# Occupancy-based fitting method
fit_occ <- function(samplesize, neutral_data, remove_ones = T){
  # Inputs #
  # samplesize: otherwise called Ns in this code, read depth (# of reads in each sample) used to normalize the abundance tables
  # neutral_data: the list output from NCM_sample_prep()
  # remove_ones: if TRUE, does not include taxa that have occupancies of 1 for NTm inference

  # Outputs #
  # A list with the following elements: 
  # $fitting_results: fitted coefficients with confidence intervals
  # $ndf_w_preds: neutral data frame with fitted, upper and lower limit predictions for expected sample occurrence frequencies 
  # $m.fit: model fit object
  
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  require(boot)
  
  ndf <- neutral_data$neutral_df  # Data frame containing data for NCM fitting 
  Ns <- samplesize                # Sample depth 
  d <- 1/Ns                       # detection limit 
  freq <- ndf$local_frequency     # local community sample occupancies/occurrence frequencies
  p <- ndf$source_abundance       # source community abundances
  
  print(paste("fitting params: sample size =", Ns, ", d =", d))
  
  # Fit NCM using the occupancy-based method
  m.fit <- nlsLM(freq ~ pbeta(d, shape1 = Nsm*p, shape2 = Nsm*(1-p), lower.tail = F), 
                 lower = 0, start = list(Nsm = 1))
  
  reg <- coef(m.fit)
  
  
  conf <- tryCatch({confint(m.fit)}, 
                   error = function(e) {c(coef(m.fit), coef(m.fit))})
  
  # Predicted sample occupancies based on fitted NTm 
  pred <- pbeta(d, shape1 = coef(m.fit)*p, shape2 = coef(m.fit)*(1-p), lower.tail = F) 
  pred_LL <- pbeta(d, shape1 = conf[1]*p, shape2 = conf[1]*(1-p), lower.tail = F)
  pred_UL <- pbeta(d, shape1 = conf[2]*p, shape2 = conf[2]*(1-p), lower.tail = F)
  ndf$predicted_freq <- pred
  ndf$pred_LL <- pred_LL
  ndf$pred_UL <- pred_UL 
  
 
  Rsq <- 1 - (sum((ndf$predicted_freq - ndf$local_frequency)^2) / 
                sum((mean(ndf$local_frequency) - ndf$local_frequency)^2))   # Model fit
  
  shared_abund <- neutral_data$summary["abundance shared"]
  
  results <- setNames(c(reg, conf, Rsq, Ns, d, shared_abund),
                      c("fit", "LL", "UL", "Rsq", "sample size", "detection limit", "shared abundance"))
  
  out <- list(results, ndf, m.fit) 
  names(out) <- c("fitting_results", "ndf_w_preds", "model_fit")
  
  return(out)
}

# Variance-based fitting method (accounts for multinomial sampling variance)
fit_var <- function(samplesize, neutral_data){
  # Inputs #
  # samplesize: read depth (# of reads in each sample) used to normalize the abundance tables, shorthand as 'Ns' in this function
  # neutral_data: the list output from NCM_sample_prep()
  
  # Outputs #
  # A list with the following elements: 
  # $fitting_results: fitted coefficients with confidence intervals
  # $ndf_w_preds: neutral dataframe with fitted, upper and lower limit predictions for expected variance in relative abundance 
  # $m.fit: model fit object
  
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  require(boot)
  
  ndf <- neutral_data$neutral_df  # Data frame containing data for NCM fitting 
  Ns <- samplesize                # Sample depth 
  d <- 1/Ns                       # detection limit 
  freq <- ndf$local_frequency     # local community sample occupancies/occurrence frequencies
  p <- ndf$source_abundance       # source community abundances
  
  print(paste("fitting params: sample size =", Ns, ", d =", d))
  
  # Fit NCM using the variance-based method
  m.fit <- nlsLM(var ~ (p*(1-p)/((Nsm + 1)) + (p*(1-p)/(Ns))),
                 lower = 0, start = list(Nsm = 1))
  
  reg <- coef(m.fit)
  conf <- confint(m.fit)
  
  # Predicted relative abundance variances based on fitted NTm 
  pred <- (p*(1-p)/((reg) + 1) + (p*(1-p)/(Ns)))
  pred_LL <- (p*(1-p)/((conf[1]) + 1) + (p*(1-p)/(Ns)))
  pred_UL <- (p*(1-p)/((conf[2]) + 1) + (p*(1-p)/(Ns)))
  ndf$predicted_var <- pred
  ndf$predvar_LL <- pred_LL
  ndf$predvar_UL <- pred_UL 
  
  
  Rsq <- 1 - (sum((ndf$predicted_var - ndf$local_variance)^2) /
                sum((mean(ndf$local_variance) - ndf$local_variance)^2))   # Model fit
  
  shared_abund <- neutral_data$summary["abundance shared"]
  
  results <- setNames(c(reg, conf, Rsq, Ns, d, shared_abund),
                      c("fit", "LL", "UL", "Rsq", "sample size", "detection limit", "shared abundance"))
  
  out <- list(results, ndf, m.fit) 
  names(out) <- c("fitting_results", "ndf_w_preds", "model_fit")
  
  return(out)
}




# Dirichlet Multinomial Loglikelihood fitting
fit_DMLL <- function(neutral_data, local, samplesize, plot = F, maxNTm = 1E12){
  # Inputs #
  # neutral_data: the list output from NCM_sample_prep()
  # local: abundance table (not relative abundance) for the local community rarefied to samplesize
  # samplesize: read depth (# of reads in each sample) used to normalize the abundance tables, otherwise called Ns in this function
  # plot: If TRUE, plots the likelihood function against NTm with the maximum likelihood estimate
  # maxNTm: Maximum value for NTm in search space, reduce if computation takes too long and if there is reason to suspect upper limit is much lower
  
  # Outputs #
  # A list with the following elements: 
  # $fitting_results: fitted coefficients with confidence intervals
  # $ndf_w_preds: neutral data frame with fitted, upper and lower limit predictions for expected sample occurrence frequencies 
  # $m.fit: model fit object
  
  # Compute log-likelihood for a single local community sample using the Dirichlet-multinomial model 
  loglik_ntm <- function(Ntm, neutral_data, local, samplesize) {
    
    Ns <- samplesize
    ndf <- neutral_data$neutral_df
    
    # Precompute constants
    constant <- log(Ns) + lgamma(Ntm) + lgamma(Ns) - lgamma(Ntm + Ns)
    
    # Function to compute log-likelihood for a single row in abundance table 
    compute_row_LL <- function(tar_row) {
      xi <- tar_row[ndf$column_index]
      pi <- ndf$source_abundance
      
      # Filter out zero entries in xi
      non_zero_mask <- xi > 0 # log-gamma of 0 is undefined, so taxon cannot be included if it has abundance of 0
      xi <- xi[non_zero_mask]
      pi <- pi[non_zero_mask]
      
      if (length(xi) == 0) return(0)  # Skip if no non-zero entries
      
      # Compute LL summation for this row
      summand <- constant - sum(log(xi) + lgamma(Ntm * pi) + lgamma(xi) - lgamma((Ntm * pi) + xi))
      
      return(summand)
    }
    
    # row-wise summation of log-likelihoods 
    LL <- sum(apply(local, 1, compute_row_LL))
    return(LL)
  }
  
  # Fitting via maximum likelihood 
  # search over LL space to determine Ntm that maximizes LL
  fit_DMLL <- optimize(loglik_ntm, interval = c(1, maxNTm),
                       neutral_data = neutral_data, local = local, samplesize = samplesize,
                       upper = maxNTm, maximum = T) 
  
  # plot likelihood function 
  if(plot){
    Ntms <- 10 ^ seq(1, log10(maxNTm), length.out = 100)
    LLs <- sapply(Ntms, loglik_ntm, local = local, neutral_data = neutral_data, N = N)
    plot(x = Ntms, y = LLs, log = 'x')
    abline(v = fit_DMLL$maximum, col = "purple3")
  }
  
  return(fit_DMLL)
}


depth_profile_ncm <- function(ps_local, ps_source, samplesizes, reps){
  # Inputs #
  # samplesizes: vector of read depths (# of reads in each sample) used to normalize the abundance tables
  # ps_local: subsetted phyloseq object with abundance table for the local/island community normalized to 'samplesize' 
  # ps_source: subsetted phyloseq object with abundance table for the source/mainland community normalized to 'samplesize' 
  # plot: If TRUE, plots the likelihood function against NTm with the maximum likelihood estimate
  # maxNTm: Maximum value for NTm in search space, reduce if computation takes too long and if there is reason to suspect upper limit is much lower
  # reps: number of repetitions for each read depth in 'samplesizes'
  
  # Outputs #
  # A data frame with the results of NTm inference sample depth profiling with estimates from the occupancy, 
  # variance and DM-LL methods at each sample depth provided in vector 'samplesizes' 
  results <- data.frame(fit = numeric(0), upper = numeric(0), lower = numeric(0),
                        shared_species = numeric(0), sample_size = numeric(0), source_samples = numeric(0), 
                        local_samples = numeric(0), method = character(0))
  
  for(R in 1:length(samplesizes)){
    for(i in 1:reps){
      local_rare <- rarefy_even_depth(ps_local, sample.size = samplesizes[R], trimOTUs = F)
      source_rare <- rarefy_even_depth(ps_source, sample.size = samplesizes[R], trimOTUs = F)
      
      neutral_data <- NCM_sample_prep(samplesize = samplesizes[R], local = local_rare, source = source_rare, correction = T,
                                      verbose = F, meta = F)
      ncm_occ <- fit_occ(samplesize = samplesizes[R], neutral_data = neutral_data, m_only = F)
      results <- rbind(results, list(fit = ncm_occ$fitting_results[1], upper = ncm_occ$fitting_results[3], 
                                     lower = ncm_occ$fitting_results[2], shared_species = neutral_data$summary[1], 
                                     sample_size = samplesizes[R], source_samples = nsamples(source_rare),
                                     local_samples = nsamples(local_rare), method = "occupancy"))
      ncm_var <- fit_var(samplesize = samplesizes[R], neutral_data = neutral_data, m_only = F)
      results <- rbind(results, list(fit = ncm_var$fitting_results[1], upper = ncm_var$fitting_results[3], 
                                     lower = ncm_var$fitting_results[2], shared_species = neutral_data$summary[1], 
                                     sample_size = samplesizes[R], source_samples = nsamples(source_rare),
                                     local_samples = nsamples(local_rare), method = "variance"))
      ncm_LL  <- fit_DMLL(neutral_data = neutral_data, local = as.matrix(otu_table(local_rare)), N = samplesizes[R], plot = F)
      results <- rbind(results, list(fit = ncm_LL$maximum, upper = NA, 
                                     lower = NA, shared_species = neutral_data$summary[1], 
                                     sample_size = samplesizes[R], source_samples = nsamples(source_rare), 
                                     local_samples = nsamples(local_rare), method = "DM-LL"))
    }
  }
  return(results)
}


# # Scenario 3 plotting function #
vis_scenario_3 <- function(s3res, local_samples, NTm_sim, legend = T){
  # Inputs: 
  # s3res: scenario 3 dataframe 
  # local_samples: visualize for this number of local community samples
  # NTm_sim: Input ground-truth NTm value for simulation 
  # legend: if TRUE, include legend for inference methods 
  
   s3 <- s3res[s3res$local_samples == local_samples & s3res$Input_NTm == NTm_sim, ]
   
   s <- data.frame(fitted_ntm = c(s3$NTm_occ, s3$NTm_var, s3$NTm_LL),
                   fitted_ntm_sd = c(s3$NTm_occ_sd, s3$NTm_var_sd, s3$NTm_LL_sd),
                   sample_size = rep(s3$sample_size, 3),
                   alpha_deviation = rep(s3$alpha_deviation, 3),
                   method = c(rep("occupancy", length(s3$NTm_occ)), 
                              rep("variance", length(s3$NTm_var)), 
                              rep("DM-LL", length(s3$NTm_LL)))
   )
   
   s$method <- factor(s$method, levels = c("occupancy", "variance", "DM-LL"))
   s$alpha_deviation <- format(s$alpha_deviation, scientific = T)
   s$category <- factor(paste0("σ = ", s$alpha_deviation),
                        levels = c("σ = 0e+00", "σ = 5e-05", "σ = 5e-04", "σ = 5e-03"),
                        labels = c("σ = 0", "σ = 5E-05", "σ = 5E-04", "σ = 5E-03"))
   
   p <- ggplot(s, aes(x = sample_size, y = fitted_ntm)) + 
     geom_line(aes(color = method), size = 0.75) + 
     geom_point(aes(color = method), size = 1) + 
     geom_errorbar(aes(ymax = fitted_ntm + fitted_ntm_sd, ymin = fitted_ntm - fitted_ntm_sd, color = method), width = 0) +
     scale_color_manual(values = c("occupancy" = "#C00000", "variance" = "#156082", "DM-LL" = "#CE4ABE")) + 
     scale_x_continuous(trans = "log10",
                        breaks = breaks_log(base = 10),
                        labels = label_log(base = 10)) + 
     scale_y_continuous(trans = "log10",
                        breaks = breaks_log(base = 10, n = 6),
                        labels = label_log(base = 10),
                        limits = c(1E1, 1E13)) +
     coord_cartesian(ylim = c(1e1, 1e6)) + 
     annotation_logticks(base = 10, outside = T, scaled = T, sides = "bl",
                         short = unit(0.05, "cm"), mid = unit(0.05, "cm"), 
                         long = unit(0.1, "cm")) +
     facet_wrap(~category, scales = 'free_y', nrow = 1) + 
     geom_hline(yintercept = NTm_sim, linetype = "dashed", color = "#000000", size = 0.75) +
     
     labs(x = "sample depth (reads per sample)", 
          y = expression(inferred~N[T]*'m'),
          color = "Inference method") +
     guides(fill = "none") + 
     theme_bw()
   
   if(!legend){
     p <- p + theme(legend.position = "none")
   }
   
   return(p)
 }


# Define function for running simulations over parameter space ####
# Permute over parameter space, simulate and fit NCM #
sample_and_fit_NCMs <- function(params, true_source, s = 1, method = "occupancy", 
                                sample_source = T, reps = 10){
  # Inputs: 
  # params: parameter set as data frame (see above), 
  # true_source: 'true' source community abundances
  # s = start position (for debugging/partial runs)
  # method (ie what method to use for fitting?): 
  # "occupancy" - sample occupancy-based fitting;
  # "variance" - beta dist variance-based fitting;
  # "loglikelihood" - likelihood based fitting
  # sample_source: T if sampling the source community, F if using the true source community for inference
  # reps: number of simulations to run for each line in 'params'
  
  # Outputs: 
  # data frame of fitted NTm values and parameters used for simulation case.
  
  
  set.seed(NULL)
  require(sads) # loaded for sads::rsad
  require(LaplacesDemon) # loaded for LaplacesDemon::rdirichlet
  
  # source community
  iww_rel <- true_source/sum(true_source)
  
  # Create results data frame
  nsims <- dim(params)[1]
  
  #nsims = 100 # override for debugging/partial runs only 
  
  if(method == "loglikelihood"){
    outputs <- c("f_NTm", "f_NTm_sd", 
                 "sample_size", "Input_NTm", "local_samples", "source_samples", "shared species")
  } else {
    outputs <- c("f_NTm", "f_NTm_LL", "f_NTm_UL", "Rsq", "f_NTm_sd", "f_NTm_LL_sd", "f_NTm_UL_sd", "Rsq_sd", 
                 "sample_size", "Input_NTm", "local_samples", "source_samples", "shared species")
  }
  
  sim_ncm_results <- data.frame(matrix(nrow = nsims, ncol = length(outputs)))
  rownames(sim_ncm_results) <- seq(1, nsims)
  colnames(sim_ncm_results) <- outputs 
  
  
  # main routine
  for(i in s:nsims){
    print(paste("simulation progress (#): ", i, "of", nsims))
    
    ## 1. Select values for params
    Ns            <- round(params$samplesize[i])
    if(sample_source) source_samps <- round(params$n_source[i]) else source_samps <- NA
    local_samps  <- round(params$n_local[i])
    #m_true_sim  <- calc_true_m(m_hat_sim, NT, Ns) # in case m is modelled instead of Ntm
    NTm_sim       <- params$sim_ntm[i]
    
    for(rep in 1:reps){
      # 2. Neutral assembly for true local community
      assemblies <- if(sample_source) max(params$n_source) else 100 # source community assemblies to make
      AS_sim  <- rdirichlet(n = assemblies, alpha = NTm_sim*iww_rel) # generate assembly
      
      ## 3. Draw samples from IWW and AS communities
      # Source
      if(sample_source){
        iww_sample <- rmultinom(n = source_samps, size = Ns, prob = iww_rel)
        
        iww_species <- 1:nrow(iww_sample)
        iww_otu     <- data.frame(matrix(nrow = source_samps, ncol = length(iww_species)))
        
        if(source_samps == 1){
          iww_otu[1, ] <- iww_sample
        } else {
          for(k in 1:source_samps){
            iww_otu[k, ] <- iww_sample[, k]
          }
        }
      } else {
        iww_otu <- data.frame(matrix(nrow = 1, ncol = length(iww_rel)))
        iww_otu[1, ] <- iww_rel
      }
      
      # local
      inds <- sample(1:assemblies, local_samps, replace = F) # randomly pick from instances of assemblies
      as_otu <- data.frame(matrix(nrow = local_samps, ncol = ncol(AS_sim)))
      
      for(j in 1:local_samps){
        ind <- inds[j] # pick assembly instance to sample from
        AS_comm <- AS_sim[ind, ] # pick one of the assemblies
        as_samp <- rmultinom(n = 1, size = Ns, prob = AS_comm) # draw abundance data
        as_otu[j, ] <- as_samp
      }
      

      ## 5. Create NCM input data frame
      neutral_data <- NCM_sample_prep_sim(samplesize = Ns, local = as_otu, source = iww_otu, 
                                          correction = F, verbose = F, meta = F)
      
      ## 6. Fit NCM 
      if(method == "occupancy"){
        ncm_res <- fit_occ(samplesize = Ns, neutral_data = neutral_data)
        fits    <- ncm_res$fitting_results[1:4]
       
      } else if(method == "variance"){
        ncm_res <- fit_var(samplesize = Ns, neutral_data = neutral_data)
        fits    <- ncm_res$fitting_results[1:4]


      } else if(method == "loglikelihood"){
        ncm_res <- fit_DMLL(neutral_data = neutral_data, local = as_otu, 
                              N = Ns, plot = F)
        fits <- ncm_res$maximum 
        }
      
      if(rep == 1){
        inferences <- data.frame(matrix(ncol = length(fits), nrow = 10))
        inferences[1, ] <- fits
      } else {
        inferences[rep, ] <- fits
      }
    }
    
    ## 7. compile results
    inference_means <- apply(inferences, MARGIN = 2, mean)
    inference_sds   <- apply(inferences, MARGIN = 2, sd)
    
    if(method == "loglikelihood"){
      res <- c(inference_means, inference_sds, Ns, NTm_sim, local_samps,
               source_samps, dim(neutral_data$neutral_df)[1]) |> unlist()
    } else {
      res <- c(unname(inference_means), unname(inference_sds), Ns, NTm_sim, local_samps,
               source_samps, dim(neutral_data$neutral_df)[1]) |> unlist()
    }
    
    ## 8. Save results 
    sim_ncm_results[i, ] <- res
    
    # uncomment below to save interim progress (currently saving every 100 simulations) # 
    # if(i%%100 == 0) saveRDS(sim_ncm_results, 
    #                         file = paste0(i, "_of_", nsims, "_sim_results.RDS"))
  }
  
  return(sim_ncm_results)
}

# Generate advantage terms for each species on the island 
make_alphas <- function(size, deviation, seed){
  set.seed(seed)
  
  # Generates noisy alpha parameters based on the normal distribution
  if(size%%2 == 0){
    a <- rnorm(size/2, 0, deviation) # this is adding random noise proportional to selection strength 
    a_sym <- c(a, -a) # makes them sum to 0
  } else {
    a <- rnorm((size-1)/2, 0, deviation)
    a_sym <- c(a, -a, 0, 0) # makes them sum to 0
  }
  
  alphas <- sample(a_sym, size, replace = T)
  
  set.seed(NULL) 
  
  return(alphas)
} 

apply_advantages <- function(neutral_assembly, alphas){
  non_neutral_assembly <- sweep(neutral_assembly, 2, alphas, "+") # add alpha terms to abundances 
  non_neutral_assembly <- apply(non_neutral_assembly, 2, function(x) ifelse(x < 0, 0, x)) # remove taxa with <0 abundance
  non_neutral_assembly <- apply(non_neutral_assembly, 1, function(x) x/sum(x)) # re-scale relative abundances
  
  return(t(non_neutral_assembly))
}

sample_and_fit_non_neutral_NCMs <- function(params, true_source, s = 1, sample_source = T, 
                                            alpha_deviation = 0, reps){
  # Inputs: 
  # params: parameter set as data frame (see above), 
  # true_source: 'true' source community abundances
  # s = start position (for debugging/partial runs)
  # sample_source: T if sampling the source community, F if using the true source community for inference
  # reps: number of simulations to run for each line in 'params'
  
  # Outputs: 
  # data frame of fitted NTm values and parameters used for simulation case.
  
  set.seed(NULL)
  require(sads) # loaded for sads::rsad
  require(LaplacesDemon) # loaded for LaplacesDemon::rdirichlet
  
  # source community
  iww_rel <- true_source/sum(true_source)
  
  # Create results data frame
  nsims <- dim(params)[1]
  
  #nsims = 100 # override for debugging/partial runs only 
  
  sim_ncm_results <- data.frame(matrix(nrow = nsims, ncol = 11))
  colnames(sim_ncm_results) <- c("NTm_occ", "NTm_occ_sd", "NTm_var", "NTm_var_sd",
                                 "NTm_LL", "NTm_LL_sd","sample_size", "Input_NTm", 
                                 "local_samples", "source_samples", "alpha_deviation")
  
  seeds <- 1:(nsims*reps)
  
  # main routine
  for(i in s:nsims){
    print(paste("simulation progress (#): ", i, "of", nsims))
    
    ## 1. Select values for params
    Ns            <- round(params$samplesize[i])
    if(sample_source) source_samps <- round(params$n_source[i]) else source_samps <- NA
    local_samps  <- round(params$n_local[i])
    NTm_sim       <- params$sim_ntm[i]
    
    for(rep in 1:reps){
      ## 2. Neutral assembly for true local community
      assemblies <-  100 # source community assemblies to make
      AS_sim  <- rdirichlet(n = assemblies, alpha = NTm_sim*iww_rel) # generate neutral assembly
      
      # create and apply advantage parameters
      alphas <- numeric(length(iww_rel))
      mm <- apply(AS_sim, 2, mean) # look at mean relative abundances on the island
      alphas[which(mm >= 1/1E10)] <- make_alphas(size = length(mm[mm >= 1/1E10]), 
                                                 deviation = alpha_deviation, 
                                                 seed = seeds[(i-1)*reps + rep])  # create advantage terms for species have at least abundance 1 given 1E10 individuals on the mainland 

      AS_sim  <- apply_advantages(AS_sim, alphas) # apply advantage terms to those species 
      
      ## 3. Draw samples from IWW and AS communities
      # Source
      if(sample_source){
        iww_sample <- rmultinom(n = source_samps, size = Ns, prob = iww_rel)
        
        iww_species <- 1:nrow(iww_sample)
        iww_otu     <- data.frame(matrix(nrow = source_samps, ncol = length(iww_species)))
        
        if(source_samps == 1){
          iww_otu[1, ] <- iww_sample
        } else {
          for(k in 1:source_samps){
            iww_otu[k, ] <- iww_sample[, k]
          }
        }
      } else {
        iww_otu <- data.frame(matrix(nrow = 1, ncol = length(iww_rel)))
        iww_otu[1, ] <- iww_rel
      }
      
      # local
      inds <- sample(1:assemblies, local_samps, replace = F) # randomly pick from instances of assemblies
      as_otu <- data.frame(matrix(nrow = local_samps, ncol = ncol(AS_sim)))
      
      for(j in 1:local_samps){
        ind <- inds[j] # pick assembly instance to sample from
        AS_comm <- AS_sim[ind, ] # pick one of the assemblies
        
        as_samp <- rmultinom(n = 1, size = Ns, prob = AS_comm) # draw abundance data
        as_otu[j, ] <- as_samp
      }
      
      ## 5. Create NCM input data frame
      neutral_data <- NCM_sample_prep_sim(samplesize = Ns, local = as_otu, source = iww_otu, 
                                          correction = T, verbose = F, meta = F)

      ## 6. Fit NCM 
      # Occupancy
      ncm_occ <- fit_occ(samplesize = Ns, neutral_data = neutral_data, m_only = F)
      fitocc   <- ncm_occ$fitting_results[1]
        
      # Variance
      ncm_var <- fit_var(samplesize = Ns, neutral_data = neutral_data, m_only = F)
      fitvar <- ncm_var$fitting_results[1]
        
      # Multinomial Dirichlet
      ncm_LL <- fit_DMLL(neutral_data = neutral_data, local = as_otu, N = Ns, plot = F)
      fitLL <- ncm_LL$maximum 
        
      if(rep == 1) i_occ <- i_var <- i_LL <- numeric(10)
        
      i_occ[rep] <- fitocc
      i_var[rep] <- fitvar
      i_LL[rep] <- fitLL
    }
    
    ## 7. compile results
    res <- c(mean(i_occ), sd(i_occ), mean(i_var), sd(i_var), mean(i_LL), sd(i_LL), 
             Ns, NTm_sim, local_samps, source_samps, alpha_deviation)
    
    
    ## 8. Save results 
    sim_ncm_results[i, ] <- res
    
  }
  
  return(sim_ncm_results)
}

