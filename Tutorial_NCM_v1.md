NCM inference tutorial V1
================
Ramis Rafay
2026-05-06

- [A guide to inferring immigration rates from abundance
  datasets](#a-guide-to-inferring-immigration-rates-from-abundance-datasets)
  - [The quick way: read depth profiling and inference result
    visualization](#the-quick-way-read-depth-profiling-and-inference-result-visualization)
  - [Manually fitting each method and inspecting plots of each
    fit](#manually-fitting-each-method-and-inspecting-plots-of-each-fit)
    - [Preparing abundance datasets for NCM fitting and
      inference](#preparing-abundance-datasets-for-ncm-fitting-and-inference)
    - [Occupancy-based inference
      method](#occupancy-based-inference-method)
    - [Variance-based inference
      method](#variance-based-inference-method)
    - [Dirichlet multinomial log-likelihood
      method](#dirichlet-multinomial-log-likelihood-method)

# A guide to inferring immigration rates from abundance datasets

NOTE: The code is generally designed to work with abundance datasets
packaged as Phyloseq objects.

## The quick way: read depth profiling and inference result visualization

The `depth_profile_ncm` function provides a convenient way to run all
three NCM inference methods across a range of rarefied read depths. This
is useful when you want to examine how sensitive the estimated
immigration parameter is to sequencing depth.

Rather than manually rarefying the data, preparing the neutral model
input, fitting each method, and repeating the process for multiple read
depths, `depth_profile_ncm` performs this workflow automatically. For
each specified read depth, the function applies the occupancy-based,
variance-based, and Dirichlet multinomial log-likelihood methods, then
returns a dataframe containing the fitted results.

In general, the function is useful for answering questions such as:

- Do the different inference methods converge toward similar estimates?
- Is there a read depth above which the estimates become relatively
  stable?

<!-- -->

    ##         fit    upper     lower shared_species sample_size source_samples
    ## 1  938.4395 1176.328  747.7239           5288       90000             19
    ## 2 1237.8080 1600.270 1009.2078           5288       90000             19
    ## 3  628.7051       NA        NA           5288       90000             19
    ## 4  879.1461 1096.247  704.2253           5296       90000             19
    ## 5 1264.2186 1645.262 1026.4790           5296       90000             19
    ## 6  625.8604       NA        NA           5296       90000             19
    ##   local_samples    method
    ## 1            15 occupancy
    ## 2            15  variance
    ## 3            15     DM-LL
    ## 4            15 occupancy
    ## 5            15  variance
    ## 6            15     DM-LL

This can be visualized as below, showing how the occupancy, variance and
DM-LL methods estimate change across read depths, and how well estimates
agree with one another.

``` r
ggplot(zheng_profile_results, aes(x = sample_size, y = fit)) +
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
```

<img src="Tutorial_NCM_v1_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

## Manually fitting each method and inspecting plots of each fit

### Preparing abundance datasets for NCM fitting and inference

The `NCM_sample_prep` function takes in three inputs:

- The sample size used for rarefaction prior to conversion to relative
  abundances
- The abundance data corresponding to the source and local communities

What it returns is list consisting of:

- \$Neutral_df: a table containing the means and variances of relative
  abundances, and their sample occupancies for the shared species
  between input source and local communities
- \$summary: summary of the sample prep process
- \$source_only: OTUs/ASVs only seen in the source community
- \$local_only: OTUs/ASVs only seen in the local communty
- \$all_data: data for all species

`Neutral_df` is required by all of the inference methods as a common
input

``` r
# Load necessary function files 
source("Function_definitions.R")

# Load abundance dataset, here we will use the Zheng et. al dataset as an example
ps_zheng <- readRDS("ps_zheng2019.RDS") # this is a phyloseq object

# Set sample size and rarefy abundance table
samplesize <- 90000
seed <- 100 # Set seed number if reproducibility is required  
ps_zheng_rare <- rarefy_even_depth(ps_zheng, sample.size = samplesize, 
                                   rngseed = seed, trimOTUs = F, verbose = F)

# Delineate source community and local community samples
zheng_source <- subset_samples(ps_zheng_rare, sample_type == "IWW")
zheng_local <- subset_samples(ps_zheng_rare, sample_type == "AS")

# Prepare for NCM fitting and inference
zheng_neutraldata <- NCM_sample_prep(samplesize = samplesize, ps_local = zheng_local, 
                                     ps_source = zheng_source, verbose = T, correction = T)
```

    ## [1] "taxa below detection limit (1/sample size) are not considered for analysis"
    ## [1] "# of shared taxa: 5287"
    ## [1] "number of shared taxa above detection limit: 3462"
    ## [1] "Abundance in shared community (%) 93.05"
    ## [1] "number of source-only taxa: 2042"
    ## [1] "number of local-only taxa: 1086"

``` r
head(zheng_neutraldata$neutral_df)
```

    ##      column_index otu_name local_frequency local_abundance local_variance
    ## ASV1            1     ASV1       0.9333333    0.0011051852   6.202996e-06
    ## ASV2            2     ASV2       1.0000000    0.0009355556   3.951326e-06
    ## ASV3            3     ASV3       1.0000000    0.0004096296   8.393063e-08
    ## ASV4            4     ASV4       1.0000000    0.0154392593   1.069456e-04
    ## ASV5            5     ASV5       1.0000000    0.0158355556   2.002515e-05
    ## ASV6            6     ASV6       1.0000000    0.0178948148   1.744902e-04
    ##      source_frequency source_abundance source_variance
    ## ASV1                1      0.027291228    6.951222e-05
    ## ASV2                1      0.023512865    5.414243e-05
    ## ASV3                1      0.018215205    2.105568e-04
    ## ASV4                1      0.003269006    7.393835e-06
    ## ASV5                1      0.002642105    4.162889e-06
    ## ASV6                1      0.003330409    1.008492e-05

### Occupancy-based inference method

The occupancy-based method looks at how frequently each taxa are
detected across samples, fitted against the NCM’s prediction.

For taxon $i$, the model assumes that the local relative abundance $x_i$
follows a beta distribution with shape parameters:

$$
\alpha_i = N_T m p_i
$$

and

$$
\beta_i = N_T m(1-p_i)
$$

where $p_i$ is the relative abundance of taxon $i$ in the source
community, $m$ is the immigration probability, and $N_T$ is the
effective community size.

The expected occurrence probability (i.e., sample occupancy) above the detection threshold $D$
is:

$$
O_i=\int_D^1\frac{x_i^{\alpha_i - 1}(1-x_i)^{\beta_i - 1}}{\mathrm{B}(\alpha_i,\beta_i)}\,dx_i
$$

This integral gives the probability that taxon $i$ occurs above the
detection limit in a local community.

``` r
# Using the Zheng et. al dataset as an example
zheng_occ_fit <- fit_occ(samplesize = samplesize, neutral_data = zheng_neutraldata,
                         remove_ones = T)

zheng_occ_fit$fitting_results
```

    ##              fit               LL               UL              Rsq 
    ##     8.977308e+02     7.176093e+02     1.121525e+03    -1.270811e+00 
    ##      sample size  detection limit shared abundance 
    ##     9.000000e+04     1.111111e-05     9.305000e+01

``` r
# This can be visualized conveniently: 
plot_occupancy <- function(ndf_w_preds, title){
  require(ggplot2)
  
  p <-  ggplot(ndf_w_preds) + aes(x = source_abundance * 100) +
    geom_point(aes(y = local_frequency * 100), colour = "#85838380") + 
    geom_line(aes(y = predicted_freq * 100), colour = "#900000") +
    geom_line(aes(y = pred_LL * 100), color = "#90000070", linetype = "dashed") +
    geom_line(aes(y = pred_UL * 100), color = "#90000070", linetype = "dashed") + 
    coord_trans(x = "log") + scale_x_continuous(breaks = c(1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 10), 
                                                labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1",
                                                           "1", "10")) + 
    xlab("Mean Relative Abundance (%)") + ylab("Sample Occupancy (%)") + ylim(c(0, 100)) +
    theme_bw() + ggtitle(title)
  
  return(p)
}

plot_occupancy(zheng_occ_fit$ndf_w_preds, title = paste("Zheng occupancy fit, R =", samplesize))
```

<img src="Tutorial_NCM_v1_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

### Variance-based inference method

The variance-based method estimates immigration by comparing the
observed variance in taxon relative abundances across local communities
with the variance expected under the Neutral Community Model.

For taxon $i$, the expected variance is:

$$
\widehat{\mathrm{Var}}[X_i]=\frac{p_i(1-p_i)}{N_T m + 1}+\frac{p_i(1-p_i)}{R}
$$

where $X_i$ is the relative abundance of taxon $(i$) in the local
community, $p_i$ is its relative abundance in the source community,
$N_T$ is the effective community size, $m$ is the immigration
probability, and $R$ is the rarefied sequencing depth.

The first term describes the variance expected from neutral community
assembly. When $N_Tm$ is large, immigration strongly links local
communities to the source community, reducing variation among local
samples. When $N_Tm$ is small, local communities are more weakly
connected to the source, and greater variation among local communities
is expected.

The second term accounts for sampling variance caused by finite
sequencing depth. This term is important because observed relative
abundances are estimated from a finite number of reads rather than
measured directly from the entire community.

``` r
# Using the Zheng et. al dataset as an example
zheng_var_fit <- fit_var(samplesize = samplesize, neutral_data = zheng_neutraldata)
```

    ## Waiting for profiling to be done...

``` r
zheng_var_fit$fitting_results
```

    ##              fit               LL               UL              Rsq 
    ##     1.276291e+03     1.037859e+03     1.656929e+03     1.663486e-02 
    ##      sample size  detection limit shared abundance 
    ##     9.000000e+04     1.111111e-05     9.305000e+01

``` r
# This can be visualized conveniently: 
plot_variance <- function(ndf_w_preds, title){
  require(ggplot2)
  
  p <-  ggplot(ndf_w_preds) + aes(x = source_abundance) +
      geom_point(aes(y = local_variance), color = "#858383CC") +
      geom_line(aes(y = predicted_var), colour = "#000090") +
      geom_line(aes(y = predvar_LL), color = "#000090BB", linetype = "dashed") +
      geom_line(aes(y = predvar_UL), color = "#000090BB", linetype = "dashed") +
      scale_y_log10(limits = c(1E-9, 1E-2), n.breaks = 8) + 
      coord_trans(x = "log") + scale_x_continuous(breaks = c(1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 10), 
                                                  labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1",
                                                             "1", "10")) + 
      xlab("Mean Relative Abundance in source (%)") + ylab("Relative Abundance Variance") + 
      theme_bw() + ggtitle(title)
  
  return(p)
}

plot_variance(zheng_var_fit$ndf_w_preds, title = paste("Zheng variance fit, R =", samplesize))
```

<img src="Tutorial_NCM_v1_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

### Dirichlet multinomial log-likelihood method

The Dirichlet multinomial log-likelihood method estimates $N_Tm$ by
maximizing the likelihood of the observed local community abundance
data. Unlike the occupancy- and variance-based methods, which use
taxon-level summary statistics, the DM-LL method uses the full local
community abundance table.

For each local community sample $j$, the observed abundance vector
$\mathbf{y}^{(j)}$ is treated as a multinomial sample drawn from an
unobserved local relative abundance vector, $\hat{\mathbf{p}}$. This
local abundance vector is assumed to vary around the source community
composition according to a Dirichlet distribution parameterized by
$N_Tm\mathbf{p}$:

$$
P(\mathbf{y})=\prod_{j=1}^{J}\int_{\hat{\mathbf{p}}}\mathrm{Mult}\left(\mathbf{y}^{(j)} \mid R,\hat{\mathbf{p}}\right)\mathrm{Dir}\left(\hat{\mathbf{p}} \mid N_T m \mathbf{p}\right)\,d\hat{\mathbf{p}}
$$

where $\mathbf{y}^{(j)}$ is the abundance vector for local sample $j$,
$R$ is the rarefied sequencing depth, $\mathbf{p}$ is the source
community relative abundance vector, and $N_Tm$ is the immigration
parameter being estimated.

After integrating over $\hat{\mathbf{p}}$, the likelihood can be
evaluated in closed form:

$$
P(\mathbf{y})=\prod_{j=1}^{J}\frac{R\,\mathrm{B}(N_T m, R)}{\prod_{k: y_k^{(j)} > 0}y_k^{(j)}\,\mathrm{B}\left(N_T m p_k,y_k^{(j)}\right)}
$$

The method then finds the value of $N_Tm$ that maximizes the
log-likelihood. The log-likelihood is used instead of the raw likelihood
for numerical stability.

Along with the `neutral_data` object returned by `NCM_sample_prep`, this
method also requires the local community abundance table as input
because it evaluates the likelihood of the observed abundance counts
directly.

``` r
zheng_DMLL_fit <- fit_DMLL(neutral_data = zheng_neutraldata, local = as.matrix(otu_table(zheng_local)), 
                           samplesize = samplesize, plot = T, maxNTm = 1E6)
```

<img src="Tutorial_NCM_v1_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />
