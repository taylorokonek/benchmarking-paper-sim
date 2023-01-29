# Author: Taylor Okonek
# Project: Benchmarking
# Simulation comparing run times for rejection sampler, MH algorithm, MCMC in Stan

# Model: Binomial w/intercept, BYM2, cluster-level random effect
# loggamma(0.1, 0.1) prior on all precisions
# flat prior on intercept in unbenchmarked model N(0,0.001^-1)
# logitbeta(0.5, 0.5) prior on phi in BYM2

library(tidyverse)
library(SUMMER)
library(INLA)
library(spdep)
library(rgdal)
library(stbench)
library(rstan)
library(posterior)

# Generate Data -----------------------------------------------------------

n <- 9 # number of areas
nsamps <- c(5, 10, 100, 1000) # number of samples in each area
z <- c(0.3, 0.29) # natl mean
var_z <- c(0.0001, 0.01)# natl variance
w_i <- rep(1/n, n) # just equal weights in each region
p_i <- 0.27 + seq(1:n)/100 # binomial probabilities (vary by area)
N <- 100 # N in each cluster (binomial totals total)

sim_df <- expand.grid(nsamps = nsamps, z = z, var_z = var_z)
sim_df$sim_id <- 1:nrow(sim_df)
sim_df$seed_start <- sim_df$sim_id*1000
sim_df$seed_stop <- sim_df$seed_start + 9

# use South Africa graph file
poly_adm1 <- readOGR(dsn = "gadm36_ZAF_shp/",
                     layer = "gadm36_ZAF_1")
admin1_mat <- poly2nb(SpatialPolygons(poly_adm1@polygons))
nb2INLA(file = "admin1_southafrica.graph",admin1_mat)
admin1_mat <- nb2mat(admin1_mat, zero.policy = TRUE)
# make row and column names NAME_1
colnames(admin1_mat) <- rownames(admin1_mat) <- factor(poly_adm1$NAME_1)
admin1_names <- data.frame(GADM = poly_adm1@data$NAME_1,
                           Internal = rownames(admin1_mat))

# which sim
for (l in 1:nrow(sim_df)) {
  which_sim <- l
  
  for (k in 1:10) {
    # fill in binomial data
    curr_seed <- sim_df$seed_start[which_sim] + (k - 1)
    set.seed(curr_seed)
    samps <- matrix(ncol = n,
                    nrow = sim_df$nsamps[which_sim])
    for (i in 1:nrow(samps)) {
      for (j in 1:ncol(samps)) {
        samps[i,j] <- rbinom(1, size = N, prob = p_i[j])
      }
    }
    
    df <- as.data.frame(samps)
    colnames(df) <- poly_adm1$NAME_1
    
    # make dataframe
    df <- df %>%
      gather(Region, Value) %>%
      mutate(N = N) %>%
      mutate(Region = factor(Region))
    
    df$cluster <- 1:nrow(df)
    
    # save dataframe
    save(df, sim_df, n, w_i, p_i, admin1_mat, file = paste0("df_",which_sim,"_seed_",curr_seed,".rds"))
  }
  
  print(paste0("Sim: ", which_sim))
}

