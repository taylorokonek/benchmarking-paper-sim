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
library(parallel)


# Catch the “task ID”
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("The task ID is:", task_id))

# load the appropriate dataset
fnames <- list.files("Data/")
# remove admin1_southafrica.graph from fnames
fnames <- fnames[fnames != "admin1_southafrica.graph"]
print(fnames[task_id])
load(file = paste0("Data/", fnames[task_id]))

# get which simulation this is based on fnames
which_sim <- fnames[task_id] %>% str_split("_") %>% unlist %>% nth(2) %>% as.numeric()
which_seed <- fnames[task_id] %>% str_split("_") %>% unlist %>% nth(4) %>% str_split("\\.") %>% unlist %>% nth(1) %>% as.numeric()

# Rejection Sampler: fit model in INLA -------------------------------------------------------

temp <- admin1_mat
rownames(temp) <- colnames(temp) <- 1:n

priors =  list(phi = list(prior="logitbeta",params=c(0.5,0.5), fixed = FALSE),
               prec = list(prior="pc.prec",params=c(1,0.01), fixed = FALSE))

formula <- Value ~ 1 + f(Region, model = "bym2", graph = temp, 
                         hyper = priors, constr = TRUE) +
  f(cluster, model = "iid", hyper = list(prec = list(prior = "loggamma",params = c(0.1,0.1), fixed = FALSE)))

start1 <- Sys.time()
result <- inla(formula = formula,
               family = "binomial",
               data = df %>% mutate(Region = as.numeric(Region)),
               Ntrials = N,
               control.compute = list(dic = TRUE, 
                                      mlik = TRUE, cpo = TRUE, openmp.strategy = "default", 
                                      config = TRUE))
stop1 <- Sys.time()

start2 <- Sys.time()
samp <- inla.posterior.sample(n = 100000, result = result)
stop2 <- Sys.time()

# gather samples into posterior draws
region_idx <- samp[[1]]$latent %>% rownames() %>% str_detect("Region") %>% which() %>% head(n)
intercept_idx <- samp[[1]]$latent %>% rownames() %>% str_detect("Intercept") %>% which()

region_mat <- matrix(nrow = length(region_idx), ncol = 100000)
intercept_mat <- matrix(nrow = 1, ncol = 100000)
cluster_tau_mat <- matrix(nrow = 1, ncol = 100000)

# fill in samples
for (i in 1:100000) {
  region_mat[,i] <- samp[[i]]$latent[region_idx]
  intercept_mat[,i] <- samp[[i]]$latent[intercept_idx]
  cluster_tau_mat[,i] <- samp[[i]]$hyperpar[3]
}
cluster_sd_mat <- 1/sqrt(cluster_tau_mat)

# lono correction
h = 16 * sqrt(3) / (15 * pi)
denom <- sqrt(1 + (h ^ 2) * (cluster_sd_mat ^ 2))

# combine samples into draws from the expit(linear predictor)
samp_mat <- region_mat + intercept_mat[rep(1,n),]
for (i in 1:ncol(samp_mat)) {
  samp_mat[,i] <- samp_mat[,i] / denom[,i]
}

samp_mat <- expit(samp_mat)

# get unbenchmarked national samples
unbenched_natl <- apply(samp_mat,2,function(x) {mean(x)})

# Rejection Sampler: Benchmark -------------------------------------------------------

start3 <- Sys.time()
benched_samps_rs <- benchmark_sampler(posterior_samps = samp_mat,
                                      nregion = n,
                                      ntime = 1,
                                      natl = sim_df[which_sim,]$z,
                                      natl_sd = sqrt(sim_df[which_sim,]$var_z),
                                      pop_weights = w_i)
stop3 <- Sys.time()

# get number of samples needed to get approximately 1000 accepted
new_nsamp <- ceiling(1000 * 1/ (benched_samps_rs$prop_accepted))

if (benched_samps_rs$prop_accepted == 0) {
  stop("NO SAMPLES ACCEPTED FOR RS - TRY AGAIN WITH LARGER INITIAL SAMPLE SIZE")
}

# Rejection Sampler: Re-do with appropriate number of samples -------------

if(new_nsamp < 1000) {
  start2 <- Sys.time()
  samp <- inla.posterior.sample(n = new_nsamp, result = result)
  stop2 <- Sys.time()
} else {
  break_vals <- unique(c(seq(0, new_nsamp, by = 1000), new_nsamp))
  
  start2 <- Sys.time()
  samp <- inla.posterior.sample(n = 1000, result = result)
  for (i in 2:(length(break_vals)-2)) {
    samp <- c(samp, inla.posterior.sample(n = 1000, result = result))
    print(i)
  }
  samp <- c(samp, inla.posterior.sample(n = new_nsamp - ((length(break_vals)-2) * 1000), result = result))
  stop2 <- Sys.time()
}


# gather samples into posterior draws
region_idx <- samp[[1]]$latent %>% rownames() %>% str_detect("Region") %>% which() %>% head(n)
intercept_idx <- samp[[1]]$latent %>% rownames() %>% str_detect("Intercept") %>% which()

region_mat <- matrix(nrow = length(region_idx), ncol = new_nsamp)
intercept_mat <- matrix(nrow = 1, ncol = new_nsamp)
cluster_tau_mat <- matrix(nrow = 1, ncol = new_nsamp)

# fill in samples
for (i in 1:new_nsamp) {
  region_mat[,i] <- samp[[i]]$latent[region_idx]
  intercept_mat[,i] <- samp[[i]]$latent[intercept_idx]
  cluster_tau_mat[,i] <- samp[[i]]$hyperpar[3]
}
cluster_sd_mat <- 1/sqrt(cluster_tau_mat)

# lono correction
h = 16 * sqrt(3) / (15 * pi)
denom <- sqrt(1 + (h ^ 2) * (cluster_sd_mat ^ 2))

# combine samples into draws from the expit(linear predictor)
samp_mat <- region_mat + intercept_mat[rep(1,n),]
for (i in 1:ncol(samp_mat)) {
  samp_mat[,i] <- samp_mat[,i] / denom[,i]
}

samp_mat <- expit(samp_mat)

# get unbenchmarked national samples
unbenched_natl <- apply(samp_mat,2,function(x) {mean(x)})

start3 <- Sys.time()
benched_samps_rs <- benchmark_sampler(posterior_samps = samp_mat,
                                      nregion = n,
                                      ntime = 1,
                                      natl = sim_df[which_sim,]$z,
                                      natl_sd = sqrt(sim_df[which_sim,]$var_z),
                                      pop_weights = w_i)
stop3 <- Sys.time()

# number of accepted samples
benched_samps_rs$natl_list[[1]] %>% length()

# save new_nsamp for return_mat
new_nsamp_rs <- ceiling(new_nsamp)

# MH: fit model in INLA ----------------------------------------------------------------

prior.fixed <- list(mean.intercept = logit(sim_df[which_sim,]$z), prec.intercept = 10)

start1_mh <- Sys.time()
result_mh <- inla(formula = formula,
                  family = "binomial",
                  data = df %>% mutate(Region = as.numeric(Region)),
                  Ntrials = N,
                  control.fixed = prior.fixed,
                  control.compute = list(dic = TRUE, 
                                         mlik = TRUE, cpo = TRUE, openmp.strategy = "default", 
                                         config = TRUE))
stop1_mh <- Sys.time()

start2_mh <- Sys.time()
samp <- inla.posterior.sample(n = 1000, result = result_mh)
for (i in 2:8) {
  samp <- c(samp, inla.posterior.sample(n = 1000, result = result_mh))
  print(i)
}
stop2_mh <- Sys.time()

# gather samples into posterior draws
region_idx <- samp[[1]]$latent %>% rownames() %>% str_detect("Region") %>% which() %>% head(n)
intercept_idx <- samp[[1]]$latent %>% rownames() %>% str_detect("Intercept") %>% which()

region_mat <- matrix(nrow = length(region_idx), ncol = 8000)
intercept_mat <- matrix(nrow = 1, ncol = 8000)
cluster_tau_mat <- matrix(nrow = 1, ncol = 8000)

# fill in samples
for (i in 1:8000) {
  region_mat[,i] <- samp[[i]]$latent[region_idx]
  intercept_mat[,i] <- samp[[i]]$latent[intercept_idx]
  cluster_tau_mat[,i] <- samp[[i]]$hyperpar[3]
}
cluster_sd_mat <- 1/sqrt(cluster_tau_mat)

# lono correction
h = 16 * sqrt(3) / (15 * pi)
denom <- sqrt(1 + (h ^ 2) * (cluster_sd_mat ^ 2))

# combine samples into draws from the expit(linear predictor)
samp_mat_mh <- region_mat + intercept_mat[rep(1,n),]
for (i in 1:ncol(samp_mat_mh)) {
  samp_mat_mh[,i] <- samp_mat_mh[,i] / denom[,i]
}

samp_mat_mh <- expit(samp_mat_mh)


# MH: Benchmark -------------------------------------------------------

# THIS SHOULD BE IN PARALLEL BUT I CAN'T FIGURE OUT HOW THE HELL TO DO THAT
start_vals <- seq(from = 0, to = ncol(samp_mat_mh), length.out = 5)
start_vals <- start_vals[1:(length(start_vals) - 1)] + 1
end_vals <- seq(from = 0, to = ncol(samp_mat_mh), length.out = 5)[2:(length(start_vals) + 1)]

benched_samps_mh_lst <- list()

start3_mh <- Sys.time()
for (i in 1:4) {
  benched_samps_mh_lst[[i]] <- benchmark_mh(posterior_samps = samp_mat_mh[,start_vals[i]:end_vals[i]],
                                            posterior_samps_fe = intercept_mat[,start_vals[i]:end_vals[i]] %>% matrix(nrow = 1),
                                            nregion = n,
                                            ntime = 1,
                                            natl = sim_df[which_sim,]$z,
                                            natl_sd = sqrt(sim_df[which_sim,]$var_z),
                                            pop_weights = w_i,
                                            fe_prior_means = prior.fixed[[1]],
                                            fe_prior_variances = prior.fixed[[2]])
}
stop3_mh <- Sys.time()

draw_mat <- cbind(benched_samps_mh_lst[[1]]$natl_list[[1]],
                  benched_samps_mh_lst[[2]]$natl_list[[1]],
                  benched_samps_mh_lst[[3]]$natl_list[[1]],
                  benched_samps_mh_lst[[4]]$natl_list[[1]])

# remove burn-in samples
draw_mat %>% dim()
draw_mat <- draw_mat[1000:1999,] 

ess_bulk_mh <- ess_bulk(draw_mat)
ess_tail_mh <- ess_tail(draw_mat)
ess_basic_mh <- ess_basic(draw_mat)

# get number of samples needed to get approximately 1000 accepted, after 1000 burn in for each chain
new_nsamp <- 1000 * (4000) / ess_bulk_mh
new_nsamp <- ceiling(new_nsamp)

# add back in number of burn-in samples we need
new_nsamp <- new_nsamp + 4000

# save new_nsamp for later
new_nsamp_mh <- new_nsamp - 4000

# MH: Re-do with appropriate number of samples -------------

if(new_nsamp < 1000) {
  start2_mh <- Sys.time()
  samp <- inla.posterior.sample(n = new_nsamp, result = result_mh)
  stop2_mh <- Sys.time()
} else {
  break_vals <- unique(c(seq(0, new_nsamp, by = 1000), new_nsamp))
  
  start2_mh <- Sys.time()
  samp <- inla.posterior.sample(n = 1000, result = result_mh)
  for (i in 2:(length(break_vals)-2)) {
    samp <- c(samp, inla.posterior.sample(n = 1000, result = result_mh))
    print(i)
  }
  samp <- c(samp, inla.posterior.sample(n = new_nsamp - ((length(break_vals)-2) * 1000), result = result_mh))
  stop2_mh <- Sys.time()
}


# gather samples into posterior draws
region_idx <- samp[[1]]$latent %>% rownames() %>% str_detect("Region") %>% which() %>% head(n)
intercept_idx <- samp[[1]]$latent %>% rownames() %>% str_detect("Intercept") %>% which()

region_mat <- matrix(nrow = length(region_idx), ncol = new_nsamp)
intercept_mat <- matrix(nrow = 1, ncol = new_nsamp)
cluster_tau_mat <- matrix(nrow = 1, ncol = new_nsamp)

# fill in samples
for (i in 1:new_nsamp) {
  region_mat[,i] <- samp[[i]]$latent[region_idx]
  intercept_mat[,i] <- samp[[i]]$latent[intercept_idx]
  cluster_tau_mat[,i] <- samp[[i]]$hyperpar[3]
}
cluster_sd_mat <- 1/sqrt(cluster_tau_mat)

# lono correction
h = 16 * sqrt(3) / (15 * pi)
denom <- sqrt(1 + (h ^ 2) * (cluster_sd_mat ^ 2))

# combine samples into draws from the expit(linear predictor)
samp_mat_mh <- region_mat + intercept_mat[rep(1,n),]
for (i in 1:ncol(samp_mat_mh)) {
  samp_mat_mh[,i] <- samp_mat_mh[,i] / denom[,i]
}

samp_mat_mh <- expit(samp_mat_mh)

# do benchmarking

start_vals <- seq(from = 0, to = ncol(samp_mat_mh), length.out = 5) %>% ceiling()
start_vals <- start_vals[1:(length(start_vals) - 1)] + 1
end_vals <- ceiling(seq(from = 0, to = ncol(samp_mat_mh), length.out = 5))[2:(length(start_vals) + 1)]

benched_samps_mh_lst <- list()

start3_mh <- Sys.time()
for (i in 1:4) {
  benched_samps_mh_lst[[i]] <- benchmark_mh(posterior_samps = samp_mat_mh[,start_vals[i]:end_vals[i]],
                                            posterior_samps_fe = intercept_mat[,start_vals[i]:end_vals[i]] %>% matrix(nrow = 1),
                                            nregion = n,
                                            ntime = 1,
                                            natl = sim_df[which_sim,]$z,
                                            natl_sd = sqrt(sim_df[which_sim,]$var_z),
                                            pop_weights = w_i,
                                            fe_prior_means = prior.fixed[[1]],
                                            fe_prior_variances = prior.fixed[[2]])
}
stop3_mh <- Sys.time()

draw_mat <- cbind(benched_samps_mh_lst[[1]]$natl_list[[1]],
                  benched_samps_mh_lst[[2]]$natl_list[[1]],
                  benched_samps_mh_lst[[3]]$natl_list[[1]],
                  benched_samps_mh_lst[[4]]$natl_list[[1]])

# remove burn-in
draw_mat <- draw_mat[1000:nrow(draw_mat),] 

ess_bulk_mh <- ess_bulk(draw_mat)
ess_tail_mh <- ess_tail(draw_mat)
ess_basic_mh <- ess_basic(draw_mat)
rhat_mh <- rhat(draw_mat)

# Zhang Bryant: fit model in STAN -------------------------------------------------------

# transform admin1_mat into an ICAR
temp <- ifelse(admin1_mat > 0, -1, admin1_mat)
diag(temp) <- rowSums(temp) * -1

# get scaling factor for ICAR
diag(temp)[1]
Q_scaled <- as.matrix(inla.scale.model(temp, constr = list(A=matrix(1,1,dim(temp)[1]), e=0)))
scaling_factor <- diag(Q_scaled)[1] / diag(temp)[1]

# create node_i, node_j vectors
sardinia_g <- read.delim("Data/admin1_southafrica.graph")
node_i <- c()
node_j <- c()
for (i in 1:nrow(sardinia_g)) {
  row_string <- sardinia_g[i,] %>% as.character() %>% str_split(" ") %>% unlist() %>% as.numeric()
  # add to node_i
  if (length(which(row_string[3:length(row_string)] > row_string[1])) > 0) {
    node_i <- c(node_i, rep(row_string[1], length(row_string[3:length(row_string)][which(row_string[3:length(row_string)] > row_string[1])])))
    # add to node_j
    node_j <- c(node_j, row_string[3:length(row_string)][which(row_string[3:length(row_string)] > row_string[1])])
  }
  
}

# create data for STAN
stan_dat_benched <- list(N = n, 
                         N_pairs = length(node_i),
                         N_ic = df$N,
                         y_ic = df$Value,
                         nclusts = length(unique(df$cluster)),
                         hiv_adj = 1,
                         region_id = as.numeric(df$Region),
                         cluster_id = df$cluster,
                         node_i = node_i, 
                         node_j = node_j, 
                         scaling_factor = unname(scaling_factor),
                         alpha = 0.01, 
                         U = 1,
                         natl = sim_df[which_sim,]$z,
                         natl_sd = sqrt(sim_df[which_sim,]$var_z),
                         pop_weights = w_i)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

start1_stan <- Sys.time()
fit_benched <- stan(file = "binary_intercepts_bym2_benched.stan",
                    data = stan_dat_benched,
                    chains = 4, # default
                    iter = 2000, # default
                    warmup = 1000, # default
                    cores = 2,
                    seed = which_seed,
                    control = list(max_treedepth = 12))

benched_admin1 <- rstan::extract(fit_benched, pars = "risk_i")
benched_admin1 <- benched_admin1$risk_i %>% t() 
end1_stan <- Sys.time()

stan_draws <- as.array(fit_benched)
# ess_bulk_stan_i <- rep(NA, n)
# ess_tail_stan_i <- rep(NA, n)
# ess_basic_stan_i <- rep(NA, n)
# for (i in 1:n) {
#   draw_mat <- posterior::extract_variable_matrix(stan_draws, variable = paste0("risk_i[",i,"]"))
#   ess_bulk_stan_i[i] <- ess_bulk(draw_mat)
#   ess_tail_stan_i[i] <- ess_tail(draw_mat)
#   ess_basic_stan_i[i] <- ess_basic(draw_mat)
# }
stan_draws_risk <- stan_draws[,,paste0("risk_i[",1:n,"]")] 
draw_mat <- apply(stan_draws_risk, c(1,2), mean)
ess_bulk_stan <- ess_bulk(draw_mat)
ess_tail_stan <- ess_tail(draw_mat)
ess_basic_stan <- ess_basic(draw_mat)
rhat_stan <- rhat(draw_mat)

# get number of samples needed to get approximately 1000 accepted, after 1000 burn in for each chain
new_nsamp <- 1000 * (4000) / ess_bulk_stan
new_nsamp <- ceiling(new_nsamp/4)

new_nsamp_stan <- new_nsamp*4

# Zhang Bryant: re-do with appropriate number of samples -------------------------------------------------------


start1_stan <- Sys.time()
fit_benched <- stan(file = "binary_intercepts_bym2_benched.stan",
                    data = stan_dat_benched,
                    chains = 4, # default
                    iter = 1000 + new_nsamp, # default
                    warmup = 1000, # default
                    cores = 2,
                    seed = which_seed,
                    control = list(max_treedepth = 12))

benched_admin1 <- rstan::extract(fit_benched, pars = "risk_i")
benched_admin1 <- benched_admin1$risk_i %>% t() 
end1_stan <- Sys.time()

stan_draws <- as.array(fit_benched)
stan_draws_risk <- stan_draws[,,paste0("risk_i[",1:n,"]")] 
draw_mat <- apply(stan_draws_risk, c(1,2), mean)

ess_bulk_stan <- ess_bulk(draw_mat)
ess_tail_stan <- ess_tail(draw_mat)
ess_basic_stan <- ess_basic(draw_mat)
rhat_stan <- rhat(draw_mat)


# Compute times -----------------------------------------------------------

# speed comparison

## Rejection Sampler
# model-fitting time
a <- (stop1 - start1)
units(a) <- "secs"
# sampling time (time it would take to get 1000 samples accepted)
b <- (stop2 - start2) 
units(b) <- "secs"
# benchmarking time
c <- (stop3 - start3)
units(c) <- "secs"

rs_time <- a + b + c
units(rs_time) <- "secs"

# MH Algorithm
# model-fitting time
a_mh <- (stop1_mh - start1_mh)
units(a_mh) <- "secs"
# sampling time (time it would take to get 1000 samples accepted)
b_mh <- (stop2_mh - start2_mh) 
units(b_mh) <- "secs"
# benchmarking time
c_mh <- (stop3_mh - start3_mh)
units(c_mh) <- "secs"

mh_time <- a_mh + b_mh + c_mh
units(mh_time) <- "secs"

# STAN 

stan_time <- difftime(end1_stan, start1_stan, units='secs') 

# Combine results into dataframe ------------------------------------------

results_df <- data.frame(sim_id = which_sim,
                         method = c("rs","mh","stan","unbenched_inla"),
                         nsamps = sim_df[which_sim,]$nsamps,
                         z = sim_df[which_sim,]$z,
                         var_z = sim_df[which_sim,]$var_z,
                         seed = which_seed,
                         fname = fnames[task_id])

results_df$model_time <- c(a %>% as.numeric, 
                           a_mh %>% as.numeric, 
                           NA, 
                           NA)

results_df$sampling_time <- c(b %>% as.numeric,
                              b_mh %>% as.numeric,
                              NA,
                              NA)
results_df$benchmarking_time <- c(c %>% as.numeric,
                                  c_mh %>% as.numeric,
                                  NA, 
                                  NA)

results_df$total_time <- c(rs_time %>% as.numeric,
                           mh_time %>% as.numeric,
                           stan_time %>% as.numeric,
                           NA)

results_df$ess_bulk <- c(NA, ess_bulk_mh, ess_bulk_stan, NA)
results_df$ess_tail <- c(NA, ess_tail_mh, ess_tail_stan, NA)
results_df$accepted_samps <- c(benched_samps_rs$natl_list[[1]] %>% length(), rep(NA, 3))
results_df$rhat <- c(NA, rhat_mh, rhat_stan, NA)

# add in posterior means, medians, variance for national level data
results_df$posterior_mean <- c(
  benched_samps_rs$natl_list[[1]] %>% mean(),
  c(benched_samps_mh_lst[[1]]$natl_list[[1]],
    benched_samps_mh_lst[[2]]$natl_list[[1]],
    benched_samps_mh_lst[[3]]$natl_list[[1]],
    benched_samps_mh_lst[[4]]$natl_list[[1]]) %>% mean(),
  draw_mat %>% as.vector() %>% mean(),
  unbenched_natl %>% mean)

results_df$posterior_var <- c(
  benched_samps_rs$natl_list[[1]] %>% var(),
  c(benched_samps_mh_lst[[1]]$natl_list[[1]],
    benched_samps_mh_lst[[2]]$natl_list[[1]],
    benched_samps_mh_lst[[3]]$natl_list[[1]],
    benched_samps_mh_lst[[4]]$natl_list[[1]]) %>% var(),
  draw_mat %>% as.vector() %>% var(),
  unbenched_natl %>% var)

results_df$posterior_median <- c(
  benched_samps_rs$natl_list[[1]] %>% median(),
  c(benched_samps_mh_lst[[1]]$natl_list[[1]],
    benched_samps_mh_lst[[2]]$natl_list[[1]],
    benched_samps_mh_lst[[3]]$natl_list[[1]],
    benched_samps_mh_lst[[4]]$natl_list[[1]]) %>% median(),
  draw_mat %>% as.vector() %>% median(),
  unbenched_natl %>% median)

# add in number of samples taken
results_df$nsamps_beforebenchmarking_withoutburnin <- c(new_nsamp_rs,
                                                        new_nsamp_mh,
                                                        new_nsamp_stan,
                                                        NA)

# Save results ------------------------------------------------------------

out_name <- fnames[task_id]
out_name <- paste0("results_", out_name)

save(results_df, file = paste0("Results/", out_name))








