# Author: Taylor Okonek
# Project: Benchmarking 
# Simulation comparing run times for rejection sampler, MH algorithm, MCMC in Stan


library(tidyverse)


# Load results ------------------------------------------------------------

fnames <- list.files("Results")

df <- NULL
for (i in 1:length(fnames)) {
  load(paste0("Results/",fnames[i]))
  df <- rbind(df, results_df)
}

# Code used to generate data (sim_df)
n <- 9 # number of areas
nsamps <- c(5, 10, 100, 1000) # number of samples in each area
z <- c(0.3, 0.29, 0.28) # natl mean
var_z <- c(0.0001, 0.01, 0.00001)# natl variance
w_i <- rep(1/n, n) # just equal weights in each region
p_i <- 0.27 + seq(1:n)/100 # binomial probabilities (vary by area)
N <- 100 # N in each cluster (binomial totals total)

sim_df <- expand.grid(nsamps = nsamps, z = z, var_z = var_z)
sim_df$sim_id <- 1:nrow(sim_df)
sim_df$seed_start <- sim_df$sim_id*1000
sim_df$seed_stop <- sim_df$seed_start + 9


# compare times -----------------------------

# Var = 0.0001
df %>%
  filter(z == 0.29, var_z == 0.0001) %>%
  filter(method != "unbenched_inla") %>%
  mutate(Method = ifelse(method == "mh", "MH",
                         ifelse(method == "rs", "RS", "STAN"))) %>%
  mutate(nsamps = factor(nsamps)) %>%
  ggplot(aes(y = total_time, x = Method, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~nsamps, scales = "free", nrow = 1) +
  xlab("Method") +
  ylab("Run-time (seconds)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) +
  scale_fill_viridis_d()
ggsave(filename = "Plots/z_29_var_0001.pdf")

df %>%
  filter(z == 0.3, var_z == 0.0001) %>%
  filter(method != "unbenched_inla") %>%
  mutate(Method = ifelse(method == "mh", "MH",
                         ifelse(method == "rs", "RS", "STAN"))) %>%
  mutate(nsamps = factor(nsamps)) %>%
  ggplot(aes(y = total_time, x = Method, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~nsamps, scales = "free", nrow = 1) +
  xlab("Method") +
  ylab("Run-time (seconds)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) +
  scale_fill_viridis_d()
ggsave(filename = "Plots/z_30_var_0001.pdf")

df %>%
  filter(z == 0.29, var_z == 0.01) %>%
  filter(method != "unbenched_inla") %>%
  mutate(Method = ifelse(method == "mh", "MH",
                         ifelse(method == "rs", "RS", "STAN"))) %>%
  mutate(nsamps = factor(nsamps)) %>%
  ggplot(aes(y = total_time, x = Method, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~nsamps, scales = "free", nrow = 1) +
  xlab("Method") +
  ylab("Run-time (seconds)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) +
  scale_fill_viridis_d()
ggsave(filename = "Plots/z_29_var_01.pdf")

# numerical summaries
df %>%
  filter(z == 0.29, var_z == 0.01, nsamps == 5) %>%
  group_by(method) %>%
  summarise(time = median(total_time))

df %>%
  filter(z == 0.29, var_z == 0.01, nsamps == 1000) %>%
  group_by(method) %>%
  summarise(time = median(total_time))

df %>%
  filter(method != "unbenched_inla") %>%
  group_by(z, var_z, nsamps, method) %>%
  summarise(time = median(total_time),
            min = min(total_time),
            max = max(total_time))

df %>%
  filter(z == 0.3, var_z == 0.01) %>%
  filter(method != "unbenched_inla") %>%
  mutate(Method = ifelse(method == "mh", "MH",
                         ifelse(method == "rs", "RS", "STAN"))) %>%
  mutate(nsamps = factor(nsamps)) %>%
  ggplot(aes(y = total_time, x = Method, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~nsamps, scales = "free", nrow = 1) +
  xlab("Method") +
  ylab("Run-time (seconds)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) +
  scale_fill_viridis_d()
ggsave(filename = "Plots/z_30_var_01.pdf")
