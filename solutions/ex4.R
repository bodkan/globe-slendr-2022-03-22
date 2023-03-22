library(slendr)
init_env()

# chimpanzee outgroup
chimp <- population("CHIMP", time = 7e6, N = 5000)

# two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)

# Neanderthal population splitting at 600 ky ago from modern humans
nea <- population("NEA", parent = afr, time = 600e3, N = 1000)

gf <- gene_flow(from = nea, to = eur, rate = 0.05, start = 55000, end = 50000)

# compile the entire model into a single object
model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  gene_flow = gf,
  generation_time = 30
)

# verify visually
plot_model(model, proportions = TRUE)
plot_model(model, sizes = FALSE, proportions = TRUE)
plot_model(model, log = TRUE, proportions = TRUE)


# this is new here --------------------------------------------------------
emh_dates <- runif(n = 40, min = 5000, max = 40000)
nea_dates <- c(70000, 40000)

present_samples <- schedule_sampling(model, times = 0,         list(chimp, 1), list(afr, 5), list(eur, 10))
nea_samples     <- schedule_sampling(model, times = nea_dates, list(nea, 1))
emh_samples     <- schedule_sampling(model, times = emh_dates, list(eur, 1))

samples <- rbind(present_samples, emh_samples, nea_samples)

ts <- msprime(model, samples = samples, sequence_length = 100e6, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)
ts_save(ts, "data/ex4.trees")

ts <- ts_load("data/ex4.trees", model)
ts


# solution to exercise #4 -------------------------------------------------

library(dplyr)
library(ggplot2)

# extract table with names and times of sampled Europeans (ancient and present day)
eur_inds <- ts_samples(ts) %>% filter(pop == "EUR")
eur_inds

# compute f4-ration statistic (this will take ~30s)
#
# from Petr et al., PNAS 2019:
#
#                f4(Altai Nea. (A), Chimp (O); X, African (C))
# alpha = ---------------------------------------------------------------
#           f4(Altai Nea. (A), Chimp (O); Vindija Nea. (B), African (C))

eur_inds$nea <- ts_f4ratio(ts, X = eur_inds$name, A = "NEA_1", B = "NEA_2", C = "AFR_1", O = "CHIMP_1")$alpha

eur_inds %>%
  ggplot(aes(time, nea)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, color = "red", linewidth = 0.5) +
    xlim(40000, 0) + coord_cartesian(ylim = c(0, 0.1)) +
    labs(x = "time [years ago]", y = "Neanderthal ancestry proportion")
