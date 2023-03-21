library(slendr)
init_env()

# chimpanzee outgroup
chimp <- population("CHIMP", time = 7e6, N = 5000)

# two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)

# Neanderthal population splitting at 600 ky ago from modern humans
# (becomes extinct by 40 ky ago)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000)

# compile the entire model into a single object
model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  generation_time = 30
)

# verify visually
plot_model(model)
plot_model(model, sizes = FALSE)
plot_model(model, log = TRUE)

# simulate 100Mb of sequence
t1 <- Sys.time()

ts <-
  msprime(model, sequence_length = 100e6, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)

t2 <- Sys.time()
t2 - t1
# Time difference of 2.799282 mins

ts_save(ts, file = "data/ex2.trees")

ts <- ts_load("data/ex2.trees", model)
ts

# extract samples as a list of names --------------------------------------

library(dplyr)
library(ggplot2)

samples <- ts_samples(ts) %>% split(., .$pop) %>% lapply(pull, "name")

str(samples)
samples$EUR %>% head(5)

# compute diversity -------------------------------------------------------

pi <- ts_diversity(ts, sample_sets = samples)

arrange(pi, diversity)

# compute divergence ------------------------------------------------------

div <- ts_divergence(ts, sample_sets = samples)

arrange(div, desc(divergence))

# F_ST --------------------------------------------------------------------

fst <- ts_fst(ts, sample_sets = samples)

arrange(fst, desc(Fst))

# f-stats -----------------------------------------------------------------

ts_f3(ts, A = "AFR_1", B = "EUR_1", C = "CHIMP_1")
ts_f4(ts, W = "AFR_1", X = "EUR_1", Y = "NEA_1", Z = "CHIMP_1")

# allele frequency spectrum -----------------------------------------------
afs_afr <- ts_afs(ts, list(sample = c("AFR_1", "AFR_2", "AFR_3", "AFR_4", "AFR_5")), polarised = TRUE) %>% .[-c(1, length(.))]
afs_nea <- ts_afs(ts, list(sample = c("NEA_1", "NEA_2", "NEA_3", "NEA_4", "NEA_5")), polarised = TRUE) %>% .[-c(1, length(.))]
afs_eur <- ts_afs(ts, list(sample = c("EUR_1", "EUR_2", "EUR_3", "EUR_4", "EUR_5")), polarised = TRUE) %>% .[-c(1, length(.))]
afs_chimp <- ts_afs(ts, list(sample = c("CHIMP_1", "CHIMP_2", "CHIMP_3", "CHIMP_4", "CHIMP_5")), polarised = TRUE) %>% .[-c(1, length(.))]

ggplot() +
  geom_line(aes(1:9, afs_afr, color = "AFR")) +
  geom_line(aes(1:9, afs_nea, color = "NEA")) +
  geom_line(aes(1:9, afs_eur, color = "EUR")) +
  geom_line(aes(1:9, afs_chimp, color = "CHIMP")) +
  scale_color_manual(values = c("green", "blue", "red", "pink")) +
  labs(y = "frequency", x = "allele count")

plot_model(model, log = TRUE)
