library(slendr)
init_env()

# chimpanzee outgroup
chimp <- population("CHIMP", time = 7e6, N = 5000)

# two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)

# Neanderthal population splitting at 600 ky ago from modern humans
nea <- population("NEA", parent = afr, time = 600e3, N = 1000)

# this is new here --------------------------------------------------------
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

# simulate 100Mb of sequence
ts <-
  msprime(model, sequence_length = 100e6, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)
ts_save(ts, file = "data/ex3.trees")

ts <- ts_load("data/ex3.trees", model)
ts

library(dplyr)
library(ggplot2)

# extract samples as a list of names

samples <- ts_samples(ts) %>% split(., .$pop) %>% lapply(pull, "name")

str(samples)
samples$EUR %>% head(5)

# solution to exercise #3 -------------------------------------------------

# we can test for the presence of introgression with:
#   f4(AFR, EUR; NEA, CHIMP)
#
#      (A)  (B)  (B)  (A)
#      (B)  (A)  (B)  (A)
#
# which compares the proportion of allele sharing patterns BABA vs ABBA

ts_f4(ts, W = "AFR_1", X = "AFR_2", Y = "NEA_1", Z = "CHIMP_1")

ts_f4(ts, W = "AFR_1", X = "EUR_1", Y = "NEA_1", Z = "CHIMP_1")

afr_samples <- samples$AFR %>% sample(25)
eur_samples <- samples$EUR %>% sample(25)

f4_afr <- lapply(afr_samples, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1")) %>% bind_rows()
f4_eur <- lapply(eur_samples, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1")) %>% bind_rows()

f4_afr$pop <- "AFR"
f4_eur$pop <- "EUR"

f4 <- rbind(f4_afr, f4_eur)

ggplot(f4, aes(pop, f4, color = pop)) +
  geom_boxplot() +
  geom_jitter() +
  geom_hline(yintercept = 0, linetype = 2)
