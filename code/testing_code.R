# Code for testing allele_inheritance functions

source("allele-inheritance_all-in.R")

# Estimate allele frequencies
frequencies(file = "../data/all_genotypes_20.csv", nLoci = 20)

# Calculate ai for observed offspring
ai(pfile = "../data/parent_genotypes_20.csv", ofile = "../data/calf_genotypes_20.csv", ffile = "freqs.csv", nLoci = 20, nTriads = 3)

# Simulate offspring and calculate their ai
sim(pfile = "../data/parent_genotypes_20.csv", ffile = "freqs.csv", nLoci = 20, nTriads = 3, iterations = 1000)

# Visualize results
obs <- read.table("observed_ai.csv", header = FALSE, sep = ",")
exp <- read.table("sim_ai.csv", header = FALSE, sep = ",")

hist(exp[, 1], xlim = c(0, 6), main = "")
abline(v = mean(obs[, 1]), lty = 4, lwd = 2, col = "red")





# Other commands to testing problems within functions
pfile = "../data/parent_genotypes_10.csv"
ofile = "../data/calf_genotypes_10.csv"
ffile = "freqs.csv"
nLoci = 10
nTriads = 3
