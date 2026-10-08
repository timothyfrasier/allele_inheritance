# Code for testing allele_inheritance functions


pfile = "../data/parents_5000.csv"
ofile = "../data/offspring_5000.csv"
ffile = "freqs.csv"
nLoci = 5000
nTriads = 3
iterations = 100


source("allele-inheritance.R")

# Estimate allele frequencies
frequencies(file = "../data/genotypes_5000.csv", nLoci = 5000)

# Calculate ai for observed offspring
ai(pfile = "../data/parents_5000.csv", ofile = "../data/offspring_5000.csv", ffile = "freqs.csv", nLoci = 5000, nTriads = 3)

# Simulate offspring and calculate their a1
sim(pfile = "../data/parents_5000.csv", ffile = "freqs.csv", nLoci = 5000, nTriads = 3, iterations = 100)

# Visualize results
obs <- read.table("observed_ai.csv", header = FALSE, sep = ",")
exp <- read.table("sim_ai.csv", header = FALSE, sep = ",")

library(ggplot2)    
ggplot(exp) +    
   theme_bw() +    
   geom_histogram(aes(x = V1), alpha = 0.6) +    
   geom_vline(xintercept = mean(obs$V1), color = "red", linewidth = 1.5, linetype = "dashed") +    
   xlab("Allele Inheritance") +    
   ylab("Frequency")  


# How Heterozygous Must Calves Be To Survive?
# crit_ai = 2.5
sim_loss(pfile = "../data/parents_5000.csv", ffile = "freqs.csv", nLoci = 5000, nTriads = 3, crit_ai = 2.5, iterations = 100)

# Visualize results
obs <- read.table("observed_ai.csv", header = FALSE, sep = ",")
exp <- read.table("sim_ai_loss.csv", header = FALSE, sep = ",")

library(ggplot2)    
ggplot(exp) +    
  theme_bw() +    
  geom_histogram(aes(x = V1), alpha = 0.6) +    
  geom_vline(xintercept = mean(obs$V1), color = "red", linewidth = 1.5, linetype = "dashed") +    
  xlab("Allele Inheritance") +    
  ylab("Frequency")  
