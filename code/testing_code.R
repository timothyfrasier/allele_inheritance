# Code for testing allele_inheritance functions

source("allele-inheritance_all-in.R")

# Estimate allele frequencies
frequencies(file = "../data/genotypes_20.csv", nLoci = 20)

# Calculate ai for observed offspring
ai(pfile = "../data/parents_20.csv", ofile = "../data/offspring_20.csv", ffile = "freqs.csv", nLoci = 20, nTriads = 3)

# Simulate offspring and calculate their ai
sim(pfile = "../data/parents_20.csv", ffile = "freqs.csv", nLoci = 20, nTriads = 3, iterations = 500)

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

hist(exp[, 1], xlim = c(0, 6), main = "")
abline(v = mean(obs[, 1]), lty = 4, lwd = 2, col = "red")

# How Heterozygous Must Calves Be To Survive?
# H = 0.1
sim_loss(pfile = "../data/parent_genotypes_20.csv", ffile = "freqs.csv", nLoci = 20, nTriads = 3, H = 0.1, iterations = 500)
exp1 <- read.table("sim_ai_loss.csv", header = FALSE, sep = ",")

`library(ggplot2)`    
`ggplot(expected) +`    
`   theme_bw() +`    
`   geom_histogram(aes(x = V1), alpha = 0.6) +`    
`   geom_vline(xintercept = mean(observed$V1), color = "red", linewidth = 1.5, linetype = "dashed") +`    
`   xlab("Allele Inheritance") +`    
`   ylab("Frequency")`  


hist(exp1[, 1], xlim = c(0, 6), main = "")
abline(v = mean(obs[, 1]), lty = 4, lwd = 2, col = "red")


# H = 0.2
sim_loss(pfile = "../data/parent_genotypes_20.csv", ffile = "freqs.csv", nLoci = 20, nTriads = 3, H = 0.2, iterations = 500)
exp2 <- read.table("sim_ai_loss.csv", header = FALSE, sep = ",")

hist(exp2[, 1], xlim = c(0, 6), main = "")
abline(v = mean(obs[, 1]), lty = 4, lwd = 2, col = "red")


# H = 0.3
sim_loss(pfile = "../data/parent_genotypes_20.csv", ffile = "freqs.csv", nLoci = 20, nTriads = 3, H = 0.3, iterations = 500)
exp <- read.table("sim_ai_loss.csv", header = FALSE, sep = ",")

hist(exp[, 1], xlim = c(0, 6), main = "")
abline(v = mean(obs[, 1]), lty = 4, lwd = 2, col = "red")


