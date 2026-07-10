#####################################################
#             allele-inheritance.R                  #
#                                                   #
# Code for testing hypotheses regarding the rate at #
# which offspring inherit maternal and paternal     #
# alleles that are different from each other.       #
#---------------------------------------------------#
# It contains 3 main functions that users can       #
# employ:                                           #
#                                                   #
#   1. frequencies: For calculating allele          #
#      frequencies within a genotype file.          #
#                                                   #
#   2. ai: For calculating allele inheritance       #
#      values for observed offspring.               #
#                                                   #
#   3. sim: For simulating offspring from parental  #
#      genotypes and calculating their ai.          #
#---------------------------------------------------#
# It also contains 2 "helper" functions, that the   #
# user does not interact with directly, but rather  #
# are called upon by other functions. These are:    #
#                                                   #
#   4. mendel: Simulates offspring from known       #
#      parents (used in the sim function).          #
#                                                   #
#   5. aisim: Calculates allele inheritance for the #
#      simulated offspring (used in the sim         #
#      function). Code is the same as in the ai     #
#      function, but some input and output options  #
#      are set to deal with the simulated data.     #
#                                                   #
#                                                   #
#                        by                         #
#                    Tim Frasier                    #
#              Last Updated: 10-JUL-2026            #
#####################################################


######################################################
#               1. frequencies                       #
#                                                    #
# This function calculates the allele frequencies    #
# (for biallelic loci) within a data set. It         #
# assumes two columns for each locus and one row     #
# per individual. The first allele should start      #
# in column 2 (the first column should have          #
# sample/individual labels). One allele should be    #
# labeled as 1, and another as 2. Missing data       #
# should be indicated as 0.                          #
#----------------------------------------------------#
# This function requires two pieces of input:        #
#                                                    #
#   1. file: The name of the file containing the     #
#      genotypes (including the relative path to     #
#      it from the working directory); and           #
#                                                    #
#   2. nLoci: The number of loci used.               #
#                                                    #
# example:                                           #
# frequencies(file = "genotypes.txt", nLoci = 10000) #
######################################################

frequencies <- function(file, nLoci) {

# Load the data data.table library
library(data.table)
  
#-------------------------------#
#   Get Allele Frequencies      #
#-------------------------------#

#--- Table for holding results ---#
freqs <- data.frame(matrix(NA, nrow = 2, ncol = nLoci))

genotypes <- data.frame(fread(file, header = FALSE, sep = ","))

# Loop through loci
for (i in 1:nLoci) {
  
  #--- Print progress to screen (every 500th step) ---#
  if (i %% 500 == 0) {
    print(paste("Calculating allele frequencies for locus: ", i))
  }
  
  #--- Get appropriate columns ---#
  temp <- genotypes[, c(i * 2, (i * 2) + 1)]

  # Count the number of times each allele appears
  one <- sum(temp == 1)
  two <- sum(temp == 2)

  # Write Frequencies
  freqs[1, i] <- one / (one + two)
  freqs[2, i] <- two / (one + two)
}

write.table(freqs, "freqs.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
print("Done! Allele frequencies written to file freqs.csv in working directory.")
}


############################################################
#                        2. ai                             #
#                                                          #
# This function calculates the allele inheritance for each #
# offspring based on whether or not they inherited a       #
# paternal allele different from that inherited from their #
# mother. It only does this for informative loci, and      #
# weights these values based on the expected               #
# heterozygosity at those same loci.                       #
#----------------------------------------------------------#
# It requires 5 pieces of information from the user:       #
#                                                          #
#   1. pfile: The name (and path to) the file containing   #
#      the parental genotypes.                             #
#                                                          #
#   2. ofile: The name (and path to) the file containing   #
#      the offspring genotypes.                            #
#                                                          #
#   3. ffile: The name (and path to) the file containing   #
#      the allele frequencies (generated using the         #
#      frequencies function of this package.               #
#                                                          #
#   4. nLoci: The number of loci used.                     #
#                                                          #
#   5. nTriads: The number of triads being analyzed.       #
#                                                          #
# example:                                                 #
# ai(pfile = "parents.csv", ofile = "offspring.csv",       #
#    ffile = "freqs.csv", nLoci = 5000, nTriads = 3)       #
############################################################

ai <- function(pfile, ofile, ffile, nLoci, nTriads) {
  
# Vector for holding ai values for each offspring
allele_inheritance <- rep(NA, times = nTriads)

# Load the data data.table library
library(data.table)

# Read in data
parents1 <- data.frame(fread(pfile, header = FALSE, sep = ","))
offspring1 <- data.frame(fread(ofile, header = FALSE, sep = ","))

for (i in 1:nTriads) {
  
  # Vector to hold Expected Heterozygosity
  het <- rep(NA, times = nLoci)
  
  # Vector to hold the number of typed loci
  typed <- rep(NA, times = nLoci)
  
  # Vector to hold the number of informative loci
  iLoci <- rep(NA, times = nLoci)
  
  # Vector to hold if inherited alleles differ
  differ <- rep(NA, times = nLoci)
  
  for (j in 1:nLoci) {
    
    #--- Print progress to screen (every 500th step) ---#
    if (j %% 500 == 0) {
      print(paste("Processing triad ", i, "locus ", j))
    }
    
    #---------------------------------------------------------#
    # Get appropriate parental, offspring, and frequency data #
    #---------------------------------------------------------#
    parents <- parents1[((i * 2) - 1):(i * 2), c(j * 2, (j * 2) + 1)]
    offspring <- offspring1[i, c(j * 2, (j * 2) + 1)]
    freqs <- data.frame(fread(ffile, header = FALSE, sep = ",", select = c(j)))

    
    #-----------------------------------#
    # Calculate Expected Heterozygosity #
    #-----------------------------------#
    if ((sum(parents == 0) + sum(offspring == 0)) > 0) {
      het[j] <- 0
      typed[j] <- 0
    } else {
      het[j] <- 1 - (as.numeric(freqs[1, 1])^2 + as.numeric(freqs[2, 1])^2)
      typed[j] <- 1
    }
    
    
    #---------------------------------#
    #  Is Locus Informative?          #
    #---------------------------------#
    if (typed[j] == 0) {
      iLoci[j] <- 0
    } else {
      if (sum(parents) == 5) {
        iLoci[j] <- 1
      } else {
        if ((sum(parents) == 6) && (as.numeric(parents[1, 1]) != as.numeric(parents[1, 2]))) {
          iLoci[j] <- 1
        } else {
          if (sum(parents) == 7) {
            iLoci[j] <- 1
          } else {
            iLoci[j] <- 0
          }
        }
      }
    }
  
    #-----------------------------------#
    #     Do Alleles Differ?            #
    #-----------------------------------#
    if (typed[j] == 0) {
      differ[j] <- 0
    } else {
      if (iLoci[j] == 0) {
        differ[j] <- 0
      } else {
        if(as.numeric(offspring[1, 1]) == as.numeric(offspring[1, 2])) {
        differ[j] <- 0
      } else {
        differ[j] <- 1
      }
      }
    }
  }
  allele_inheritance[i] <- (sum(differ) / sum(iLoci)) / (sum(het) / sum(typed))
}
write.table(allele_inheritance, "observed_ai.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
print("Done! Results written to file observed_ai.csv.")
}

############################################################
#                        3. sim                            #
#                                                          #
# This function runs two other functions (mendel and aisim)#
# to simulate offspring from each mating pair by randomly  #
# selecting one allele from each parent, and then          #
# calculates the average allele inheritance across the     #
# simulated data set. Thus, the output is one value of     #
# allele inheritance (the mean for all offspring) for each #
# iteration.                                               #
#----------------------------------------------------------#
# It requires 5 pieces of information from the user:       #
#                                                          #
#   1. pfile: The name (and path to) the file containing   #
#      the parental genotypes.                             #
#                                                          #
#   2. ffile: The name (and path to) the file containing   #
#      the allele frequencies (generated using the         #
#      frequencies function of this package.               #
#                                                          #
#   3. nLoci: The number of loci used.                     #
#                                                          #
#   4. nTriads: The number of triads being analyzed.       #
#                                                          #
#   5. iterations: The number of iterations to conduct     #
#      (i.e., the number of simulated data sets to create  #
#      and then analyse).                                  #
#                                                          #
############################################################

sim <- function(pfile, ffile, nLoci, nTriads, iterations) {
 
  #-------------------------------#
  #    Load Necessary Packages    #
  #-------------------------------#
  library(data.table)
  
  # Vector for holding results
  sim_ai <- rep(NA, times = iterations)
  
  for (n in 1:iterations) {
    
    # Generate simulated offspring
    mendel(pfile = pfile, nTriads = nTriads, nLoci = nLoci, n = n)
    
    # Calculated their allele-inheritance
    sim_ai[n] <- aisim(pfile = pfile, ffile = ffile, nLoci = nLoci, nTriads = nTriads, n = n)
  }
  write.table(sim_ai, "sim_ai.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
  print("Done! Results written to file sim_ai.csv.")
}  
  
  
#----------------------------------------------------------------------------#
#                             HELPER FUNCTIONS                               #
#                                                                            #
# The functions below are not used directly by the user. Instead, they are   #
# called upon by other functions.                                            #
#----------------------------------------------------------------------------#

############################################
#                 mendel                   #
#                                          #
# This function takes parental genotypes   #
# and generates simulated offspring by     #
# randomly selecting one allele from each  #
# parent. It then calculates the allele-   #
# inheritance of those offspring and       #
# returns the average ai for those         #
# simulated offspring.                     #
#------------------------------------------#
############################################

mendel <- function(pfile, nTriads, nLoci, n) {

# Create file for holding simulated offspring
offspring <- matrix(NA, nrow = nTriads, ncol = ((2 * nLoci) + 1))

# Get parental data
parents1 <- data.frame(fread(pfile, header = FALSE, sep = ","))

for (i in 1:nTriads) {
  
  # Number simulated offspring
  offspring[i, 1] <- i
  
  for (j in 1:nLoci) {
    
    #--- Print progress to screen (every 500th step) ---#
    if (j %% 500 == 0) {
      print(paste("Iteration", n, "Simulating offspring", i, "locus", j))
    }
    
    #-------------------------------#
    # Get appropriate parental data #
    #-------------------------------#
    parents <- parents1[((i * 2) - 1):(i * 2), c(j * 2, (j * 2) + 1)]
    
    #-----------------------#
    #  Generate Offspring   #
    #-----------------------#
    if (sum(parents == 0) > 0) {
      offspring[i, j * 2] <- 0
      offspring[i, (j * 2) + 1] <- 0
    } else {
      offspring[i, (j * 2)] <- sample(as.numeric(parents[1, ]), size = 1, replace = FALSE)
      offspring[i, (j * 2) + 1] <- sample(as.numeric(parents[2, ]), size = 1, replace = FALSE)
    }
  }
}
write.table(offspring, "simOff.csv", col.names = FALSE, row.names = FALSE, sep = ",")
}


############################################################
#                        aisim                             #
#                                                          #
# This is just slighly different than the original ai      #
# function, designed for easier automation for the         #
# analysis of simulated offspring.                         #
############################################################

aisim <- function(pfile, ffile, nLoci, nTriads, n) {
  
  # Vector for holding ai values for each offspring
  allele_inheritance <- rep(NA, times = nTriads)
  
  # Read in offspring file
  offspring1 <- data.frame(fread("simOff.csv", header = FALSE, sep = ","))
  
  # Read in parental file
  parents1 <- data.frame(fread(pfile, header = FALSE, sep = ","))
  
  for (i in 1:nTriads) {
    
    # Vector to hold Expected Heterozygosity
    het <- rep(NA, times = nLoci)
    
    # Vector to hold the number of typed loci
    typed <- rep(NA, times = nLoci)
    
    # Vector to hold the number of informative loci
    iLoci <- rep(NA, times = nLoci)
    
    # Vector to hold if inherited alleles differ
    differ <- rep(NA, times = nLoci)
    
    for (j in 1:nLoci) {
      
      #--- Print progress to screen (every 500th step) ---#
      if (j %% 500 == 0) {
        print(paste("Iteration", n, "Calculating ai for triad", i, "locus", j))
      }
      
      #---------------------------------------------------------#
      # Get appropriate parental, offspring, and frequency data #
      #---------------------------------------------------------#
      parents <- parents1[((i * 2) - 1):(i * 2), c(j * 2, (j * 2) + 1)]
      offspring <- offspring1[i, c(j * 2, (j * 2) + 1)]
      
      freqs <- data.frame(fread(ffile, header = FALSE, sep = ",", select = c(j)))
      
      
      #-----------------------------------#
      # Calculate Expected Heterozygosity #
      #-----------------------------------#
      if ((sum(parents == 0) + sum(offspring == 0)) > 0) {
        het[j] <- 0
        typed[j] <- 0
      } else {
        het[j] <- 1 - (as.numeric(freqs[1, 1])^2 + as.numeric(freqs[2, 1])^2)
        typed[j] <- 1
      }
      
      
      #---------------------------------#
      #  Is Locus Informative?          #
      #---------------------------------#
      if (typed[j] == 0) {
        iLoci[j] <- 0
      } else {
        if (sum(parents) == 5) {
          iLoci[j] <- 1
        } else {
          if ((sum(parents) == 6) && (as.numeric(parents[1, 1]) != as.numeric(parents[1, 2]))) {
            iLoci[j] <- 1
          } else {
            if (sum(parents) == 7) {
              iLoci[j] <- 1
            } else {
              iLoci[j] <- 0
            }
          }
        }
      }
      
      #-----------------------------------#
      #     Do Alleles Differ?            #
      #-----------------------------------#
      if (typed[j] == 0) {
        differ[j] <- 0
      } else {
        if (iLoci[j] == 0) {
          differ[j] <- 0
        } else {
          if(as.numeric(offspring[1, 1]) == as.numeric(offspring[1, 2])) {
            differ[j] <- 0
          } else {
            differ[j] <- 1
          }
        }
      }
    }
    allele_inheritance[i] <- (sum(differ) / sum(iLoci)) / (sum(het) / sum(typed))
  }
  return(mean(allele_inheritance))
}
