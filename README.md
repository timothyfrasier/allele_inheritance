# INSTRUCTIONS FOR USING THE ALLELE INHERITENCE R SCRIPTS

These scripts are an R adaptation of the C code described in Frasier (2008), which is implemented in the STORM program for testing hypotheses regarding allele inheritance. Specifically, they allow users to test the hypothesis that offspring inherit paternal alleles differing from the maternal alleles more/less often than expected given the parents genotypes and Mendelian inheritance.    

There are two versions of the scripts: one that reads entire data files into R, and another that only reads part of the data files into R at a time. There are two example data sets, each containing three files: a file containing offspring genotypes, a file containing the parental genotypes, and a file containing all genotypes combined for use in calculating allele frequencies.    

These instructions assume that the files are laid out as follows:
    - The data files are in a "data" directory,
    - The scripts are in a "code" directory, and
    - You have set R's working directory to the
      "code" directory    

|--| code/    
|-----|- allele-inheritance_all-in.R    
|-----|- allele-inheritance_partial.R    
|    
|--| data/    
|-----|- genotypes_full.csv    
|-----|- genotypes_5000.csv    
|-----|- offspring_5000.csv    
|-----|- offspring_full.csv    
|-----|- parents_5000.csv    
|-----|- parents_full.csv    

Instructions and notes are also present within the code itself that may be helpful as a guide.    

------

## File Descriptions:     

    1. allele-inheritance_all_in.R
        - Contains functions for conducting analyses reading in the entire data files into R.

    2. allele-inheritance_partial.R
        - Contains functions for conducting analyses reading in just part of the file(s) at a time.

    3. offspring_5000.csv
        - Data for 3 offspring genotyped at 5000 SNPs

    4. offspring_full.csv
        - Data for 3 offspring genotyped at 57,704 SNPs

    5. parents_5000.csv
        - Genotypes for the 3 parental pairs at 5000 SNPs
        - Parents pairs must be in the same order as the offspring. For example, rows 1 and 2 must be the genotypes of the parents of the offspring in row 1 of the offspring file. Rows 3 and 4 must be the genotypes of the parents of the offspring in row 2 of the offspring file. And so on.

    6. parents_full.csv
        - Genotypes for the 3 parental pairs at 57,704 SNPs
        - Parents pairs must be in the same order as the offspring. For example, rows 1 and 2 must be the genotypes of the parents of the offspring in row 1 of the offspring file. Rows 3 and 4 must be the genotypes of the parents of the offspring in row 2 of the offspring file. And so on.

    7. genotypes_5000.csv
        - Combined offspring and parental genotypes at 5000 SNPs

    8. genotypes_full.csv
        - Combined offspring and parental genotypes at 57,704 SNPs

-----

## Running analyses reading full data files into R (example using the 5000 SNP data set, but can also use with larger data set).    

1. Open RStudio and make the "code" directory R's working directory    

2. Load the script into R, so that the functions are available to you    
`source("allele-inheritance_all-in.R")`

3. Make sure that you have the `data.table` package installed    

4. Use the `frequencies` function to calculate the allele frequencies. This will write the results to a file called `freqs.csv` in R's working directory. Requires two arguments:    
    a. The name and relative path to the genotypes file ("file"). File *must* be comma-delimited    
    b. The number of loci ("nLoci")    
`frequencies(file = "../data/genotypes_5000.csv", nLoci = 5000)`

5. Use the `ai` function ("ai" for allele inheritance, *not* artificial intelligence) to calculate observed allele inheritance. This will write the results (a single column with one value for each triad) into a file called `observed_ai.csv` within R's working directory. Requires five arguments:    
    a. The name and relative path to the parental genotype file ("pfile"). File *must* be comma-delimited    
    b. The name and relative path to the offspring genotype file ("ofile"). File *must* be comma-delimited    
    c. The name and relative path to the allele frequency file ("ffile")    
    d. The number of loci ("nLoci")    
    e. The number of mother-father-triads ("nTriads")    
`ai(pfile = "../data/parents_5000.csv", ofile = "../data/offspring_5000.csv", ffile = "freqs.csv", nLoci = 5000, nTriads = 3)`

6. Use the `sim` function to generate simulated offspring from each parental pair. **This will take a while to run!** Progress updates will be written to your screen so that you can keep track of progress. But you will need to be patient and have other things to do in the meantime (or perhaps go to sleep and check it in the morning). The output will be a file called `sim_ai.csv` that contains *one value* for each iteration (the average ai across *all* simulated  offspring for that iteration). The result will be a distribution of expected allele inheritance values given just Mendelian inheritance. Requires five arguments:    
    a. The name and relative path to the parental genotype file ("pfile"). File *must* be comma-delimited    
    b. The name and relative path to the allele frequency file ("ffile")    
    c. The number of loci ("nLoci")    
    d. The number of triads ("nTriads")    
    e. The number of iterations that you want to conduct (i.e., how many simulated offspring do you want to generate for each parental pair?). I just use 20 here, but to get a good distribution of "expected" values, you should run more (e.g., 500 or 1000).    
`sim(pfile = "../data/parents_5000.csv", ffile = "freqs.csv", nLoci = 5000, nTriads = 3, iterations = 10)`

7. You can then compare and visualize the difference between observed and expected values using standard R functions. One example is to plot a histogram of the expected values, and then a red dashed line where the observed value is. Example commands for this are below. *Note that the calculation is how often a paternal allele is inherited that is different from the maternal allele, so observed values larger than expected mean that homozygotes are missing from your data set, and vice versa.*     
`observed <- read.table("observed_ai.csv", header = FALSE, sep = ",")`    
`expected <- read.table("sim_ai.csv", header = FALSE, sep = ",")`    
`library(ggplot2)`    
`ggplot(expected) +`    
`   theme_bw() +`    
`   geom_histogram(aes(x = V1), alpha = 0.6) +`    
`   geom_vline(xintercept = mean(observed$V1), color = "red", linewidth = 1.5, linetype = "dashed") +`    
`   xlab("Allele Inheritance") +`    
`   ylab("Frequency")`    
-----

## Running analyses reading only a portion of the data files into R at a time (example using the 5000 SNP data set, but can also use with larger data set).

The commands for this are the same as above, except for which script you read in.        
`source("allele-inheritance_partial.R")`

## References
Frasier TR (2008) STORM: Software for testing hypotheses of relatedness and mating patterns. *Molecular Ecology Resources* **8**: 1264-1266.
