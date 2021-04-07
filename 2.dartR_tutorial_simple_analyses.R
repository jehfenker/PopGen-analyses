####      Jessica Fenker 2020        ####
#### Australian National University  ####
####    jehfenker@gmail.com          ####


## This script is a continuation of "1.dartR_tutorial_basics.R", and
## include some basic popgen analyses


# First, open the saved project,loading you data and calling all the packages
library(SNPRelate)
library(dartR)
library(adegenet)
library(ggplot2)
library(data.table)
library(irelr)
library(strataG)

###########################
#### BASIC STATISTICS #####
###########################

#Example with just one species
#visualize your data
gl.map.interactive(bilin_ind, provider = "Esri.NatGeoWorldMap")

#Basic stats (Hs, Ho, Fs, etc) 
bilin_basic<-gl.basic.stats(bilin_ind)
bilin_basic

#There are many "gl.report" functions in the package. Here is an example
bilin_heteroz<-gl.report.heterozygosity(bilin_ind)
bilin_heteroz


#For Tajima's D, use this function, updated by Renee Catullo
get_tajima_D <- function(x){
  require(dartR) # possibly not needed for a function in an R package?
  
  # Find allele frequencies (p1 and p2) for every locus in every population
  allele_freqs <- gl.percent.freq(x)
  names(allele_freqs)[names(allele_freqs) == "frequency"] <- "p1"
  allele_freqs$p1 <- allele_freqs$p1 / 100
  allele_freqs$p2 <- 1 - allele_freqs$p1
  
  # Get the names of all the populations
  pops <- unique(allele_freqs$popn)
  
  #split each population
  allele_freqs_by_pop <- split(allele_freqs, allele_freqs$popn)
  
  # Internal function to calculate pi
  calc_pi <- function(allele_freqs) {
    n = allele_freqs$nobs * 2  # vector of n values
    pi_sqr <- allele_freqs$p1 ^ 2 + allele_freqs$p2 ^ 2
    h = (n / (n - 1)) * (1 - pi_sqr) # vector of values of h
    sum(h) # return pi, which is the sum of h across loci
  }
  
  get_tajima_D_for_one_pop <- function(allele_freqs_by_pop) {
    pi <- calc_pi(allele_freqs_by_pop)
    
    #Calculate number of segregating sites, ignoring missing data (missing data will not appear in teh allele freq calcualtions)
    #S <- sum(!(allele_freqs_by_pop$p1 == 0 | allele_freqs_by_pop$p1 == 1))
    S <- sum(allele_freqs_by_pop$p1 >0 & allele_freqs_by_pop$p1 <1)
    if(S == 0) {
      warning("No segregating sites")
      data.frame(pi = NaN, 
                 S = NaN, 
                 D = NaN, 
                 Pval.normal = NaN, 
                 Pval.beta = NaN)
    }
    
    n <- mean(allele_freqs_by_pop$nobs * 2 )
    
    tmp <- 1:(n - 1)
    a1 <- sum(1/tmp)
    a2 <- sum(1/tmp^2)
    b1 <- (n + 1)/(3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
    c1 <- b1 - 1/a1
    c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    
    
    #calculate D and do beta testing
    D <- (pi - S/a1) / sqrt(e1 * S + e2 * S * (S - 1))
    Dmin <- (2/n - 1/a1)/sqrt(e2)
    Dmax <- ((n/(2*(n - 1))) - 1/a1)/sqrt(e2)
    tmp1 <- 1 + Dmin * Dmax
    tmp2 <- Dmax - Dmin
    a <- -tmp1 * Dmax/tmp2
    b <- tmp1 * Dmin/tmp2
    p <- pbeta((D - Dmin)/tmp2, b, a)
    p <- ifelse(p < 0.5, 2 * p, 2 * (1 - p))
    
    data.frame(pi = pi, 
               S = S, 
               D = D, 
               Pval.normal = 2 * pnorm(-abs(D)), 
               Pval.beta = p)
  }
  
  output <- do.call("rbind", lapply(allele_freqs_by_pop, 
                                    get_tajima_D_for_one_pop))
  data.frame(population = rownames(output), output, row.names = NULL)
}
​
#note that for Tajima's D it is important to have a good number of individuals​
bilin_tajima<-get_tajima_D(bilin.filt)
bilin_tajima


###############
## You can also obtain some statistics (including Nei's Fst) using the packages
## hierfstat and tidyverse

library(hierfstat)
library(tidyverse)

#convert genlight to genind
gi_bilin <- gl2gi(gl = bilin.filt)
gi_bilin_ind<-gi_bilin
pop(gi_bilin_ind)<-indNames(gi_bilin_ind)

#convert genind to hierfstat object
hfs_bilin <- genind2hierfstat(gi_bilin)
hfs_bilin_ind <- genind2hierfstat(gi_bilin_ind)
pop.vector <- levels(as.factor(gi_bilin$pop))
​
basic_stats <- basic.stats(hfs_bilin_ind)
​
Ho <- basic_stats$Ho %>%
  as.data.frame() %>%
  dplyr::summarise_all(list(mean)) %>% t() %>%
  as.data.frame() %>%
  dplyr::rename(Ho = "V1") %>%
  rownames_to_column(var = "Populations") %>%
  dplyr::mutate(Ho = round(Ho, 2))
​
He <- basic_stats$Hs %>%
  as.data.frame() %>%
  dplyr::summarise_all(list(mean), na.rm = TRUE) %>% t() %>%
  as.data.frame() %>%
  dplyr::rename(He = "V1") %>%
  rownames_to_column(var = "Populations") %>%
  dplyr::mutate(He = round(He, 2))
​
Fis <- basic_stats$Fis %>%
  as.data.frame() %>%
  dplyr::summarise_all(list(mean), na.rm = TRUE) %>% t() %>%
  as.data.frame() %>%
  dplyr::rename(Fis = "V1") %>%
  rownames_to_column(var = "Populations") %>%
  dplyr::mutate(Fis = round(Fis, 2)) %>% arrange(Fis)
​
stats_pop <- Ar %>% left_join(Ho, by = "Populations") %>%
  left_join(He, by = "Populations") %>%
  left_join(Fis, by = "Populations")
​
FST_bilin <- pairwise.fst(gi_bilin_ind, pop = NULL, res.type = c("dist", "matrix"))
GD_Cp_bilin <- genet.dist(hfs_bilin_ind, method = "Cp")
GD_Nei_bilin <- genet.dist(hfs_bilin_ind, method = "Da")
GD_TnN_bilin <- genet.dist(hfs_bilin_ind, method = "Dch")
​
Fst_bilin<-as.matrix(FST_bilin)
write.csv(Fst_bilin,"bilin_fst.csv")



#############################
#### SPATIAL STATISTICS #####
#############################

#Euclidean distance is the geographic distance between two points
Edis_bilin <- as.matrix(dist(bilin_ind@other$latlong))
dim(Edis_bilin)
write.csv(Edis_bilin,"bilin_Edis2.csv")


#There are numerous ways to measure genetic distance
#one is the proportion of shared alleles in our dataset
bilin_gi<-gl2gi(bilin_ind)
Gdis_bilin <- 1 - propShared(bilin_gi)


#to visualise the matrix, you can use
image(Gdis_bilin)
table.value(Gdis_bilin, csize = 0.4)
write.csv(Gdis_bilin,"bilin_Gdis2.csv")

#an example of a IBD graph
ibd_bilin<-gl.ibd(bilin_ind)


#Done!
​