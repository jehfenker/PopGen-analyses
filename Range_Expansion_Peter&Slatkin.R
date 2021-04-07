####      Jessica Fenker 2019        ####
#### Australian National University  ####
####    jehfenker@gmail.com          ####


## Range expansion analysis, to infer the origin and to infer the strength of the founder effect
## Choose the models that you would like to test
## the outgroup is usually the sister taxa
## Details of the methods are found in Peter & Slatkin (2013), Evolution and Peter & Slatkin (2015)
## https://github.com/BenjaminPeter/rangeexpansion  

source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
library(snpStats)
devtools::install_github("BenjaminPeter/rangeExpansion", ref="package")

library(rangeExpansion)


## First step is to modifi your SNAPP file: without headings and with comma as a separator

#open this new file, including the OUTGROUP
snp.file_bilin_out <- "bilin_with_outgroup.csv" 
coord.file_bilin_out <- "bilin_out_coords.csv" #file with individual id (same as SNP file), coords, outgroup and pop as optional
ploidy <- 2 #diploid individuals


#Regions are sets of populations that can be analyzed independently, i.e. they correspond
#to clusters that a priory are thought to have a different origin.
#Each list entry corresponds to an analysis, i. e. the following command
#will analyze the populations `REGION_1`, `REGION_2` and `REGION_3` individually, but will
#also jointly analyze `REGION_1` and `REGION_2`. NULL corresponds to all
#individuals analyzed

region.bilin_out<- list(NULL)

raw.data_bilin_out<- load.data.snapp(snp.file_bilin_out, coord.file_bilin_out, sep=',', na.strings = "?", ploidy=ploidy)


pop_bilin_out <- make.pop(raw.data_bilin_out, ploidy)
psi_bilin_out <- get.all.psi(pop_bilin_out)
write.csv(psi_bilin_out,file="psi_table_bilin_withOUT.csv")


res.bilin_out <- run.regions(region = region.bilin_out, pop=pop_bilin_out, psi=psi_bilin_out)

summary(res.bilin_out)
quartz()
plot(res.bilin_out)
pdf("pop_plots_bilin_withOUT.pdf",width = 6, height = 6) 
par(mfrow=c(1,1))
