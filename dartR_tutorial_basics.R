####      Jessica Fenker 2020        ####
#### Australian National University  ####
####    jehfenker@gmail.com          ####


## This is a simple and initial script to work with the dartR
## package, performing analyses with SNP data, specially the 
## ones obtained from Diversity Arrays Tech (https://www.diversityarrays.com/)
## Note that they frequently update the package, so some functions 
## can be a bit different and there could be new cool analysis available
## be sure that you have the last version installed


#install necessary packages
install.packages("devtools")
library(devtools)
install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "qvalue"))
packages_to_install <-c("dartR","adegenet","ggplot2", "data.table")
install.packages(packages_to_install)
install_github("andersgs/irelr")
install_github('ericarcher/strataG')

#load libraries
library(SNPRelate)
library(dartR)
library(adegenet)
library(ggplot2)
library(data.table)
library(irelr)
library(strataG)


#set working directory
setwd("~/jessica/SNP/Dipes")


######
## SILICO DART
#First thing first. Look at your silicoDArT data (presence absence of snps to determine
#major groups) to check for the quality of your data. One idea is to sum up the total 
#number of SNPs for each individual and compare them all

################################
######## LOAD YOUR DATA ########
################################
#Call the SNP csv document and an metafile that contains additional information on
#individuals. In this case, it jut contains the population that each individual belongs
#to and geographic coordinates. The name of the ind most be the same in both tables

dipes<-gl.read.dart("Dipes.csv",ind.metafile ="Dipes_metadata.csv",probar=TRUE)
dipes #check if your latlong information uploaded correctly 

#visualize in the map where the population is
gl.map.interactive(dipes, provider = "Esri.NatGeoWorldMap")

#Report the number of loci, individuals and populations
nLoc(dipes)
nPop(dipes)
table(pop(dipes))
names (dipes@other$loc.metrics)


################################
###### INVESTIGATE/REPORT ######
################################

#If you want, generate a smear plot of genotypes across samples and loci 
#shows the loci that were not genotyped for a particular individual (white pixels)
glPlot(dipes)

par(mfcol = c(1,2))
hist(dipes$other$loc.metrics$CallRate, xlab = "Call Rate", main = NULL)
hist(dipes$other$loc.metrics$RepAvg, xlab = "Repeatability", main = NULL)

#smear plot of the genotypes across samples and loci can be an informative starting point
#shows the loci that were not genotyped for a particular individual (white pixels)
par(mfcol = c(1,1))
glPlot(dipes)

#check all the gl.report functions on the dartR package
#some examples:
#gives no. loci with percentages of missing values, which the call rate exceeds a range of thresholds
gl.report.callrate(dipes)

#gives no. loci which the reproducibility (averaged over the two allelic states) exceeds a range of thresholds
gl.report.RepAvg(dipes)

#SNPs with a low read depth may be more likley to contain sequence errors
#The SNPs identified in these loci may not be real but an artefact of different but similar genes being grouped together (paralogs)
gl.report.rdepth(dipes)



################################
######## LOAD YOUR DATA ########
################################

# dartR has some good default locus filtering options. Play with your data, test different filterings to see how they work.
# can be a good idea to name your filtered dataset something new so you can backtrack if needed
# check the package CRAN for extra filtering


#Call Rate 95% -to filter out inds that sequenced poorly (amount of missing data per individual)
#note that you can also filter the call rate by individual (method="ind")
dipes.filt1 <- gl.filter.callrate(dipes, method = "loc", threshold = 0.95, v = 3)
#if you don't want to use the call rate PER INDIVIDUAL, you can still remove a specific ind from your analysis
dipes.filt2 <- gl.drop.ind(dipes.filt1, "species_name",recalc = TRUE, mono.rm = FALSE)

#Drop loci with minor allele freq  <2% - threshold depend on the number of individuals you have
#a general rule of 1/2n is applied
dipes.filt3<- gl.filter.maf(dipes.filt2, threshold = 0.003)

#RepAvg measures repeatability - values towards 1 is high reproducability 
#Create an index for repeatability to filter data by (a meassurement of quality per loci)
dipes.filt4<-gl.filter.RepAvg(dipes.filt3,t=0.99)

#Read Depth 5x - filterloci with less than 5x read depth (5 overlapping sequences to assure its a good locus)
dipes.filt5 <- gl.filter.rdepth(dipes.filt4, lower=5, upper=100)

#Remove monomorphic SNPs (as they do not provide information for population structure)
#avoid to remove if you want to do anaylsis(Fst and related) that are related to the heterozigosity
dipes.filt6<-gl.filter.monomorphs(dipes.filt5, v=0)

#Secondaries - get only one SNP per locus, removing duplicated SNPs in the same fragment
dipes.filt7 <- gl.filter.secondaries(dipes.filt6)

#Recalculate all the metrics and report now, after filtering, the number of loci, individuals and populations
#Also, check gl.report functions for each filter if wanted
dipes.filt <- gl.recalc.metrics(dipes.filt7, verbose=3)
dipes.filt
nLoc(dipes.filt)
nInd(dipes.filt)
levels(pop(dipes.filt))

#compare the results
par(mfcol = c(1,2))
glPlot(dipes)
glPlot(dipes.filt)
par(mfcol = c(1,1))

################################
#### VISUALIZING YOUR DATA #####
################################

#### PCOA
#  PCoA is a simple way to visualize the structure of your data
pcoa <- gl.pcoa(dipes.filt)
gl.pcoa.scree(pcoa)
gl.pcoa.plot(pcoa, dipes.filt)
gl.pcoa.plot(pcoa, dipes.filt,labels="ind")
gl.pcoa.plot(pcoa, dipes.filt,xaxis=1,yaxis=3) #Axis 1 and 3

#3d plot is nice too
#The orientation of the plot can be move by clicking and draging the cursor
gl.pcoa.plot.3d(pcoa, nanagr.filt, title = "PCoA", xaxis = 1, yaxis = 2, zaxis = 3, shape = "sphere", radius = 2, legend = "topright")

#I usually visualize the PCoA here and then plot the results in ggplot for a nicer visualization
#including the name of all the points in the PCoA
pcoa_ggplot <- cbind(pcoa$scores[,1], pcoa$scores[,2], pcoa$scores[,3])
colnames(pcoa_ggplot) <- c("PC1", "PC2", "PC3")
pcoa_ggplot <- as.data.frame(pcoa_ggplot)
names <- dipes.filt@pop
ggplot(pcoa_ggplot, aes(x = PC1, y = PC2)) +  #PC1 X PC2
  geom_point(aes(PC1, PC2), size = 5, color = 'grey') +
  geom_label_repel(aes(PC1, PC2, fill = factor(names), label = rownames(pcoa_ggplot)),
                   fontface = 'bold', color = 'white',
                   box.padding = unit(1, "lines"),
                   point.padding = unit(1, "lines"),
                   segment.color = 'grey50',size=3 ) +
  labs(x = "PCoA Axis 1 (26.8 %)", y = "PCoA Axis 2 (14.2 %)") + ## change the values MANUALY here!!!
  scale_fill_discrete(name = "Populations") +
  theme_classic(base_size = 16)


##############################
#### EXPORTING YOUR DATA #####
##############################

#Export data to run phylogenetic trees
#Before that, define that you want it at the individual level, not pop
dipes.filt_ind<-dipes.filt
pop(dipes.filt_ind)<-indNames(dipes.filt_ind)


#save your SNP data in a format that can be read on SVD-QUARTETS to run the tree
gl2svdquartets(dipes.filt_ind, outfile = "dipes_svd.nex",
  method = 2,
  verbose = 3,
  outpath = getwd())

#For IQ-TREE, first you need to save your SNP data in fasta
gl2fasta(dipes.filt_ind,outfile="dipes.fasta",verbose=3, outpath = getwd())

#then, solve ambiguity codes
#this part of the script was provided from Carlos Pavon
library(ape)
library(seqinr)
​seq.tmp <- read.fasta("dipes.fasta",seqtype="DNA",forceDNAtolower=F)
​
for (i in 1:length(seq.tmp)){
  seq.tmp[[i]] <- seq.tmp[[i]][-which(seq.tmp[[i]]==" ")]
}
amb.codes <- list(M=c("A","C"),R=c("A","G"),W=c("A","T"),S=c("C","G"),Y=c("C","T"),K=c("G","T"))
for (i in 1:length(seq.tmp)){
  seq.trash <- seq.tmp[[i]]
  rep.trash <- which(seq.trash%in%names(amb.codes))
  for (j in rep.trash){
    seq.trash[j] <- sample(amb.codes[[which(names(amb.codes)==seq.trash[j])]],1)
  }
  seq.tmp[[i]] <- seq.trash
}
del.pos <- c()
for (i in 1:length(seq.tmp[[1]])){
  pos.trash <- unlist(map(seq.tmp,i))
  length.trash <- length(unique(pos.trash))
  if(length.trash==1){
    del.pos <- c(del.pos,i)
  } else {
    if(length.trash==2&"N"%in%pos.trash){
      del.pos <- c(del.pos,i)
    }
  }
}
for(i in 1:length(seq.tmp)){
  seq.tmp[[i]] <- seq.tmp[[i]][-del.pos]
}
​
write.fasta(seq.tmp,names=names(seq.tmp),nbchar=length(seq.tmp[[1]]),file.out="dipes_solved.fasta")

#Now, use this fasta doc to run your tree



############################
#### SPLITTING SPECIES #####
############################

#Kind of a repetition of the steps above, but per species

table(pop(dipes))


#there are a few commands to select samples that can be useful,as keep or skip a population 
#or even drop individuals. A few examples:
bilin<-gl.keep.pop(dipes, pop.list=c("bilin"), recalc = TRUE, mono.rm = FALSE)
magna<-gl.keep.pop(dipes,pop.list=c("magnaETE","magnaWTE","magnaKIM"), recalc = TRUE, mono.rm = FALSE)
magna2<-gl.drop.ind(magna,ind.list=c("magnaETE_CCM0001","magnaETE_R0002"), recalc = TRUE, mono.rm = FALSE)
dipes_not_perp<-gl.drop.pop(dipes,pop.list=c("perplexa"), recalc = TRUE, mono.rm = FALSE)

# there is also gl.edit.recode.ind and gl.edit.recode.pop, where you can edit and reassign
# individuals or population, using a recode table. Check website if necessary.
# as a genlight object, you have the full capabilities of the adegenet package at your fingertips 


###########
## FILTERING, now per species
#using just one species as example
gl.map.interactive(bilin, provider = "Esri.NatGeoWorldMap")

bilin.filt1 <- gl.filter.callrate(bilin, method='ind', threshold=0.4,recalc = T)
bilin.filt2 <- gl.filter.callrate(bilin.filt1, threshold = 0.9,recalc = T)
bilin.filt4 <- gl.filter.RepAvg(bilin.filt3, threshold = 0.995)
bilin.filt6 <- gl.filter.monomorphs(bilin.filt5, verbose=3)
bilin.filt <- gl.recalc.metrics(bilin.filt7, verbose=3)
bilin.filt #final file

##########
##Plot a Principal Component Analysis
pcoa_bilin <- gl.pcoa (bilin.filt)
gl.pcoa.plot(pcoa_bilin, bilin.filt) #2d
#3d plot - this will cause R to open another window to show the plot.
#The orientation of the plot can be move by clicking and draging the cursor.
gl.pcoa.plot.3d(pcoa_bilin,bilin.filt)


### ggplot - PCoA
pcoa_bilin2 <- cbind(pcoa_bilin$scores[,1], pcoa_bilin$scores[,2], pcoa_bilin$scores[,3])
colnames(pcoa_bilin2) <- c("PC1", "PC2", "PC3")
pcoa_bilin2 <- as.data.frame(pcoa_bilin2)
names <- bilin.filt@pop
#PC1 X PC2
ggplot(pcoa_bilin2, aes(x = PC1, y = PC2)) + #PC1 and PC2
  geom_point(aes(PC1, PC2), size = 5, color = 'grey') +
  geom_label_repel(aes(PC1, PC2, fill = factor(names), label = rownames(pcoa6.s)),
                   fontface = 'bold', color = 'white',
                   box.padding = unit(1, "lines"),
                   point.padding = unit(1, "lines"),
                   segment.color = 'grey50',size=3 ) +
  labs(x = "PCoA Axis 1 (26.8 %)", y = "PCoA Axis 2 (14.2 %)") + ## change the values MANUALY here!!!
  scale_fill_discrete(name = "Populations") +
  theme_classic(base_size = 16)



###########################
#### SAVING YOUR DATA #####
###########################

# Save your data script, command and objects to access it later
# A good practice is to always create a project as save it when finished
# use ?save for help

save.image("~/Desktop/dipes_SNP_data.RData")
saveRDS(object = dipes_SNP, file = "dipes_SNP_analyses.rds") #for single objects

#Done!