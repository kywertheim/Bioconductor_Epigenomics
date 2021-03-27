#Load the libraries needed for the following tasks.
library(AnnotationHub)
library(GenomicRanges)
library(rtracklayer)

#Create an annotation hub.
ah <- AnnotationHub()

#Question 1.
#Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.
#How many islands exists on the autosomes?

#Search for data on CpG islands in the human genome.
ah_human <- subset(ah, species == 'Homo sapiens')
ah_humanCpG <- query(ah_human, 'CpG Islands')

#Download the data associated with one of the available human genomes.
ah_humanCpG$genome
ah_humanCpG[1]
ah_hg19CpG <- ah_humanCpG[['AH5086']]

#Check that there are many chromosomes in this dataset.
seqlevels(ah_hg19CpG)

#Split the dataset into groups defined by chromosome.
ah_hg19CpG_split <- split(ah_hg19CpG, seqnames(ah_hg19CpG))

#Create a filter for the autosomes.
autosome <- paste('chr', 1:22, sep='')

#Apply the filter to subset the CpG islands lying on the autosomes.
autosome_CpG <- unlist(ah_hg19CpG_split[autosome])

#Question 2.
#How many CpG Islands exists on chromosome 4?

autosome_CpG[4]

#Question 3.
#Obtain the data for the H3K4me3 histone modification for the H1 cell line from Epigenomics Roadmap, using AnnotationHub. Subset these regions to only keep regions mapped to the autosomes (chromosomes 1 to 22).
#How many bases does these regions cover?

#Search for the Epigenomics Roadmap data on H3K4me3 histone modification for the H1 cell line.
ah_h1H3K4me3 <- query(ah, c('H3K4me3', 'H1 cell', 'roadmap'))

#Download the narrow-peak data.
ah_h1H3K4me3_narrow <- ah_h1H3K4me3[['AH29884']]

#Check that there are many chromosomes in this dataset.
seqlevels(ah_h1H3K4me3_narrow)

#Split the dataset into groups defined by chromosome.
ah_h1H3K4me3_narrow_split <- split(ah_h1H3K4me3_narrow, seqnames(ah_h1H3K4me3_narrow))

#Apply the filter to subset the regions lying on the autosomes.
autosome_h1H3K4me3 <- unlist(ah_h1H3K4me3_narrow_split[autosome])

#Count the number of bases in these regions.
sum(width(autosome_h1H3K4me3))

#Question 4.
#Obtain the data for the H3K27me3 histone modification for the H1 cell line from Epigenomics Roadmap, using the AnnotationHub package. Subset these regions to only keep regions mapped to the autosomes. In the return data, each region has an associated "signalValue".
#What is the mean signalValue across all regions on the standard chromosomes?

#Search for the Epigenomics Roadmap data on H3K27me3 histone modification for the H1 cell line.
ah_h1H3K27me3 <- query(ah, c('H3K27me3', 'H1 cell', 'roadmap'))

#Download the narrow-peak data.
ah_h1H3K27me3_narrow <- ah_h1H3K27me3[['AH29892']]

#Check that there are many chromosomes in this dataset.
seqlevels(ah_h1H3K27me3_narrow)

#Split the dataset into groups defined by chromosome.
ah_h1H3K27me3_narrow_split <- split(ah_h1H3K27me3_narrow, seqnames(ah_h1H3K27me3_narrow))

#Apply the filter to subset the regions lying on the autosomes
autosome_h1H3K27me3 <- unlist(ah_h1H3K27me3_narrow_split[autosome])

#Average the signal values of these regions.
mean(autosome_h1H3K27me3$signalValue)

#Question 5.
#Bivalent regions are bound by both H3K4me3 and H3K27me3.
#Using the regions we have obtained above, how many bases on the standard chromosomes are bivalently marked?

#Find the bivalent regions by intersecting the two sets of regions.
bivalent <- intersect(autosome_h1H3K4me3, autosome_h1H3K27me3)

#Count the number of bases in the intersection.
sum(width(bivalent))

#Question 6.
#We will examine the extent to which bivalent regions overlap CpG Islands.
#How big a fraction (expressed as a number between 0 and 1) of the bivalent regions, overlap one or more CpG Islands?

#Use the bivalent regions as the query and the CpG islands as the subject.
bivalentONCpG <- findOverlaps(bivalent, autosome_CpG)

#Find the unique set of query hits (a bivalent region may overlap with multiple CpG islands) and divide its size by the total number of bivalent regions.
length(unique(queryHits(bivalentONCpG)))/length(bivalent)

#Question 7.
#How big a fraction (expressed as a number between 0 and 1) of the bases which are part of CpG Islands, are also bivalent marked?

#Find the regions containing the bivalently marked CpG islands.
bivalentINTERSECTCpG <- intersect(bivalent, autosome_CpG)

#After reducing the ranges to prevent counting a base more than once, count the number of bases in these ranges and divide it by the number of bases in the CpG islands.
sum(width(reduce(bivalentINTERSECTCpG)))/sum(width(autosome_CpG))

#Question 8.
#How many bases are bivalently marked within 10kb of CpG Islands?

#Add flanking regions (10kb each) to the CpG islands.
autosome_CpG_resized <- resize(autosome_CpG, width=20000+width(autosome_CpG), fix='center')

#Find the regions of bivalently marked bases in the extended CpG islands and then count the number of bases there. 
sum(width(intersect(autosome_CpG_resized, bivalent)))

#Question 9.
#How big a fraction (expressed as a number between 0 and 1) of the human genome is contained in a CpG Island?

#Download the hg19 genome.
query(ah, 'RefSeq')
ah_hg19 <- ah[['AH5040']]

#Find out the sequence lengths of the chromosomes, select the autosomes only, and then sum the base counts.
ah_hg19_size <- sum(as.numeric(seqlengths(ah_hg19)[1:22]))

#Find out the number of bases in the CpG islands and then divide it by the genome size.
sum(width(autosome_CpG))/ah_hg19_size

#Question 10.
#Compute an odds-ratio for the overlap of bivalent marks with CpG islands.

#Create a contingency table.
ConTable <- matrix(0, nrow = 2, ncol = 2)
row.names(ConTable) <- c('bivalent in', 'bivalent out')
colnames(ConTable) <- c('CpG in', 'CpG out')

#Fill out the contingency table with the base count associated with each category.
ConTable[1,1] <- sum(width(bivalentINTERSECTCpG))
ConTable[1,2] <- sum(width(bivalent))-sum(width(bivalentINTERSECTCpG))
ConTable[2,1] <- sum(width(autosome_CpG))-sum(width(bivalentINTERSECTCpG))
ConTable[2,2] <- ah_hg19_size-sum(ConTable)

#Calculate the odds ratio.
ConTable[1,1]*ConTable[2,2]/ConTable[1,2]/ConTable[2,1]