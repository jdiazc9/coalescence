rm(list=ls())

# load data
metadata <- read.table("metaData.csv",sep=",",header=T)
taxonomy <- read.table("taxonomy_silva.csv",sep=",",header=T)
otus <- read.table("OTU_table.csv",sep=",",header=T)

data <- merge(taxonomy,otus,by="X")

# are there any species (rows) that have abundance 0 in all samples?\
n <- which(rowSums(data[,grepl("S.*\\_",colnames(data))])==0)

# keep only the metadata of experiments in the OTU table
metadata <- metadata[,]

# normalize abundances (by just sum of columns? )
data[,grepl("S.*\\_",colnames(data))] <-
  data[,grepl("S.*\\_",colnames(data))]/colSums(data[,grepl("S.*\\_",colnames(data))])
