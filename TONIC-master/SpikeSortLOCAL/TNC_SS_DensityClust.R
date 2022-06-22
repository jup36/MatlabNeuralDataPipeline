#!/usr/local/R-3.0.2/bin/Rscript

#
# Copyright (C) 2016 by Howard Hughes Medical Institute.
#

# Take as input TSV file of features

library("densityClust")
library("Rcpp")

traceback()

cmd_args = commandArgs();
cat("length(cmd_args)=", length(cmd_args), "\n")
nargs = length(cmd_args)
cat("nargs=", nargs, "\n")
if (nargs < 8)
{
  cat("\nUsage: TNC_SS_SortSpikesDensityClust.R table segment_id shank_id \n")
  quit("yes")
}

table  <- toString(cmd_args[6]) # table(s) with training data (and new data)
cat("table=", table, "\n")

seg    <- as.integer(cmd_args[7]) # animal_ids with good traces
shank  <- as.integer(cmd_args[8]) # animal_ids with bad  traces

all_data <- data.frame(read.csv(table, header=T, sep="\t"))

data_seg       <- subset(all_data, subset=(as.integer(all_data[[1]]) == as.integer(seg)))
data_seg_shank <- subset(data_seg, subset=(as.integer(data_seg[[2]]) == as.integer(shank)))

data_matrix    <- as.matrix(data_seg_shank)

dist_matrix    <- dist(data_matrix[,4:length(data_matrix[1,])])

dataClust      <- densityClust(dist_matrix, gaussian=TRUE)

X11()
plot(dataClust)

dataClust      <- findClusters(dataClust)

plotMDS(dataClust, rho=500000, delta=12)

message("Press Return To Continue")
invisible(readLines("stdin", n=1))

