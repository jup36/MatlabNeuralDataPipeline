#!/home/denisovg/R/R-3.1.2/bin/Rscript

#
# Copyright (C) 2016 by Howard Hughes Medical Institute.
#

# Take as input TSV file of features

library("ADPclust")
library("Rcpp")

traceback()

cmd_args = commandArgs();
cat("length(cmd_args)=", length(cmd_args), "\n")
nargs = length(cmd_args)
cat("nargs=", nargs, "\n")
if (nargs < 8)
{
  cat("\nUsage: TNC_SS_ADPclust.R table segment_id shank_id \n")
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

ans            <- adpclust(data_matrix, nclust = 2:10, centroids = "auto", ac = 1, f.cut = 0.1, fdelta = "mnorm",
                           mycols = NULL, dmethod = "euclidean", verbose = FALSE, draw = TRUE, findSil = TRUE)
summary(ans)

X11()

plot(ans)

message("Press Return To Continue")
invisible(readLines("stdin", n=1))

