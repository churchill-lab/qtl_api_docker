#
# use the plumber library
#
library(plumber)


#
# get the command line argumnets
#
args <-commandArgs(TRUE)


if (length(args) < 2) {
    stop("Please supply the script and the Rdata file")
}

#
# script to run
#
script <- args[1]


#
# Rdata to load
#
rdatafile <- args[2]
print(paste0("Loading ", rdatafile))
load(rdatafile)


#
# create the plumber object 
#
r <- plumb(script)


#
# Run!
#
r$run()
