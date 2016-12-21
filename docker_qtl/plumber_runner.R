#
# use the plumber library
#
library(plumber)


#
# get the command line argumnets
#
args <-commandArgs(TRUE)


#
# script to run
#
script <- args[1]

#
# create the plumber object 
#
r <- plumb(script)

#
# Run!
#
r$run()
