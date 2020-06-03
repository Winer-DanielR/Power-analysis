

#load libraries and functions
#!!!!check to make sure all packages are used
library("methylKit")
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(emdbook)

#using function from Wreczycka et al. 2017 to simulate data (dataSim2 function)
source("./functions/dataSim2.R")

##Simulating Datasets##
#modified from Wreczycka et al. 2017

#for different effect sizes
effects = c(5, 10, 15, 20, 25)
cores=20
replicates = c(2,4,6,8)
models.res=list()
set.seed(111)

for(replicate in replicates){
  treatments = as.list
  
  for(effect in effects){
    
    # Effect by the treatment
    print(effect)
    
    # Generate simulated data using methylKit library
    sim.methylBase = dataSim2(replicates=replicate,
                              sites=878645,
                              treatment=c(1,1,1,0,0,0),
                              percentage=1,
                              effect=effect,
                              add.info=TRUE)
    # Run models 
    models.res[[as.character(effect)]] = run.models(sim.methylBase, cores=cores,
                                                    difference=5, qvalue=0.01)
    
  }}