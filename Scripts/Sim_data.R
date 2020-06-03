

# load libraries and functions
# !!!!check to make sure all packages are used
library("methylKit")
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(emdbook)

# using function from Wreczycka et al. 2017 
# to simulate data (dataSim2 function)
source("./Scripts/functions/dataSim2.R")

# and a modified function from Wrecyzka et al. 2017 
# to run the methylKit model on simulated data
source("./Scripts/functions/runModels.R")

## Running the Simulation##

#set different effect sizes
effects = c(5, 10, 15, 20, 25)
cores=20
#set different replicates
replicates = c(2,4,6,8)
#set empty list for results
models.res=list()
set.seed(111)

#iterate through each replicate
for(replicate in replicates){
  # set treatment groups so that 1/2 are 1 (exposed) 
  # and 1/2 are 0 (control)
  treatments = rep(0, replicate/2)
  treatments - append(treatments, rep(1,replicate/2))
  print(replicate)
  #iterate through each effect size
  for(effect in effects){
    
    # Effect by the treatment
    print(effect)
    
    # Generate simulated data using methylKit library and dataSim2 function
    sim.methylBase = dataSim2(replicates=replicate,
                              sites=878645,
                              treatment=treatments,
                              percentage=1,
                              effect=effect,
                              add.info=TRUE)
    
    # Run models using runModels function and return the matrix
    # of true positives (TP), false positives (FP), true negatives (TN)
    # and false negatives (FN) 
    models.res[[as.character(effect)]] = run.models(sim.methylBase, cores=cores,
                                                    difference=5, qvalue=0.01)
    
  }}

models.res.diff.orig = lapply(models.res, function(x) x$diff.list)
models.res.orig = lapply(models.res, function(x) x$rates)
names(models.res.orig) = effects
names(models.res.diff.orig) = effects

#maybe remove this and adjust? 
models.res=models.res.orig
models.res.diff=models.res.diff.orig

#convert to matrix
models.res.ma = do.call("rbind", models.res)
#convert to dataframe 
models.res.df = data.frame(models.res.ma)

models.res.df = cbind(models.res.df, 
                      tool=rownames(models.res.ma),
                      effect=as.factor(as.numeric(sapply(effects, 
                                                         function(x) 
                                                           rep(x, nrow(models.res[[1]])  )))))
