# load libraries and functions
library("methylKit")
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(emdbook)
library(tidyverse)
library(dplyr)
library(LaCroixColoR)
library(lemon)

#set working directory
setwd("~/Documents/GitHub/Power-analysis")

# using function from Wreczycka et al. 2017 
# to simulate data (dataSim2 function)
source("./Scripts/functions/dataSim2.R")

# and a modified function from Wrecyzka et al. 2017 
# to run the methylKit model on simulated data
source("./Scripts/functions/runModels.R")

# add in functions that I wrote to extract and convert datasets
source("./Scripts/functions/extractConvert.R")
## Running the Simulation##

#set different effect sizes
effects = c(5, 10, 15, 20, 25)
cores=20
#set different replicates
replicates = c(2,4,6,8)
set.seed(111)

#now put each list into a list for every replicate
model.res <-setNames(replicate(length(replicates), list()), replicates)

#set an index to keep track of which replicate loop is on
index=0
#iterate through each replicate
for(replicate in replicates){
  # set treatment groups so that 1/2 are 1 (exposed) 
  # and 1/2 are 0 (control)
  treatments = rep(0, replicate/2)
  treatments = append(treatments, rep(1,replicate/2))
  print(replicate)
  #add one to the index for every turn of this loop
  index <- index + 1
  #iterate through each effect size
  for(effect in effects){
    
    # Effect by the treatment
    print(effect)
    
    # Generate simulated data using methylKit library and dataSim2 function
    sim.methylBase = dataSim2(replicates=replicate,
                              sites=50000,
                              treatment=treatments,
                              percentage=5,
                              effect=effect,
                              add.info=TRUE)
    
    
    # Run models using runModels function and return the matrix
    # of true positives (TP), false positives (FP), true negatives (TN)
    # and false negatives (FN) along with other evaluations
    # add the result of the model to the model.res list
    
    model.res[[index]][as.character(effect)] <- run.models(sim.methylBase, cores=cores,
                                                           difference=5, qvalue=0.01)
    
  }}

#set up separate datasets for each type of evaluation using function
true_positives <- extract.data(model.res, replicates, effects, "TP")
sensitivity <- extract.data(model.res, replicates, effects, "sens")
specificity <- extract.data(model.res, replicates, effects, "spec")
f_score <- extract.data(model.res, replicates, effects, "f_score")
accuracy <- extract.data(model.res, replicates, effects, "acc")

#create function for plotting data
# takes data and y label for inputs
# gives plots as outputs
bar.plot <- function(data, ylabel){
  ggplot(data, aes(replicate, eval,fill=effect)) + geom_bar( 
                        stat="identity", colour="black", position="dodge") +
    xlab("Sample Size") + ylab(ylabel) + ylim(0,1) + theme_bw() +
    theme(axis.text=element_text(size=11), legend.position="none", 
          axis.title=element_text(size=13)) + 
    scale_fill_manual(values = lacroix_palette("PassionFruit"))  + labs(fill="Effect Size")
}

#create plots for all evaluation parameters
sensitivity_bar_plot <- bar.plot(sensitivity,"Sensitivity")  
specificity_bar_plot <- bar.plot(specificity, "Specificity")
f_score_bar_plot <- bar.plot(f_score, "F-Score")
accuracy_bar_plot <- bar.plot(accuracy, "Accuracy")

#arrange plots in a grid with a shared legend
grid_arrange_shared_legend(sensitivity_bar_plot, specificity_bar_plot, f_score_bar_plot, 
   accuracy_bar_plot, ncol=2, nrow=2)


