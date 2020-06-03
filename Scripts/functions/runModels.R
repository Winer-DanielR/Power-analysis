# This function is adapted from Wreczycka et al. 2017
# Author: Katarzyna Wrecycka
# Date: 2017
# Title: compareUsingSimulated.md
# Availablity: https://github.com/BIMSBbioinfo/Strategies_for_analyzing_BS-seq

#' Calculate rates of model compared to simulation adapted from Wreczyka et al. 2017
#' basically determines all of the true positives, true negatives, false positives, and false negatives
#' 
#' @param simOutput this is the output of dataSim2 function
#' @param sub.methylDiff this is the q-value filtered methylDiff object
#'                       output of getMethylDiff()
#' @return returns a vector of accuracy metrics, TP, FP, Sensivity, etc
calc.rates<-function(simOutput, sub.methylDiff){
  
  all=paste(simOutput[[1]][[1]],simOutput[[1]][[2]],
            simOutput[[1]][[3]])
  
  true.dm=all[simOutput[[2]]]
  true.ndm=all[-simOutput[[2]]]
  
  pred.dm=paste(sub.methylDiff[[1]],sub.methylDiff[[2]],
                sub.methylDiff[[3]])
  pred.ndm=all[! all %in% pred.dm]
  
  TP=sum(true.dm %in% pred.dm)
  
  FN=sum(pred.ndm %in% true.dm)
  
  FP=sum(pred.dm %in% true.ndm)
  
  TN=sum(pred.ndm %in% true.ndm)
  
  p = TP / (TP + FP)
  r = TP / (TP+FN)
  f_score = 2*((p*r)/(p+r))
  
  return(c(TP=TP,FN=FN,FP=FP,TN=TN,
           acc=(TP+TN)/length(all),
           spec=TN/(TN+FP) ,
           sens=TP/(TP+FN),
           f_score= f_score,
           precision=as.numeric(TP / (TP + FP)),
           recall=r,
           NPV= as.numeric(TN / (TN + FN))
  ) )
}

# Call differentially methylated cytosines using methylKit
# It calculate true positive positives (TP), false negatives (FN), false positives (FP),
# accuracy (acc), specificity (spec), sensiticity (sens) and F-score (f_score).
#
#' @param sim.methylBase a methylBase object from the methylKit library
#' @param cores a number of cores
#' @param difference cutoff for absolute value of methylation percentage change
#'                   between test and control (default:5)
#' @param qvalue cutoff for qvalue of differential methylation statistic
#'               (default:0.01)
#' @return returns a matrix with TP, FN, FP, TN, acc, spec, sens, f_score (columns)
#'         using tools that calculate differentially methylated regions (rows)

run.models = function(sim.methylBase, cores=1,
                      difference=5, qvalue=0.01){
  
  require(methylKit)
  
  
  ## run methylkit
  combined = data.frame(test=c("F", "Chisq","F", "Chisq"),
                        adjust="qvalue",
                        overd=c("none","none", "MN", "MN"),
                        name=c("methylKit.F.qvalue.none",
                               "methylKit.Chisq.qvalue.none",
                               "methylKit.F.qvalue.MN",
                               "methylKit.Chisq.qvalue.MN"), 
                        stringsAsFactors = FALSE)
  diff.list = list()
  methylKit.list=list()
  for(i in 1:nrow(combined)){
    co = combined[i,]
    methylkit.obj <- calculateDiffMeth(sim.methylBase[[1]], 
                                       overdispersion=co$overd,
                                       adjust = co$adjust,
                                       test=co$test,
                                       mc.cores=cores)
    methylkit.obj.diff = getMethylDiff(methylkit.obj, 
                                       difference=difference,qvalue=qvalue)
    diff.list[[i]] <- methylkit.obj.diff
    methylKit.list[[i]]=calc.rates(sim.methylBase,
                                   methylkit.obj.diff)
    
  }
  names(methylKit.list) <- combined$name
  names(diff.list) <- combined$name
  
  list(
    diff.list=diff.list,
    rates=do.call("rbind",methylKit.list)
  )
}