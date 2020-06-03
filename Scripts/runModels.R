# Run methylKit model adapted from Wreczyka et al. 2017
#
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