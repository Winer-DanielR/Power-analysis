#' Calculate rates of models compared to simulation adapted from Wreczyka et al. 2017
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
