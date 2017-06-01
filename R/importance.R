extractGatesFromTree<-function(tree){
 xe<-new.env()
 xe$listOfGates<-list()
 extract<-function(tree){
  if(length(tree)>1){
   xe$listOfGates[[length(xe$listOfGates) + 1]]<-list(gate=tree$gate,score=tree$score,parameter=tree$parameter)
   extract(tree$left_subtree)
   extract(tree$right_subtree)
  }
 }
 extract(tree)
 return(xe$listOfGates)
}

#' Extract all gates used to build flowForest
#'
#' @param flowForest Object of a class \code{flowForest}, as that created by the function flowForest.
#' @param nOfTube Number of tube in tubeset on which flowForest was build; only gates on this tube will be extracted.
#' @export
extractGatesFromflowForest<-function(flowForest,nOfTube){
 listOfGates<-list()
 for(t in flowForest$forest)
  listOfGates<-c(listOfGates, extractGatesFromTree(t))
 listOfGates<-listOfGates[sapply(listOfGates, FUN=function(gate, nOfTube) gate$gate[1]+1 == nOfTube, nOfTube)]
 return(listOfGates)
}

scoreCells<-function(tubeset,nOfTube=1,cellsToCheck,listOfGates,sizePenalty=FALSE){
  listOfBareGates<-lapply(listOfGates, FUN=function(g) g$gate)
  ginisOfGates<-sapply(listOfGates, FUN=function(g) g$score)  #giniIndex; lower is better
  max<-max(ginisOfGates)
  min<-min(ginisOfGates)
  if(min < max){
 scoresOfGates<-sapply(ginisOfGates,
   FUN=function(s) (max-s)/(max-min))
   #linear transformation to [0,1]; now higher is better;
  } else {
 scoresOfGates<-rep(1, length(listOfBareGates))
  }

  scores<-rep(0, length(cellsToCheck))
  for(i in 1:length(listOfBareGates)){
   ing<-inGate(tubeset, listOfBareGates[[i]], cellsToCheck)
   scores<-scores + (ing * ifelse(is.null(scoresOfGates), 1, scoresOfGates[[i]]) / ifelse(sizePenalty, max(10, sum(sqrt(ing))), 1))
  }
  return(scores)
}

detectOutliers<-function(dataFrame, q=.005){
  upperQuantiles<-apply(dataFrame, 2, stats::quantile, 1-q)
  lowerQuantiles<-apply(dataFrame, 2, stats::quantile, q)
  outliers<-rep(FALSE, nrow(dataFrame))
  for(i in 1:ncol(dataFrame)){
 outliers<-outliers | dataFrame[,i] > upperQuantiles[i]
 outliers<-outliers | dataFrame[,i] < lowerQuantiles[i]
  }
  return(outliers)
}

#' Number of events in gate on given tubeset
#'
#' @param gate Object of a class \code{gate}, as that extracted from flowForest by function \code{extractGatesFromflowForest}.
#' @param tubeset Flow cytometry measurements, in a form of a tubeset, as that created by the function \code{loadTubeset}.
#' @param cellsToCheck Vector of integers indicating, to which cells in tube, on which the gate is made, narrow counting. By default all cells in this tube.
#' @export
gateSize<-function(gate, tubeset, cellsToCheck){
  sum(inGate(tubeset, gate, cellsToCheck))
}

#' Identify predictive cell populations
#'
#' @param tubeset Flow cytometry measurements, in a form of a tubeset, as that created by the function \code{loadTubeset}.
#' @param flowForest Object of a class \code{flowForest}, trained by the function flowForest on \code{tubeset}.
#' @param nOfTube Number of tube in tubeset on which flowForest was build; only gates on this tube will be extracted.
#' @param threshold A vector of indices of objects in tubeset which are to be assumed to reside in the test set. By default those not in \code{trainMask}.
#' @param transparency Transparency of events on the scatterplot produced by the function.
#' @param sizePenalty Shall gates comprising many events be considered less important?
#' @param additionalPlots If \code{FALSE}, only histogram of scores of events and scatterplot with important events will be printed. Otherwise, additionally barplot with frequency of usage of each criterion on gates in splits will be printed too.
#' @return Data frame with coordinates of events from the tubeset (all columns except last two), it's scores of importance (the penultimate column) and logical vector indicating, if they were considered as important or not (the last column).
#' @export
importance<-function(tubeset, flowForest, nOfTube=1, threshold=0.95, transparency=0.1, sizePenalty=TRUE, additionalPlots=FALSE){
  if(!suppressPackageStartupMessages(requireNamespace("GGally", quietly=TRUE)) | !suppressPackageStartupMessages(requireNamespace("ggplot2", quietly=TRUE)))
   stop("Viewing importance requires ggplot2 and GGally packages.")
  listOfGates<-extractGatesFromflowForest(flowForest, nOfTube)
  listOfBareGates<-lapply(listOfGates, FUN=function(g) g$gate)
  nCellsInTube<-tubeset[4+nOfTube]
  nOfParameters<-tubeset[3]
  cellsToCheck<-sample(1:nCellsInTube, min(100000/nOfParameters, nCellsInTube)) #for plotting it is sane to reduce number of events;

  message("Number of gates on tube ", nOfTube, " extracted from flowForest: ", length(listOfGates))
  gateParameters<-sapply(listOfGates, FUN=function(x) x$parameter)

  sizesOfGates<-sapply(listOfBareGates, gateSize, tubeset, cellsToCheck)
if(additionalPlots){
  h<-ggplot2::qplot(sizesOfGates, geom='histogram',xlab="Size of gate", ylab="Frequency", bins=15) +
  ggplot2::ggtitle("Histogram of sizes of gates")

   print(h)
   cat ("Press [enter] to see plot of parameters frequency")
   line<-readline()
   h<-ggplot2::qplot(as.factor(gateParameters-1), xlab="Number of parameter", ylab="Frequency") +
  ggplot2::ggtitle("Plot of parameters frequency")
   print(h)

   cat ("Press [enter] to see histogram of scores of events")
   line<-readline()
  }
  X<-extractEvents(tubeset,nOfTube,cellsToCheck)
  outliers<-detectOutliers(X, .005)
  scores<-scoreCells(tubeset, nOfTube, cellsToCheck, listOfGates, sizePenalty)
  q<-stats::quantile(scores, threshold)
  important<-(scores > q)
  X$scores<-scores
  X$important<-important
  X<-X[!outliers,]
  message("Threshold: ", round(q, digits=5), ", events with score above threshold: ", sum(X$important), ", all events:", nrow(X), ", removed outliers: ", sum(outliers))
  p<-ggplot2::qplot(scores, geom='histogram',xlab="Score", ylab="Frequency", bins=20, fill=important) +
   ggplot2::ggtitle("Histogram of scores of events") + ggplot2::scale_fill_manual(values=c("black", "red"))
  print(p)
  cat ("Press [enter] to see scatterplot with importance")
  line<-readline()
  gg_lower<-function(data, mapping, ..., transparency=.1, mask=important) {
 ggplot2::ggplot(data=data, mapping=mapping) +
   ggplot2::geom_point(data=data[!data$important,], alpha=transparency, size=.5, colour="black") +
   ggplot2::geom_point(data=data[data$important,], alpha=transparency, size=1, colour="red")
  }
  gg_diag<-function(data, mapping, ...) {
 ggplot2::ggplot(data=data, mapping=mapping) +
   ggplot2::geom_density(data=data[!data$important,], fill="black", alpha=.5) +
   ggplot2::geom_density(data=data[data$important,], fill="red", alpha=.5)
  }
  p<-GGally::ggpairs(X, columns=1:(ncol(X)-2) , lower=list(continuous=GGally::wrap(gg_lower, transparency=transparency)), diag=list(continuous=gg_diag), upper=list(continuous=GGally::wrap("cor", size=4)))
  print(p)
  return(invisible(X))
}
