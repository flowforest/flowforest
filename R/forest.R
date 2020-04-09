predictTree<-function(set,inputMask=1:set[2],tree) {
 if(length(tree)==1)
  return(rep(tree$decision,length(inputMask)))
 # Which means: we are at leaf, return decision
 boolToLeft<-(dotGate(set,tree$gate,inputMask)[,tree$parameter] > tree$threshold)
 boolToRight<-!boolToLeft
 ans<-rep(NA,length(inputMask))
 ans[boolToLeft]<-predictTree(set,inputMask[boolToLeft],tree$left_subtree)
 ans[boolToRight]<-predictTree(set,inputMask[boolToRight],tree$right_subtree)

 return(ans)
}

#' Predict method for flowForest objects
#'
#' This function predicts classes of new objects with given \code{flowForest} object.
#' @method predict flowForest
#' @param object \code{flowForest} object (as that created by the function \code{flowForest}), that will be used for prediction.
#' @param x Test data, flow cytometry measurements in a form of a \code{tubeset}, as that created by the function \code{loadTubeset}.
#' @param testMask A vector of indices of objects in tubeset which are to be assumed to reside in the test set. By default the whole set.
#' @param type \code{"score"} or \code{"class"} - types of returned value. In first case score from trees voting for each object is returned; otherwise - the most frequent class in voting.
#' @param ... Ignored
#' @return If \code{type == "class"} - predicted class of each object in tubeset is returned; if \code{type == "score"} - score is returned - a real number from between 0 and 1, obtained from voting.
#' @export
predict.flowForest<-function(object,x,testMask=1:x[2],type="response",...) {
 stopifnot(inherits(object,"flowForest"))
 stopifnot(inherits(x,"tubeset"))
 forest<-object$forest
 votes<-rep(0,length(testMask))
 for(tree in forest){
  votes<-votes+predictTree(x,1:x[2],tree)[testMask]
 }
 if(identical(type,"score")) return(votes/length(forest))
 return((votes/length(forest))>.5)
}

#' Train a flowForest model
#'
#' @param x Flow cytometry measurements, in a form of an object of the class \code{tubeset}, as that created by the function \code{loadTubeset}.
#' @param y Decision; must be a logical vector of a same length as the number of objects in \code{x}. Unknown decisions may be marked with \code{NA}; they can't be in the training set, though
#' @param trainMask A vector of indices of objects in tubeset which are to be assumed to reside in the trainig set. By default the whole set
#' @param testMask A vector of indices of objects in tubeset which are to be assumed to reside in the test set. By default those not in \code{trainMask}
#' @param nOfTrees Number of trees to grow
#' @param gateTries Number of gates to try at each split; this does not count gates rejected due to insufficient event fraction
#' @param gateMode Mode of the gate generation; only \code{pseudosphere} is currently supported. See \code{\link{makeRandomGate}} for more information.
#' @param asymmetric If \code{TRUE}, only gates which hold on average more events for the \code{TRUE} class will be taken into account. This is useful to force flowForest to focus on event populations characteristic for the \code{TRUE} class.
#' @param keepForest Shall the forest structure be saved? This is required for importance and predicting new data, but requires more memory.
#' @param verbose If \code{TRUE}, performance scores of \code{flowForest} counted on test set and OOB will be printed after growing each tree. Otherwise, only after building whole \code{flowForest}.
#' @param relaxedGateAccept Accept gates with some event, not only at least one event for all objects.
#' @return Object of a class \code{flowForest}
#' @export
flowForest<-function(
 x,y,
 trainMask=1:x[2],testMask=(1:x[2])[!(1:x[2])%in%trainMask],
 nOfTrees=100,gateTries=5,gateMode='pseudosphere',asymmetric=TRUE,
 keepForest=TRUE,verbose=TRUE,relaxedGateAccept=FALSE
){
 # Decision must be binary
 stopifnot(is.logical(y))
 # Data must a tubeset
 stopifnot(inherits(x,"tubeset"))
 # Decision and data must match
 stopifnot(length(y)==x[2])
 # Train and test have to be vectors of integers/doubles
 stopifnot(is.numeric(trainMask)&&is.numeric(testMask))
 # Train and test mustn't overlap
 stopifnot(all(!(testMask%in%trainMask)))
 stopifnot(is.numeric(nOfTrees)&&is.numeric(gateTries))
 stopifnot(nOfTrees>0 &&gateTries>0)
 stopifnot(gateMode%in%c('pseudosphere'))
 stopifnot(is.logical(keepForest)&&(length(keepForest)==1))

 crit<-if(relaxedGateAccept){function(x) max(x)>1} else {function(x) min(x)>1}

 makeBestSplit<-function(mask=1:x[2]){
  bestScore<-Inf # Better score is smaller score
  bestGate<-NULL

  tried<-good<-impr<-0
  while(TRUE){
   if((gateTries == 0)&&!is.null(bestGate))
    break
   gate<-makeRandomGate(x,gateMode=gateMode)
   tried<-tried+1
   dotResult<-dotGate(x,gate,mask)

   # Gate has to have averagely more events in flowframes of objects with positive class and cannot be too small
   if(
    crit(dotResult[,1]) && #There is something in the gate
    (!asymmetric||(mean(dotResult[y[mask],1])>mean(dotResult[!y[mask],1]))) #More events for TRUE class
   ){
    # Gate is OK, score a try
    good<-good+1
    gateTries<-gateTries-1
    # ... and evaluate it
    evaluation<-colGI(dotResult[,],y[mask])
    iBest<-which.min(evaluation[1,])

    if(evaluation[1,iBest]<bestScore){
      impr<-impr+1
      bestScore<-evaluation[1,iBest]
      bestThre<-evaluation[2,iBest]
      bestSplit<-dotResult[,iBest]>bestThre
      bestParam<-iBest
      bestGate<-gate
      if(bestScore==0) break #we cannot find better gate, so move on!
    }
   }
  }
  if(verbose)
   message(sprintf("Best gate found among %d generated, %d correct, %d improving score.",
    tried,good,impr))
  return(list(bestGate=bestGate,bestScore=bestScore,bestParam=bestParam,bestThre=bestThre,bestSplit=bestSplit))
 }

 maxDepth<-30 #maximum depth of tree; can handle at least 10^9 objects
 tree<-function(mask=1:x[2],depth=0){
  # Class going downstream
  in_cls<-y[mask]
  if(depth>=maxDepth||sum(in_cls)==0||sum(!in_cls)==0){
   # We have a leaf
   leafDecision<-mean(in_cls)>.5
   return(list(
    decision=leafDecision
   ))
  }else{
   newSplit<-makeBestSplit(mask)
   maskLeft<-mask[newSplit$bestSplit]
   maskRight<-mask[!newSplit$bestSplit]
   return(list(
    gate=newSplit$bestGate,
    parameter=newSplit$bestParam,
    threshold=newSplit$bestThre,
    score=newSplit$bestScore,
    left_subtree=tree(maskLeft,depth+1),
    right_subtree=tree(maskRight,depth+1)
   ))
  }
 }

 forest<-list()
 oob<-NULL
 test<-NULL
 oobVotes<-rep(0,length(trainMask))
 oobTries<-rep(0,length(trainMask))
 testVotes<-rep(0,length(testMask))

 for(i in 1:nOfTrees){
  if(verbose) message(sprintf("Growing tree %s/%s...",i,nOfTrees))

  mask<-bootstrap(trainMask)
  newTree<-tree(mask,depth=0)

  # Predict the whole set
  fullPreds<-predictTree(x,1:x[2],newTree)

  # Extract OOB errors
  mask_oob<-!(trainMask %in% mask)
  oobVotes[mask_oob]<-oobVotes[mask_oob]+fullPreds[trainMask][mask_oob]
  oobTries[mask_oob]<-oobTries[mask_oob]+1
  oobSc<-scores((oobVotes/oobTries),y[trainMask])
  oob<-rbind(oob,oobSc)

  if(length(testMask)>0){
   # Extract test errors
   testVotes<-testVotes + fullPreds[testMask]
   tstSc<-scores(testVotes/i,y[testMask])
   test<-rbind(test,tstSc)
   if(verbose||i==nOfTrees) message(sprintf(" OOB %s\nTEST %s",fscores(oobSc),fscores(tstSc)))
  }else{
   if(verbose||i==nOfTrees) message(sprintf(" OOB %s",fscores(oobSc)))
  }

  if(keepForest)
   forest[[length(forest)+1]]<-newTree
 }
 ans<-list(
  forest=forest,
  oob=oob,
  test=test,
  oobVotes=oobVotes,
  testVotes=testVotes,
  oobTries=oobTries,
  params=list(
   nOfTrees=nOfTrees,
   gateTries=gateTries,
   trainMask=trainMask,
   testMask=testMask,
   gateMode=gateMode,
   asymmetric=asymmetric,
   xMetrics=stats::setNames(as.numeric(x)[1:4],c("Tubes","Objects","Parameters","AveEventCount"))
  )
 )
 class(ans)<-"flowForest"
 return(ans)
}
