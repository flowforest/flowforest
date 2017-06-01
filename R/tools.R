bootstrap<-function(mask){
 if(length(mask)==1)
  return(mask)
 unique(sample(mask,length(mask),replace=TRUE))
}

AUROC<-function(score,cls)
 mean(rank(score)[cls]-1:sum(cls))/sum(!cls)

scores<-function(scores,true){
 if(length(scores)==0||sum(is.nan(scores))!=0)
  return(c(AUROC=NA,accuracy=NA,fScore=NA,precision=NA,recall=NA))
 pred<-scores>0.5
 prec<-sum(pred & true)/sum(pred)
 reca<-sum(pred & true)/sum(true)
 c(
  AUROC=ifelse(length(unique(true))==2,AUROC(scores,true),NA),
  accuracy=mean(pred==true),
  fScore=2/(1/prec+1/reca),
  precision=prec,
  recall=reca)
}

fscores<-function(x) paste(sprintf("%s=%0.3f",names(x),x),collapse=", ")

#' Print gate
#'
#' Prints the gate object.
#' @param x Object of a class \code{gate}
#' @param ... Ignored
#' @return invisibly \code{x}
#' @export
print.gate<-function(x,...) {
  cat(sprintf(" Gate on tube #%s:\n",x[1]+1))
  z<-matrix(x[-1],3)
  cat(paste(
   sprintf(
    "  parameter #%d %s %0.3f",
    as.integer(z[1,]+1),
    ifelse(z[2,]>0,">","<"),
    z[3,]
   ),collapse=" &\n"))
  cat("\n")
  return(invisible(NULL))
}

#' Print flowForest
#'
#' Prints the flowForest object.
#' @param x Object of a class \code{flowForest}
#' @param ... Ignored
#' @return invisibly \code{x}
#' @export
print.flowForest<-function(x,...) {
 cat(" flowForest model\n\n")
 cat(sprintf("      Number of trees: %s\n",x$params$nOfTrees))
 cat(sprintf("            Gate mode: %s\n",x$params$gateMode))
 cat(sprintf(" Number of gate tries: %s\n",x$params$gateTries))
 cat("\n")
 cat("  OOB estimate of error: \n")
 print(utils::tail(x$oob,1))
 cat("\n")
 if(length(x$params$testMask) > 0){
  cat(" TEST estimate of error: \n")
  print(utils::tail(x$test,1))
 }
 return(invisible(x))
}
