convertFCS<-function(fNames,fNameOut,reduce=1,removeColumns=NULL,verbose=TRUE) {
 output<-file(fNameOut,"wb")
 writeBin(as.integer(length(fNames)),output,size=4)
 if(verbose)
  message(sprintf("Reading FCS %s for metadata...",fNames[1]))
 Q<-flowCore::read.FCS(fNames[1])
 if(!is.null(removeColumns))
  Q@exprs<-Q@exprs[,-removeColumns]
 M<-ncol(Q@exprs)
 writeBin(as.integer(M),output,size=4)
 btm<-rep(Inf,M*2)
 top<-rep(-Inf,M*2)
 for(f in fNames){
  if(verbose)
   message(sprintf("Reading FCS %s...",f))
  Q<-flowCore::read.FCS(f)
  if(!is.null(removeColumns))
   Q@exprs<-Q@exprs[,-removeColumns]
  if(reduce<1)
   Q<-Q[sample(nrow(Q),nrow(Q)*reduce)]
  writeBin(as.integer(nrow(Q)),output,size=4)
  J<-as.matrix(Q@exprs)
  apply(J,2,stats::quantile,c(0,.05,.95,1))->qtl
  btm<-ifelse(btm<qtl[1:2,],btm,qtl[1:2,])
  top<-ifelse(top>qtl[3:4,],top,qtl[3:4,])
  # Should be written by-row
  writeBin(as.numeric(t(J)),output,size=4)
  if(verbose)
   message(sprintf("Appended %s.",f))
 }
 minVals<-matrix(btm,nrow=2)[1,]
 maxVals<-matrix(top,nrow=2)[2,]
 loQtls<-matrix(btm,nrow=2)[2,]
 hiQtls<-matrix(top,nrow=2)[1,]
 writeBin(as.numeric(c(minVals,maxVals,loQtls,hiQtls)),output,size=4)
 writeBin(1337L,output)
 close(output)
 if(verbose)
  message(sprintf("File %s completed.",fNameOut))
}

#' Convert a set of FCS files into a tubeset
#'
#' This function allows one to create a >>tubeset<< (a specific form of cytometry data on which flowForest can operate) from a set of FCS files
#'
#' @param fNamesInTubes A list of character vectors; n-th element corresponds to n-th tube, its m-th element is a path to an FCS file describing m-th object
#' @param dirOut Path to write output files; may not exist, function will try to recursively create it
#' @param removeColumns A numeric vector of parameter indices which should be dropped from the converted set; useful to remove nonsense data like event ID or time
#' @param verbose Set to \code{TRUE} to enable process progress reporting
#' @param reduce Only extract random sample of \code{reduce * 100\%} of events. Set to 1 to disable any manipulation
#' @return Invisible \code{NULL}; the important side effect is a set of files \code{<dirOut>/tube<N>.bin} holding the converted data, which can be loaded with \code{loadTubeset} function.
#' @note You don't have to create separate tubesets for training and test data; this can be specified as a mask given to the flowForest function.
#'
#' Converter will overwrite existing \code{tube<N>.bin} files, yet it is adviced to use empty paths (for instance, N-tube set written over N+1-tube one will leave one old tube file, leading to either load failure or unexpected results).
#' @export
createTubeset<-function(fNamesInTubes,dirOut,reduce=1,removeColumns=NULL,verbose=FALSE) {
 # Load flowCore
 if(!requireNamespace("flowCore",quietly=TRUE))
  stop("Converting functionality requires Bioconductor flowCore package.")

 # Check identical number of records for each tube
 if(!all(diff(sapply(fNamesInTubes,length)) == 0))
  stop("Each tube must have the same number of objects.")

 # Check reduce parameter
 stopifnot(length(reduce) == 1 && is.numeric(reduce) && reduce > 0 && reduce <= 1)

 dir.create(dirOut,recursive=TRUE,showWarnings=FALSE)
 for (e in 1:length(fNamesInTubes)) {
  convertFCS(fNamesInTubes[[e]],sprintf("%s/tube%d.bin",dirOut,e),reduce=1,removeColumns=removeColumns,verbose=verbose)
 }
 return(invisible(NULL))
}

#' Convert a set of data.frames into a tubeset
#'
#' This function allows one to create a >>tubeset<< (a specific form of cytometry data on which flowForest can operate) from a set of \code{data.frame}s.
#'
#' @param dfList A list of data.frames or a list of lists of data.frames. In the first case, it is assumed that there is a single tube, otherwise that each top-levl list corresponds to a tube. Within a tube, each data.frame corresponds to an object, and its rows to events recorded for this object. Colums of a data.frame represent the recorded event properties, thus their count must be equal for all data.frames from a tube.
#' @param dirOut Path to write output files; may not exist, function will try to recursively create it
#' @param verbose Set to \code{TRUE} to enable process progress reporting
#' @return Invisible \code{NULL}; the important side effect is a set of files \code{<dirOut>/tube<N>.bin} holding the converted data, which can be loaded with \code{\link{loadTubeset}} function.
#' @note You don't have to create separate tubesets for training and test data; this can be specified as a mask given to the flowForest function.
#'
#' Converter will overwrite existing \code{tube<N>.bin} files, yet it is adviced to use empty paths (for instance, N-tube set written over N+1-tube one will leave one old tube file, leading to either load failure or unexpected results).
#' @export
createTubsetRaw<-function(dfList,dirOut,verbose=FALSE){
 if(inherits(dfList[[1]],"data.frame")){
  #Assuming one tube;
  dfList<-list(dfList);
 }

 dir.create(dirOut,recursive= TRUE,showWarnings=FALSE);

 #Check sizes for equality
 if(!all(diff(sapply(dfList,length))==0))
  stop("All lists should have equal size.");

 #Check input structure
 if(!all(sapply(dfList,sapply,inherits,"data.frame")))
  stop("Input should be a list of lists of data.frames!");

 #Check df sizes
 if(!all(sapply(dfList,function(x){
  all(diff(sapply(x,ncol))==0)
 }))) stop("All data.frames within one list should have equal number of columns.");

 for(tubeN in 1:length(dfList)){
  sprintf("%s/tube%d.bin",dirOut,tubeN)->fNameOut;
  if(verbose)
   message(sprintf("Writing %s...",fNameOut));

  tube<-dfList[[tubeN]];
  output<-file(fNameOut,"wb");
  writeBin(as.integer(length(tube)),output,size=4);
  M<-ncol(tube[[1]]);
  writeBin(as.integer(M),output,size=4);
  btm<-rep(Inf,M*2);
  top<-rep(-Inf,M*2);
  for(frame in tube){
   writeBin(as.integer(nrow(frame)),output,size=4);
   J<-apply(frame,2,as.numeric);
   if(any(is.na(J)))
 stop("NAs or non-numeric values in input!");
   apply(J,2,stats::quantile,c(0,.05,.95,1))->qtl;
   btm<-ifelse(btm<qtl[1:2,],btm,qtl[1:2,]);
   top<-ifelse(top>qtl[3:4,],top,qtl[3:4,]);
   writeBin(as.numeric(t(J)),output,size=4); #As floats
  }
  minVals<-matrix(btm,nrow=2)[1,];
  maxVals<-matrix(top,nrow=2)[2,];
  loQtls<-matrix(btm,nrow=2)[2,];
  hiQtls<-matrix(top,nrow=2)[1,];
  writeBin(as.numeric(c(minVals,maxVals,loQtls,hiQtls)),output,size=4);
  writeBin(1337L,output);
  close(output);
  if(verbose) message("File completed.");
 }
 return(invisible(NULL));
}
