#' Load tubset
#'
#' Load a conveted set of FCS files
#' @param path A path to a directory with tube<N>.bin files representing tubeset (i.e. output of \code{createTubeset})
#' @return Object representing data
#' @note Tubeset will be freed automatically by R; loading is done by
#' memory mapping, hence loading the same data from other process on the same
#' machine will not consume additional memory.
#' Also, flowForest may load datasets larger than available RAM; still this
#' may impose a substantial slowdown.
#' @export
#' @useDynLib flowForest C_loadTubeset
loadTubeset<-function(path){
 path<-normalizePath(path)
 toCheck<-list.files(path=path, pattern="tube[0123456789]+\\.bin")
 idx<-as.numeric(gsub("tube([0123456789]+)\\.bin","\\1",toCheck))
 if((length(idx)==0)||any(is.na(idx))||(min(idx)!=1)||!all(diff(sort(idx))==1))
  stop("Tubset is wrong, expected <path>/tube1.bin, tube2.bin, ..., tubeN.bin.")
 files<-sprintf("%s/tube%d.bin",path,1:max(idx))
 stopifnot(length(files)>0)
 if(!all(file.exists(files)))
  stop("Some files are not reachable!")
 tubeset<-.Call(C_loadTubeset,files)
 class(tubeset)<-"tubeset"
 tubeset
}

#' Print tubset
#'
#' Prints the tubeset object.
#' @param x Object of a class \code{tubeset}
#' @param ... Ignored
#' @return invisibly \code{x}
#' @export
print.tubeset<-function(x,...){
 cat(sprintf(" Tubeset with %s tubes, %s objects and %s parameters.\n Average flow frame has about %s events.\n\n",
  x[1],x[2],x[3],x[4]))
 return(invisible(x))
}

#' Extract events
#'
#' Read a given subset of events from a single tube. Useful for visualisation.
#' @param set Object representing data
#' @param tube Tube index (must be a single number)
#' @param events Indices of the selected events from a flowset composed of all flow frames in a given tube
#' @return Data frame, rows are events in a requested order, columns are parameters as in the original flow frame.
#' @export
#' @useDynLib flowForest C_extractEvents
extractEvents<- function(set,tube,events){
 tube<-as.integer(tube)
 stopifnot(length(tube)==1)
 stopifnot(tube>0)
 stopifnot(tube<=set[1])

 events<-as.integer(events)
 stopifnot(all(!is.na(events)))
 stopifnot(min(events)>0)
 stopifnot(max(events)<=set[4+tube])
 .Call("C_extractEvents",set,tube,events)->z
 data.frame(matrix(z,ncol=set[3],byrow=TRUE))
}

#' Generate a random gate
#'
#' This function generates a random gate on a randomly selected tube; flowForest builds many such gates to find ones which can be useful for separating classes. In principle, there can be many ways in which such gate is generated, but the current version of flowForest only supports one, \code{pseudosphere}.
#'
#' In the \code{pseudosphere} mode, the code first selects a random tube. Then, a random event from the tube, which will be the geometrical centre of the final gate. Finally, the gate is expanded in all dimensions, by a random fraction of the range between 5-th and 95-th percentile in a certain dimension.
#' @param set Flow cytometry measurements, in a form of an object of the class \code{tubeset}, as that created by the function \code{loadTubeset}.
#' @param gateMode Mode of the gate generation; only \code{pseudosphere} currently supported.
#' @return Object of a class \code{gate}, random gate of a mode d.
#' @export
#' @useDynLib flowForest C_randomGate
makeRandomGate<-function(set,gateMode){
 if(identical(gateMode,'pseudosphere')){
  ans<-.Call(C_randomGate,set,0L,0L)
 }else{
  stop(sprintf("Unsupported gateMode %s!",gateMode));
 }
 class(ans)<-"gate"
 ans
}

#' Apply a gate on a tubeset
#'
#' @param set Tubeset
#' @param gate Gate
#' @param mask A numeric mask specifying on which objects gate shall be applied.
#' @return Matrix of a dim number of objects x number of parameters+1.
#' Rows correspond to objects, first column to a number of events inside gate,
#' n+1 column to average value of a parameter n over events inside gate.
#' @export
#' @useDynLib flowForest C_dotGateMasked
##!! mask must be an integer mask, not logical!
dotGate<-function(set,gate,mask=1:set[2])
 t(.Call(C_dotGateMasked, set, as.integer(mask), gate))

#' Check whether certain events are enclosed by a gate
#'
#' @param set Tubeset
#' @param gate Gate
#' @param idx Indices of the selected events, similar to what \code{\link{extractEvents}} expects
#' @return Boolean vector of a same size and orfer as \code{idx}. Holds \code{TRUE}s for events inside the given gate.
#' @export
#' @useDynLib flowForest C_inGate
inGate<-function(set,gate,idx) .Call(C_inGate,set,gate,idx)
