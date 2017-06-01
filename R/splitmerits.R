#' Gini index and optimal threshold for columns of scores
#'
#' @description Information gain and optimal threshold for columns of scores
#' @param scores A matrix of scores, each column is one variant
#' @param y A logical decision
#' @return A matrix of results; cols correspond to variants, rows to respectively GI and optimal threshold
#' @export
#' @useDynLib flowForest C_colGI
#' @examples
#' colGI(cbind(iris[,-5],rep(77,150)),iris$Species=='versicolor')

colGI<-function(scores,y) {
 y<-as.logical(y)
 if (any(is.na(y)))
  stop("No NAs in y allowed.")
 if (any(is.na(scores)))
  stop("No NAs in scores allowed.")
 scores<-as.matrix(scores)
 if (!is.double(scores))
  stop("Scores can't be converted into a double matrix.")
 .Call(C_colGI,scores,y)
}
