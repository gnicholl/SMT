


#' Viterbi alignments for IBM1 and IBM2
#'
#' Finds best alignment for two sentences given IBM1 or IBM2 model estimates.
#' @export
viterbi_align = function(object,target_sen,source_sen) {
  stopifnot(class(object) %in% c("IBM1","IBM2"))

  le = length(target_sen)
  lf = length(source_sen)
  a = rep(NA,le)
  for (j in 1:le) {
    if (class(object)=="IBM1") a[j] = which.max(sapply(X=1:lf, FUN=function(i) object$tmatrix[[target_sen[j]]][[source_sen[i]]] ))
    if (class(object)=="IBM2") a[j] = which.max(sapply(X=1:lf, FUN=function(i) object$tmatrix[[target_sen[j]]][[source_sen[i]]] * object$amatrix[[le]][[lf]][j,i] ))
  }
  return(a)
}
