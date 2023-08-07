

#' @export
object.size.IBM1 = function(object) {
  totsize = 0
  for (a in ls(envir=object$tmatrix)) {
    for (b in ls(envir=object$tmatrix[[a]]) ) {
      totsize = totsize + object.size(object$tmatrix[[a]][[b]])
    }
  }
  return(list("tmatrix"=totsize))
}

#' @export
object.size.IBM2 = function(object) {
  totsize = 0
  for (a in ls(envir=object$tmatrix)) {
    for (b in ls(envir=object$tmatrix[[a]]) ) {
      totsize = totsize + object.size(object$tmatrix[[a]][[b]])
    }
  }
  return(
    list(
      "tmatrix"=totsize,
      "amatrix"=object.size(object$amatrix)
      )
    )
}

#' @export
object.size.IBM3 = function(object) {
  totsize = 0
  for (a in ls(envir=object$tmatrix)) {
    for (b in ls(envir=object$tmatrix[[a]]) ) {
      totsize = totsize + object.size(object$tmatrix[[a]][[b]])
    }
  }

  fsize = 0
  for (a in ls(envir=object$fmatrix)) {
    fsize = fsize + object.size(object$fmatrix[[a]])
  }

  return(
    list(
      "tmatrix"=totsize,
      "dmatrix"=object.size(object$dmatrix),
      "fmatrix"=fsize
    )
  )
}


