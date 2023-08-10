

#' @export
object.size = function(x) {
  if (class(x)=="IBM1") {
    object.size.IBM1(x)
  } else if (class(x)=="IBM2") {
    object.size.IBM2(x)
  } else if (class(x)=="IBM3") {
    object.size.IBM3(x)
  } else {
    utils::object.size(x)
  }
}

object.size.IBM1 = function(object) {
  totsize = 0
  for (a in ls(envir=object$tmatrix)) {
    for (b in ls(envir=object$tmatrix[[a]]) ) {
      totsize = totsize + object.size(object$tmatrix[[a]][[b]])
    }
  }
  return(list("tmatrix"=totsize))
}


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


