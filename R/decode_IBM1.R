
#' (IBM1) Stack decoder
#'
#' Uses the same algorithm as IBM2, and simply calls `decode.IBM2` with argument `IBM1=TRUE`.
#' Only difference is that in IBM1 all alignments
#' are equally likely. See `?decode.IBM2` for syntax.
#' @import collections
#' @export
decode.IBM1 = function(...) {

  return(decode.IBM2(..., IBM1=TRUE))

} # decode.IBM1
