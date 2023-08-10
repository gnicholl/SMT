
#' Generic decoding function.
#'
#' This generic function translates ("decodes") a given sentence in one language
#' to a sentence in the output language given the model estimates from `object`.
#' Currently supports decoding based on `IBM1` and `IBM2`. See `?decode.IBM1`
#' and `?decode.IBM2`.
#'
#' @export
decode = function(object,...) UseMethod("decode")

decode.default = function(object,...) {
  warning("no method defined for provided class")
}
