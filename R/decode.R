
#' @export
decode = function(object,...) UseMethod("decode")

decode.default = function(object,...) {
  warning("no method defined for provided class")
}
