
#' (IBM2) Compute translation probabilities given a foreign sentence.
#'
#' Takes a sentence in a foreign language f and produces a list of possible translations
#' in language e and their respective probabilities based on the IBM1 model. Only returns
#' translations with probabilities > 0. Currently works by considering all permutations
#' (unlike IBM1, order matters, making the list even longer)
#' of most likely (based on threshold parameter) words from language e. Thus, caution
#' should be taken in predicting very long sentences from language f.
#' @param object result from IBM2()
#' @param fsentence sentence in f language you'd like to translate to e language (vector or space-delimited string)
#' @param threshold reduce the number of e language words to consider by eliminating ones with Prob<=threshold (default=0.001).
#' @param maxlength only consider e translations which have maxlength words or fewer. By default, will look for translations with the same number of words (or fewer) as fsentence.
#' @return A data.frame with two columns, sorted in descending order of pr:
#'    \item{esentence}{Possible translation of fsentence.}
#'    \item{pr}{Translation probability (up to a constant).}
#' @examples
#' # download english-french sentence pairs
#' temp = tempfile()
#' download.file("http://www.manythings.org/anki/fra-eng.zip",temp);
#' ENFR = readr::read_tsv(file=unz(temp,"fra.txt"),col_names=c("en","fr","details"));
#' unlink(temp);
#'
#' # a bit of pre-processing
#' e = tolower(stringr::str_squish(tm::removePunctuation(ENFR$en[1:200])));
#' f = tolower(stringr::str_squish(tm::removePunctuation(ENFR$fr[1:200])));
#'
#' # estimate model
#' out = IBM2(e,f,maxiter=50,eps=0.01,init.IBM1=5);
#'
#' # possible english translations and their probabilities
#' decode(out, fsentence="une bière sil vous plaît")
#' @import Matrix
#' @export
decode.IBM2 = function(object, fsentence, threshold=0.001,
                        maxlength=length(unlist(stringr::str_split(fsentence, pattern=" ")))  ) {

  # extract relevant values from tmatrix
  fsentence = unlist(stringr::str_split(fsentence, pattern=" "))
  lf = length(fsentence)
  tmp = object$tmatrix[,fsentence,drop=FALSE]
  tmp = tmp[rowSums(tmp)>threshold,,drop=FALSE]

  # compute probabilities for all possible sentences of length <= maxlength
  returnme = data.frame(esentence=character(0), pr=numeric(0))
  for (le in 1:maxlength) {

    candidates = gtools::permutations(n=nrow(tmp),r=le,v=rownames(tmp),repeats.allowed=TRUE) # use all permutations now!!!!
    returnme_le = data.frame(esentence=rep("",nrow(candidates)),
                             pr=rep(-1,nrow(candidates)))
    for (i in 1:nrow(candidates)) {
      returnme_le$esentence[i] = stringr::str_flatten(candidates[i,], collapse=" ")

      if (length(out$amatrix) < le) {
        returnme_le$pr[i] = 0
      } else if (length(out$amatrix[[le]]) < lf) {
        returnme_le$pr[i] = 0
      } else if (is.null(out$amatrix[[le]][[lf]])) {
        returnme_le$pr[i] = 0
      } else {
        returnme_le$pr[i] = prod(colSums(  tmp[candidates[i,],,drop=FALSE]*
                                             out$amatrix[[le]][[lf]]))
      }

    } # end for

    returnme = rbind(returnme,returnme_le)

  } # end for

  # sort dataframe, drop 0-probabilities, and return
  returnme = returnme[order(returnme$pr,decreasing=TRUE),]
  rownames(returnme) = NULL
  returnme = returnme[returnme$pr>0,,drop=FALSE]
  return(returnme)

} # decode.IBM2
