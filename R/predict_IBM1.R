
#' (IBM1) Compute translation probabilities given a foreign sentence.
#'
#' Takes a sentence in a foreign language and produces a list of possible translations and their probabilities based on the IBM1 model.
#' @param fsentence sentence in language you'd like to translate (vector or space-delimited string)
#' @param tmatrix as output from IBM1()
#' @param threshold reduce the number of words to consider by eliminating very unlikely ones with Prob<=threshold (default=0.001).
#' @param maxlength only consider translations which have maxlength words or fewer. By default, will look for translations with the same number of words (or fewer) as fsentence.
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
#' out = IBM1(e,f,maxiter=50,eps=0.01);
#'
#' # possible english translations and their probabilities
#' predict.IBM1(fsentence="une bière sil vous plaît",tmatrix=out$tmatrix)
#' @import Matrix
#' @export
predict_IBM1 = function(fsentence, tmatrix, threshold=0.001,
                       maxlength=length(unlist(stringr::str_split(fsentence, pattern=" ")))  ) {

  # extract relevant values from tmatrix
  fsentence = unlist(stringr::str_split(fsentence, pattern=" "))
  lf = length(fsentence)
  tmp = tmatrix[,fsentence]
  tmp = tmp[rowSums(tmp)>threshold,]

  # compute probabilities for all possible sentences of length <= maxlength
  returnme = data.frame(esentence=character(0), pr=numeric(0))
  for (le in 1:maxlength) {

    candidates = gtools::combinations(n=nrow(tmp),r=le,v=rownames(tmp),repeats.allowed=TRUE)
    returnme_le = data.frame(esentence=rep("",nrow(candidates)),
                             pr=rep(-1,nrow(candidates)))
    for (i in 1:nrow(candidates)) {
      returnme_le$esentence[i] = stringr::str_flatten(candidates[i,], collapse=" ")
      returnme_le$pr[i] = prod(colSums(tmp[candidates[i,],,drop=FALSE])) / (lf^le)
    } # end for

    returnme = rbind(returnme,returnme_le)

  } # end for

  # sort dataframe and return
  returnme = returnme[order(returnme$pr,decreasing=TRUE),]; rownames(returnme) = NULL
  return(returnme)

} # predict.IBM1
