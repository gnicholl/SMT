
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
#' @param useIBM1 (default FALSE) If TRUE, use IBM1 decoding method first to find most likely translations based on combinations, then use permutations only of this subset of translations.
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
#' translations = decode(out, fsentence="une bière sil vous plaît")
#'
#' # 10 most likely translations
#' translations = translations[1:10,,drop=FALSE]
#'
#' # might be a bit faster to screen translations with IBM1 method first:
#' translations = decode(out, fsentence="une bière sil vous plaît", useIBM1=TRUE)
#' @import Matrix
#' @export
decode.IBM2 = function(object, fsentence, threshold=0.001, useIBM1=FALSE,
                        maxlength=length(unlist(stringr::str_split(fsentence, pattern=" ")))  ) {

  # extract relevant values from tmatrix
  fsentence = unlist(stringr::str_split(fsentence, pattern=" "))
  lf = length(fsentence)
  tmp = object$tmatrix[,fsentence,drop=FALSE]
  tmp = tmp[rowSums(tmp)>threshold,,drop=FALSE]

  if (!useIBM1) {

    # compute probabilities for all possible sentences of length <= maxlength
    returnme = data.frame(esentence=character(0), pr=numeric(0))
    for (le in 1:maxlength) {

      candidates = gtools::permutations(n=nrow(tmp),r=le,v=rownames(tmp),repeats.allowed=TRUE) # use all permutations now!!!!
      returnme_le = data.frame(esentence=rep("",nrow(candidates)),
                               pr=rep(-1,nrow(candidates)))
      for (i in 1:nrow(candidates)) {
        returnme_le$esentence[i] = stringr::str_flatten(candidates[i,], collapse=" ")

        if (length(object$amatrix) < le) {
          returnme_le$pr[i] = 0
        } else if (length(object$amatrix[[le]]) < lf) {
          returnme_le$pr[i] = 0
        } else if (is.null(object$amatrix[[le]][[lf]])) {
          returnme_le$pr[i] = 0
        } else {
          returnme_le$pr[i] = prod(colSums(  tmp[candidates[i,],,drop=FALSE]*
                                               object$amatrix[[le]][[lf]]))
        }

      } # end for

      returnme = rbind(returnme,returnme_le)

    } # end for

  } else { # useIBM=TRUE

    # get most likely word combos from IBM1 method
    class(object) = "IBM1"
    IBM1predictions = decode(object=object, fsentence=fsentence, threshold=threshold, maxlength=maxlength)
    class(object) = "IBM2"

    # for each IBM1 combo, try all permutations and probabilities using IBM2 method
    returnme = data.frame(esentence=character(0), pr=numeric(0))
    for (j in 1:nrow(IBM1predictions)) {
      esentence   = unlist(stringr::str_split(IBM1predictions$esentence[j], pattern=" "))
      le          = length(esentence)
      candidates  = unique(gtools::permutations(n=le,r=le,v=esentence,set=FALSE,repeats.allowed=FALSE))
      ncand       = nrow(candidates)
      returnme_le = data.frame(esentence=rep("",ncand), pr=rep(-1,ncand))

      for (i in 1:ncand) {
        returnme_le$esentence[i] = stringr::str_flatten(candidates[i,], collapse=" ")

        if (length(object$amatrix) < le) {
          returnme_le$pr[i] = 0
        } else if (length(object$amatrix[[le]]) < lf) {
          returnme_le$pr[i] = 0
        } else if (is.null(object$amatrix[[le]][[lf]])) {
          returnme_le$pr[i] = 0
        } else {
          returnme_le$pr[i] = prod(colSums(  tmp[candidates[i,],,drop=FALSE]*
                                               object$amatrix[[le]][[lf]]))
        }

      } # for i

      returnme = rbind(returnme,returnme_le)
    } # for j

  }

  # sort dataframe, drop 0-probabilities, and return
  returnme = returnme[order(returnme$pr,decreasing=TRUE),]
  rownames(returnme) = NULL
  returnme = returnme[returnme$pr>0,,drop=FALSE]
  return(returnme)

} # decode.IBM2
