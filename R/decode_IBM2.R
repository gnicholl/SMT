
#' (IBM2) Stack decoder.
#'
#' Implements the IBM2 decoder Algorithm 1 provided in `Wang and Waibel (1998)`.
#' Relies on three models: IBM2 translation model, a sentence length model,
#' and a kgrams language model. If the latter two aren't provided, they are estimated
#' using a similar approach to `Wang and Waibel (1998)`: sentence length is estimated
#' using a poisson regression, and the language model is estimated using 3rd-order kgrams.
#' Heuristics (`max.length`,`threshold`,`max.nsen`) are used to reduce computational time,
#' at the expense of an increased likelihood of not finding the optimal translation.
#'
#' Decoding is based on bayes rule: P(e|f) ~ P(f|e) P(e). (Here ~ means "proportional to").
#' So in order to translate from f to e, we actually estimate IBM2 from e to f to obtain P(f|e).
#' P(e) is the language model and estimated using the `kgrams` package. The `kgrams`
#' documentation provides more details on the different estimation methods available.
#'
#' Actually, the above equation implicitly conditioned on the length of sentence e.
#' Let m be the length of sentence e and n be the length of sentence f.
#' Then we really have P(e|f,m,n) ~  P(f|e,m,n) P(e|m,n).
#' When we decode, we do not know m beforehand. Thus we really need to model P(e,m|f,n).
#' To do this, we assume P(e,m|f,n) = P(e|f,m,n) P(m|n) where P(m|n) is called
#' the sentence length model. The factor of P(m|n) ensures that our decoder
#' selects a translation of reasonable length.
#'
#' @param object result from IBM2()
#' @param target.sentence sentence in f language we'd like to decode (character string)
#' @param max.length the maximum length of sentence in e language to consider. If left NULL (default), it searches the training corpus for the greatest e sentence length associated with the length of the given f sentence.
#' @param threshold Vector of two threshold probabilities. First threshold is a cut off for which subset of e words to consider in our algorithm. The second is a cut-off for whether a hypothesis translation will be added to the stack. Smaller values mean greater chance of finding the right translation, but more computational time.
#' @param max.nsen For each iteration, this is the max number of sentences which are allowed to be added to the stack. Larger values mean greater chance of finding the right translation, but more computational time.
#' @param senlength.model Matrix such that `senlength.model[m,n]` is the probability that e sentence is length m given f sentence is length n. If not provided, it is estimated form the training corpus using a simple poisson regression: `glm(e_lengths ~ f_lengths, family="poisson")`.
#' @param language.model The result of a `kgrams::language_model` estimation. If not provided, we estimate a degree-3 kgrams model from the training corpus using the default "ml" (maximum likelihood) approach. See `kgrams` package for more details.
#' @param IBM1 Default is FALSE. If TRUE, ignores the alignment probabilities of the IBM2 estimation, treating it like an IBM1 model.
#' @return List of best translations found using the stack decoder.
#' @examples
#' # download english-french sentence pairs
#' temp = tempfile()
#' download.file("http://www.manythings.org/anki/fra-eng.zip",temp);
#' ENFR = readr::read_tsv(file=unz(temp,"fra.txt"),col_names=c("en","fr","details"));
#' unlink(temp);
#'
#' # a bit of pre-processing
#' e = removePunctuation(ENFR$en[40001:44000])
#'   e = str_squish(e)
#'   e = tolower(e)
#' f = removePunctuation(ENFR$fr[40001:44000])
#'   f = str_squish(f)
#'   f = tolower(f)
#'
#' # estimate model
#' model = IBM2(target=f,source=e, maxiter=30, init.IBM1=30)
#'   # notice e is source and f is target, even though ultimately
#'   # we want to translate from f to e.
#'
#' # possible english translations and their probabilities
#' best_translation = decode(model, target.sentence="il est un peu rouillé")
#'   # returns "hes a little rusty", as expected
#'
#' @import collections
#' @export
decode.IBM2 = function(object, target.sentence,
                       max.length=NULL, threshold=c(1e-5,1e-10,1e-12), max.nsen=100,
                       senlength.model=NULL, language.model=NULL,
                       IBM1=FALSE) {
  # params
  target.sentence = unlist(stringr::str_split(target.sentence," "))
  ltarget = length(target.sentence)
  if (is.null(max.length)) {
    max.length = max(object$corpus$source_lengths[object$corpus$target_lengths==ltarget]) - 1
  }
  if (length(threshold)==1) threshold = rep(threshold,3)

  # if no sentence length model, use basic poisson regression
  if (is.null(senlength.model)) {
    f_lengths = object$corpus$target_lengths; max.f.length = max(f_lengths)
    e_lengths = object$corpus$source_lengths; max.e.length = max(e_lengths)
    slmodel = glm(e_lengths ~ f_lengths, family="poisson")
    lambda = predict(slmodel, newdata=data.frame(f_lengths=1:max.f.length), type="response")
    senlength.model = matrix(0,nrow=max.e.length,ncol=max.f.length)
    for (i in 1:max.e.length) senlength.model[i,] = dpois(i, lambda=lambda)
  }

  # if no language model, use mle from kgrams
  if (is.null(language.model)) {
    freqs = kgrams::kgram_freqs(object$corpus$source, N = 3, verbose = FALSE)
    language.model = kgrams::language_model(freqs)
  }

  # list of "promising words" -> those with +ve probabilities
  promising_words = NULL
  for (targetword in target.sentence) {
    for ( sourceword in ls(object$tmatrix[[targetword]]) ) {
      if (object$tmatrix[[targetword]][[sourceword]] > threshold[1]) {
        promising_words = c(promising_words,sourceword)
      }
    }
  }
  promising_words = unique(promising_words)
  promising_words = promising_words[promising_words!="<NULL>"]

  # get word frequencies
  e_sentences = lapply(X=object$corpus$source,FUN=function(s) unlist(stringr::str_split(s, " ")))
  e_sentences = unique(unlist(e_sentences))
  wfrq = table(e_sentences)
  wfrq = wfrq[promising_words]
  wfrq = wfrq / sum(wfrq)

  # average translation probability when don't know future words
  predict_cost = function(gj) {
    res = 0
    for (wk in promising_words) {
      if (!is.null(object$tmatrix[[gj]][[wk]])) {
        res = res + wfrq[wk]*object$tmatrix[[gj]][[wk]]
      }
    }
    return(as.numeric(res))
  }
  ex_costs = sapply(X=target.sentence, FUN=predict_cost)

  # scoring function
  tau_kl = function(j,i,prefix,lsource) {
    if (!IBM1) { #IBM2: use alignment probs
      if (length(object$amatrix) < ltarget) return(0)
      if (length(object$amatrix[[ltarget]]) < lsource+1) return(0)
      if (is.null(object$amatrix[[ltarget]][[lsource+1]])) return(0)
      res = object$amatrix[[ltarget]][[lsource+1]][j,i+1]
    } else { #IBM1: all alignments equaly likely
      res = (1/lsource)
    }

    k = length(prefix)
    if (0<=i & i<=k) {
      if (is.null(object$tmatrix[[ target.sentence[j] ]][[ prefix[i] ]])) return(0)
      return(res * object$tmatrix[[ target.sentence[j] ]][[ prefix[i] ]] )
    } else {
      return(res * ex_costs[j])
    }
  }

  tau_l = function(prefix,lsource) {
    return(
      prod(sapply(
        X=1:ltarget,
        FUN = function(j) sum( sapply(
          X=1:lsource,
          FUN=function(i) tau_kl(j,i,prefix,lsource) ) )
      )))
  }

  tau = function(prefix) {
    return(
      sum(
        sapply(
          X=length(prefix):max.length,
          FUN=function(lsource) senlength.model[lsource+1,ltarget]*tau_l(prefix,lsource)
        )
      )
    )
  }

  scoreIBM2 = function(prefix) {
    res = tau(prefix)
    for (i in 1:length(prefix)) {
      if (i==1) {
        res = res*kgrams::probability(prefix[1] %|% "<NULL>", model=language.model)
      } else {
        res = res*kgrams::probability(prefix[i] %|% paste0(prefix[1:(i-1)], collapse=" "), model=language.model)
      }
    }
    return(res)
  }

  # set up queues
  all_Q = vector(mode = "list", length = max.length+1)
  for (k in 1:(max.length+1)) {
    all_Q[[k]] = priority_queue()
  }
  keepers = priority_queue()

  # start with null hypothesis
  H0 = list("length"=0,"prefix"=NULL, "score"=threshold[2])
  all_Q[[1]]$push(H0, priority=threshold[2])

  # build up the stacks
  thresh = threshold[2]
  for (k in 1:(max.length+1)) {
    if (k == 3) thresh = threshold[3]
    nsen = 0
    while( (all_Q[[k]]$size() > 0)  &  (nsen<max.nsen)  ) {
      H = all_Q[[k]]$pop()
      for (w in promising_words) {
        Hnew = list("length"=H$length+1,
                    "prefix"=c(H$prefix,w),
                    "score"=scoreIBM2( c(H$prefix,w) ) )
        if (Hnew$score >= thresh) {
          all_Q[[Hnew$length+1]]$push(Hnew, priority=Hnew$score); nsen = nsen+1

          Hnew$finalscore = senlength.model[Hnew$length + 1,ltarget]*
            tau_l(Hnew$prefix,Hnew$length)*
            kgrams::probability(paste0(c("<NULL>", Hnew$prefix), collapse=" "), model=language.model)
          keepers$push(Hnew,priority=Hnew$finalscore)
        }
        if (nsen>=max.nsen) break
      }

    }
  }

  #return(paste0(keepers$pop()$prefix,collapse=" "))
  return(keepers$as_list())

} # decode.IBM2


