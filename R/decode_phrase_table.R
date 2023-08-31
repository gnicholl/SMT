


#' Phrase-based translation stack decoder.
#'
#' Implements the stack decoder described in `Koehn (2009)` with help from the NLTK (python package) stack decoder code.
#' One key difference is that I don't currently implement future cost estimation,
#' which can allow better ranking of partial hypotheses.
#' Relies on three models: phrase translation table, sentence length model,
#' and a kgrams language model. If the latter two aren't provided, they are estimated
#' using a similar approach to `Wang and Waibel (1998)`: sentence length is estimated
#' using a poisson regression, and the language model is estimated using 3rd-order kgrams.
#' Heuristics (`max_stack_size`) are used to reduce computational time,
#' at the expense of an increased possibility of not finding the best translation.
#'
#' @param object result from build_phrase_table()
#' @param target.sentence sentence in f language we'd like to decode (character string)
#' @param max_stack_size (default 100) the maximum number of partial translations to keep from each stack. Lower number means faster algorithm, but more likely best translation is missed.
#' @param senlength.model Matrix such that `senlength.model[m,n]` is the probability that e sentence is length m given f sentence is length n. If not provided, it is estimated form the training corpus using a simple poisson regression: `glm(e_lengths ~ f_lengths, family="poisson")`.
#' @param language.model The result of a `kgrams::language_model` estimation. If not provided, we estimate a degree-3 kgrams model from the training corpus using the default "ml" (maximum likelihood) approach. See `kgrams` package for more details.
#' @return List of best translations found using the stack decoder.
#' @examples
#' # download english-french sentence pairs
#' temp = tempfile()
#' download.file("http://www.manythings.org/anki/fra-eng.zip",temp);
#' ENFR = readr::read_tsv(file=unz(temp,"fra.txt"),col_names=c("en","fr","details"));
#' unlink(temp);
#'
#' # a bit of pre-processing
#' e = removePunctuation(ENFR$en[1:10000])
#'   e = str_squish(e)
#'   e = tolower(e)
#' f = removePunctuation(ENFR$fr[1:10000])
#'   f = str_squish(f)
#'   f = tolower(f)
#'
#' # estimate models in both directions
#' model3_etof = IBM3(target=e,source=f,maxiter=5, init.IBM1=15, init.IBM2=15,  verbose=100)
#' model3_ftoe = IBM3(target=f,source=e,maxiter=5, init.IBM1=15, init.IBM2=15,  verbose=100)
#'
#' # extract phrases from the IBM models and compute phrase probabilities
#' phtable = build_phrase_table(model3_etof,model3_ftoe)

#' # possible english translations and their probabilities
#' best_translation = decode(phtable,target.sentence = phtable$corpus$target[[6000]]) # "laissemoi réfléchir"
#'   # returns "let me think"
#'
#' @import collections
#' @export
decode.phrase_table = function(object, target.sentence, max_stack_size=100,
                               senlength.model=NULL, language.model=NULL)
{

  # params
  target.sentence = unlist(stringr::str_split(target.sentence," "))
  ltarget = length(target.sentence)
  max.length = max(object$corpus$source_lengths[object$corpus$target_lengths==ltarget])

  # reverse order of phrase table
  pht_reverse = new.env()
  for (ph1 in ls(object$phrase_table)) {
    for (ph2 in ls(object$phrase_table[[ph1]])) {
      if (is.null(pht_reverse[[ph2]])) pht_reverse[[ph2]] = new.env()
      pht_reverse[[ph2]][[ph1]] = object$phrase_table[[ph1]][[ph2]]
    }
  }

  # sentence length model
  if (is.null(senlength.model)) {
    f_lengths = object$corpus$target_lengths; max.f.length = max(f_lengths)
    e_lengths = object$corpus$source_lengths; max.e.length = max(e_lengths)
    slmodel = glm(e_lengths ~ f_lengths, family="poisson")
    lambda = predict(slmodel, newdata=data.frame(f_lengths=1:max.f.length), type="response")
    senlength.model = matrix(0,nrow=max.e.length,ncol=max.f.length)
    for (i in 1:max.e.length) senlength.model[i,] = dpois(i, lambda=lambda)
  }

  # language model
  if (is.null(language.model)) {
    freqs = kgrams::kgram_freqs(object$corpus$source, N = 3, verbose = FALSE)
    language.model = kgrams::language_model(freqs)
  }

  # extract all possible phrases from target sentence
  possible_expansions = function(hypothesis) {
    alltrgphrases = list()
    for (i in 1:ltarget) {
      for (j in 1:ltarget) {
        if (j+i-1 <= ltarget) {
          phrase = paste0(target.sentence[j:(j+i-1)],collapse=" " )
          if (!is.null(pht_reverse[[phrase]])) {
            if (  !any( j:(j+i-1) %in% hypothesis$f_translated)  ) {
              alltrgphrases[[length(alltrgphrases)+1]] = c(j,j+i-1)
            }
          }
        }
      }
    }
    return(alltrgphrases)
  }

  # list of "promising phrases" -> those with +ve probabilities
  translations_for = function(f_phrase) {
    return(ls(pht_reverse[[f_phrase]]))
  }

  # predict future cost
  # future_score_table = vector(mode = "list", length = ltarget)
  # for (lseq in 1:ltarget) {
  #   for (start in 1:(ltarget-lseq+1)) {
  #     end = start + lseq
  #     phrase = paste0(target.sentence[start:end],collapse=" ")
  #     if ( !is.null( pht_reverse[[phrase]] ) ) {
  #       max_score = 0
  #       for (candidate in ls(pht_reverse[[phrase]])) {
  #         max_score = max(max_score,pht_reverse[[phrase]][[candidate]])
  #       }
  #       future_score_table[[start]][[end]] = max_score
  #     } else {
  #       future_score_table[[start]][[end]] = 0
  #     }
  #
  #     if (end-start > 1) {
  #       for ( mid in (start+1):(end-1) ) {
  #         newscore = future_score_table[[start]][[mid]] + future_score_table[[mid]][[end]]
  #         if (newscore > future_score_table[[start]][[end]]) future_score_table[[start]][[end]] = newscore
  #       }
  #     }
  #   }
  # }

  # scoring function
  tau = function(hypothesis) {
    score = 1

    nphrases = length(hypothesis$e_spans)
    for (i in 1:nphrases) {
      # phrase prob
      e_phrase = hypothesis$e_prefix[hypothesis$e_spans[[i]][1]:hypothesis$e_spans[[i]][2]]
        e_phrase = paste0(e_phrase,collapse=" ")
      f_phrase = target.sentence[hypothesis$f_spans[[i]][1]:hypothesis$f_spans[[i]][2]]
        f_phrase = paste0(f_phrase,collapse=" ")
      if (is.null(pht_reverse[[f_phrase]][[e_phrase]])) return(0)
      score = score*pht_reverse[[f_phrase]][[e_phrase]]

      # distortion prob
      if (i == 1) {
        score = score * 0.5^abs(hypothesis$f_spans[[i]][1] - 0 - 1)
      } else {
        score = score * 0.5^abs(hypothesis$f_spans[[i]][1] - hypothesis$f_spans[[i-1]][2] - 1)
      }
    }

    # language model
    for (i in 1:length(hypothesis$e_prefix)) {
      if (i==1) {
        langscore = kgrams::probability(hypothesis$e_prefix[1] %|% BOS(), model=language.model)
      } else {
        langscore = kgrams::probability(hypothesis$e_prefix[i] %|% paste0(c(BOS(),hypothesis$e_prefix[1:(i-1)]), collapse=" "), model=language.model)
      }
      if (is.na(langscore)) return(0)
      score = score*langscore
    }

    return(score)
  }



  # set up queues
  all_Q = vector(mode = "list", length =ltarget+1)
  for (k in 1:(ltarget+1)) {
    all_Q[[k]] = priority_queue()
  }
  finalQ = priority_queue()

  # start with empty hypothesis
  H0 = list(
    e_prefix = NULL,
    f_translated = NULL,
    e_spans = list(),
    f_spans = list(),
    score = 0
  )
  all_Q[[1]]$push(H0, priority=H0$score)

  # build up the stacks
  for (k in 1:(ltarget+1)) {
    iter = 1
    while( (all_Q[[k]]$size() > 0) & (iter <= max_stack_size) ) {
      H = all_Q[[k]]$pop(); iter = iter+1
      for (fspan in possible_expansions(H)) {
        f_phrase2 = target.sentence[fspan[1]:fspan[2]]
        f_phrase = paste0(f_phrase2 , collapse=" ")
        for (e_phrase in translations_for(f_phrase)) {
          e_phrase2 = unlist(stringr::str_split(e_phrase," "))

          Hnew = list()
          Hnew$e_prefix = c(H$e_prefix,e_phrase2)
          Hnew$f_translated = c(H$f_translated, fspan[1]:fspan[2])
          Hnew$e_spans = H$e_spans
            Hnew$e_spans[[length(Hnew$e_spans)+1]] = c(length(H$e_prefix)+1, length(Hnew$e_prefix))
          Hnew$f_spans = H$f_spans
            Hnew$f_spans[[length(Hnew$f_spans)+1]] = fspan
          Hnew$score = tau(Hnew)

          if(  (length(Hnew$e_prefix) <= max.e.length) & (Hnew$score > 0)  ) {
            all_Q[[length(Hnew$f_translated)+1]]$push(Hnew, priority=Hnew$score)
            Hnew$finalscore = Hnew$score*senlength.model[length(Hnew$e_prefix),ltarget]
            Hnew$finalscore = Hnew$finalscore*kgrams::probability(EOS() %|% paste0(c(BOS(),Hnew$e_prefix), collapse=" "), model=language.model)
            if (Hnew$finalscore>0) finalQ$push(Hnew, priority=Hnew$finalscore)
          }
        }
      }
    }
  }

  return(finalQ$as_list())

} # decod.phrase_table








