

#' IBM1 Model
#'
#' The simplest SMT model from Brown et al. (1993)
#' @param target vector of sentences in language we want to translate to. Function assumes sentences are space-delimited.
#' @param source vector of sentences in language we want to translate from. Function assumes sentences are space-delimited.
#' @param maxiter max number of EM iterations allowed
#' @param eps convergence criteria for perplexity (i.e. negative log-likelihood)
#' @param add.null.token If TRUE (default), adds <NULL> to beginning of each source sentence. Allows target words to be aligned with "nothing".
#' @param init.tmatrix tmatrix from a previous estimation. If not provided, algorithm starts with uniform probabilities.
#' @param verbose If >=1, shows progress bar which updates every `verbose` steps, plus a summary when each iteration is complete. If 0.5 (default), only shows the summary without progress bars. If 0, shows nothing.
#' @return
#'    \item{tmatrix}{Environment object containing translation probabilities for target-source word pairs. E.g. tmatrix$go$va (equivalently, tmatrix[["go"]][["va"]]) gives the probability of target="go" given source="va".}
#'    \item{numiter}{Number of iterations}
#'    \item{maxiter}{The `maxiter` argument supplied by the user.}
#'    \item{eps}{The `eps` argument supplied by the user.}
#'    \item{converged}{TRUE if algorithm stopped once eps criteria met. FALSE otherwise.}
#'    \item{perplexity}{Final likelihood/perplexity value.}
#'    \item{time_elapsed}{Time in minutes the algorithm ran for.}
#'    \item{corpus}{data frame containing the target and source sentences and their lengths}
#' @examples
#' # download english-french sentence pairs
#' temp = tempfile()
#' download.file("http://www.manythings.org/anki/fra-eng.zip",temp);
#' ENFR = readr::read_tsv(file=unz(temp,"fra.txt"),col_names=c("en","fr","details"));
#' unlink(temp);
#'
#' # a bit of pre-processing
#' e = tolower(stringr::str_squish(tm::removePunctuation(ENFR$en[1:1000])));
#' f = tolower(stringr::str_squish(tm::removePunctuation(ENFR$fr[1:1000])));
#'
#' # models
#' initmodel  = SMT::IBM1(e,f,maxiter=5,eps=0.01);
#' finalmodel = SMT::IBM1(e,f,maxiter=40,eps=0.01,init.tmatrix=initmodel$tmatrix);
#' @import progress
#' @export
IBM1 = function(target,source,maxiter=30,eps=0.01,add.null.token=TRUE,init.tmatrix=NULL,verbose=0.5) {

  # error checking
  stopifnot(maxiter>=0,eps>0,add.null.token %in% c(TRUE,FALSE), (verbose %in% c(0,0.5)) | (verbose>=1),
    is.null(init.tmatrix) | is.environment(init.tmatrix))


  start_time = Sys.time()

  # add NULL token
  if (add.null.token) source = paste0("<NULL> ",source)

  # some things we'll need repeatedly
  # split sentences into words
  e_sentences = lapply(X=target,FUN=function(s) unlist(stringr::str_split(s, " ")))
  f_sentences = lapply(X=source,FUN=function(s) unlist(stringr::str_split(s, " ")))

  # list of all unique words
  e_allwords = unique(unlist(stringr::str_split(target, pattern=" ")))
  f_allwords = unique(unlist(stringr::str_split(source, pattern=" ")))
  n = length(e_sentences); n_eword = length(e_allwords); n_fword = length(f_allwords)
  # perplexity calculation stuff
  target_lengths = sapply(X=1:n, FUN=function(i) length(e_sentences[[i]]) )
  source_lengths = sapply(X=1:n, FUN=function(i) length(f_sentences[[i]]) )
  perp_const = sum(target_lengths*log(source_lengths))

  # initialize matrices
  t_e_f = new.env()
  if (!is.null(init.tmatrix)) {
    for (a in ls(init.tmatrix)) {
      t_e_f[[a]] = new.env()
      for (b in ls(init.tmatrix[[a]])) {
        t_e_f[[a]][[b]] = init.tmatrix[[a]][[b]]
      }
    }
  }
  c_e_f = new.env()
  s_total = new.env()
  total_f = new.env()
  for ( fword in f_allwords ) {
    total_f[[fword]] = 0
  }
  for ( eword in e_allwords ) {
    if (is.null(init.tmatrix)) t_e_f[[eword]] = new.env()
    c_e_f[[eword]] = new.env()
  }
  for (i in 1:n) {
    for (eword in e_sentences[[i]]) {
      for (fword in f_sentences[[i]]) {
        if (is.null(init.tmatrix)) t_e_f[[eword]][[fword]] = 1/n_eword
        c_e_f[[eword]][[fword]] = 0
      }
    }
  }

  # initial setup message
  time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
  if(verbose>=0.5) print(paste0("initial setup ;;; time elapsed: ",time_elapsed,"min"))

  # EM algorithm
  iter = 1
  prev_perplex = 0
  total_perplex = Inf
  while (iter<=maxiter & abs(total_perplex - prev_perplex)>eps) {

    # E Step
    if (verbose>=1) {
      pb = progress_bar$new(total=n,clear=TRUE,format=paste0("iter: ",iter," (E-step) [:bar] :current/:total (eta: :eta)")  )
      pb$tick(0)
    }
    for (s in 1:n) {

      # normalization
      for (eword in e_sentences[[s]]) {
        s_total[[eword]] = 0
        for (fword in f_sentences[[s]]) {
          s_total[[eword]] = s_total[[eword]] + t_e_f[[eword]][[fword]]
        }
      }

      # counts
      for (eword in e_sentences[[s]]) {
        for (fword in f_sentences[[s]]) {
          c_e_f[[eword]][[fword]] = c_e_f[[eword]][[fword]] + t_e_f[[eword]][[fword]]/s_total[[eword]]
          total_f[[fword]] = total_f[[fword]] + t_e_f[[eword]][[fword]]/s_total[[eword]]
        }
      }

      # progress
      if(verbose>=1 & s%%verbose==0) pb$tick(verbose)

    } # end for s in 1:n
    if (verbose>=1) pb$terminate()

    # M Step, and reset counts to 0
    if (verbose>=1) {
      pb = progress_bar$new(total=n_eword+n_fword,clear=TRUE,format=paste0("iter: ",iter," (M-step) [:bar] :current/:total (eta: :eta)")  )
      pb$tick(0); k = 1
    }
    for (eword in e_allwords) {
      for (fword in ls(t_e_f[[eword]])) {
        t_e_f[[eword]][[fword]] = c_e_f[[eword]][[fword]] / total_f[[fword]]
        c_e_f[[eword]][[fword]] = 0
      }

      if(verbose>=1){if(k%%verbose==0) pb$tick(verbose)}
      if(verbose>=1) k = k+1
    }
    for (fword in f_allwords) {
      total_f[[fword]] = 0

      if(verbose>=1){if(k%%verbose==0) pb$tick(verbose)}
      if(verbose>=1) k = k+1
    }
    if (verbose>=1) pb$terminate()

    # compute perplexity
    if (verbose>=1) {
      pb = progress_bar$new(total=n,clear=TRUE,format=paste0("iter: ",iter," (perplexity) [:bar] :current/:total (eta: :eta)")  )
      pb$tick(0)
    }
    prev_perplex = total_perplex
    total_perplex = perp_const
    for (s in 1:n) {
      for (eword in e_sentences[[s]]) {
        tmp = 0
        for (fword in f_sentences[[s]]) {
          tmp = tmp + t_e_f[[eword]][[fword]]
        }
        total_perplex = total_perplex - log(tmp)
      }

      if(verbose>=1 & s%%verbose==0) pb$tick(verbose)
    }
    if (verbose>=1) pb$terminate()

    time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
    if(verbose>=0.5) print(paste0("iter: ",iter,"; perplexity value: ",total_perplex, "; time elapsed: ",time_elapsed,"min"))

    iter = iter + 1

  } # end while
  corpus = data.frame(target, source, target_lengths, source_lengths)
  retobj = list(
    "tmatrix"=t_e_f,
    "numiter"=iter-1,
    "maxiter"=maxiter,
    "eps"=eps,
    "converged"=abs(total_perplex - prev_perplex)<=eps,
    "perplexity"=total_perplex,
    "time_elapsed"=time_elapsed,
    "corpus"=corpus
  )
  class(retobj) = "IBM1"
  return(retobj)

} # end function IBM1








