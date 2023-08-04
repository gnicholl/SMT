

#' IBM1 Model
#'
#' The simplest SMT model from Brown et al. (1993)
#' @param e vector of sentences in language we want to translate to
#' @param f vector of sentences in language we want to translate from
#' @param maxiter max number of EM iterations allowed
#' @param eps convergence criteria for perplexity (i.e. negative log-likelihood)
#' @param add.null.token If TRUE (default), adds <NULL> to beginning of each f sentence. Allows e words to be aligned with "nothing".
#' @param init.tmatrix tmatrix from a previous estimation. If not provided, algorithm starts with uniform probabilities.
#' @return
#'    \item{tmatrix}{Environment object containing translation probabilities for e-f word pairs. E.g. tmatrix$go$va (equivalently, tmatrix[["go"]][["va"]]) gives the probability of e="go" given f="va".}
#'    \item{numiter}{Number of iterations}
#'    \item{maxiter}{The `maxiter` argument supplied by the user.}
#'    \item{eps}{The `eps` argument supplied by the user.}
#'    \item{converged}{TRUE if algorithm stopped once eps criteria met. FALSE otherwise.}
#'    \item{perplexity}{Final likelihood/perplexity value.}
#'    \item{time_elapsed}{Time in minutes the algorithm ran for.}
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
#' # models
#' initmodel  = SMT::IBM1(e,f,maxiter=5,eps=0.01);
#' finalmodel = SMT::IBM1(e,f,maxiter=40,eps=0.01,init.tmatrix=initmodel$tmatrix);
#'
#' # try with and without sparseMatrix (fast match)
#' test3 = SMT::IBM1(e,f,maxiter=200,eps=0.01,sparse=FALSE,fmatch=TRUE);
#' test4 = SMT::IBM1(e,f,maxiter=200,eps=0.01,sparse=TRUE,fmatch=TRUE);
#' @export
IBM1 = function(e,f,maxiter=30,eps=0.01,add.null.token=TRUE,init.tmatrix=NULL) {

  start_time = Sys.time()

  # add NULL token
  if (add.null.token) f = paste0("<NULL> ",f)

  # some things we'll need repeatedly
  # split sentences into words
  e_sentences = lapply(X=e,FUN=function(s) unlist(stringr::str_split(s, " ")))
  f_sentences = lapply(X=f,FUN=function(s) unlist(stringr::str_split(s, " ")))

  # list of all unique words
  e_allwords = unique(unlist(stringr::str_split(e, pattern=" ")))
  f_allwords = unique(unlist(stringr::str_split(f, pattern=" ")))
  n = length(e_sentences); n_eword = length(e_allwords); n_fword = length(f_allwords)
  # perplexity calculation stuff
  e_lengths = sapply(X=1:n, FUN=function(i) length(e_sentences[[i]]) )
  f_loglengths = sapply(X=1:n, FUN=function(i) log(length(f_sentences[[i]])) )
  perp_const = sum(e_lengths*f_loglengths)

  # initialize matrices
  if (is.null(init.tmatrix))  t_e_f = new.env()
  if (!is.null(init.tmatrix)) t_e_f = init.tmatrix
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
  print(paste0("initial setup ;;; time elapsed: ",time_elapsed,"min"))

  # EM algorithm
  iter = 1
  prev_perplex = 0
  total_perplex = Inf
  while (iter<=maxiter & abs(total_perplex - prev_perplex)>eps) {

    # E Step
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

    } # end for s in 1:n

    # M Step, and reset counts to 0
    for (eword in e_allwords) {
      for (fword in ls(t_e_f[[eword]])) {
        t_e_f[[eword]][[fword]] = c_e_f[[eword]][[fword]] / total_f[[fword]]
        c_e_f[[eword]][[fword]] = 0
      }
    }
    for (fword in f_allwords) {
      total_f[[fword]] = 0
    }

    # compute perplexity
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
    }

    time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
    print(paste0("iter: ",iter,"; perplexity value: ",total_perplex, "; time elapsed: ",time_elapsed,"min"))

    iter = iter + 1

  } # end while

  retobj = list(
    "tmatrix"=t_e_f,
    "numiter"=iter-1,
    "maxiter"=maxiter,
    "eps"=eps,
    "converged"=abs(total_perplex - prev_perplex)<=eps,
    "perplexity"=total_perplex,
    "time_elapsed"=time_elapsed
  )
  class(retobj) = "IBM1"
  return(retobj)

} # end function IBM1








