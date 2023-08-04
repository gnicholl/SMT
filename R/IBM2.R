
#' IBM2 Model
#'
#' The second SMT model from Brown et al. (1993)
#' @param e vector of sentences in language we want to translate to
#' @param f vector of sentences in language we want to translate from
#' @param maxiter max number of EM iterations allowed
#' @param eps convergence criteria for perplexity (i.e. negative log-likelihood)
#' @param init.IBM1 number of iterations (integer>=0) of IBM1 to perform to initialize IBM2 algorithm
#' @param add.null.token If TRUE (default), adds <NULL> to beginning of each f sentence. Allows e words to be aligned with "nothing".
#' @param init.tmatrix tmatrix from a previous estimation. If not provided, algorithm starts with uniform probabilities.
#' @param init.amatrix amatrix from a previous estimation. If not provided, algorithm starts with uniform probabilities.
#' @return
#'    \item{tmatrix}{Environment object containing translation probabilities for e-f word pairs. E.g. tmatrix$go$va (equivalently, tmatrix[["go"]][["va"]]) gives the probability of e="go" given f="va".}
#'    \item{amatrix}{A list of alignment probability matrices. E.g. amatrix[[3]][[4]] gives the alignment probability matrix for e sentences of length 3 and f sentences of length 4.}
#'    \item{numiter}{Number of iterations}
#'    \item{maxiter}{As above}
#'    \item{eps}{As above}
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
#' e = tolower(stringr::str_squish(tm::removePunctuation(ENFR$en[1:1000])));
#' f = tolower(stringr::str_squish(tm::removePunctuation(ENFR$fr[1:1000])));
#'
#' # IBM1 guarantees global optimum, while IBM2 doesn't => initialize with IBM1.
#' # IBM1 also faster, may make IBM2 faster overall.
#' # Experiment with number of initial IBM1 iterations:
#' test0 = SMT::IBM2(e,f,maxiter=50,eps=0.01,init.IBM1=0);
#' test1 = SMT::IBM2(e,f,maxiter=50,eps=0.01,init.IBM1=1);
#' test2 = SMT::IBM2(e,f,maxiter=50,eps=0.01,init.IBM1=2);
#' test3 = SMT::IBM2(e,f,maxiter=50,eps=0.01,init.IBM1=3);
#' @export
IBM2 = function(e,f,maxiter=30,eps=0.01,init.IBM1=10,add.null.token=TRUE,init.tmatrix=NULL,init.amatrix=NULL) {

  start_time = Sys.time()
  print(paste0("------running ",init.IBM1," iterations of IBM1-------"))
  out_IBM1 = IBM1_v2(e=e,f=f,maxiter=init.IBM1,eps=eps,add.null.token=add.null.token,init.tmatrix=init.tmatrix)
  print(paste0("------running ",maxiter," iterations of IBM2-------"))

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

  # initialize t matrix and counts
  t_e_f = out_IBM1$tmatrix
  c_e_f = new.env()
  s_total = new.env()
  total_f = new.env()
  for ( fword in f_allwords ) {
    total_f[[fword]] = 0
  }
  for ( eword in e_allwords ) {
    c_e_f[[eword]] = new.env()
  }
  for (i in 1:n) {
    for (eword in e_sentences[[i]]) {
      for (fword in f_sentences[[i]]) {
        c_e_f[[eword]][[fword]] = 0
      }
    }
  }

  # initialize a matrix and counts
  s_lengths = data.frame(
    e_lengths = sapply(X=1:n, FUN=function(i) length(e_sentences[[i]]) ),
    f_lengths = sapply(X=1:n, FUN=function(i) length(f_sentences[[i]]) )
  )
  combos = unique(s_lengths)
  if (is.null(init.amatrix)) {
    aprob = list()
    for (k in unique(combos$e_lengths)) {
      aprob[[k]] = list()
    }
  } else {
    aprob = init.amatrix
  }
  acount = aprob
  for (k in 1:nrow(combos)) {
    le = combos$e_lengths[k]
    lf = combos$f_lengths[k]
    if (is.null(init.amatrix)) aprob[[le]][[lf]] = matrix(1/lf,nrow=le,ncol=lf)
    acount[[le]][[lf]] = matrix(0,nrow=le,ncol=lf)
  }

  # initial setup message
  time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
  print(paste0("initial setup ;;; time elapsed: ",time_elapsed,"min"))

  # EM algorithm
  iter = 1
  prev_perplex = 0
  total_perplex = Inf
  while (iter<=maxiter & abs(total_perplex - prev_perplex)>eps) {

    for (s in 1:n) {

      e_sen = e_sentences[[s]]; le = length(e_sen)
      f_sen = f_sentences[[s]]; lf = length(f_sen)

      # normalization
      for (j in 1:le) {
        eword = e_sen[j]
        s_total[[eword]] = 0
        for (i in 1:lf) {
          fword = f_sen[i]
          s_total[[eword]] = s_total[[eword]] + t_e_f[[eword]][[fword]] * aprob[[le]][[lf]][j,i]
        }
      }

      # update counts
      for (j in 1:le) {
        eword = e_sen[j]
        for (i in 1:lf) {
          fword = f_sen[i]
          c = t_e_f[[eword]][[fword]]*aprob[[le]][[lf]][j,i] / s_total[[eword]]
          c_e_f[[eword]][[fword]] = c_e_f[[eword]][[fword]] + c
          total_f[[fword]] = total_f[[fword]] + c
          acount[[le]][[lf]][j,i] = acount[[le]][[lf]][j,i] + c
        }
      }

    } # end for i in 1:n


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

    # update a probs, reset counts
    for (k in 1:nrow(combos)) {
      le = combos$e_lengths[k]
      lf = combos$f_lengths[k]
      aprob[[le]][[lf]][] = acount[[le]][[lf]] / Rfast::rowsums(acount[[le]][[lf]])
      acount[[le]][[lf]][] = 0
    }

    # compute perplexity
    prev_perplex = total_perplex
    total_perplex = 0
    for (s in 1:n) {
      e_sen = e_sentences[[s]]; le = length(e_sen)
      f_sen = f_sentences[[s]]; lf = length(f_sen)

      for (j in 1:le) {
        eword = e_sen[j]
        tmp = 0
        for (i in 1:lf) {
          fword = f_sen[i]
          tmp = tmp + t_e_f[[eword]][[fword]]*aprob[[le]][[lf]][j,i]
        }
        total_perplex = total_perplex - log(tmp)
      }
    }


    time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
    print(paste0("iter: ",iter,"; perplexity value: ",total_perplex, "; time elapsed: ",time_elapsed,"min"))

    iter = iter + 1

  } # end while

  # return list of results
  retobj = list(
    "tmatrix"=t_e_f,
    "amatrix"=aprob,
    "numiter"=iter-1,
    "maxiter"=maxiter,
    "eps"=eps,
    "converged"=abs(total_perplex - prev_perplex)<=eps,
    "perplexity"=total_perplex,
    "time_elapsed"=time_elapsed
  )
  class(retobj) = "IBM2"
  return(retobj)

} # end function IBM2









