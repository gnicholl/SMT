
#' IBM2 Model
#'
#' The second SMT model from Brown et al. (1993)
#' @param e vector of sentences in language we want to translate to
#' @param f vector of sentences in language we want to translate from
#' @param maxiter max number of EM iterations allowed
#' @param eps convergence criteria for perplexity
#' @param init.IBM1 number of iterations of IBM1 to perform to initialize IBM2
#' @param sparse If FALSE (default), use base R matrices. If TRUE, use sparseMatrix from the Matrix package.
#' @param fmatch If TRUE, use fmatch from fastmatch package for faster lookup. Otherwise use base R lookup.
#' @return
#'    \item{tmatrix}{Matrix of translation probabilities (cols are words from e, rows are words from f). If sparse=TRUE, tmatrix will be a sparseMatrix from the Matrix package, and will generally take up substantially less memory.}
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
#' e = tolower(stringr::str_squish(tm::removePunctuation(ENFR$en[1:200])));
#' f = tolower(stringr::str_squish(tm::removePunctuation(ENFR$fr[1:200])));
#'
#' # try with and without sparseMatrix (slow match)
#' test1 = SMT::IBM1(e,f,maxiter=200,eps=0.01,sparse=FALSE,fmatch=FALSE);
#' test2 = SMT::IBM1(e,f,maxiter=200,eps=0.01,sparse=TRUE,fmatch=FALSE);
#'
#' # try with and without sparseMatrix (fast match)
#' test3 = SMT::IBM1(e,f,maxiter=200,eps=0.01,sparse=FALSE,fmatch=TRUE);
#' test4 = SMT::IBM1(e,f,maxiter=200,eps=0.01,sparse=TRUE,fmatch=TRUE);
#' @importFrom fastmatch fmatch
#' @import Matrix
#' @export
IBM2 = function(e,f,maxiter=30,eps=0.01,init.IBM1=2,sparse=FALSE,fmatch=FALSE,cl=NULL) {

  start_time = Sys.time()
  print(paste0("------running ",init.IBM1," iterations of IBM1-------"))
  out_IBM1 = IBM1(e=e,f=f,maxiter=init.IBM1,eps=eps,sparse=sparse,fmatch=fmatch,cl=cl)
  print("------now run IBM2-------")

  # some things we'll need repeatedly
  # split sentences into words
  if (is.null(cl)) {
    e_sentences = lapply(X=e,FUN=function(s) unlist(stringr::str_split(s, " ")))
    f_sentences = lapply(X=f,FUN=function(s) unlist(stringr::str_split(s, " ")))
  } else {
    e_sentences = parallel::parLapply(X=e,fun=function(s) unlist(stringr::str_split(s, " ")), cl=cl)
    f_sentences = parallel::parLapply(X=f,fun=function(s) unlist(stringr::str_split(s, " ")), cl=cl)
  }
  # list of all unique words
  e_allwords = unique(unlist(stringr::str_split(e, pattern=" ")))
  f_allwords = unique(unlist(stringr::str_split(f, pattern=" ")))
  n = length(e_sentences); n_eword = length(e_allwords); n_fword = length(f_allwords)
  # placeholder for perplexity
  perplex = rep(0,n)

  # initialize t matrix and counts
  if (sparse) {

    t_e_f = sparseMatrix(i=NULL,j=NULL,dims=c(n_eword,n_fword),x=1)
    c_e_f = sparseMatrix(i=NULL,j=NULL,dims=c(n_eword,n_fword),x=1)
    rownames(t_e_f) = e_allwords
    colnames(t_e_f) = f_allwords
    rownames(c_e_f) = e_allwords
    colnames(c_e_f) = f_allwords

    for (i in 1:n) {
      t_e_f[e_sentences[i][[1]],f_sentences[i][[1]]] = 1/n_eword
    }

  } else {

    t_e_f = out_IBM1$tmatrix
    c_e_f = matrix(0,nrow=n_eword,ncol=n_fword)
    rownames(c_e_f) = e_allwords
    colnames(c_e_f) = f_allwords

  } # if (sparse)

  # initialize a matrix and counts
  s_lengths = data.frame(
    e_lengths = sapply(X=1:n, FUN=function(i) length(e_sentences[i][[1]]) ),
    f_lengths = sapply(X=1:n, FUN=function(i) length(f_sentences[i][[1]]) )
  )
  combos = unique(s_lengths)
  aprob = list()
  for (k in unique(combos$e_lengths)) {
    aprob[[k]] = list()
  }
  acount = aprob
  for (k in 1:nrow(combos)) {
    le = combos$e_lengths[k]
    lf = combos$f_lengths[k]
    aprob[[le]][[lf]] = matrix(1/lf,nrow=le,ncol=lf)
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

    for (i in 1:n) {

      e_sen = e_sentences[i][[1]]; le = length(e_sen); u_e_words = unique(e_sen)
      f_sen = f_sentences[i][[1]]; lf = length(f_sen); u_f_words = unique(f_sen)

      # update count matrices
      if (!fmatch) {
        tmp = t_e_f[e_sen,f_sen,drop=FALSE]
        tmp = tmp*aprob[[le]][[lf]]
        tmp = tmp / rowSums(  tmp  )
        acount[[le]][[lf]] = acount[[le]][[lf]] + tmp
        c_e_f[u_e_words,u_f_words] = c_e_f[u_e_words,u_f_words] + biaggregate(tmp)[u_e_words,u_f_words]

      } else {
        match_e = fmatch(u_e_words,e_allwords)
        match_f = fmatch(u_f_words,f_allwords)
        tmp = t_e_f[match_e,match_f]
        d = diag(f_wordfreq)
        tmp = tmp/rowSums(  tmp %*% d  )
        tmp = (tmp*as.vector(e_wordfreq)) %*% d
        c_e_f[match_e,match_f] = c_e_f[match_e,match_f] + tmp
      }
    } # end for i in 1:n

    # update probs, reset counts
    if (sparse) {
      t_e_f[] = c_e_f %*% Diagonal(x=1/colSums(c_e_f))
      c_e_f[] = sparseMatrix(i=NULL,j=NULL,dims=c(n_eword,n_fword),x=1)
    } else {
      # compute tprob and reset count
      t_e_f[] = c_e_f %*% diag(1/colSums(c_e_f))
      c_e_f[] = matrix(0,nrow=n_eword,ncol=n_fword)

      # compute aprob and reset count
      for (k in 1:nrow(combos)) {
        le = combos$e_lengths[k]
        lf = combos$f_lengths[k]
        aprob[[le]][[lf]] = acount[[le]][[lf]] / rowSums(acount[[le]][[lf]])
        acount[[le]][[lf]] = matrix(0,nrow=le,ncol=lf)
      }
    }

    # compute perplexity
    if (!fmatch) {
      if (is.null(cl)) {
        for (i in 1:n) {
          e_sen = e_sentences[i][[1]]; le = length(e_sen)
          f_sen = f_sentences[i][[1]]; lf = length(f_sen)
          perplex[i] = sum(log(rowSums(t_e_f[e_sen,f_sen]*aprob[[le]][[lf]])))
        }
      } else {
        clusterExport(cl=cl,c('t_e_f','e_sentences','f_sentences'))
        perplex = parallel::parSapply(
          X=1:n,
          FUN=function(i) sum(t_e_f[e_sentences[i][[1]], f_sentences[i][[1]]] )
          , cl=cl)
      } # if is.null
    } else {
      if (is.null(cl)) {
        for (i in 1:n) {
          perplex[i] = sum(t_e_f[fmatch(e_sentences[i][[1]],e_allwords),
                                 fmatch(f_sentences[i][[1]],f_allwords)] )
        }
      } else {
        clusterExport(cl=cl,c('t_e_f','e_sentences','f_sentences','e_allwords','f_allwords'))
        perplex = parallel::parSapply(
          X=1:n,
          FUN=function(i) sum(t_e_f[fmatch(e_sentences[i][[1]],e_allwords),
                                    fmatch(f_sentences[i][[1]],f_allwords)] )
          , cl=cl )
      } # if is.null
    } # if !fmatch
    prev_perplex = total_perplex
    total_perplex = -sum(perplex)
    time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
    print(paste0("iter: ",iter,"; perplexity value: ",total_perplex, "; time elapsed: ",time_elapsed,"min"))

    iter = iter + 1

  } # end while

  return(list(
    "tmatrix"=t_e_f,
    "amatrix"=aprob,
    "numiter"=iter-1,
    "maxiter"=maxiter,
    "eps"=eps,
    "converged"=abs(total_perplex - prev_perplex)<=eps,
    "perplexity"=total_perplex,
    "time_elapsed"=time_elapsed
  ))

} # end function IBM2









