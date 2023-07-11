

#' IBM1 Model
#'
#' The simplest SMT model from Brown et al. (1993)
#' @param e vector of sentences in language we want to translate to
#' @param f vector of sentences in language we want to translate from
#' @param maxiter max number of EM iterations allowed
#' @param eps convergence criteria for perplexity
#' @param sparse If FALSE (default), use base R matrices. If TRUE, use sparseMatrix from the Matrix package.
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
#' # try with and without sparseMatrix
#' test1 = SMT::IBM1(e,f,maxiter=200,eps=0.01,sparse=FALSE);
#' test2 = SMT::IBM1(e,f,maxiter=200,eps=0.01,sparse=TRUE);
#' @import Matrix
#' @export
IBM1 = function(e,f,maxiter=30,eps=0.01,sparse=FALSE) {

  start_time = Sys.time()

  # keep list of words in each sentence
  e_sentences = lapply(X=e,FUN=function(s) unlist(stringr::str_split(s, " ")))
  f_sentences = lapply(X=f,FUN=function(s) unlist(stringr::str_split(s, " ")))
  e_allwords = unique(unlist(stringr::str_split(e, pattern=" ")))
  f_allwords = unique(unlist(stringr::str_split(f, pattern=" ")))
  n = length(e_sentences); n_eword = length(e_allwords); n_fword = length(f_allwords)

  # initialize matrices
  if (sparse) {

    t_e_f = Matrix::sparseMatrix(i=NULL,j=NULL,dims=c(n_eword,n_fword),x=1)
    c_e_f = Matrix::sparseMatrix(i=NULL,j=NULL,dims=c(n_eword,n_fword),x=1)
    rownames(t_e_f) = e_allwords
    colnames(t_e_f) = f_allwords
    rownames(c_e_f) = e_allwords
    colnames(c_e_f) = f_allwords

    for (s in 1:n) {
      t_e_f[e_sentences[s][[1]],f_sentences[s][[1]]] = 1/n_eword
    }

  } else {

    t_e_f = matrix(1/n_eword,nrow=n_eword,ncol=n_fword)
    c_e_f = matrix(0,nrow=n_eword,ncol=n_fword)
    rownames(t_e_f) = e_allwords
    colnames(t_e_f) = f_allwords
    rownames(c_e_f) = e_allwords
    colnames(c_e_f) = f_allwords

  } # if (sparse)

  # initial setup message
  time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
  print(paste0("initial setup ;;; time elapsed: ",time_elapsed,"min"))

  # EM algorithm
  iter = 1
  prev_perplex = 0
  total_perplex = Inf
  while (iter<=maxiter & abs(total_perplex - prev_perplex)>eps) {

    for (i in 1:n) {
      e_wordfreq = table(e_sentences[i][[1]]); u_e_words = names(e_wordfreq)
      f_wordfreq = table(f_sentences[i][[1]]); u_f_words = names(f_wordfreq)

      # update count matrices
      tmp = t_e_f[u_e_words,u_f_words]
      tmp = tmp/rowSums(  tmp %*% diag(f_wordfreq)  )
      tmp = (tmp*as.vector(e_wordfreq)) %*% diag(f_wordfreq)
      c_e_f[u_e_words,u_f_words] = c_e_f[u_e_words,u_f_words] + tmp
    }

    # t probs
    if (sparse) {
      t_e_f[e_allwords,f_allwords] = c_e_f[e_allwords,f_allwords] %*% diag(1/colSums(c_e_f))
      attributes(t_e_f)$Dimnames = attributes(c_e_f)$Dimnames
    } else {
      t_e_f = c_e_f %*% diag(1/colSums(c_e_f))
      attributes(t_e_f) = attributes(c_e_f)
    }

    # compute perplexity
    perplex = rep(0,n)
    for (i in 1:n) {
      perplex[i] = sum(t_e_f[e_sentences[i][[1]],f_sentences[i][[1]]]) - length(e_sentences[i][[1]])*log(length(f_sentences[i][[1]]))
    }
    prev_perplex = total_perplex
    total_perplex = -sum(perplex)
    time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
    print(paste0("iter: ",iter,"; perplexity value: ",total_perplex, "; time elapsed: ",time_elapsed,"min"))

    iter = iter + 1

  } # end while

  return(list(
    "tmatrix"=t_e_f,
    "numiter"=iter-1,
    "maxiter"=maxiter,
    "eps"=eps,
    "converged"=abs(total_perplex - prev_perplex)<=eps,
    "perplexity"=total_perplex,
    "time_elapsed"=time_elapsed
    ))

} # end function IBM1









