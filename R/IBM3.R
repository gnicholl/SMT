

#' IBM3 Model
#'
#' The third SMT model from Brown et al. (1993)
#' @param e vector of sentences in language we want to translate to
#' @param f vector of sentences in language we want to translate from
#' @param maxiter max number of EM iterations allowed
#' @param eps convergence criteria for perplexity (i.e. negative log-likelihood)
#' @param heuristic If TRUE (default) use a heuristic hill-climbing algorithm to find most likely alignments. If FALSE, search over all alignments (not recommended unless only looking at small sentences.)
#' @param maxfert Maximum number of e words ("fertility") which an f word is allowed to be mapped to. The default is 5. In practice fertility tends to be small, so this number should normally be sufficient.
#' @param init.IBM1 number of iterations (integer>=0) of IBM1 to perform to initialize IBM2 algorithm
#' @param init.IBM2 number of iterations (integer>=0) of IBM2 to perform to initialize IBM3 algorithm
#' @param sparse If TRUE, uses sparse matrices from Matrix package (default FALSE).
#' @return
#'    \item{tmatrix}{Matrix of translation probabilities (cols are words from e, rows are words from f). If sparse=TRUE, tmatrix will be a sparseMatrix from the Matrix package, and will generally take up substantially less memory.}
#'    \item{dmatrix}{A list of distortion probability matrices. E.g. dmatrix[[3]][[4]] gives the distortion probability matrix for e sentences of length 3 and f sentences of length 4.}
#'    \item{fertmatrix}{Matrix of fertility probabilities.}
#'    \item{p_null}{Probability of NULL token insertion.}
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
#' # a model
#' model3 = IBM3(e=e,f=f,maxiter=10, init.IBM1=10, init.IBM2=20)
#'
#' @importFrom Rfast colsums
#' @importFrom Rfast rowsums
#' @import index0
#' @import Matrix
#' @importFrom fastmatch fmatch
#' @export
IBM3 = function(e, f, maxiter=30, eps=0.01, heuristic=TRUE, maxfert=5, init.IBM1=5, init.IBM2=10, sparse=FALSE, fmatch=FALSE) {

  # initialize with IBM1 and IBM2
  start_time = Sys.time()
  out_IBM2 = IBM2(e=e,f=paste0("<NULL> ",f),maxiter=init.IBM2,init.IBM1=init.IBM1,sparse=sparse, fmatch=fmatch)
  print(paste0("------running ",maxiter," iterations of IBM3-------"))

  # set what functions to use
  if (!sparse) rowSums = function(...) rowsums(...)
  if (!sparse) colSums = function(...) colsums(...)

  ### HELPER FUNCTIONS #########################################################
    # collects set of neighbouring alignments
    neighbouring = function(a,jpegged) {
      N = matrix(0, nrow=(le-1)*(lf+1) + (le-1)*(le-2), ncol=le)
      Nitr = 1
      for (j in (1:le)[-jpegged]) {
        for (i in 0:lf) {
          new_a = a
          new_a[j] = i
          N[Nitr,] = new_a
          Nitr = Nitr + 1
        } # for i
      } # for j

      for (j1 in (1:le)[-jpegged]) {
        for (j2 in (1:le)[-c(jpegged,j1)]) {
          new_a = a
          new_a[j1] = new_a[j2]
          N[Nitr,] = new_a
          Nitr = Nitr + 1
        } # for j2
      } # for j1

      return(unique(N))
    }

    # finds best choice among all neighbours
    hillclimb = function(a,jpegged) {

      aold = -1
      while (!all(a == aold)) {
        aold = a
        neighbours = neighbouring(a,jpegged)
        for (r in 1:nrow(neighbours)) {
          if ( alignmentprob(neighbours[r,]) > alignmentprob(a) ) a = neighbours[r,]
        }# for
      }# while

      return(a)
    }# hillclimb

    # sampleIBM3: returns subset of alignments to iterate over
    sampleIBM3 = function(){
      f_sen_null = as.index0(c("<NULL>",f_sen))
      A = NULL
      a = rep(0,le)
      for (j in 1:le) {
        for (i in 0:lf) {
          a[j] = i # pegging one alignment point
          for (j_prime in (1:le)[-j]) {
            pr_prev = 0
            for (i_prime in 0:lf) {
              if (fmatch) pr_align_IBM2 = out_IBM2$tmatrix[fmatch(e_sen[j_prime],tmat_rnames),fmatch(f_sen_null[i_prime],tmat_cnames)]*
                                           out_IBM2$amatrix[[le]][[lf+1]][j_prime,i_prime+1]
              if (!fmatch)  pr_align_IBM2 = out_IBM2$tmatrix[e_sen[j_prime],f_sen_null[i_prime]]*out_IBM2$amatrix[[le]][[lf+1]][j_prime,i_prime+1]
              if (pr_align_IBM2 > pr_prev) {
                a[j_prime] = i_prime
                pr_prev = pr_align_IBM2
              }
            }
          }
          a = hillclimb(a,j)
          A = rbind(A,a,neighbouring(a,j))
        } # for i
      } # for j

      return(unique(A))
    } # sampleIBM3

    # alignmentprob(a) - computes probability of given alignment a
    #   inherits: e_sen, f_sen, le, lf, p_null, fertmatrix, t_e_f, dprob
    alignmentprob = function(a) {
      f_sen_null = as.index0(c("<NULL>",f_sen))
      f_aj = as.index1(f_sen_null[a]) # f word to which each e word aligned
      phi = as.index0( sapply(X=0:lf, FUN=function(i) sum(a==i)) ) # fertilities

      # NULL insertion
      prob_null = choose(le-phi[0],phi[0]) * ( p_null^phi[0] ) * ( max(1 - p_null,0.0000001)^(le - 2*phi[0]) )

      # fertility
      if (!fmatch) fert = sapply(X=1:lf, FUN=function(i) fertmatrix[ f_sen[i], min(phi[i],maxfert)+1 ])
      if (fmatch)  fert = sapply(X=1:lf, FUN=function(i) fertmatrix[ fmatch(f_sen[i],f_allwords), min(phi[i],maxfert)+1 ])
      prob_fert = prod( factorial(phi[1:lf]) * fert )

      # translation
      if (!fmatch) prob_t = prod(sapply(X=1:le, FUN=function(j) t_e_f[ e_sen[j],f_aj[j] ] ))
      if (fmatch)  prob_t = prod(sapply(X=1:le, FUN=function(j) t_e_f[ fmatch(e_sen[j],tmat_rnames),  fmatch(f_aj[j],tmat_cnames) ] ))

      # distortion
      prob_d = prod(sapply(X=1:le, FUN=function(j) dprob[[le]][[lf]][j,a[j]+1] ))

      # output
      return(as.numeric(prob_null*prob_fert*prob_t*prob_d))
    }
  ##############################################################################



  ### INITIALIZE MATRICES ######################################################
    # split into words
    e_sentences = lapply(X=e,FUN=function(s) unlist(stringr::str_split(s, " ")))
    f_sentences = lapply(X=f,FUN=function(s) unlist(stringr::str_split(s, " ")))

    # list of all unique words
    e_allwords = unique(unlist(stringr::str_split(e, pattern=" ")))
    f_allwords = c("<NULL>", unique(unlist(stringr::str_split(f, pattern=" "))) )
    n = length(e_sentences); n_eword = length(e_allwords); n_fword = length(f_allwords)

    # init tmatrix
    t_e_f = out_IBM2$tmatrix # translation probabilities
    if (!sparse) c_e_f = matrix(0,nrow=n_eword,ncol=n_fword)         # num. times f word translated as e word
    if (sparse)  c_e_f = sparseMatrix(i=NULL,j=NULL,dims=c(n_eword,n_fword),x=1)
    rownames(c_e_f) = e_allwords; colnames(c_e_f) = f_allwords
    tmat_rnames = rownames(t_e_f); tmat_cnames = colnames(t_e_f)

    # init dmatrix
    s_lengths = data.frame(
      e_lengths = sapply(X=1:n, FUN=function(i) length(e_sentences[i][[1]]) ),
      f_lengths = sapply(X=1:n, FUN=function(i) length(f_sentences[i][[1]]) )
    )
    combos = unique(s_lengths)
    dprob = list()
    for (k in unique(combos$e_lengths)) {
      dprob[[k]] = list()
    }
    dcount = dprob
    for (k in 1:nrow(combos)) {
      le = combos$e_lengths[k]
      lf = combos$f_lengths[k]
      dprob[[le]][[lf]] = matrix(1/le,nrow=le,ncol=lf+1)
      dcount[[le]][[lf]] = matrix(0,nrow=le,ncol=lf+1)
    }

    # init NULL prob
    p_null = 0.5
    p1count = 0
    p0count = 0

    # init fertility matrix
    fertmatrix = matrix(1/(maxfert+1), nrow=n_fword, ncol=maxfert+1)
    fertcount = matrix(0, nrow=n_fword, ncol=maxfert+1)
    rownames(fertmatrix) = f_allwords
    rownames(fertcount) = f_allwords
  ##############################################################################


  time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
  print(paste0("initial setup ;;; time elapsed: ",time_elapsed,"min"))

  # EM algorithm
  iter = 1
  perplex_vec = rep(0,n)
  prev_perplex = 0
  total_perplex = Inf
  while (iter<=maxiter & abs(total_perplex - prev_perplex)>eps) { # until convergence

    ### E STEP #################################################################
      pb = progress::progress_bar$new(total=n,clear=FALSE,
        format=paste0("iteration ",iter," (:what) [:bar] :current/:total (eta: :eta)")  )
      for (k in 1:n) { # for all sentence pairs
        # extract kth sentence pair
        e_sen = e_sentences[k][[1]]; le = length(e_sen)
        f_sen = f_sentences[k][[1]]; lf = length(f_sen)

        # too many possible alignments -> get list of most likely ones from IBM2
        pb$tick(tokens=list(what="0 step; sampling alignments"))
        if (heuristic & le>1)  {
          A = sampleIBM3()
        } else {
          A = gtools::permutations(n=lf+1,r=le,v=0:lf,repeats.allowed=TRUE)
        }

        # compute probability of each alignment
        pb$tick(tokens=list(what="E step; compute align probs"))
        ctotal = sapply(X=1:nrow(A), FUN=function(r) alignmentprob(A[r,]))
        perplex_vec[k] = mean(ctotal)
        ctotal = ctotal / sum(ctotal)

        # update expected counts given probabilities
        pb$tick(tokens=list(what="E step; expected counts    "))
        for (r in 1:nrow(A)) {
          a = A[r,]
          f_sen_null = as.index0(c("<NULL>",f_sen))
          f_aj = as.index1(f_sen_null[a]) # f word to which each e word aligned
          phi = as.index0( sapply(X=0:lf, FUN=function(i) sum(a==i)) ) # fertilities

          for (j in 1:le) {
            if (!fmatch) c_e_f[e_sen[j],f_aj[j]]      = c_e_f[e_sen[j],f_aj[j]]      + ctotal[r]
            if (fmatch)  {
              rlookup = fmatch(e_sen[j],e_allwords)
              clookup = fmatch(f_aj[j], f_allwords)
              c_e_f[rlookup,clookup]      = c_e_f[rlookup,clookup]      + ctotal[r]
            }
            dcount[[le]][[lf]][j,a[j]+1] = dcount[[le]][[lf]][j,a[j]+1] + ctotal[r]
          }

          p1count = p1count + phi[0]*ctotal[r]
          p0count = p0count + (le - 2*phi[0])*ctotal[r]

          for (i in 0:lf) {
            if (!fmatch) fertcount[ f_sen_null[i], min(phi[i],maxfert)+1 ] = fertcount[ f_sen_null[i], min(phi[i],maxfert)+1 ] + ctotal[r]
            if (fmatch) {
              rlookup = fmatch(f_sen_null[i], f_allwords)
              fertcount[ rlookup, min(phi[i],maxfert)+1 ] = fertcount[ rlookup, min(phi[i],maxfert)+1 ] + ctotal[r]
            }
          }

        } # for a in A

        pb$tick()

      } # for k (all sentences)
    ############################################################################


    ### M STEP #################################################################
      pb$tick(tokens=list(what="M step; update counts      "))
      # translation probs
      tmp = colSums(c_e_f)
      t_e_f[,tmp>0] = c_e_f[,tmp>0] %*% diag(1/tmp[tmp>0])

      # distortion probs
      for (k in 1:nrow(combos)) {
        le = combos$e_lengths[k]
        lf = combos$f_lengths[k]
        if (le==1) {
          dprob[[le]][[lf]][] = 1
        } else {
          tmp = colSums(dcount[[le]][[lf]])
          dprob[[le]][[lf]][,tmp>0] = dcount[[le]][[lf]][,tmp>0] %*% diag(1/tmp[tmp>0])
        }
      }

      # fertility probs
      tmp = rowSums(fertcount)
      fertmatrix[tmp>0,] = fertcount[tmp>0,] / tmp[tmp>0]

      # NULL token probs
      p_null = p1count / (p1count + p0count)
    ############################################################################


    ### HOUSEKEEPING ###########################################################
      # fix NaNs
      t_e_f[is.infinite(t_e_f)] = 1; t_e_f[is.nan(t_e_f)] = 0

      # reinitialize to 0
      c_e_f[] = 0
      for (k in 1:nrow(combos)) {
        le = combos$e_lengths[k]
        lf = combos$f_lengths[k]
        dcount[[le]][[lf]][] = 0
      }
      fertcount[] = 0
      p1count = 0; p0count = 0

      # perplexity

      # iterate
      prev_perplex = total_perplex
      total_perplex = -sum(log(perplex_vec))
      time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
      print(paste0("iter: ",iter,"; perplexity value: ",total_perplex, "; time elapsed: ",time_elapsed,"min"))
      iter = iter + 1
    ############################################################################

  } # while not converged

  retobj = list(
    "tmatrix"=t_e_f,
    "dmatrix"=dprob,
    "fertmatrix"=fertmatrix,
    "p_null"=as.numeric(p_null),
    "numiter"=iter-1,
    "maxiter"=maxiter,
    "eps"=eps,
    "converged"=abs(total_perplex - prev_perplex)<=eps,
    "perplexity"=total_perplex,
    "time_elapsed"=time_elapsed
  )
  class(retobj) = "IBM3"
  return(retobj)

}# IBM3



