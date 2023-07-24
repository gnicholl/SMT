
#' IBM4 Model
#'
#' The fourth SMT model from Brown et al. (1993)
#' @param e vector of sentences in language we want to translate to
#' @param e_wordclass named vector where names are each unique e word and values are the "class" of each e word numbered from 1 to n_e_class where n_e_class is the total number of word classes.
#' @param f vector of sentences in language we want to translate from
#' @param f_wordclass named vector where names are each unique f word and values are the "class" of each f word numbered from 1 to n_f_class where n_f_class is the total number of word classes.
#' @param maxiter max number of EM iterations allowed (default 5)
#' @param eps convergence criteria for perplexity (i.e. negative log-likelihood)
#' @param heuristic If TRUE (default) use a heuristic hill-climbing algorithm to find most likely alignments. If FALSE, search over all alignments (not recommended unless only looking at small sentences.)
#' @param maxfert Maximum number of e words ("fertility") which an f word is allowed to be mapped to. The default is 5. In practice fertility tends to be small, so this number should normally be sufficient.
#' @param init.IBM1 number of iterations (integer>=0) of IBM1 to perform to initialize IBM2 algorithm (default 5)
#' @param init.IBM2 number of iterations (integer>=0) of IBM2 to perform to initialize IBM3 algorithm (default 5)
#' @param init.IBM3 number of iterations (integer>=0) of IBM3 to perform to initialize IBM4 algorithm (default 3)
#' @param sparse If TRUE, uses sparse matrices from Matrix package (default FALSE).
#' @return
#'    \item{tmatrix}{Matrix of translation probabilities (cols are words from e, rows are words from f). If sparse=TRUE, tmatrix will be a sparseMatrix from the Matrix package, and will generally take up substantially less memory.}
#'    \item{d1}{Array of relative distortion probabilities for the first word in a cept.}
#'    \item{dg1}{Array of relative distortion probabilities for subsequent words in a cept.}
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
#' e_all = unique(unlist(stringr::str_split(e, pattern=" ")))
#' f_all = unique(unlist(stringr::str_split(f, pattern=" ")))
#'
#' # get POS tagging models from udpipe package
#' m_eng = udpipe::udpipe_download_model(language = "english-ewt")
#' m_eng = udpipe::udpipe_load_model(m_eng)
#' m_fr = udpipe::udpipe_download_model(language = "french-gsd")
#' m_fr = udpipe::udpipe_load_model(m_fr)
#'
#' # annotate e and f words
#' en_annotate = udpipe::udpipe_annotate(m_eng, x = e_all) %>%
#'   as.data.frame() %>%
#'   dplyr::select(-sentence)
#' fr_annotate = udpipe::udpipe_annotate(m_fr, x = f_all) %>%
#'   as.data.frame() %>%
#'   dplyr::select(-sentence)
#'
#'# code nouns as 1, everything else (including missings) as 2
#' ecl = as.numeric(en_annotate$upos=="NOUN")+1
#'   names(ecl) = e_all
#'   ecl[is.na(ecl)] = 2
#' fcl = as.numeric(fr_annotate$upos=="NOUN")+1
#'   names(fcl) = f_all
#'   fcl[is.na(fcl)] = 2
#'
#' # run model
#' model4 = IBM4(
#'   e=e, e_wordclass=ecl, f=f, f_wordclass=fcl,
#'   maxiter=10, init.IBM1=10, init.IBM2=20, init.IBM3=10
#' )
#'
#' @importFrom Rfast colsums
#' @importFrom Rfast rowsums
#' @import index0
#' @import Matrix
#' @export
IBM4 = function(e, e_wordclass, f, f_wordclass, maxiter=5, eps=0.01, heuristic=TRUE, maxfert=5, init.IBM1=5, init.IBM2=5, init.IBM3=3, sparse=FALSE) {

  # check inputs
  if (any(is.na(e_wordclass))) stop("NAs present in e_wordclass")
  if (any(is.na(f_wordclass))) stop("NAs present in f_wordclass")
  if (any(as.integer(e_wordclass)!=e_wordclass)) stop("e_wordclass should be named vector of positive integers")
  if (any(as.integer(f_wordclass)!=f_wordclass)) stop("f_wordclass should be named vector of positive integers")
  if (any(e_wordclass<1)) stop("e_wordclass should be named vector of *positive* integers")
  if (any(f_wordclass<1)) stop("f_wordclass should be named vector of *positive* integers")
  if (heuristic==FALSE) warning("heuristic=FALSE means all possible alignments are considered. This may be slow when corpus contains long sentences.",immediate. = TRUE)

  # split into words
  e_sentences = lapply(X=e,FUN=function(s) unlist(stringr::str_split(s, " ")))
  n = length(e_sentences)
  max_le = 0
  for (s in 1:n) max_le = max(max_le,length(e_sentences[[s]]))
  f_sentences = lapply(X=f,FUN=function(s) unlist(stringr::str_split(s, " ")))

  # list of all unique words
  e_allwords = unique(unlist(stringr::str_split(e, pattern=" ")))
  f_allwords = c("<NULL>", unique(unlist(stringr::str_split(f, pattern=" "))) )
  n_eword = length(e_allwords); n_fword = length(f_allwords)
  if (sparse==FALSE) warning(paste0("sparse=FALSE means the translation matrix with dimensions ",n_eword," by ",n_fword," will be stored in memory despite being mostly 0s. Set sparse=TRUE if memory is an issue."),immediate. = TRUE)

  # initialize with IBM1, IBM2, IBM3
  start_time = Sys.time()
  out_IBM3 = IBM3(e=e,f=f,maxiter=init.IBM3,init.IBM2=init.IBM2,init.IBM1=init.IBM1,sparse=sparse,fmatch=TRUE)
  print(paste0("------running ",maxiter," iterations of IBM4-------"))

  # set what functions to use
  if (!sparse) rowSums = function(...) rowsums(...)
  if (!sparse) colSums = function(...) colsums(...)

  ### HELPER FUNCTIONS #########################################################
  # collects set of neighbouring alignments
  neighbouring = function(a,jpegged) {
    N = NULL
    for (j in (1:le)[-jpegged]) {
      for (i in 0:lf) {
        new_a = a
        new_a[j] = i
        N = rbind(N,new_a)
      } # for i
    } # for j

    for (j1 in (1:le)[-jpegged]) {
      for (j2 in (1:le)[-c(jpegged,j1)]) {
        new_a = a
        new_a[j1] = new_a[j2]
        N = rbind(N,new_a)
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
        if (is.nan(alignmentprob(a))) {
          return("error")
        }
        if (is.nan(alignmentprob(neighbours[r,]))) {
          return("error")
        }
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
            pr_align_IBM3 = alignmentprobIBM3(a)
            if (pr_align_IBM3 > pr_prev) {
              a[j_prime] = i_prime
              pr_prev = pr_align_IBM3
            }
          }
        }
        a = hillclimb(a,j)
        if (all(a=="error")) return("error")
        A = rbind(A,neighbouring(a,j))
      } # for i
    } # for j

    return(unique(A))
  } # sampleIBM3

  # computes alignment probability based on IBM3 method (see IBM3() for documentation)
  alignmentprobIBM3 = function(a) {
    f_sen_null = as.index0(c("<NULL>",f_sen))
    f_aj = as.index1(f_sen_null[a]) # f word to which each e word aligned
    phi = as.index0( sapply(X=0:lf, FUN=function(i) sum(a==i)) ) # fertilities
    prob_null = choose(le-phi[0],phi[0]) * ( out_IBM3$p_null^phi[0] ) * ( max(1 - out_IBM3$p_null,0.0000001)^(le - 2*phi[0]) )
    fert = sapply(X=1:lf, FUN=function(i) out_IBM3$fertmatrix[ f_sen[i], min(phi[i],maxfert)+1 ])
    prob_fert = prod( factorial(phi[1:lf]) * fert )
    prob_t = prod(sapply(X=1:le, FUN=function(j) out_IBM3$tmatrix[ e_sen[j],f_aj[j] ] ))
    prob_d = prod(sapply(X=1:le, FUN=function(j) out_IBM3$dmatrix[[le]][[lf]][j,a[j]+1] ))

    return(as.numeric(prob_null*prob_fert*prob_t*prob_d))
  }

  # alignmentprob(a) - computes probability of given alignment a
  #   inherits: e_sen, f_sen, le, lf, p_null, fertmatrix, t_e_f, dprob
  alignmentprob = function(a) {
    f_sen_null = as.index0(c("<NULL>",f_sen))
    f_aj = as.index1(f_sen_null[a]) # f word to which each e word aligned
    phi = as.index0( sapply(X=0:lf, FUN=function(i) sum(a==i)) ) # fertilities

    # NULL insertion
    prob_null = choose(le-phi[0],phi[0]) * ( p_null^phi[0] ) * ( max(1 - p_null,0.0000001)^(le - 2*phi[0]) )

    # fertility
    fert = sapply(X=1:lf, FUN=function(i) fertmatrix[ f_sen[i], min(phi[i],maxfert)+1 ])
    prob_fert = prod( factorial(phi[1:lf]) * fert )

    # translation
    prob_t = prod(sapply(X=1:le, FUN=function(j) t_e_f[ e_sen[j],f_aj[j] ] ))

    # distortion
      # determine cept numbering
      cepts = NULL
      for (i in 1:lf) {
        cepts = rbind(cepts,a==i)
      }
      wcount = rowSums(cepts)
      ceptcount = 1
      ceptnumber = rep(0,lf)
      for (i in 1:lf) {
        if(wcount[i]>0) {
          ceptnumber[i] = ceptcount
          ceptcount = ceptcount+1
        }
      }

      # position distance for each e word
      distortion = function(j) {
        if (a[j]==0) return(1) # ignore NULL token

        ceptnumber = as.index0(c(0,ceptnumber))
        ceptj = ceptnumber[a[j]]
        if (j==which(ceptnumber[a]==ceptj)[1]) { # first position in cept
          prevcept = ceptj - 1
          if (prevcept==0) {
            return(d1array[ f_wordclass[nclass_f+1] , e_wordclass[e_sen[j]] ,  which(dmap==(j - 0))  ])
          } else {
            prevcept_positions = which(ceptnumber[a]==prevcept)
            prevcenter = ceiling(mean(prevcept_positions))
            return(d1array[ f_wordclass[f_sen[which(ceptnumber==prevcept)-1]]   , e_wordclass[e_sen[j]] ,  which(dmap==(j - prevcenter))  ])
          }
        } else { # position 2 or greater in cept
          for (m in 2:length(which(ceptnumber[a]==ceptj))) {
            if (j==which(ceptnumber[a]==ceptj)[m]) {
              prevposition = which(ceptnumber[a]==ceptj)[m-1]
              return(dg1array[ e_wordclass[e_sen[j]] , which(dmap==(j - prevposition))])
            }
          }
        } # endif
      } # distortion()

      # probability
      prob_d = prod(sapply(X=1:le, FUN=distortion) )

    # output
    return(as.numeric(prob_null*prob_fert*prob_t*prob_d))
  }
  ##############################################################################



  ### INITIALIZE MATRICES ######################################################
  # init tmatrix
  t_e_f = out_IBM3$tmatrix # translation probabilities
  if (!sparse) c_e_f = matrix(0,nrow=n_eword,ncol=n_fword)         # num. times f word translated as e word
  if (sparse)  c_e_f = sparseMatrix(i=NULL,j=NULL,dims=c(n_eword,n_fword),x=1)
  rownames(c_e_f) = e_allwords; colnames(c_e_f) = f_allwords

  # init distortion arrays
  nclass_f = max(f_wordclass)
  nclass_e = max(e_wordclass)
  d1array  = array(1/(2*(max_le)+1), dim=c(nclass_f+1, nclass_e, 2*(max_le)+1   )) # nclass_f+1 is NULL token!
  dg1array = array(1/(2*(max_le)+1), dim=c(nclass_e, 2*(max_le)+1   ))
  c_d1array  = array(0, dim=c(nclass_f+1, nclass_e, 2*(max_le)+1   ))
  c_dg1array = array(0, dim=c(nclass_e, 2*(max_le)+1   ))
  dmap = c((-max_le):-1 ,0,  1:(max_le))

  # init NULL prob
  p_null = max(out_IBM3$p_null,0.25)
  p1count = 0
  p0count = 0

  # init fertility matrix
  fertmatrix = out_IBM3$fertmatrix
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
    for (k in 1:n) { # for all sentence pairs

      # extract kth sentence pair
      e_sen = e_sentences[k][[1]]; le = length(e_sen)
      f_sen = f_sentences[k][[1]]; lf = length(f_sen)

      # too many possible alignments -> get list of most likely ones from IBM2
      if (heuristic & le>1)  {
        A = sampleIBM3()
        if (all(A=="error")) {
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
        }
      } else {
        A = gtools::permutations(n=lf+1,r=le,v=0:lf,repeats.allowed=TRUE)
      }

      # compute probability of each alignment
      ctotal = sapply(X=1:nrow(A), FUN=function(r) alignmentprob(A[r,]))
      perplex_vec[k] = mean(ctotal)
      ctotal = ctotal / sum(ctotal)

      # update expected counts given probabilities
      for (r in 1:nrow(A)) {
        a = A[r,]
        f_sen_null = as.index0(c("<NULL>",f_sen))
        f_aj = as.index1(f_sen_null[a]) # f word to which each e word aligned
        phi = as.index0( sapply(X=0:lf, FUN=function(i) sum(a==i)) ) # fertilities

        # translation counts
        for (j in 1:le) {
          c_e_f[e_sen[j],f_aj[j]] = c_e_f[e_sen[j],f_aj[j]] + ctotal[r]
        }

        # distortion counts
          # determine cept numbering
          cepts = NULL
          for (i in 1:lf) {
            cepts = rbind(cepts,a==i)
          }
          wcount = rowSums(cepts)
          ceptcount = 1
          ceptnumber = rep(0,lf)
          for (i in 1:lf) {
            if(wcount[i]>0) {
              ceptnumber[i] = ceptcount
              ceptcount = ceptcount+1
            }
          }

          # position distance for each e word
          ceptnumber = as.index0(c(0,ceptnumber))
          for (j in 1:le) {
            if (a[j] > 0) {
              ceptj = ceptnumber[a[j]]
              if (j==which(ceptnumber[a]==ceptj)[1]) { # if first word in cept
                prevcept = ceptj - 1
                if (prevcept==0) { #previous cept is <NULL> (i.e. this is first cept)
                  c_d1array[ f_wordclass[nclass_f+1] , e_wordclass[e_sen[j]] ,  which(dmap==(j - 0))  ] =
                    c_d1array[ f_wordclass[nclass_f+1] , e_wordclass[e_sen[j]] ,  which(dmap==(j - 0))  ] + ctotal[r]
                } else {
                  prevcept_positions = which(ceptnumber[a]==prevcept)
                  prevcenter = ceiling(mean(prevcept_positions))
                  fclass = f_wordclass[f_sen[which(ceptnumber==prevcept)-1]]
                  eclass = e_wordclass[e_sen[j]]
                  d = which(dmap==(j - prevcenter))
                  c_d1array[fclass,eclass,d] = c_d1array[fclass,eclass,d] + ctotal[r]
                }

              } else { # second or more word in cept
                for (m in 2:length(which(ceptnumber[a]==ceptj))) {
                  if (j==which(ceptnumber[a]==ceptj)[m]) {
                    prevposition = which(ceptnumber[a]==ceptj)[m-1]
                    eclass =  e_wordclass[e_sen[j]]
                    d = which(dmap==(j - prevposition))
                    c_dg1array[eclass,d] = c_dg1array[eclass,d] + ctotal[r]
                  }
                }
              } # if
            } # if a[j]>0
          } # end for

        # null token counts
        p1count = p1count + phi[0]*ctotal[r]
        p0count = p0count + (le - 2*phi[0])*ctotal[r]

        # fertility counts
        for (i in 0:lf) {
          fertcount[ f_sen_null[i], min(phi[i],maxfert)+1 ] = fertcount[ f_sen_null[i], min(phi[i],maxfert)+1 ] + ctotal[r]
        }

      } # for a in A

    } # for k (all sentences)
    ############################################################################


    ### M STEP #################################################################
    # translation probs
    tmp = colSums(c_e_f)
    t_e_f[,tmp>0] = c_e_f[,tmp>0] %*% diag(1/tmp[tmp>0])

    # distortion probs
    for (h in 1:nclass_e) {
      if (sum(c_dg1array[h,])>0) dg1array[h,] = c_dg1array[h,] / sum(c_dg1array[h,])
      for (m in 1:(nclass_f+1)) {
        if (sum(c_d1array[m,h,])>0) d1array[m,h,] = c_d1array[m,h,] / sum(c_d1array[m,h,])
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
    c_d1array[]  = 0
    c_dg1array[] = 0
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
    "d1"=d1array,"dg1"=dg1array,
    "fertmatrix"=fertmatrix,
    "p_null"=as.numeric(p_null),
    "numiter"=iter-1,
    "maxiter"=maxiter,
    "eps"=eps,
    "converged"=abs(total_perplex - prev_perplex)<=eps,
    "perplexity"=total_perplex,
    "time_elapsed"=time_elapsed
  )
  class(retobj) = "IBM4"
  return(retobj)

}# IBM4



