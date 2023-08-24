

#' IBM4 Model
#'
#' The fourth SMT model from Brown et al. (1993). Note that, unlike IBM1 and IBM2 functions,
#' NULL insertion is always assumed for IBM3 and 4, and cannot be turned off. Note also that I use
#' the basic hill climbing algorithm, and do not use pegging to increase the number
#' of alignments considered.
#' @param target vector of sentences in language we want to translate to. Function assumes sentences are space-delimited.
#' @param target_wordclass named vector where names are each unique target word and values are the "class" of each target word numbered from 1 to n_target_class where n_target_class is the total number of target word classes.
#' @param source vector of sentences in language we want to translate from. Function assumes sentences are space-delimited.
#' @param source_wordclass named vector where names are each unique source word and values are the "class" of each source word numbered from 1 to n_source_class where n_source_class is the total number of source word classes. (do NOT include <NULL> token)
#' @param maxiter max number of EM iterations allowed
#' @param eps convergence criteria for perplexity (i.e. negative log-likelihood)
#' @param heuristic If TRUE (default) use a heuristic hill-climbing algorithm to find most likely alignments. If FALSE, search over all alignments (not recommended unless only looking at short sentences). Sentences that are length 3 or smaller with always search over all alignments, even if heuristic=TRUE.
#' @param maxfert Maximum number of e words ("fertility") to which an f word can be aligned. The default is 5.
#' @param init.IBM1 number of iterations (integer>=0) of IBM1 to perform to initialize IBM2 algorithm
#' @param init.IBM2 number of iterations (integer>=0) of IBM2 to perform to initialize IBM3 algorithm
#' @param init.IBM3 number of iterations (integer>=0) of IBM3 to perform to initialize IBM4 algorithm
#' @param init.tmatrix tmatrix from a previous estimation. Used to initialize IBM1 (if `init.IBM1`>0), IBM2 (if `init.IBM2` >0), or IBM3 otherwise. If not provided, algorithm starts with uniform probabilities.
#' @param init.amatrix amatrix from a previous estimation. Used to initialize IBM2 if `init.IBM2` >0. If not provided, algorithm starts with uniform probabilities.
#' @param init.fmatrix fmatrix from a previous estimation. Used to initialize IBM3. If not provided, algorithm starts with uniform probabilities.
#' @param init.dmatrix dmatrix from a previous estimation. Used to initialize IBM3. If not provided, algorithm starts with uniform probabilities.
#' @param init.p_null  p_null  from a previous estimation. Used to initialize IBM3. If not provided, algorithm starts with `p_null=0.5`.
#' @param init.d1array d1array from a previous estimation. Used to initialize IBM4. If not provided, algorithm starts with uniform probabilities.
#' @param init.dg1array dg1array from a previous estimation. Used to initialize IBM4. If not provided, algorithm starts with uniform probabilities.
#' @param verbose If >=1, shows progress bar which updates every `verbose` steps, plus a summary when each iteration is complete. If 0.5 (default), only shows the summary without progress bars. If 0, shows nothing.
#' @return
#'    \item{tmatrix}{Environment object containing translation probabilities for target-source word pairs. E.g. tmatrix$go$va (equivalently, tmatrix[["go"]][["va"]]) gives the probability of target="go" given source="va".}
#'    \item{d1array}{An array}
#'    \item{dg1array}{An array}
#'    \item{fmatrix}{Environment object containing vectors of fertility probabilities for each source word. E.g. fmatrix[["va"]][1] is probability that "va" is aligned with 0 words. fmatrix[["va"]][3] is probability "va" is aligned with 2 words. Each vector is length `maxfert+1`.}
#'    \item{p_null}{Probability of NULL token insertion.}
#'    \item{numiter}{Number of iterations}
#'    \item{maxiter}{As above}
#'    \item{eps}{As above}
#'    \item{converged}{TRUE if algorithm stopped once eps criteria met. FALSE otherwise.}
#'    \item{perplexity}{Final likelihood/perplexity value.}
#'    \item{time_elapsed}{Time in minutes the algorithm ran for.}
#'    \item{corpus}{data frame containing the target and source sentences and their lengths}
#'    \item{prev_best_aligns}{list containing best alignments (i.e. "viterbi" alignments) for each target sentence}
#' @examples
#' # download english-french sentence pairs
#' temp = tempfile()
#' download.file("http://www.manythings.org/anki/fra-eng.zip",temp);
#' ENFR = readr::read_tsv(file=unz(temp,"fra.txt"),col_names=c("en","fr","details"));
#' unlink(temp);
#'
#' # a bit of pre-processing
#' e = tolower(stringr::str_squish(tm::removePunctuation(ENFR$en[1:10000])));
#' f = tolower(stringr::str_squish(tm::removePunctuation(ENFR$fr[1:10000])));
#'
#' # don't use any word classes, i.e. word class is 1 for all words
#' all_e_words = unique(unlist(lapply(X=e,FUN=function(s) unlist(stringr::str_split(s, " ")))))
#' all_f_words = unique(unlist(lapply(X=f,FUN=function(s) unlist(stringr::str_split(s, " ")))))
#' class1 = rep(1,length(all_e_words)); names(class1)=all_e_words
#' class2 = rep(1,length(all_f_words)); names(class2)=all_f_words
#'
#' # a model
#' model4 = IBM4(target=e,target_wordclass=class1,
#'                 source=f,source_wordclass=class2,
#'                 maxiter=5, init.IBM1=10, init.IBM2=10, init.IBM3=2,
#'                 verbose=100)
#'
#' @import progress
#' @export
IBM4 = function(target, source, target_wordclass, source_wordclass, maxiter=30, eps=0.01, heuristic=TRUE, maxfert=5,
                init.IBM1=10, init.IBM2=10, init.IBM3=10,
                init.tmatrix=NULL, init.amatrix=NULL, init.fmatrix=NULL, init.dmatrix=NULL,
                init.d1array=NULL, init.dg1array=NULL, init.p_null=NULL,
                verbose=0.5)
{

  # error checking
  samplemethod=1

  # initialize with IBM1 and IBM2 and IBM3
  start_time = Sys.time()
  out_IBM3 = IBM3(target=target,source=source,eps=eps,maxiter=init.IBM3,
                  init.IBM2=init.IBM2,init.IBM1=init.IBM1,
                  maxfert=maxfert, heuristic=heuristic,
                  init.tmatrix=init.tmatrix,init.amatrix=init.amatrix,
                  init.fmatrix=init.fmatrix,init.dmatrix=init.dmatrix,init.p_null=init.p_null,
                  verbose=verbose,samplemethod=samplemethod)
  if(verbose>=0.5) print(paste0("------running ",maxiter," iterations of IBM4-------"))

  # add NULL token
  e = target
  f = paste0("<NULL> ",source)

  ### HELPER FUNCTIONS #########################################################
  move = function(a, i, j) {
    a[j]=i
    return(a)
  }
  swap = function(a, j1, j2) {
    old_a = a
    a[j1] = old_a[j2]
    a[j2] = old_a[j1]
    return(a)
  }
  neighbourhood = function(a, le, lf) {
    N = matrix(0, nrow=le*lf + le*(le-1), ncol=le)
    Nit = 1
    for (j in 1:le) {
      for (i in 1:lf) {
        N[Nit,] = move(a, i, j)
        Nit = Nit+1
      }
    }
    for (j1 in 1:le) {
      for (j2 in (1:le)[-j1]) {
        N[Nit,] = swap(a,j1,j2)
        Nit = Nit+1
      }
    }
    return(unique(N))
  }
  viterbi_IBM2 = function(e_sen,f_sen) {
    le = length(e_sen)
    lf = length(f_sen)
    a = rep(NA,le)
    for (j in 1:le) {
      a[j] = which.max(sapply(X=1:lf, FUN=function(i) out_IBM2$tmatrix[[e_sen[j]]][[f_sen[i]]] * out_IBM2$amatrix[[le]][[lf]][j,i] ))
    }
    return(a)
  }
  pr_IBM3 = function(a, e_sen, f_sen) {
    le = length(e_sen); lf = length(f_sen)
    f_aj = f_sen[a] # f word to which each e word aligned
    phi = sapply(X=1:lf, FUN=function(i) sum(a==i))  # fertilities

    # NULL insertion
    prob_null = choose(le-phi[1],phi[1]) * ( p_null^phi[1] ) * ( max(1 - p_null,0.0000001)^(le - 2*phi[1]) )
    if (prob_null==0) return(0)

    # fertility
    fert = sapply(  X=2:lf, FUN=function(i) fertmatrix[[f_sen[i]]][min(phi[i],maxfert)+1] )
    prob_fert = prod( factorial(phi[2:lf]) * fert )
    if (prob_fert==0) return(0)

    # translation
    prob_t = prod(sapply(X=1:le, FUN=function(j) t_e_f[[ e_sen[j] ]][[ f_aj[j] ]] ))
    if (prob_t==0) return(0)

    # distortion
    prob_d = prod(sapply(X=1:le, FUN=function(j) dprob[[le]][[lf]][j,a[j]] ))
    if (prob_d==0) return(0)

    # output
    return(as.numeric(prob_null*prob_fert*prob_t*prob_d))
  }
  pr_IBM4 = function(a, e_sen, f_sen) {
    le = length(e_sen); lf = length(f_sen)
    f_aj = f_sen[a] # f word to which each e word aligned
    phi = sapply(X=1:lf, FUN=function(i) sum(a==i))  # fertilities

    # NULL insertion
    prob_null = choose(le-phi[1],phi[1]) * ( p_null^phi[1] ) * ( max(1 - p_null,0.0000001)^(le - 2*phi[1]) )
    if (prob_null==0) return(0)

    # fertility
    fert = sapply(  X=2:lf, FUN=function(i) fertmatrix[[f_sen[i]]][min(phi[i],maxfert)+1] )
    prob_fert = prod( factorial(phi[2:lf]) * fert )
    if (prob_fert==0) return(0)

    # translation
    prob_t = prod(sapply(X=1:le, FUN=function(j) t_e_f[[ e_sen[j] ]][[ f_aj[j] ]] ))
    if (prob_t==0) return(0)

    # distortion #############################
      # determine cept numbering
      ceptnumber = rep(0,lf)
      ceptcount = 1
      for (i in 1:lf) {
        if(phi[i]>0) {
          ceptnumber[i] = ceptcount
          ceptcount = ceptcount+1
        }
      }

      # position distance for each e word
      distortion = function(j) {
        if (a[j]==1) return(1) # ignore words aligned with NULL token

        t_class = target_wordclass[fastmatch::fmatch(e_sen[j],lookup_tclass)]
        ceptj = ceptnumber[a[j]]

        if (j==which(ceptnumber[a]==ceptj)[1]) { # first position in cept

          prevcept = ceptj - 1
          if(ceptj==1) { # if in the first cept, treat previous cept as <NULL> cept without distortion
            return(d1array[ nclass_s+1   , t_class ,  which(dmap==(j - 0))  ])
          } else {
            fwordaligned = f_sen[which(ceptnumber==prevcept)]
            if (fwordaligned=="<NULL>") { # previous cept is <NULL>
              s_class = nclass_s+2
            } else { # previous cept is non-NULL word
              s_class = source_wordclass[ fastmatch::fmatch(fwordaligned, lookup_sclass) ]
            }
            prevcept_positions = which(ceptnumber[a]==prevcept)
            prevcenter = ceiling(mean(prevcept_positions))
            return(d1array[ s_class   , t_class ,  which(dmap==(j - prevcenter))  ])
          }


        } else { # position 2 or greater in cept

          for (m in 2:phi[a[j]] ) {
            if (j==which(ceptnumber[a]==ceptj)[m]) {
              prevposition = which(ceptnumber[a]==ceptj)[m-1]
              return(dg1array[ t_class , which(dmap==(j - prevposition))])
            }
          }
        } # endif
      } # distortion()

      # compute probability
      prob_d = prod(sapply(X=1:le, FUN=distortion) )
      if (prob_d==0) return(0)
    ########################################

    # output
    return(as.numeric(prob_null*prob_fert*prob_t*prob_d))
  }
  hillclimb = function(a, e_sen, f_sen) {
    le = length(e_sen)
    lf = length(f_sen)
    aold = -1
    while (!all(a == aold)) {
      aold = a
      neighbours = neighbourhood(a,le,lf)
      oldprob = alignment_cache[[deparse(a)]]
      if (is.null(oldprob)) {
        oldprob = pr_IBM4(a, e_sen, f_sen)
        alignment_cache[[deparse(a)]] = oldprob
      }

      for (r in 1:nrow(neighbours)) {
        newprob = alignment_cache[[deparse(neighbours[r,])]]
        if (is.null(newprob)) {
          newprob = pr_IBM4(neighbours[r,], e_sen, f_sen)
          alignment_cache[[deparse(neighbours[r,])]] = newprob
        }
        if ( newprob > oldprob ) {
          a = neighbours[r,]
          oldprob=newprob
        }
      }# for
    }# while

    return(a)
  }
  sample_IBM4 = function(e_sen, f_sen, s) {
    if (samplemethod==1) {
      return(sample_IBM4_method1(e_sen, f_sen, s))
    } else {
      return(sample_IBM4_method2(e_sen, f_sen, s))
    }
  }
  prev_best_aligns = out_IBM3$prev_best_aligns
  sample_IBM4_method1 = function(e_sen, f_sen, s) {
    a_0 = prev_best_aligns[[s]]
    a_n = hillclimb(a_0,e_sen,f_sen)
    prev_best_aligns[[s]] <<- a_n
    return(neighbourhood(a_n,length(e_sen),length(f_sen)))
  }
  ##############################################################################



  ### INITIALIZE MATRICES ######################################################
  # split into words
  e_sentences = lapply(X=e,FUN=function(s) unlist(stringr::str_split(s, " ")))
  f_sentences = lapply(X=f,FUN=function(s) unlist(stringr::str_split(s, " ")))

  # list of all unique words
  e_allwords = unique(unlist(stringr::str_split(e, pattern=" ")))
  f_allwords = unique(unlist(stringr::str_split(f, pattern=" ")))
  n = length(e_sentences); n_eword = length(e_allwords); n_fword = length(f_allwords)

  # initialize t matrix and counts
  t_e_f = out_IBM3$tmatrix
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


  # init distortion arrays and word class lookup
  max_le = max(out_IBM3$corpus$target_lengths)
  nclass_t = max(target_wordclass)
  nclass_s = max(source_wordclass)
  d1array  = array(1/(2*max_le+1), dim=c(nclass_s+2, nclass_t, 2*(max_le)+1   )) # nclass_s+1 and nclass_s+2 deal with NULL token
  dg1array = array(1/(2*max_le+1), dim=c(nclass_t, 2*max_le+1   ))
  c_d1array  = array(0, dim=c(nclass_s+2, nclass_t, 2*max_le+1   ))
  c_dg1array = array(0, dim=c(nclass_t, 2*max_le+1   ))
  dmap = c((-max_le):-1 ,0,  1:max_le)
  lookup_tclass = names(target_wordclass)
  lookup_sclass = names(source_wordclass)


  # init NULL prob
  p_null = out_IBM3$p_null
  p1count = 0
  p0count = 0

  # init fertility matrix
  fertmatrix = out_IBM3$fmatrix
  fertcount  = new.env()
  for ( fword in f_allwords ) {
    fertcount[[fword]]  = rep(0, maxfert+1)
  }
  ##############################################################################


  time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
  if(verbose>=0.5) print(paste0("initial setup ;;; time elapsed: ",time_elapsed,"min"))

  # EM algorithm
  iter = 1
  perplex_vec = rep(0,n)
  prev_perplex = 0
  total_perplex = Inf
  while (iter<=maxiter & abs(total_perplex - prev_perplex)>eps) { # until convergence

    ### E STEP #################################################################
    if (verbose>=1) {
      pb = progress_bar$new(total=n,clear=TRUE,format=paste0("iter: ",iter," (E-step) [:bar] :current/:total (eta: :eta)")  )
      pb$tick(0)
    }
    for (s in 1:n) { # for all sentence pairs
      # extract kth sentence pair
      e_sen = e_sentences[[s]]; le = length(e_sen)
      f_sen = f_sentences[[s]]; lf = length(f_sen)

      # too many possible alignments -> get list of most likely ones from IBM2
      if (samplemethod==1) alignment_cache = new.env()
      if (le==1 | (le<=3 & lf<=3) ) { # always find all permutations for small sentences
        A = gtools::permutations(n=lf,r=le,v=1:lf,repeats.allowed=TRUE)
      } else if (heuristic)  {
        A = sample_IBM4(e_sen,f_sen,s)
      } else { # find all permutations even for large sentences if heuristic=FALSE (not recommended)
        A = gtools::permutations(n=lf,r=le,v=1:lf,repeats.allowed=TRUE)
      }

      # compute probability of each alignment
      ctotal = rep(NA, nrow(A))
      if (samplemethod==1) {
        for (r in 1:nrow(A)) {
          tmp = alignment_cache[[deparse(A[r,])]]
          if (is.null(tmp)) tmp = pr_IBM4(A[r,],e_sen,f_sen)
          ctotal[r] = tmp
        }
      } else {
        for (r in 1:nrow(A)) {
          ctotal[r] = pr_IBM4(A[r,],e_sen,f_sen)
        }
      }
      prev_best_aligns[[s]] = A[which.max(ctotal),]
      perplex_vec[s] = mean(ctotal)
      if (sum(ctotal)==0) warning(paste0("sentence ",s,": all sampled alignments had probability 0"),immediate. = TRUE)
      if (sum(ctotal)>0) ctotal = ctotal / sum(ctotal)

      # update expected counts given probabilities
      for (r in 1:nrow(A)) {
        if (ctotal[r]>0) { # only care about possible alignments
          a = A[r,]
          f_aj = f_sen[a] # f word to which each e word aligned
          phi = sapply(X=1:lf, FUN=function(i) sum(a==i))  # fertilities
          # determine cept numbering
          ceptnumber = rep(0,lf)
          ceptcount = 1
          for (i in 1:lf) {
            if(phi[i]>0) {
              ceptnumber[i] = ceptcount
              ceptcount = ceptcount+1
            }
          }

          for (j in 1:le) {
            c_e_f[[ e_sen[j] ]][[ f_aj[j] ]] = c_e_f[[ e_sen[j] ]][[ f_aj[j] ]] + ctotal[r]
            total_f[[f_aj[j]]] = total_f[[f_aj[j]]] + ctotal[r]

            ####################################################################
            # DISTORTION
              if (a[j]>1) { # ignore words aligned with NULL token

                t_class = target_wordclass[fastmatch::fmatch(e_sen[j],lookup_tclass)]
                ceptj = ceptnumber[a[j]]

                if (j==which(ceptnumber[a]==ceptj)[1]) { # first position in cept

                  prevcept = ceptj - 1
                  if(ceptj==1) { # if in the first cept, treat previous cept as <NULL> cept without distortion
                    c_d1array[ nclass_s+1   , t_class ,  which(dmap==(j - 0))  ] =
                      c_d1array[ nclass_s+1   , t_class ,  which(dmap==(j - 0))  ] + ctotal[r]
                  } else {
                    fwordaligned = f_sen[which(ceptnumber==prevcept)]
                    if (fwordaligned=="<NULL>") { # previous cept is <NULL>
                      s_class = nclass_s+2
                    } else { # previous cept is non-NULL word
                      s_class = source_wordclass[ fastmatch::fmatch(fwordaligned, lookup_sclass) ]
                    }
                    prevcept_positions = which(ceptnumber[a]==prevcept)
                    prevcenter = ceiling(mean(prevcept_positions))
                    c_d1array[ s_class   , t_class ,  which(dmap==(j - prevcenter))  ] =
                      c_d1array[ s_class   , t_class ,  which(dmap==(j - prevcenter))  ] + ctotal[r]
                  }


                } else { # position 2 or greater in cept

                  for (m in 2:phi[a[j]] ) {
                    if (j==which(ceptnumber[a]==ceptj)[m]) {
                      prevposition = which(ceptnumber[a]==ceptj)[m-1]
                      c_dg1array[ t_class , which(dmap==(j - prevposition))] =
                        c_dg1array[ t_class , which(dmap==(j - prevposition))] + ctotal[r]
                    }
                  }

                } # endif

              } # if (a[j]>1)
            ##################################################################
          }

          p1count = p1count + phi[1]*ctotal[r]
          p0count = p0count + (le - 2*phi[1])*ctotal[r]

          for (i in 2:lf) {
            fertcount[[ f_sen[i] ]][ min(phi[i],maxfert)+1 ] = fertcount[[ f_sen[i] ]][ min(phi[i],maxfert)+1 ] + ctotal[r]
          }
        } # if ctotal[r]>0
      } # for a in A

      if(verbose>=1 & s%%verbose==0) pb$tick(verbose)

    } # for k (all sentences)
    if(verbose>=1) pb$terminate()
    ############################################################################


    ### M STEP #################################################################
    if (verbose>=1) {
      pb = progress_bar$new(total=n_eword+(nclass_t*(nclass_s+2))+n_fword+1,clear=TRUE,format=paste0("iter: ",iter," (M-step) [:bar] :current/:total (eta: :eta)")  )
      pb$tick(0); k = 1
    }
    # translation probs
    for (eword in e_allwords) {
      for (fword in ls(t_e_f[[eword]])) {
        if (total_f[[fword]] == 0) {
          warning(paste0("total_f[[",fword,"]] = 0"))
        } else {
          t_e_f[[eword]][[fword]] = c_e_f[[eword]][[fword]] / total_f[[fword]]
        }
        c_e_f[[eword]][[fword]] = 0
      }

      if(verbose>=1){if(k%%verbose==0) pb$tick(verbose)}
      if(verbose>=1) k = k+1
    }

    # distortion probs
    # distortion probs
    for (h in 1:nclass_t) {
      if (sum(c_dg1array[h,])>0) dg1array[h,] = c_dg1array[h,] / sum(c_dg1array[h,])
      c_dg1array[h,] = 0
      for (m in 1:(nclass_s+2)) {
        if (sum(c_d1array[m,h,])>0) d1array[m,h,] = c_d1array[m,h,] / sum(c_d1array[m,h,])
        c_d1array[m,h,] = 0

        if(verbose>=1){if(k%%verbose==0) pb$tick(verbose)}
        if(verbose>=1) k = k+1
      }
    }

    # fertility probs
    for ( fword in f_allwords ) {
      if (sum(fertcount[[fword]]) == 0){
        if(!fword=="<NULL>") warning(paste0("fertcount[[",fword,"]] = 0"),immediate. = TRUE)
      } else {
        fertmatrix[[fword]] = fertcount[[fword]] / sum(fertcount[[fword]])
      }
      total_f[[fword]] = 0
      fertcount[[fword]][]  = 0

      if(verbose>=1){if(k%%verbose==0) pb$tick(verbose)}
      if(verbose>=1) k = k+1
    }

    # NULL token probs
    p_null = p1count / (p1count + p0count)
    p1count = 0
    p0count = 0
    if (verbose>=1) pb$terminate()
    ############################################################################


    ### HOUSEKEEPING ###########################################################
    # iterate
    prev_perplex = total_perplex
    total_perplex = -sum(log(perplex_vec[perplex_vec>0]))
    time_elapsed = round(difftime(Sys.time(),start_time,units='min'),3)
    print(paste0("iter: ",iter,"; perplexity value: ",total_perplex, "; time elapsed: ",time_elapsed,"min"))
    iter = iter + 1
    ############################################################################

  } # while not converged

  retobj = list(
    "tmatrix"=t_e_f,
    "d1array"=d1array,
    "dg1array"=dg1array,
    "fmatrix"=fertmatrix,
    "p_null"=as.numeric(p_null),
    "numiter"=iter-1,
    "maxiter"=maxiter,
    "eps"=eps,
    "converged"=abs(total_perplex - prev_perplex)<=eps,
    "perplexity"=total_perplex,
    "time_elapsed"=time_elapsed,
    "corpus"=out_IBM3$corpus,
    "prev_best_aligns"=prev_best_aligns
  )
  class(retobj) = "IBM4"
  return(retobj)

}# IBM4



