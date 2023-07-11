


IBM1 = function(e,f,maxiter=30,eps=0.01,sparse=FALSE) {

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



  # EM algorithm
  iter = 1
  prev_perplex = 0
  total_perplex = Inf
  start_time = Sys.time()
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

  return(t_e_f)

} # end function IBM1

