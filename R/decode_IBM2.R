
#' (IBM2) Compute translation probabilities given a foreign sentence.
#'
#' Takes a sentence in a foreign language f and produces the best translation
#' it can find based on IBM2 model. It starts with a word-by-word translation
#' into e language. Then it iterates between 4 heuristics: dropping a word,
#' reordering words, trying alternate translations of words, or adding
#' alternate translations of words. These heuristics are necessary to make
#' decoding computationally reasonable, but may result in a local optimum
#' rather than a global optimum.
#' @param object result from IBM2()
#' @param fsentence sentence in f language you'd like to translate to e language (vector or space-delimited string)
#' @param max_swap maximum distance away two words can be to consider swapping them during decoding (default 2).
#' @param numiter number of times heuristics should be tried (default=4). Usually only a couple iterations are needed to find local optimum.
#' @param verbose If TRUE (default), then print results of each heuristic step to the screen.
#' @return The best translation in e language found after numiter rounds of heuristics.
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
#' # estimate model
#' out = IBM2(e,f,maxiter=50,eps=0.01,init.IBM1=5);
#'
#' # possible english translations and their probabilities
#' besttranslation = decode(out, fsentence="une bière sil vous plaît")
#'
#' @import Matrix
#' @export
decode.IBM2 = function(object, fsentence, max_swap=2, numiter=4, verbose=TRUE) {

  f_eg = unlist(stringr::str_split(fsentence, pattern=" ")); lf=length(f_eg)

  # translation matrix relevant to given f sentence
  tmp = object$tmatrix[,f_eg]
  tmp = tmp[rowSums(tmp)>=0.1,,drop=FALSE]

  # find most likely translations for each f word
  besttrans = list()
  for (j in 1:lf) {
    allcand = tmp[,f_eg[j]]
    besttrans[[j]] = allcand[allcand>=0.1]
  }

  # starting point: best translation for each word
  currentbest = rep("",lf)
  maporig = 1:lf
  for (j in 1:lf) {
    if (length(besttrans[[j]]) > 0) {
      ewords = names(besttrans[[j]])
      currentbest[j] = ewords[which.max(besttrans[[j]])]
    }
  }
  maporig = maporig[currentbest!=""]
  currentbest = currentbest[currentbest!=""]
  le = length(currentbest)
  maporig = 1:le
  if (length(object$amatrix[[le]]) >= lf) {
    if (!is.null(object$amatrix[[le]][[lf]])) {
      currentprob = prod(rowSums(tmp[currentbest,f_eg,drop=FALSE]*object$amatrix[[le]][[lf]]))
    } else {
      currentprob = 0
    }
  } else {
    currentprob = 0
  }
  if(verbose) print(paste0(c("START:", currentbest,"-------- Probability:",currentprob), collapse=" "))


  iter = 0
  lastprob = -1
  while((iter < numiter) & (currentprob!=lastprob)) {

    lastprob=currentprob

    # heuristic 1: delete a word
    if (length(object$amatrix[[le-1]]) >= lf ) {
      if ( !is.null( object$amatrix[[le-1]][[lf]] ) ) {
        if(verbose) print("----Heuristic 1 (Delete)----")
        new_i = NULL
        for (i in 1:le) {
          new_esen = currentbest[-i]
          new_prob = prod(rowSums(tmp[new_esen,f_eg,drop=FALSE]*object$amatrix[[le-1]][[lf]]))
          if(verbose) print(paste0(c(new_esen,"-------- Probability:",new_prob),collapse=" "))
          if (new_prob > currentprob) {
            new_i = i
            currentprob = new_prob
          }
        }
        if (!is.null(new_i)) {
          maporig = maporig[-new_i]
          currentbest = currentbest[-new_i]; le = length(currentbest)
          if(verbose) print(paste0(c("UPDATE:", currentbest,"-------- Probability:",currentprob), collapse=" "))
        }
      }
    }

    # if e sentence still too long, drop "worst" e words
    while (length(object$amatrix[[le]]) < lf ) {
      dropme = which.min(sapply(X=besttrans,FUN=max)[maporig])
      currentbest = currentbest[-dropme]
      maporig = maporig[-dropme]
      le = length(currentbest)
    }

    # heuristic 2: local reordering
    if(verbose) print("----Heuristic 2 (Local reordering)----")
    new_i = NULL
    new_i_swap = NULL
    for (i in 1:le) {
      for (dist in 1:max_swap) {
        if (i + dist <= le) {
          new_esen = currentbest
          new_esen[i] = currentbest[i + dist]
          new_esen[i + dist] = currentbest[i]
          new_prob = prod(rowSums(tmp[new_esen,f_eg,drop=FALSE]*object$amatrix[[le]][[lf]]))
          if(verbose) print(paste0(c(new_esen,"-------- Probability:",new_prob),collapse=" "))
          if (new_prob > currentprob) {
            new_i = i
            new_i_swap = i+dist
            currentprob = new_prob
          }
        }
      }
    }
    if (!is.null(new_i)) {
      maporig_tmp = maporig
      currentbest_temp = currentbest

      currentbest[new_i] = currentbest_temp[new_i_swap]
      currentbest[new_i_swap] = currentbest_temp[new_i]
      maporig[new_i] = maporig_tmp[new_i_swap]
      maporig[new_i_swap] = maporig_tmp[new_i]

      if(verbose) print(paste0(c("UPDATE:", currentbest,"-------- Probability:",currentprob), collapse=" "))
    }

    # heuristic 3: swapping a word with alternative translation
    if(verbose) print("----Heuristic 3 (Swap translations)----")
    new_i = NULL
    new_word = ""
    for (i in 1:le) {
      for (etry in names(besttrans[[maporig[i]]]) ) {
        new_esen = currentbest
        new_esen[i] = etry
        new_prob = prod(rowSums(tmp[new_esen,f_eg,drop=FALSE]*object$amatrix[[le]][[lf]]))
        if(verbose) print(paste0(c(new_esen,"-------- Probability:",new_prob),collapse=" "))
        if (new_prob > currentprob) {
          new_i = i
          new_word = etry
          currentprob = new_prob
        }
      }
    }
    if (!is.null(new_i)) {
      currentbest[new_i] = new_word
      if(verbose) print(paste0(c("UPDATE:", currentbest,"-------- Probability:",currentprob), collapse=" "))
    }

    # heuristic 4: add a translation word
    if ( length(object$amatrix[[le+1]]) >= lf )
      if ( !is.null( object$amatrix[[le+1]][[lf]] ) ) {
        if(verbose) print("----Heuristic 4 (Add translation)----")
        new_i = NULL
        new_word = ""
        after=0
        for (i in 1:le) {
          for (etry in names(besttrans[[maporig[i]]]) ) {
            new_esen = currentbest
            new_esen = append(new_esen,etry,after=i)
            new_prob = prod(rowSums(tmp[new_esen,f_eg,drop=FALSE]*object$amatrix[[le+1]][[lf]]))
            if(verbose) print(paste0(c(new_esen,"-------- Probability:",new_prob),collapse=" "))
            if (new_prob > currentprob) {
              new_i = i
              new_word = etry
              after=0
              currentprob = new_prob
            }

            if (length(names(besttrans[[maporig[i]]])) > 1) {
              new_esen = currentbest
              new_esen = append(new_esen,etry,after=i-1)
              new_prob = prod(rowSums(tmp[new_esen,f_eg,drop=FALSE]*object$amatrix[[le+1]][[lf]]))
              if(verbose) print(paste0(c(new_esen,"-------- Probability:",new_prob),collapse=" "))
              if (new_prob > currentprob) {
                new_i = i
                new_word = etry
                after = -1
                currentprob = new_prob
              }
            }
          }
        }
        if (!is.null(new_i)) {
          maporig     = append(maporig, maporig[new_i], after=new_i+after)
          currentbest = append(currentbest, new_word, after=new_i+after)
          le = length(currentbest)
          if(verbose) print(paste0(c("UPDATE:", currentbest,"-------- Probability:",currentprob), collapse=" "))
        }
      }

    iter = iter + 1

  } # while

  return( paste0(currentbest, collapse = " ") )

} # decode.IBM2
