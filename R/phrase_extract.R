



#' Symmetrization of Word Alignments
#'
#' Takes two vectors representing alignments of two sentences going in opposite directions.
#' Uses the heuristic algorithm described in Koehn (2009) to combine the alignments
#' into a single set of "symmetrized" alignments.
#' The result is a list of length-2 vectors representing the final set of alignment points.
#' @export
symalign = function(a1,a2) {

  # setup
  neighbours = list(c(-1,0),c(0,-1),c(1,0),c(0,1),
                    c(-1,-1),c(-1,1),c(1,-1),c(1,1))
  unaligned1 = 1:length(a1)
  unaligned2 = 1:length(a2)

  # intersection
    a_intersect = list()
    for (i in 1:length(a1)) {
      if (a1[i]!=0) {
        if (i == a2[a1[i]]) {
          a_intersect[[length(a_intersect)+1]] = c(i,a1[i])
          unaligned1 = unaligned1[-which(unaligned1==i)]
          unaligned2 = unaligned2[-which(unaligned2==a1[i])]
        }
      }
    }

  # union
    a_union = NULL
    for (i in 1:length(a1)) {
      if (a1[i]>0) a_union[[length(a_union)+1]] = c(i,a1[i])
    }
    for (j in 1:length(a2)) {
      if (a2[j]>0) a_union[[length(a_union)+1]] = c(a2[j],j)
    }
    a_union = unique(a_union)

  # heuristic
    # start with intersection
    a_final = a_intersect
    # add neighbours if unaligned
    for (i in 1:length(a1)) {
      for (j in 1:length(a2)) {
        if (list(c(i,j)) %in% a_final) {
          for (n in neighbours) {
            newpoint = c(i,j) + n
            if ((min(newpoint)>0) & (newpoint[1]<=length(a1)) & (newpoint[2]<=length(a2))) {
              if ((newpoint[1] %in% unaligned1) | (newpoint[2] %in% unaligned2)) {
                if (list(newpoint) %in% a_union) {
                  a_final[[length(a_final)+1]] = newpoint
                  unaligned1 = unaligned1[-which(unaligned1==newpoint[1])]
                  unaligned2 = unaligned2[-which(unaligned2==newpoint[2])]
                }
              }
            }
          }
        }
      }
    }
    # add more from union if still unaligned words
    for (i in 1:length(a1)) {
      for (j in 1:length(a2)) {
        if ((i %in% unaligned1) | (j %in% unaligned2)) {
          if (list(c(i,j)) %in% a_union) {
            a_final[[length(a_final)+1]] = c(i,j)
            unaligned1 = unaligned1[-which(unaligned1==i)]
            unaligned2 = unaligned2[-which(unaligned2==j)]
          }
        }
      }
    }

    return(a_final)
}


#' Extract phrases from alignment points
#'
#' Takes a list of alignment points (e.g. as produced by symalign function)
#' and returns a list of phrase pairs extracted from the alignments.
#' Based on algorithm provided in Koehn (2009).
#' @export
phrase_extract = function(align_points,l1,l2) {

  # globals
  phrase_pairs = list()
  apmatrix = matrix(unlist(align_points),ncol=2,byrow=TRUE)

  # helper function
  help_extract = function(e_start,e_end,f_start,f_end) {
    if (f_end==0) return(NULL)
    for (a in align_points) {
      if ( ((a[1]<e_start) | (a[1]>e_end)) & ((f_start <= a[2]) & (a[2] <= f_end)) ) return(NULL)
    }

    fs = f_start
    repeat {
      fe = f_end
      repeat {
        phrase_pairs[[length(phrase_pairs)+1]] <<- list(target=c(e_start,e_end),source=c(fs,fe))
        fe = fe+1
        if (fe > l2) break
        if (fe %in% apmatrix[,2]) break
      }
      fs = fs-1
      if (fs < 1) break
      if (fs %in% apmatrix[,2]) break
    }
  } # help_extract

  # extract phrase pairs
  for (e_start in 1:l1) {
    for (e_end in e_start:l1) {
      f_start = l2; f_end = 0
      for (a in align_points) {
        if ((e_start <= a[1]) & (a[1] <= e_end)) {
          f_start = min(f_start,a[2])
          f_end = max(f_end,a[2])
        }
      }
      help_extract(e_start,e_end,f_start,f_end)
    }
  }

  return(phrase_pairs)

} # phrase_extract




#' Build phrase table of probabilities
#'
#' Takes two models (any of IBM1 through IBM4) translating in opposite directions.
#' For each sentence, symmetrizes the best (viterbi) alignments
#' and extracts phrases based on the alignments. Produces the phrase translation
#' probability table as an R environment object.
#' @return
#'    \item{phrase_table}{Environment object containing phrase translation probabilities.}
#'    \item{corpus}{Data frame containing the corpus used to train the models.}
#' @examples
#'  # data
#'  temp = tempfile()
#'  download.file("http://www.manythings.org/anki/fra-eng.zip",temp);
#'  ENFR = readr::read_tsv(file=unz(temp,"fra.txt"),col_names=c("en","fr","details"));
#'  unlink(temp);
#'
#'  # a bit of pre-processing
#'  e = removePunctuation(ENFR$en[1:10000])
#'  e = str_squish(e)
#'  e = tolower(e)
#'  f = removePunctuation(ENFR$fr[1:10000])
#'  f = str_squish(f)
#'  f = tolower(f)
#'
#'  # estimate models in both directions
#'  model1 = IBM1(target=e,source=f, maxiter=30)
#'  model2 = IBM1(target=f,source=e, maxiter=30)
#'
#'  # phrase table (prob. of f phrase given e phrase)
#'  phtable = build_phrase_table(model1,model2)
#'
#'  # e.g. phrase probabilities for e phrase "calm down"
#'  phtable$phrase_table[["calm down"]][["calmetoi"]]
#'  phtable$phrase_table[["calm down"]][["calmezvous"]]
#'  phtable$phrase_table[["calm down"]][["du calme"]]
#'  phtable$phrase_table[["calm down"]][["tranquille"]]
#'
#'  # we can get phrase table in opposite direction by reversing model order:
#'  phtable2 = build_phrase_table(model2,model1)
#'  phtable$phrase_table[["calmetoi"]][["calm down"]]
#'  phtable$phrase_table[["calmetoi"]][["chill out"]]
#'  phtable$phrase_table[["calmetoi"]][["relax"]]
#' @export
build_phrase_table = function(model1, model2, ondisk=FALSE) {

  if (!ondisk) {
    trg_sentences = lapply(X=model1$corpus$target,FUN=function(s) unlist(stringr::str_split(s, " ")))
    src_sentences = lapply(X=model2$corpus$target,FUN=function(s) unlist(stringr::str_split(s, " ")))
    n = length(trg_sentences)

    phrase_table = new.env()
    for (i in 1:n) {
      b1=model1$best_alignments[[i]]-1
      b2=model2$best_alignments[[i]]-1
      myphrases = phrase_extract(symalign(b1,b2), length(b1), length(b2))
      for (ph in myphrases) {
        trgphrase = paste0(trg_sentences[[i]][ph$target[1]:ph$target[2]],collapse=" ")
        srcphrase = paste0(src_sentences[[i]][ph$source[1]:ph$source[2]],collapse=" ")

        if (is.null(phrase_table[[trgphrase]]))              phrase_table[[trgphrase]] = new.env()
        if (is.null(phrase_table[[trgphrase]][[srcphrase]])) phrase_table[[trgphrase]][[srcphrase]] = 0
        if (is.null(phrase_table[[trgphrase]][["<TOTAL>"]])) phrase_table[[trgphrase]][["<TOTAL>"]] = 0

        phrase_table[[trgphrase]][[srcphrase]] = phrase_table[[trgphrase]][[srcphrase]] + 1
        phrase_table[[trgphrase]][["<TOTAL>"]] = phrase_table[[trgphrase]][["<TOTAL>"]] + 1
      }
    }

    # compute probabilities
    for (e in ls(phrase_table)) {
      for ( f in ls(phrase_table[[e]]) ) {
        if (f!="<TOTAL>") phrase_table[[e]][[f]] = phrase_table[[e]][[f]] / phrase_table[[e]][["<TOTAL>"]]
      }
      rm("<TOTAL>",envir=phrase_table[[e]])
    }


    retobj = list()
    retobj$phrase_table = phrase_table
    retobj$corpus = data.frame(target=model2$corpus$target,
                               source=model1$corpus$target,
                               target_lengths=model2$corpus$target_lengths,
                               source_lengths=model1$corpus$target_lengths)
    class(retobj) = "phrase_table"
    return(retobj)
  }
}
























