
#' Evaluate Translation
#'
#' Calculates a series of evaluation metrics for a hypothesis translation based
#' on a given reference translation. Based on Koehn's "Statistical Machine Translation"
#' Chapter 8. Word Error Rate (WER) adapted from https://www.thepythoncode.com/article/calculate-word-error-rate-in-python
#' @param hypothesis the translation we wish to evaluate
#' @param reference the "ground truth" translation against which we compare our hypothesis
#' @param type character vector containing any of "precision","recall","f-measure","PER","WER","BLEU" corresponding to each type of evaluation metric. Default is all of them.
#' @param BLEU.maxorder If "BLEU" specified, this is the maximum ngram precision to consider in the calculation. (integer >= 1)
#' @param BLEU.weights If "BLEU" specified, this is a numeric vector of weights applied to the ngram precisions. Thus, its length should be equal to `BLEU.maxorder`. If sum of weights don't equal 1, they are rescaled.
#' @return A list containing the values of each evaluation metric. "BLEU" contains a sublist with information about the assumptions for the BLEU calculation.
#' @examples
#' systemA = "Israeli officials responsibility of airport safety"
#' systemB = "airport security Israeli officials are responsible"
#' reference   = "Israeli officials are responsible for airport security"
#' evaluate(systemA,reference)
#' evaluate(systemB,reference)
#' @import ngram
#' @export
evaluate = function(
    hypothesis, reference,
    type=c("precision","recall","f-measure","PER","WER","BLEU"),
    BLEU.maxorder=4,
    BLEU.weights=rep(1/BLEU.maxorder,BLEU.maxorder)
    )
{

  # convert to word vectors
  hypothesis_split = unlist(stringr::str_split(hypothesis, pattern=" "))
  reference_split   = unlist(stringr::str_split(reference, pattern=" "))

  # empty list to be filled
  retlist = list()

  # compute number of correct
  t_counts = table(hypothesis_split)
  r_counts = table(reference_split)
  n_correct = 0
  for (s in names(t_counts)) {
    if (s %in% reference_split) n_correct = n_correct + min(t_counts[s], r_counts[s])
  }

  # evaluation metrics
  if ("precision" %in%  type) retlist[["precision"]] = n_correct / length(hypothesis_split)
  if ("recall" %in%  type)    retlist[["recall"]]    = n_correct / length(reference_split)
  if ("f-measure" %in% type)  retlist[["fmeasure"]]  = 2*n_correct / (length(hypothesis_split) + length(reference_split))
  if ("PER" %in% type)        retlist[["PER"]]       = 1 - (n_correct - max(0,length(hypothesis_split) - length(reference_split)))/length(reference_split)

  if ("WER" %in% type) {
    d = matrix(0, nrow=length(reference_split)+1, ncol=length(hypothesis_split)+1)
    d[,1] = 1:nrow(d) - 1
    d[1,] = 1:ncol(d) - 1
    for (i in 1:length(reference_split)) {
      for(j in 1:length(hypothesis_split)) {
        if (reference_split[i] == hypothesis_split[j]) {
          d[i+1,j+1] = d[i,j]
        } else {
          subst  = d[i  ,j  ] + 1
          insert = d[i+1,j  ] + 1
          delete = d[i  ,j+1] + 1
          d[i+1,j+1] = min(subst,insert,delete)
        }
      }#for j
    }#for i

    retlist[["WER"]] = d[nrow(d),ncol(d)] / length(reference_split)
  }#WER

  if ("BLEU" %in% type) {
    BLEU.weights = BLEU.weights / sum(BLEU.weights)
    ng_precisions = rep(NA,BLEU.maxorder)
    for (n in 1:BLEU.maxorder) {
      ng_hyp = get.ngrams(ngram(hypothesis, n=n))
      ng_ref = get.ngrams(ngram(reference,n=n))
      t_counts = table(ng_hyp)
      r_counts = table(ng_ref)
      n_correct = 0
      for (s in names(t_counts)) {
        if (s %in% ng_ref) n_correct = n_correct + min(t_counts[s], r_counts[s])
      }
      ng_precisions[n] = n_correct / length(ng_hyp)
    }

    retlist[["BLEU"]] = list(
      "value"=min(1,length(hypothesis)/length(reference)) * prod(ng_precisions^BLEU.weights),
      "ng_precisions"=ng_precisions,
      "max.order"=BLEU.maxorder,
      "weights"=BLEU.weights)
  }

  return(retlist)
}





