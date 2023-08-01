




#' Part-of-Speech (POS) Tagger
#'
#' This a simple wrapper for udpipe POS tagging. It takes a single sentence (can be longer text as well) and returns the same sentence segmented and with POS tags in square brackets attached to each word.
#' @param sen A sentence (character string)
#' @param udpipe.model a pre-trained tagger model returned from a udpipe::udpipe_load_model call
#' @returns Return a character string representing a segmented and tagged version of sen. See example below.
#' @examples
#' m_eng = udpipe::udpipe_download_model(language = "english-ewt")
#' m_eng = udpipe::udpipe_load_model(m_eng)
#' pos.tagger("This is an example",m_eng) # produces "this[PRON] is[AUX] an[DET] example[NOUN]"
#' @export
pos.tagger = function(sen, udpipe.model) {
  udpipeoutput = udpipe::udpipe_annotate(udpipe.model, x = sen) %>% as.data.frame()
  words = tolower(udpipeoutput$token)
  tags = udpipeoutput$upos
  sen1_pos = paste0(words, "[", udpipeoutput$upos ,"]")
  return(paste0(sen1_pos,collapse=" "))
}
