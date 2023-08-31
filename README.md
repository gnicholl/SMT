# SMT: Statistical Machine Translation in R

This R package implements the translation models IBM1, IBM2, IBM3, and IBM4 from `Brown et al (1993)` as well as phrase-based translation. While SMT methods have largely been overtaken by neural network-based methods, the IBM models remain historically important and are still used in conjunction with neural methods thanks to their word-alignment capabilities.

## Noisy-Channel Model

Our objective is to find the best translation $\mathbf{e} = (e_1,\dots,e_n)$ for the sentence $\mathbf{f} = (f_1,\dots,f_m)$. Rather than model $\Pr(\mathbf{e} | \mathbf{f})$ directly, we model the reverse direction $\Pr(\mathbf{f} | \mathbf{e})$ and use the noisy-channel approach:

$$
\mathbf{e} = \text{argmax}_{\mathbf{e}} \Pr(\mathbf{e} | \mathbf{f}) = \text{argmax} _{\mathbf{e}} \Pr(\mathbf{f} | \mathbf{e}) \Pr(\mathbf{e})
$$

where $\Pr(\mathbf{f} | \mathbf{e})$ is the translation model (which will be estimated using the IBM models or phrase tables) and $\Pr(\mathbf{e})$ is the language model (here estimated using the `kgrams` R package).

Actually, there is one piece missing. The translation and language models we will consider implicitly treat $n$, the length of sentence $\mathbf{e}$, as given. But, we don't know in advance what length the sentence should be when we are translating--i.e. $n$ is a random variable! We thus add the sentence length model $\Pr(n | m)$, assuming that the length of $\mathbf{e}$ depends only on the length of $\mathbf{f}$. Then:

$$
\mathbf{e} = \text{argmax}_{\mathbf{e}} \Pr(\mathbf{e} | \mathbf{f}, m) = \text{argmax} _{\mathbf{e}} \Pr(\mathbf{e} | \mathbf{f}, n, m) \Pr(n | m) = \text{argmax} _{\mathbf{e}} \Pr(\mathbf{f} | \mathbf{e}, n, m) \Pr(\mathbf{e} | n) \Pr(n | m)
$$

In the decoding algorithm of `Wang & Waibel (1998)`, they state that they model $\Pr(n | m)$ using a poisson regression, which is what we do here by default.

To summarise:

- $\Pr(\mathbf{f} | \mathbf{e}, n, m)$: translation probability estimated using IBM models or phrase tables; implemented in this package
- $\Pr(\mathbf{e} | n)$: language model estimated using ngrams; implemented in the `kgrams` R package
- $\Pr(n | m)$: sentence length model; by default estimated using the built-in `glm` function to run a simple poisson regression

## Estimating the IBM Models

Algorithms for IBM1 and IBM2 largely follow those provided in `Koehn (2009)`. 
For IBM3 and IBM4, I use the heuristic algorithm _without_ pegging as described in the appendix of `Och & Ney (2003)`. 
Given a vector of output sentences `target` and input sentences `source`, each function
follows the basic syntax of `IBMX(target,source,maxiter,eps)` where `X` is 1, 2, or 3, `maxiter`
is the maximum number of iterations allowed, and `eps` is the stopping criteria
for convergence. 

The sentences in `target` and `source` are assumed to be strings containing space-delimited
words. Notice then that "Wow!" would be considered one word while "Wow !" would
be considered two words. You will likely want to do some pre-processing on the sentences
before throwing them into the models. Some examples of handy pre-processing functions are
`tm::removePunctuation`, `stringr::str_squish`, and `base::tolower`.
The only pre-processing that is done by the `IBMX` functions is to add "\<NULL>" to the beginning of all sentences in `source`. 
E.g. "do you have time ?" becomes "\<NULL> do you have time ?" This is done to allow words in the `target` sentences to be aligned
with "nothing".

As explained in the cited papers, lower-numbered models are meant to be fed as input
into higher-numbered models. For this reason, each function calls the previous one
by default in order to initialize the algorithm. 
You can specify the number of iterations for the initialization, e.g.,
`IBM3(e,f,maxiter,eps,init.IBM1=30,init.IBM2=15)`. You can also set the `init`
arguments to 0 and use results from previous model estimations instead. 
See help pages for more details.

## Estimating Phrase Tables

Phrase-based translation traditionally takes the following steps for estimating phrase probability tables (see `Koehn (2009)`):

1. Estimate IBM models in both directions (e.g. English-to-French _and_ French-to-English)
2. Combine alignments from the two models to obtain a "symmetrized" set of sentence alignments
3. Extract phrases using the symmetrized alignments
4. Estimate phrase probabilities based on the frequency with which each phrase occurs

Steps 2 to 4 are handled by the `build_phrase_table` function. For example:

```
# estimate models in both directions
model3_ftoe = IBM3(target=e,source=f,maxiter=5, init.IBM1=15, init.IBM2=15,  verbose=100)
model3_etof = IBM3(target=f,source=e,maxiter=5, init.IBM1=15, init.IBM2=15,  verbose=100)

# extract phrases from the IBM models and compute phrase probabilities
phtable = build_phrase_table(model3_ftoe, model3_etof)
```


## Data Types

These methods rely on being able to take a word or words and look up the corresponding probability in a table.
To do this fast I use R environment objects, which are the closest equivalent to python
dictionaries available in base R. These use hashing, making the table lookups
O(1) time. For example, each IBM model produces a `tmatrix`, which gives the translation
probability of an output word given an input word. E.g., `tmatrix[["fish"]][["poisson"]]`
(or equivalently `tmatrix$fish$poisson`) gives the probability that "fish" in English is the proper translation of "poisson" in French.

Environments are a bit less intuitive than matrices. Firstly, the built-in R function `object.size`
does not actually compute the total size of all objects in the environment.
Because of this I provide functions `object.size.IBM1`, `object.size.IBM2`, `object.size.IBM3`
which iterate through each level of the environments and adds up the sizes.
Secondly, unlike matrices, `tmatrix2 = tmatrix` does not make a copy of the enviornment `tmatrix`. Instead,
it copies only the _pointer_ to `tmatrix`. This means that making changes to `tmatrix2`
will also affect `tmatrix1`. 

## Decoding

To generate most likely translations for a given sentence in language f to language e, 
use the generic method `decode(object,target.sentence,senlength.model,language.model)`
where `object` is an object returned from one `IBM1()`, `IBM2()`, or `build_phrase_table()`;
`target.sentence` is the sentence in the language we want to translate _from_,
and `senlength.model`, `language.model` are the sentence-length and language models described
previously in the noisy channel section. See `?decode.IBM2` for how they should be specified.
If `senlength.model` or `language.model` aren't specified, they have defaults:

- `language.model` is a 3-gram model based on maximum likelihood without any smoothing
- `senlength.model` is a matrix of predicted probabilities based on the regression `glm(object$corpus$target_lengths ~ object$corpus$source_lengths, family="poisson")`

It is computationally infeasible to consider all
possible translations while decoding, so in practice we rely on heuristic algorithms.
For IBM1 and IBM2 I implement Algorithm 1 of `Wang & Waibel (1998)`, which is a stack decoding approach. 
For phrase translation I implement the methods described in `Koehn (2009)` chapter 6.

## Evaluate

Given a decoded translation, we then need to evaluate how good the translation is!
Ideally one would use human evaluators, but our next-best option is 
some automated evaluation measures. I provide the `evaluate` function for this purpose, 
which implements precision, recall, f-measure, PER, WER, and BLEU as decribed in `Koehn (2009)`. 
See `?evaluate` for more details.

## Help

See `?IBM1`, `?IBM2`, `?IBM3`, `?IBM4`, `?build_phrase_table`, `?decode.IBM1`, `?decode.IBM2`, `?decode.phrase_table`, `?evaluate` for more details.

## Sources

`Brown, P. F., Della Pietra, S. A., Della Pietra, V. J., & Mercer, R. L. (1993). The mathematics of statistical machine translation: Parameter estimation.`

`Koehn, P. (2009). Statistical machine translation. Cambridge University Press.`

`Och, F. J., & Ney, H. (2003). A systematic comparison of various statistical alignment models. Computational linguistics, 29(1), 19-51.`

`Wang, Y. Y., & Waibel, A. (1998). Fast decoding for statistical machine translation. In Fifth International Conference on Spoken Language Processing.`
