# SMT
Statistical Machine Translation (SMT) in R

## Introduction

This R package implements the translation models IBM1, IBM2, and IBM3 from `Brown et al (1993)`. While SMT methods have largely been overtaken by neural network-based methods, the IBM models remain historically important and are still used in conjunction with neural methods thanks to their word-alignment capabilities.

## Data

These methods rely on being able to take a word or words and look up the corresponding probability in a table.
To do this fast I use R environment objects, which are the closest equivalent to python
dictionaries available in base R. These use hashing, making the table lookups
(on average) O(1) time. For example, each model has `tmatrix`, which gives the translation
probability of an output word given an input word. E.g., `tmatrix[["fish"]][["poisson"]]`
gives the probability that "fish" in English is the proper translation of "poisson" in French.
The downside is that you cannot extract "columns" in the same way you would a matrix,
and the size of the object is not immediately obvious. Because of this
I provide functions `object.size.IBM1`, `object.size.IBM2`, `object.size.IBM3`.

## Estimation

Algorithms largely follow `Koehn (2009)`. For IBM3, I use the simpler heuristic
without pegging as described in the appendix of `Och & Ney (2003)`. Given vector of
output sentences `e` and input sentences `f`, each function
follows the basic syntax of `IBMX(e,f,maxiter,eps)` where `X` is 1, 2, or 3, `maxiter`
is the maximum number of iterations allowed, and `eps` is the stopping criteria
for convergence. The sentences in `e` and `f` are assumed to contain space-delimited
words. Notice then that "Wow!" would be considered one word while "Wow !" would
be considered two words. You will likely want to do some pre-processing on the sentences
before throwing them in the models.

As explained in the cited papers, lower-numbered models are meant to be fed as input
into higher-numbered models. For this reason, each function calls the previous one
by default in order to initialize the algorithm. 
You can specify the number of iterations for the initialization, e.g.,
`IBM3(e,f,maxiter,eps,init.IBM1=30,init.IBM2=15)`. You can also set the `init`
arguments to 0 and use results from previous model estimations instead. 
See help pages for more details.

## Decoding

To generate most likely translations for a given sentence in language f to language e, 
use the generic method `decode`. It is computationally infeasible to consider all
possible translations, so in practice we rely on heuristic algorithms.
Here I implement Algorithm 1 of `Wang & Waibel (1998)`, which is a stack decoding
approach for IBM1 and IBM2. See `?decode.IBM1`, `?decode.IBM2` for more details.

## Evaluate

Given a list of decoded translations, we then need to evaluate how good they are.
I provide the `evaluate` function for this purpose, which provides several measures
decribed in `Koehn (2009)`. See `?evaluate`.

## Help

See `?IBM1`, `?IBM2`, `?IBM3`, `?decode.IBM1`, `?decode.IBM2`, `evaluate` for more details.

## Sources

`Brown, P. F., Della Pietra, S. A., Della Pietra, V. J., & Mercer, R. L. (1993). The mathematics of statistical machine translation: Parameter estimation.`

`Koehn, P. (2009). Statistical machine translation. Cambridge University Press.`

`Och, F. J., & Ney, H. (2003). A systematic comparison of various statistical alignment models. Computational linguistics, 29(1), 19-51.`

`Wang, Y. Y., & Waibel, A. (1998). Fast decoding for statistical machine translation. In Fifth International Conference on Spoken Language Processing.`
