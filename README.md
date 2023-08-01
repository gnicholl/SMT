# SMT
Statistical Machine Translation (SMT) in R

## Estimation

This R package implements models IBM1 and IBM2 from `Brown, P. F., Della Pietra, S. A., Della Pietra, V. J., & Mercer, R. L. (1993). The mathematics of statistical machine translation: Parameter estimation.` While SMT methods have largely been overtaken by neural network-based methods, the IBM models remain historically important and are still used in conjunction with neural methods thanks to their word-alignment capabilities.

By default, the `IBM1()`, `IBM2()` functions use base R matrix methods. As the corpus size grows larger, these matrices grow vast in size. I thus experiment with using sparse matrices from the `Matrix` package for more efficient memory usage. This can be invoked with `sparse=TRUE`.

I also experiment with two methods of speedup: using `fmatch` from the `fastmatch` package for faster string lookup, and parallelization for parts of the code using the `parallel` package.

## Decoding

To generate most likely translations for a given sentence in language f to language e, 
use the generic method `decode`. For IBM1 models, it considers all combinations of
possible output words, as IBM1 is unable to distinguish different alignments of words.

For IBM2, in theory one would have to consider all permutations of words.
Instead I use a simple heuristic involving adding, deleting, and swapping words.

## Vignette

See [here for examples.](https://htmlpreview.github.io/?https://github.com/gnicholl/SMT/blob/main/vignettes/translateAndDecode.html)

## Help

See `?IBM1`, `?IBM2`, `?decode.IBM1`, `?decode.IBM2` for more details.

## Sources

`Brown, P. F., Della Pietra, S. A., Della Pietra, V. J., & Mercer, R. L. (1993). The mathematics of statistical machine translation: Parameter estimation.`

`Koehn, P. (2009). Statistical machine translation. Cambridge University Press.`
