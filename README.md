# SMT
Statistical Machine Translation (SMT) in R

## Estimation

This R package implements models IBM1 and IBM2 from `Brown, P. F., Della Pietra, S. A., Della Pietra, V. J., & Mercer, R. L. (1993). The mathematics of statistical machine translation: Parameter estimation.` While SMT methods have largely been overtaken by neural network-based methods, the IBM models remain historically important and are still used in conjunction with neural methods thanks to their word-alignment capabilities. This package opens up the IBM models to a broader audience (the methods seem only to be implemented in the C++ program `GIZA++` with not-great documentation).

By default, the `IBM1()` and `IBM2()` functions use base R matrix methods. As the corpus size grows larger, these matrices grow vast in size. I thus experiment with using sparse matrices from the `Matrix` package for more efficient memory usage. This can be invoked with `sparse=TRUE`.

I also experiment with two other methods of speedup: using `fmatch` from the `fastmatch` package for faster string lookup, and parallelization for parts of the code using the `parallel` package.

## Decoding

To generate most likely translations for a given sentence in language f to language e, 
use the generic method `decode`. Currently this method considers all combinations (for IBM1) or permutations (for IBM2) of e words. To reduce search time, it considers only e words above a threshold of probability. However,
this can still take a lot of time and memory for even moderately long f sentences,
especially for IBM2
(for example, there are over 800,000 possible 5-word sentences based on permutations,
less than 2,000 5-word sentences based on combinations). You can improve search time and memory usage
for IBM2 further with the option `useIBM1=TRUE`, which first finds most likely combinations
using the IBM1 method, and then finds the most likely permutations from this shorter subset.

## Help

See `?IBM1`, `?IBM2`, `?decode.IBM1`, `?decode.IBM2` for more details.
