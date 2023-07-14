# SMT
Statistical Machine Translation (SMT) in R

This R package implements models IBM1 and IBM2 from `Brown, P. F., Della Pietra, S. A., Della Pietra, V. J., & Mercer, R. L. (1993). The mathematics of statistical machine translation: Parameter estimation.` While SMT methods have largely been overtaken by neural network-based methods, the IBM models remain historically important and are still used in conjunction with neural methods thanks to their word-alignment capabilities. This package opens up the IBM models to a broader audience (the methods seem only to be implemented in the C++ program `GIZA++` with not-great documentation).

By default, the `IBM1()` and `IBM2()` functions use base R matrix methods. As the corpus size grows larger, these matrices grow vast in size. I thus experiment with using sparse matrices from the `Matrix` package for more efficient memory usage. This can be invoked with `sparse=TRUE`.

I also experiment with two other methods of speedup: using `fmatch` from the `fastmatch` package for faster string lookup, and parallelization for parts of the code using the `parallel` package.

See `?IBM1`, `?IBM2`, `?predict.IBM1` for more details.
