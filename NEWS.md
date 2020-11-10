# monoClust 1.0.1

* Fix some typos and clarify some documentation
* Min version of dependency `tibble()` is 3.0.0 because `tibble::add_row()` is
  used with the new behavior.

# monoClust 1.0.0

* Package is now fully working with all features intended.
* Rewrote a lot of internal functions to clean up unnecessary codes.
* Added ggplot2 versions of CV plot.
* PCP plot for circular data is now named `ggpcp` and uses ggplot2.
* New circular add/subtract operators: `%cd+%`, `%cd-%` (in degree), `%cr+%`, 
  and `%cr-%` (in radian)
* Added `wind_sensit_2008` data set.
* Added a vignette using R Markdown.
* Updated documentation.

## Changes to functions

* `MonoClust()` removed `perm.test` and `alpha` arguments. Users now need to run 
  `perm.test()` separately to perform permutation test on a MonoClust.
  * Removed `labels`, `corders`, `ran` from the output.
  * `Membership` and `Dist` outputs of `MonoClust()` are now `membership` and 
    `dist`, respectively.
* `plot.MonoClust()` added `uniform`, `branch`, `minbranch`, `stats`, 
  `cols.type`, and `show.pval`.
* `cv.test()` now returns unified output for both LOOCV and k-fold. It includes 
  MSE and SE table and a note with the type of cross-validation. 
* `abbrev` argument in `plot.MonoClust()` and `print.MonoClust()` now accepts 
  descriptive options `"no"`, `"short"`, `"abbreviate"`.
* `method` argument in `perm.test()` now uses descriptive options `"sw"`, `"rl"` 
  and `"rn"`. Users can now decide whether to apply Bonferroni correction with 
  `bon.adj` argument.
* `predict.MonoClust()` removed `na.action`. `type` argument now accepts more 
  meaningful options `"centroid"` or `"medoids"`.
* Added parallel processing capability to `perm.test()` and `cv.test()`.

# monoClust 0.4.0

* Removed the categorical variable feature to make sure it works well for all 
other cases.

# monoClust 0.3.0

* Added support for `foreach` when searching for the best split.
* Added a `NEWS.md` file to track changes to the package.
