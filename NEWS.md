# monoClust (development version)

## Changes to functions
* `plot.MonoClust()` allows the `cols` argument to have the length greater
  than the number of leaves. In that case, function will not throw error and 
  only a subset of it will be used.
  
## Fix bugs
* Fix some typos.

# monoClust 1.1.0

## Changes to functions
* Remove `bipartvar` column from `frame` object of MonoClust object. `var` 
  column should be sufficient for showing splitting variable names.
* `alt` column in `frame` object of MonoClust is now a nested tibble containing
  alternate split details. However, the package does not support specifying an
  alternate splitting route so users may have to run step-by-step by indicating
  `nclusters = 2` in `MonoClust()` and then on each branch.
* Remove `arctic_2019` data set. It was not used anywhere in the examples.

## Improvements
* Add documentation for `MonoClust.object` to explain its structure.

## Fix bugs
* Min version of dependency `tibble()` is 3.0.0 because `tibble::add_row()` is
  used with the new behavior.
* Fix some typos and clarify some documentation.

# monoClust 1.0.0

* Package is now fully working with all features intended.
* Rewrote a lot of internal functions to clean up unnecessary codes.
* Added ggplot2 versions of CV plot.
* PCP plot for circular data is now named `ggpcp` and uses ggplot2.
* New circular add/subtract operators: `%cd+%`, `%cd-%` (in degree), `%cr+%`, 
  and `%cr-%` (in radian).
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
