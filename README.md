# sisters

This is an R package for doing sister group and other diversity comparisons. The main functions are to take a trait, discretize it into a binary trait if needed, and compute all sister group comparisons: groups where one clade is monomorphic for one state, and its sister is monomorphic for a different state.

The basic workflow is to get all sister pairs of taxa (which can be slow) then use this object to compute the diversity in different clades separated by traits.

To install:

```
devtools::install_github("bomeara/sisters")
```

Look at `?sis_test` for help on using the package. 

Please note that the 'sisters' project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.
