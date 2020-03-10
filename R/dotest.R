#' Try but returns NA rather than error
#' @param code Code to run
#' @param silent Print error if TRUE
tryNA <- function(code, silent = FALSE) {
  tryCatch(code, error = function(c) {
    msg <- conditionMessage(c)
    if (!silent) message(c)
    invisible(NA)
  })
}


#' Compute multiple tests based on sister group comparisons
#' @param pairs Data.frame with one row per sister group comparison, with one column for number of taxa in state 0, and one column for the number of taxa in state 1.
#' @param drop_matches Drop sister group comparisons with equal numbers of taxa
#' @param warn Some tests will fail with warnings (too few sister groups or other reasons). Setting this to FALSE will suppress those
#' @return A vector with the results of many tests, as well as summary data for the comparisons
#' @export
#' @examples
#' data(geospiza, package="geiger")
#' cleaned <- sis_clean(geospiza$phy, geospiza$dat)
#' phy <- cleaned$phy
#' traits <- cleaned$traits
#' trait <- sis_discretize(traits[,1])
#' sisters <- sis_get_sisters(phy)
#' sisters_comparison <- sis_format_comparison(sisters, trait, phy)
#' pairs <- sis_format_simpified(sisters_comparison)
#' sis_test(pairs)
sis_test <- function(pairs, drop_matches=TRUE, warn=TRUE) {
  if(!warn) {
    oldw <- getOption("warn")
    options(warn = -1)
  }
  equal_count <- length(which(pairs$ntax.trait1==pairs$ntax.trait0))
  if(drop_matches) {
    pairs <- pairs[which(pairs$ntax.trait0!=pairs$ntax.trait1),]
  }
  kafermousset_1_derived <- data.frame(m=pairs$ntax.trait0+pairs$ntax.trait1, d=pairs$ntax.trait1)
  kafermousset_0_derived <- data.frame(m=pairs$ntax.trait0+pairs$ntax.trait1, d=pairs$ntax.trait0)

  result <- c(
    number.comparisons.trait0.bigger = length(which(pairs$ntax.trait0>pairs$ntax.trait1)),
    number.comparisons.trait1.bigger = length(which(pairs$ntax.trait1>pairs$ntax.trait0)),
    number.comparisons.trait0.equal.trait1 = equal_count,
    median.proportion.in.state.zero = stats::median(pairs$ntax.trait0/(pairs$ntax.trait0 + pairs$ntax.trait1)),
    median.ntax.diff.zero.minus.one = stats::median(pairs$ntax.trait0 - pairs$ntax.trait1),
    pvalue.sign.test = stats::binom.test(x=length(which(pairs$ntax.trait0<pairs$ntax.trait1)), n=length(which(pairs$ntax.trait0!=pairs$ntax.trait1)), alternative="two")$p.value,
    pvalue.diversity.contrast.ratiolog = tryNA(ape::diversity.contrast.test(pairs, method = "ratiolog")),
    pvalue.diversity.contrast.proportion = tryNA(ape::diversity.contrast.test(pairs, method = "proportion")),
    pvalue.diversity.contrast.difference = tryNA(ape::diversity.contrast.test(pairs, method = "difference")),
    pvalue.diversity.contrast.logratio = tryNA(ape::diversity.contrast.test(pairs, method = "logratio")),
    pvalue.slowinskiguyer.test = tryNA(ape::slowinskiguyer.test(pairs, detail=FALSE)$P.val),
    pvalue.mcconwaysims.test = tryNA(ape::mcconwaysims.test(pairs)$P.val),
    pvalue.richness.yule.test = tryNA(ape::richness.yule.test(pairs, rep(1e3, nrow(pairs)))$P.val), #following Paradis 2011, setting t to arbitrarily large size
    pvalue.kafermousset.1.derived = tryNA(min(1,scc.test(kafermousset_1_derived)$p.value)),
    pvalue.kafermousset.0.derived = tryNA(min(1,scc.test(kafermousset_0_derived)$p.value))

  )
  if(!warn) {
    options(warn = oldw )
  }
  return(result)
}

#' Do a test with a single cutoff value
#' @param cutoff Value to use as cutoff. If percentile, 0.3 = 30th percentile, etc.
#' @param x Vector of continuous trait values
#' @param use_percentile If TRUE, use cutoff as percentile
#' @param phy A phylo object
#' @param sisters Data.frame from sis_get_sisters()
#' @param warn Some tests will fail with warnings (too few sister groups or other reasons). Setting this to FALSE will suppress those
#' @return vector of outpout from sis_test()
sis_iterate_single_run <- function(cutoff, x, use_percentile=TRUE, phy, sisters=sis_get_sisters(phy), warn=FALSE) {
  trait <- sis_discretize(x, cutoff=cutoff, use_percentile=use_percentile)
  comparison <- sis_format_simpified(sis_format_comparison(sisters, trait, phy))
  test <- sis_test(comparison, warn=warn)
  test["ntax0"] <- as.numeric(length(which(trait==0)))
  test["ntax1"] <- as.numeric(length(which(trait==1)))
  if(use_percentile) {
    absolute_cutoff <- stats::quantile(x, probs=cutoff, na.rm=TRUE)
  } else {
    absolute_cutoff <- cutoff
  }
  test["absolute.cutoff"] <- as.numeric(unname(absolute_cutoff))
  return(test)
}

#' Iterate tests trying a variety of cutoff values
#'
#' This is a way of looking at the effect of using different cutoff values on the sister group comparisons. Do clades with a higher value have more species than their sister, and is this robust to what cutoff value is used? At the extremes (the min and max value) this is almost certainly not the case, unless you have many taxa with the same maximum or minimum values.
#'
#' This is a very dangerous function to use. Someone could use this to find the perfect cutoff value to find a significant result. This is one of the many forms of p-hacking. \strong{So, if you use this function and then report on significance using some cutoff, you MUST mention somewhere in your manuscript that you've tried a variety of cutoff values, and include a discussion of why you used a particular cutoff.} Ideally, you should have some biological intuition about what cutoff value is reasonable before using this function, as well.
#' @param x Vector of continuous trait values
#' @param nsteps Number of thresholds to try
#' @param phy A phylo object
#' @param sisters Data.frame from sis_get_sisters()
#' @return A data.frame, where each column is for a different cutoff percentile and every row is a number returned from sis_test()
#' @export
#' @examples
#' data(geospiza, package="geiger")
#' cleaned <- sis_clean(geospiza$phy, geospiza$dat)
#' phy <- cleaned$phy
#' trait.x <- cleaned$traits[,1]
#' sis_iterate(trait.x, phy=phy)
sis_iterate <- function(x, nsteps=11, phy, sisters=sis_get_sisters(phy)) {
  cutoffs <- seq(from=0, to=1, length.out=nsteps+2)
  cutoffs <- cutoffs[-1] #get rid of extreme
  cutoffs <- cutoffs[-length(cutoffs)] #get rid of extreme
  results <- sapply(cutoffs, sis_iterate_single_run, x=x, phy=phy, sisters=sisters)
  colnames(results) <- cutoffs
  return(results)
}
