#' Compute multiple tests based on sister group comparisons
#' @param pairs Data.frame with one row per sister group comparison, with one column for number of taxa in state 0, and one column for the number of taxa in state 1.
#' @param drop_matches Drop sister group comparisons with equal numbers of taxa
#' @return A vector with the results of many tests, as well as summary data for the comparisons
#' @export
#' @examples
#' data(geospiza, package="geiger")
#' cleaned <- sis_clean(geospiza$phy, geospiza$dat)
#' phy <- cleaned$phy
#' traits <- cleaned$traits
#' trait <- sis_discretize(traits[,1])
#' sisters <- sis_get_sisters(phy)
#' sisters_comparison <- sis_format_comparison(sisters, trait)
#' pairs <- sis_format_simpified(sisters_comparison)
#' sis_test(pairs)
sis_test <- function(pairs, drop_matches=TRUE) {
  if(drop_matches) {
    pairs <- pairs[which(pairs$ntax.trait0!=pairs$ntax.trait1),]
  }
  kafermousset_1_derived <- data.frame(m=pairs$ntax.trait0+pairs$ntax.trait1, d=pairs$ntax.trait1)
  kafermousset_0_derived <- data.frame(m=pairs$ntax.trait0+pairs$ntax.trait1, d=pairs$ntax.trait0)

  result <- c(
    #median.proportion.in.state.zero = median(pairs$ntax.trait0/(pairs$ntax.trait0 + pairs$ntax.trait1)),
    #median.ntax.diff.zero.minus.one = median(pairs$ntax.trait0 - pairs$ntax.trait1),
    sign.test = min(1,2*(1-pbinom(length(which(pairs$ntax.trait0<pairs$ntax.trait1)), nrow(pairs), 0.5))),
    diversity.contrast.ratiolog = ape::diversity.contrast.test(pairs, method = "ratiolog"),
    diversity.contrast.proportion = ape::diversity.contrast.test(pairs, method = "proportion"),
    diversity.contrast.difference = ape::diversity.contrast.test(pairs, method = "difference"),
    diversity.contrast.logratio = ape::diversity.contrast.test(pairs, method = "logratio"),
    slowinskiguyer.test = ape::slowinskiguyer.test(pairs, detail=FALSE)$P.val,
    mcconwaysims.test = ape::mcconwaysims.test(pairs)$P.val,
    richness.yule.test = ape::richness.yule.test(pairs, rep(1e3, nrow(pairs)))$P.val, #following Paradis 2011, setting t to arbitrarily large size
    kafermousset.1.derived = min(1,scc.test(kafermousset_1_derived)$p.value),
    kafermousset.0.derived = min(1,scc.test(kafermousset_0_derived)$p.value)

  )
  return(result)
}
