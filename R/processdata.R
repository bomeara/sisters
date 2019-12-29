#' Clean up trait and tree
#'
#' Does basic formatting and cleanup: makes sure the taxa are the same order in both, makes sure row names of the data are taxa, etc. Relies on geiger's treedata function. The first_col_names is for software like hisse, where the first column is often taxon names.
#'
#' @param phy A phylo object
#' @param traits A data.frame of traits
#' @param first_col_names Boolean on whether the first column has names.
#' @export
#' @return a list with phy and traits elements
sis_clean <- function(phy, traits, first_col_names = FALSE) {
  if(first_col_names) {
    rownames(traits) <- traits[,1]
    traits <- traits[,-1]
  }
  pruned <- geiger::treedata(phy, traits, sort=TRUE, warnings=FALSE)
  return(list(phy=pruned$phy, traits=pruned$traits))
}

#' Discretize continuous trait data
#'
#' Converts a vector of numbers into a vector of 0 and 1 based on whether they are below or above some value.
#' There are two ways to do this: based on percentile or based on a numeric cutoff.
#' By default, it will separate it based on the 50th percentile (cutoff of 0.5), but you can change the cutoff value and whether it is used as percentile or trait value.
#'
#' @param x Vector of values
#' @param cutoff Value to use as cutoff. If percentile, 0.5 = 50th percentile, etc.
#' @param use_percentile If TRUE, use cutoff as percentile
#' @export
#' @return a vector of 0 and 1 (and NAs)
sis_discretize <- function(x, cutoff=0.5, use_percentile=TRUE) {
  if(use_percentile) {
    absolute_cutoff <- quantile(x, probs=cutoff, na.rm=TRUE)
  } else {
    absolute_cutoff <- cutoff
  }
  result <- as.numeric(x>absolute_cutoff)
  return(result)
}
