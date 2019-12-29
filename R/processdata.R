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
  return(list(phy=pruned$phy, traits=pruned$data))
}

#' Discretize continuous trait data
#'
#' Converts a vector of numbers into a vector of 0 and 1 based on whether they are below or above some value.
#' There are two ways to do this: based on percentile or based on a numeric cutoff.
#' By default, it will separate it based on the 50th percentile (cutoff of 0.5), but you can change the cutoff value and whether it is used as percentile or trait value.
#'
#' @param x Vector of trait values
#' @param cutoff Value to use as cutoff. If percentile, 0.3 = 30th percentile, etc.
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
  if(length(unique(result))==1) {
    warning("No variation in output vector")
  }
  names(result) <- names(x)
  return(result)
}

#' Get sister groups for a node
#'
#' For a node, gives the taxa on each side. Note that the output is a data.frame with lists
#'
#' @param node Node number
#' @param phy A phylo object
#' @return a data.frame with the node numbers and columns with the tip labels of the two descendant clades
#' @export
sis_get_sister_pair <- function(node, phy) {
  descendants <- phangorn::Descendants(phy, node, type="children")
  if(length(descendants)==2) {
    return(data.frame(node=node, left=I((phangorn::Descendants(phy, descendants[1], type="tips"))), right=I((phangorn::Descendants(phy, descendants[2], type="tips")))))
  } else {
    return(c(NA, NA, NA))
  }
}

#' Get sister groups for all internal nodes
#'
#' For each node, return the vector of tip numbers for taxa on each side. It is sorted so that sister groups with fewer taxa are arranged at the top.
#'
#' @param phy A phylo object
#' @param ncores How many cores to use to run this in parallel
#' @return a data.frame with the node numbers and columns with the tip labels of the two descendant clades, plus additional info on the sister groups
#' @export
sis_get_sisters <- function(phy, ncores=parallel::detectCores()) {
  result <- do.call(rbind.data.frame,parallel::mclapply(unique(phy$edge[,1]), sis_get_sister_pair, phy=phy, mc.cores=ncores))
  result$ntax.left <- sapply(result$left, length)
  result$ntax.right <- sapply(result$right, length)
  result$ntax.total <- result$ntax.left + result$ntax.right
  result <- result[order(result$ntax.total, decreasing=FALSE),]
  return(result)
}


#' Is the taxon in one of the sister groups
#'
#' Utility function for tossing out taxa already used
#'
#' @param taxon Node number of taxon
#' @param sisters Data.frame from sis_get_sisters()
#' @return data.frame of whether the taxon is in the left or right sister group, or any
sis_find_taxon <- function(taxon, sisters) {
  result <- data.frame(left=sapply(sapply(sisters$left, is.element, taxon), any), right=sapply(sapply(sisters$right, is.element, taxon), any))
  result$any <- apply(result, 1, any)
  return(result)
}

#' Get comparison format
#'
#' Convert a data.frame of all sister groups (from sis_get_sisters) and a vector of 0 and 1 (with names equal to taxon names) to a data.frame with the sister groups that differ in traits.
#'
#' @param sisters Data.frame from sis_get_sisters()
#' @param trait vector of 0/1 data
#' @return data.frame where each row is a sister group comparison.
#' @export
#' @examples
#' data(geospiza, package="geiger")
#' cleaned <- sis_clean(geospiza$phy, geospiza$dat)
#' phy <- cleaned$phy
#' traits <- cleaned$traits
#' trait <- sis_discretize(traits[,1])
#' sisters <- sis_get_sisters(phy)
#' sisters_comparison <- sis_format_comparison(sisters, trait)
#' print(sisters_comparison)
sis_format_comparison <- function(sisters, trait) {
  sisters$left.trait <- lapply(sisters$left, sis_get_trait_values, phy=phy, trait=trait)
  sisters$right.trait <- lapply(sisters$right, sis_get_trait_values, phy=phy, trait=trait)
  sisters$left.unique <- lapply(sisters$left.trait, sis_get_monomorphic)
  sisters$right.unique <- lapply(sisters$right.trait, sis_get_monomorphic)
  sisters$ntax.trait0 <- NA
  sisters$ntax.trait1 <- NA
  for (i in sequence(nrow(sisters))) {
    if(!is.na(sisters$left.unique[i]) & !is.na(sisters$right.unique[i])) {
      if(sisters$left.unique[i]==0 & sisters$right.unique[i]==1) {
        sisters$ntax.trait0[i] <- sisters$ntax.left[i]
        sisters$ntax.trait1[i] <- sisters$ntax.right[i]
      }
      if(sisters$left.unique[i]==1 & sisters$right.unique[i]==0) {
        sisters$ntax.trait1[i] <- sisters$ntax.left[i]
        sisters$ntax.trait0[i] <- sisters$ntax.right[i]
      }
    }
  }
  return(sisters)
}

#' Get trait values for tip numbers
#'
#' @param nodes vector of node numbers (tip numbers, actually)
#' @param phy A phylo object
#' @param trait A trait vector with names equal to taxon names
sis_get_trait_values <- function(nodes, phy, trait) {
  return(I(trait[phy$tip.label[nodes]]))
}

#' Get monomorphic trait
#'
#' @param trait Vector of trait values
#' @return The state all taxa have if monomorphic; NA otherwise
sis_get_monomorphic <- function(trait) {
  unique_state <- unique(trait)
  if(length(unique_state)==1) {
    return(unique_state)
  } else {
    return(NA)
  }
}
