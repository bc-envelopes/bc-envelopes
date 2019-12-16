#' Calculate overlap percentage between hypervolumes
#'
#' @param x - Hypervolume set
#'
#' @return Numeric
#' @export
#'
#' @examples
hypervolume_overlap_percent <- function(x)
{
  stopifnot(inherits(x, "HypervolumeList"))
  # calculate overlap percentage
  union_vol <- x@HVList$Union@Volume 
  intersect_vol <- x@HVList$Intersection@Volume
  
  pct_overlap <- (intersect_vol /union_vol) * 100
  
  #   print(glue::glue("
  # There are {x@HVList$HV1@Dimensionality} axes in the set.
  # Percent of overlap without weighing is {round(pct_overlap, 3)},
  #                    and with weighing is {round(pct_overlap / x@HVList$HV1@Dimensionality, 3)}"))
  return(pct_overlap)
}

