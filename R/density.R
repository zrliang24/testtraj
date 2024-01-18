#' Create frequency plot
#'
#' @description Calculates conditional frequency above or below a given fraction
#' @param traj Data frame of trajectory points
#' @param coordinates Vector of latlon coordinates associated with trajectories
#' @param gritty Blank grid of the extent of trajectories
#' @param siteDF Data frame containing the emission concentrations and
#'   coordinates
#' @param fraction Number between 0 and 1 to calculate conditional frequency
#'   above or below
#' @param type Whether to calculate above or below the fraction
#' @param remove.outliers Removes all points above or below the fraction.
#'   Default set to true.
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom sensR findcr
#' @importFrom stats sd
#' @importFrom stats quantile
#' @importFrom stats median
#'
#' @export
#' 
conditional_frequency_fraction <- function(
    traj, coordinates, gritty, siteDF, fraction = 0.16, type = "above",
    remove.outliers = T) {
  
  # Get start and end dates from trajectories
  start_date <- min(traj$trajTime)
  end_date <- max(traj$trajTime)
  
  
  if (remove.outliers) {
    
    # Remove low end if no MDLs
    cutpoint <- median(siteDF$PM) - sd(
      subset(siteDF, PM > quantile(siteDF$PM, 0.16) &
               PM < quantile(siteDF$PM, 0.84))$PM)
    
    if (type == "above") {
      siteDF <- filter(siteDF, PM >= cutpoint)
    } else if (type == "below") {
      siteDF <- filter(siteDF, PM <= cutpoint)
    } else {
      stop("Allowable types are 'above' and 'below'")
    }
    
  }
  
  
  # Subset traj df to only dates above cutoff
  lim_dates <- siteDF$date
  # trajAbr <- subset(traj, trajTime %in% round(statDF$date, units = 'hours'))
  lim_traj <- subset(traj, trajTime %in% round(lim_dates, units = 'hours'))
  lim_traj_coords <- cbind(lim_traj$lon, lim_traj$lat)
  
  # Get unrestricted and special-case timestamp density
  unrestricted <- rasterize(coords, gritty, fun = 'count')
  limited <- rasterize(lim_traj_coords, gritty, fun = 'count')
  
  # Calculate conditional frequency (the ratio of limited to
  # unrestricted)
  conditional <- limited / unrestricted
  
  # Replace NA with 0 for calculating critical value
  unrestricted[is.na(unrestricted[])] <- 0
  

  # Find the critical value at each grid cell
  critical <- app(unrestricted, fun = function(a) {
    findcr2(a, frac = fraction)})
  
  # Values of the limited set that can be considered significant
  # Just testing > for curiosity here - shouldn't it be >=
  mask <- limited > critical
  
  # Finally, return the conditional frequency gated by the
  # significance test
  cond <- conditional * mask
  cond[cond[] == 0] <- NA
  return(cond)
  
}


significant_conditional_frequency <- function(
    traj, coordinates, gritty, site, hitRate = 0.5) {
  
  # Get start and end dates from trajectories
  start_date <- min(traj$trajTime)
  end_date <- max(traj$trajTime)
  
  # Determine dates and trajectories for the fraction of interest
  # arbitrary cutpoint for purple air data
  cutpoint <- median(site$PM) - sd(
    subset(site, PM > quantile(site$PM, 0.16) &
             PM < quantile(site$PM, 0.84))$PM)
  dfAbr <- subset(site, PM >= cutpoint)
  
  # Subset traj df to only dates above cutoff
  lim_dates <- site$date
  lim_traj <- subset(traj, trajTime %in% lim_dates)
  lim_traj_coords <- cbind(lim_traj$lon, lim_traj$lat)
  
  # Get unrestricted and special-case timestamp density
  unrestricted <- rasterize(coords, gritty, fun = 'count')
  limited <- rasterize(lim_traj_coords, gritty, fun = 'count')
  
  # Calculate conditional frequency (the ratio of limited to unrestricted)
  conditional <- limited / unrestricted
  
  # Replace NA with 0 for calculating critical value
  unrestricted[is.na(unrestricted[])] <- 0
  
  # Find the critical value at each grid cell
  # critical <- critical_value(unrestricted, hitRate)
  critical <- critical_value2(unrestricted, nrow(dfAbr)-2)
  
  # We want to encode three possibilities - significant/high,
  # not significant/high, not high
  high <- conditional >= hitRate
  
  # Not sure if this should technically be >=, using > reduces noise
  significant <- limited > critical
  
  significant_high <- high * significant
  
  # significant/high = 2; not-significant/high = 1; not high = NA
  result <- high + significant_high
  result[result[] == 0] <- NA
  
  return(result)
  
}

# Helper Function
findcr2 <- function(x, frac) {
  
  crValues <- unlist(lapply(x, function(a) {
    
    if (a == 0) {
      
      return(0)
    } else {
      sensR::findcr(a, alpha = 0.05, p0 = frac, pd0 = 0,
                    test = "difference")
    }
  }))
  
  return(crValues)
  
}
