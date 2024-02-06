#' Creates a conditional frequency fraction map from a given region of Vallejo
#' and statistic
#'
#' @param region Region of Vallejo. Options are: North, South, Downtown,
#'   Foothills, Mare Island
#' @param statistic Options are: Ninetieth, Fiftieth, Tenth, Mean
#' @param basemap ggmap basemap
#' @param buffer Fraction out of 100 to exclude.
#'
#' @export
#' 
make_map <- function(region, 
                     statistic, 
                     basemap,
                     buffer = 10) {
  
  stats <- c('ninetieth',
             'fiftieth',
             'tenth',
             'mean')
  
  statFraction <- case_when(
    statistic == stats[1] ~ 0.9,
    statistic == stats[2] ~ 0.5,
    statistic == stats[3] ~ 0.1,
    statistic == stats[4] ~ 0,
    T ~ NA)
  
  # First, subset the PM data frame
  abr <- subset(siteDF, Region == region)
  
  if (is.na(statFraction)) {
    
    print(paste('Please choose from',
                paste(stats, collapse = ', ')))
    
    return(NULL)
    
  } else if (statFraction == 0) {
    
    statDF <- abr %>%
      mutate(Target = mean(PM),
             Upper = Target * (1 + buffer/100),
             Lower = Target * (1 - buffer/100)) %>%
      filter(between(PM, Lower, Upper)) %>%
      ungroup()
    
  } else {
    
    statDF <- abr %>%
      mutate(Target = unname(quantile(PM, statFraction)),
             Upper = Target * (1 + buffer/100),
             Lower = Target * (1 - buffer/100)) %>%
      filter(between(PM, Lower, Upper)) %>%
      ungroup()
    
  }
  
  trajAbr <- subset(traj, trajTime %in% round(statDF$date, units = 'hours'))
  
  # this is more straightforward and works too
  if (nrow(trajAbr) > 100000) {
    gritty <- rast(trajExt,
                   resolution = 0.001,
                   crs = '+proj=longlat +datum=WGS84')
  } else {
    gritty <- rast(trajExt,
                   resolution = 0.01,
                   crs = '+proj=longlat +datum=WGS84')
  }
  
  cff <- conditional_frequency_fraction(trajAbr,
                                        coords,
                                        gritty,
                                        statDF,
                                        remove.outliers = F)
  
  outMap <- ggmap(basemap) +
    tidyterra::geom_spatraster(mapping = aes(), 
                               data = cff, 
                               alpha = 0.5) +
    scale_fill_viridis_c(option = 'plasma') +
    labs(title = paste(region),
         subtitle = paste('statistic:', statistic, '|',
                          'buffer:', buffer, '|',
                          'n =', nrow(trajAbr)))
  # if (justmap) {
  #   print(outMap)
  # }
  
  return(outMap)
  
}