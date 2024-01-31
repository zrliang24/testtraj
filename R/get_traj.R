#' Conducts interpolated HYSPLIT trajectory runs
#'
#' @description The function modifies the hysplit_trajectory() function from the
#'   SplitR package to return a data frame of interpolated HYSPLIT trajectories.
#'   This allows trajectory points to be interpolated down to 5 minute
#'   intervals, rather than
#'   every hour.
#' @param trajDate character string of dates in YYYY-MM-DD format
#' @param OutDir a file path to where the trajectory output files will be stored
#' @param dur the duration of each trajectory model run in hours
#' @param Lat the starting latitude (in decimal degrees)
#' @param Lon the starting longitude (in decimal degrees)
#' @param Dir the direction of the trajectory model run. Options are: 'forward'
#'   or 'backward'
#' @param Met an option to select meteorological data files. The options are
#'   gdas1 (Global Data Assimilation System 1-degree resolution data),
#'   reanalysis (NCAR/NCEP global reanalysis data), and narr (North American
#'   Regional Reanalysis). Visit the NOAA meteorological data archives for more
#'   information (https://www.ready.noaa.gov/archives.php).
#' @param MetDir an optional file path to store or find existing meteorological
#'   data files
#' @param ExecDir an optional file path to the hyts_std executable file, if the
#'   hyts_std file is not in the working directory.
#' @param ht the starting height for the model run (in meters above ground
#'   level)
#' @importFrom magrittr "%>%"
#' @importFrom purrr map_dfr
#' @returns A dataframe of interpolated HYSPLIT trajectory points.
#' @export

get_traj <- function(trajDate,
                     dur = 8,
                     OutDir,
                     Lat = 38.104841,
                     Lon = -122.245442,
                     Dir = 'backward',
                     Met = 'reanalysis',
                     ExecDir = NULL,
                     MetDir = NULL,
                     ht = 1000) {
  
  
  print(trajDate)
  # Added time out to prevent hysplit from crashing
  Sys.sleep(1)
  
  # Run hysplit
  newtraj <- hysplit_trajectory(
    lat = Lat,
    lon = Lon,
    out_dir = OutDir,
    height = ht,
    duration = dur,
    run_period = trajDate,
    daily_hours = seq(0,23),
    direction = Dir,
    met_type = Met,
    extended_met = F,
    met_dir = MetDir,
    exec_dir = ExecDir,
    model_height = 20000)
  
  # Interpolate 5 min intervals using cubic.r and slopes.r
  # Here, "date" is the date/time stamp of each unique trajectory, 24
  # trajectories per day
  dataf <- unique(newtraj$date) %>%
    map_dfr(~ run_peter(., newtraj))
  
  return(dataf)

}

