% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_traj.R
\name{get_traj}
\alias{get_traj}
\title{Conducts interpolated HYSPLIT trajectory runs}
\usage{
get_traj(
  trajDate,
  dur = 8,
  OutDir,
  Lat = 38.104841,
  Lon = -122.245442,
  Dir = "backward",
  Met = "reanalysis",
  ExecDir = NULL,
  MetDir = NULL,
  ht = 1000
)
}
\arguments{
\item{trajDate}{character string of dates in YYYY-MM-DD format}

\item{dur}{the duration of each trajectory model run in hours}

\item{OutDir}{a file path to where the trajectory output files will be stored}

\item{Lat}{the starting latitude (in decimal degrees)}

\item{Lon}{the starting longitude (in decimal degrees)}

\item{Dir}{the direction of the trajectory model run. Options are: 'forward'
or 'backward'}

\item{Met}{an option to select meteorological data files. The options are
gdas1 (Global Data Assimilation System 1-degree resolution data),
reanalysis (NCAR/NCEP global reanalysis data), and narr (North American
Regional Reanalysis). Visit the NOAA meteorological data archives for more
information (https://www.ready.noaa.gov/archives.php).}

\item{ExecDir}{an optional file path to the hyts_std executable file, if the
hyts_std file is not in the working directory.}

\item{MetDir}{an optional file path to store or find existing meteorological
data files}

\item{ht}{the starting height for the model run (in meters above ground
level)}
}
\value{
A dataframe of interpolated HYSPLIT trajectory points.
}
\description{
The function modifies the hysplit_trajectory() function from the
SplitR package to return a data frame of interpolated HYSPLIT trajectories.
This allows trajectory points to be interpolated down to 5 minute
intervals, rather than every hour. intervals, rather than every hour.
}
