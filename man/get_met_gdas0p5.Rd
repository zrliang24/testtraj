% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_met_gdas0p5.R
\name{get_met_gdas0p5}
\alias{get_met_gdas0p5}
\title{Get GDAS0.5 meteorology data files}
\usage{
get_met_gdas0p5(
  files = NULL,
  years = NULL,
  months = NULL,
  days = NULL,
  path_met_files
)
}
\arguments{
\item{files}{a vector list of exact filenames for
the GDAS0.5 files.}

\item{years}{a vector list of years for which
GDAS0.5 files are to be obtained via FTP.}

\item{months}{a vector list of months for which
GDAS0.5 files are to be obtained via FTP.}

\item{days}{a vector list of days for which GDAS0.5
files are to be obtained via FTP.}

\item{path_met_files}{a full path should be provided
for the location of the meteorological data files;
downloaded files will be saved in this location.}
}
\description{
This function downloads GDAS0.5
meteorology data files from the NOAA FTP server and
saves them to a specified folder. Files can be
downloaded either by specifying a list of filenames
(in the form of "RP\link[=MM]{YYYY}.gbl") or through bulk
download of a year of files.
}
