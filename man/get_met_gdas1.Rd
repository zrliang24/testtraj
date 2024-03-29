% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_met_gdas1.R
\name{get_met_gdas1}
\alias{get_met_gdas1}
\title{Get GDAS1 meteorology data files}
\usage{
get_met_gdas1(files = NULL, years = NULL, months = NULL, path_met_files)
}
\arguments{
\item{files}{a vector list of exact filenames for
the GDAS1 files.}

\item{years}{a vector list of years for which GDAS1
files are to be obtained via FTP.}

\item{months}{a vector list of months for which
GDAS1 files are to be obtained via FTP.}

\item{path_met_files}{a full path for the download
location of the meteorological data files.}
}
\description{
This function downloads GDAS1
meteorology data files from the NOAA FTP server and
saves them to a specified folder. Files can be
downloaded either by specifying a list of filenames
(in the form of "RP\link[=MM]{YYYY}.gbl") or through bulk
download of a year of files.
}
