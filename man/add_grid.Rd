% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_grid.R
\name{add_grid}
\alias{add_grid}
\title{Add lat/lon grid to a model}
\usage{
add_grid(
  model,
  lat = NULL,
  lon = NULL,
  range = c(5, 5),
  division = c(0.5, 0.5),
  start_day = NULL,
  start_hour = NULL,
  end_day = NULL,
  end_hour = NULL,
  duration = NULL,
  heights = NULL,
  samp_type = "avg",
  samp_interval = 24,
  name = NULL
)
}
\arguments{
\item{model}{a SplitR modeling object}

\item{lat}{a latitude value in decimal degrees
for the point of reference on the grid.}

\item{lon}{a longitude value in decimal degrees
for the point of reference on the grid.}

\item{range}{the latitude and longitude range about
the \code{grid_ref}.}

\item{division}{the division distances across the
latitude and longitude ranges.}

\item{start_day}{the day that the grid will become
active and measuring particle concentrations. This
should take the form of a single-length vector for a
day (\code{"YYYY-MM-DD"}).}

\item{start_hour}{the associated hour for the
\code{start_day} variable, taking the form of a
single integer hour (from \code{0} to \code{23}).}

\item{end_day}{the day that the grid will cease to
be active and no longer measuring particle
concentrations. This should take the form of a
single-length vector for a day (\code{"YYYY-MM-DD"}).}

\item{end_hour}{the associated hour for the
\code{end_day} variable, taking the form of a
single integer hour (from \code{0} to \code{23}).}

\item{duration}{a length of time in hours that the
grid will remain active from the start date-time.}

\item{heights}{a vector of heights for which there
will be horizontal sampling grids.}

\item{samp_type}{the method of reporting for the
sampling grid. The default is \code{avg} for
reporting average concentrations at every sampling
interval. Other options are \code{snapshot} and
\code{max} for concentrations at the time of
sampling and for maximum concentrations, respectively.}

\item{samp_interval}{the sampling interval in units
of hours.}

\item{name}{an identifier for this set of grid
parameters.}
}
\description{
Create and add a grid of latitude and
longitude points to a model object
}
