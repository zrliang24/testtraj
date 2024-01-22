library(tidyverse)

# This is a helper function for finding the slopes of the data at each point
# that will ensure that the output has a continuous second derivative. It uses a
# linear equation derived separately
make_slopes <- function(data){ # Generate slopes at each data point for interpolation.
  l <- length(data$tme) - 2
  mat <- matrix(rep(0, l^2), ncol=l)
  res <- numeric(l)
  res <- 3 * (data$A/(data$dt^2))
  res <- res + lag(res)
  res <- res[seq.int(2, length(res) - 1)]
  for (n in seq.int(2, l-1)){
    mat[n, n-1] <- 1 / data$dt[n]
    mat[n, n+1] <- 1 / data$dt[n+1]
    mat[n, n] <- 2 * (mat[n, n-1] + mat[n, n+1])
  }
  mat[1, 2] <- 1 / data$dt[2]
  mat[1, 1] <- 2 * mat[1, 2] + 1 / data$dt[1]
  res[1] <- 3 * data$A[2] / data$dt[2] + data$A[1] / data$dt[1]
  mat[l, l-1] <- 1 / data$dt[l]
  mat[l, l] <- 2 * mat[l, l-1] + 1 / data$dt[l+1]
  res[l] <- 3 * data$A[l] / data$dt[l] + data$A[l+1] / data$dt[l+1]
  out <- solve(mat, res) # Slopes generated using linear equation
  return(out)
}
# Interpolates data in the form (X: Y: tme:) and returns a data set in the same form
#   with more data points.
# Uses cubic interpolation: between each pair of original points, the path is
#   parameterized by a pair of cubics in the variable (tme).
# At each initial point (where two cubics meet), the derivatives dX/dt and dY/dt
#   and the second derivatives are continuous.
interpolate <- function(data, dens){
  # Get increment sizes of x, y, and time.
  data <- mutate(data, AX = lead(X) - X, AY = lead(Y) - Y,
                 dt = lead(tme) - tme)

  # Use increment sizes to make slopes
  slopeX <- c(NA, make_slopes(select(data, dt, tme, A=AX)), NA)
  slopeY <- c(NA, make_slopes(select(data, dt, tme, A=AY)), NA)

  # Add slopes to data
  data <- mutate(data, slopeX=slopeX, slopeY=slopeY)

  # Get parametric cubic equations for x(t) and y(t) in each interval
  cdX <- cubicData(select(data, tme=tme, Y=X, S=slopeX))
  cdY <- cubicData(select(data, tme=tme, Y=Y, S=slopeY))

  # Initialize result data
  out <- tibble(X=c(), Y=c(), tme=c())
  start <- data$tme[1]
  for (i in seq.int(1, length(data$tme) - 1)){

    # Get time interval for cubic equation with index i
    tme <- seq(start, data$tme[i+1], dens)

    # Create data for this interval
    outData <- tibble(tme=tme,
                      X=call(c(cdX$cub[i], cdX$squ[i], cdX$sl[i], cdX$y[i]), tme),
                      Y=call(c(cdY$cub[i], cdY$squ[i], cdY$sl[i], cdY$y[i]), tme))

    # Add data to output set
    out <- bind_rows(out, outData)

    # Define start time of next interval (end of this interval + one increment)
    start <- last(tme) + dens
  }
  return(out)
}

run_peter <- function(trajTime, trajDF) {

  # Here, the "date" is actually a date/time associated
  # with when the trajectory is initiated from the receptor
  trajAbr <- subset(trajDF, date == trajTime)
  # Run peter's stuff to get the output data frame
  output <- trajAbr %>%
    select(tme = hour.inc, X = lon, Y = lat) %>%
    # interpolate fx does not like negative time values
    mutate(tme = -tme) %>%
    # interpolate by 5min intervals
    interpolate(1/12) %>%
    # change time back to POSIX format and add time zone
    mutate(time = as.POSIXct(-3600 * tme,
                             origin = trajTime,
                             tz = 'UTC'),
           trajTime = trajTime) %>%
    select(lon = X, lat = Y, time, trajTime)
  # The trajTime parameter will let us associate the PurpleAir
  # receptor concentration to the trajectory data frame

  return(output)
}
