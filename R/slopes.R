library(tidyverse)

make_slopes <- function(data){
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
  out <- solve(mat, res)
  return(out)
}

interpolate <- function(data, dens){
  data <- mutate(data, AX = lead(X) - X, AY = lead(Y) - Y,
                 dt = lead(tme) - tme)
  slopeX <- c(NA, make_slopes(select(data, dt, tme, A=AX)), NA)
  slopeY <- c(NA, make_slopes(select(data, dt, tme, A=AY)), NA)
  data <- mutate(data, slopeX=slopeX, slopeY=slopeY)
  cdX <- cubicData(select(data, tme=tme, Y=X, S=slopeX))
  cdY <- cubicData(select(data, tme=tme, Y=Y, S=slopeY))
  out <- tibble(X=c(), Y=c(), tme=c())
  start <- data$tme[1]
  for (i in seq.int(1, length(data$tme) - 1)){
    tme <- seq(start, data$tme[i+1], dens)
    outData <- tibble(tme=tme,
                      X=call(c(cdX$cub[i], cdX$squ[i], cdX$sl[i], cdX$y[i]), tme),
                      Y=call(c(cdY$cub[i], cdY$squ[i], cdY$sl[i], cdY$y[i]), tme))
    out <- bind_rows(out, outData)
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