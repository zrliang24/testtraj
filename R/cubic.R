
# Get (coefficients of) quadratic function based on the first point, the second
#   point, and the slope at the second point.
quadratic1 <- function(p1, p2, s2){
  pdiff <- p2 - p1
  sp <- pdiff[2] / pdiff[1]
  squ <- (s2 - sp) / pdiff[1]
  sl <- s2 - 2 * squ * p2[1]
  y <- p1[2] - squ * (p1[1] ^ 2) - sl * p1[1]
  return(c(squ, sl, y))
}

# Get quadratic equation as before, but with the slope of the first point instead.
quadratic2 <- function(p1, s1, p2){
  pdiff <- p2 - p1
  sp <- pdiff[2] / pdiff[1]
  squ <- (sp - s1) / pdiff[1]
  sl <- s1 - 2 * squ * p1[1]
  y <- p1[2] - squ * (p1[1] ^ 2) - sl * p1[1]
  return(c(squ, sl, y))
}

# Get tibble containing the coefficients of the cubic equation for each interval
#   in the data.
cubicData <- function(data){
  # Initialize values from data
  pXi <- data$tme
  pYi <- data$Y
  si <- data$S
  p1X <- pXi
  p1Y <- pYi
  s1 <- si
  p2X <- lead(p1X)
  p2Y <- lead(p1Y)
  s2 <- lead(s1)
  # First iteration: find coefficient of x^3
  pdiffX <- p2X - p1X
  pdiffY <- p2Y - p1Y
  sp <- pdiffY / pdiffX
  sq1 <- (sp - s1) / pdiffX
  sq2 <- (s2 - sp) / pdiffX
  cub <- (sq2 - sq1) / pdiffX
  # Subtract cubic from points and slope
  p1Y <- p1Y - cub * (p1X ^ 3)
  p2Y <- p2Y - cub * (p2X ^ 3)
  s1 <- s1 - 3 * cub * (p1X ^ 2)
  # Find quadratic that traces p1Y, p2Y and s1. This should be the correct quadratic
  #  component since the cubic component has been removed.
  pdiffY <- p2Y - p1Y
  sp <- pdiffY / pdiffX
  squ <- (sp - s1) / pdiffX
  sl <- s1 - 2 * squ * p1X # Get x coefficient by subtracting the slope of the x^2 term
  # Get constant term by subtracting y-values of other terms.
  y <- p1Y - squ * (p1X ^ 2) - sl * p1X
  cub[1] <- 0 # First and last intervals are quadratic. No cubic term.
  # Use quadratic generators to get the functions in the first and last intervals
  quad <- quadratic1(c(pXi[1], pYi[1]), c(pXi[2], pYi[2]), si[2])
  squ[1] <- quad[1]
  sl[1] <- quad[2]
  y[1] <- quad[3]
  l <- length(cub) - 1
  cub[l] <- 0
  quad <- quadratic2(c(pXi[l], pYi[l]), si[l], c(pXi[l+1], pYi[l+1]))
  squ[l] <- quad[1]
  sl[l] <- quad[2]
  y[l] <- quad[3]
  return(tibble(cub=cub, squ=squ, sl=sl, y=y))
}

call <- function(cub, x){
  return(cub[1] * (x^3) + cub[2] * (x^2) + cub[3] * x + cub[4])
}


