
## helper functions

## function to generate normal truncated distribution
## given M, SD, min, and max

rtnorm <- function(n, mean = 0, sd = 1, min = 0, max = 1) {
  bounds <- pnorm(c(min, max), mean, sd)
  u <- runif(n, bounds[1], bounds[2])
  qnorm(u, mean, sd)
}

## function to generate distribution given a priory correlation with y vector
corr.vector <- function(y, rho, x = NULL) {
  # x is optional
  if (is.null(x)) x <- rnorm(length(y))
  y.perp <- residuals(lm(x ~ y))
  z <- rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
  return(z)
}


## function to derive target statistics

get.target.stats <- function(nD, nR, target.gden,
                             target.nodematch,
                             target.nodefactor.nR,
                             target.concurrent) {
  n <- nD + nR
  edges <- target.gden * n/2
  nodefactor <- target.nodefactor.nR * nR
  nodematch <- target.nodematch * edges
  concurrent <- target.concurrent * n
  c(edges = edges, nodematch = nodematch, nodefactor = nodefactor, concurrent = concurrent)
}

## function to create repeated rows and columns from a given value(vector)
rep.row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}

rep.col <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}
