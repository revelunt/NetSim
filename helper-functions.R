
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


## function to get no. of "exposed" from netsim object

get.exposed.sim <- function(sim) {
  if(!(class(sim) == "netsim")) {
    stop("requires 'netsim' object from EpiModel with custom modules. Check your data...")
  }

  ninit <- length(sim$epi)
  epi.by <- sim$control$epi.by

  ## add number of those being exposed
  vars.to.add <- c("e.num")

  ## if there's "epi.by" arguments
  if(!(is.null(epi.by) == TRUE)) {
    category <- names(table(get.vertex.attribute.active(sim$network$sim1,epi.by, at = 1)))
    ncat <- length(category)
    vars.to.add <- c(paste(paste(vars.to.add, epi.by, sep = "."), category, sep = ""))
    base.vars <- names(sim$epi)[grepl("num\\.[[:alpha:]]+", names(sim2$epi))]

    ## add placeholder
    for (var in vars.to.add) {
      sim$epi[[var]] <- as.data.frame(matrix(nrow = dim(sim$epi[[1]])[1], ncol = dim(sim$epi[[1]])[2]))
    }

    ## total n = suspected + exposed + infected
    ## exposed = total n - (suspected + infected)
    for (k in seq_len(ncat)) {
      sim$epi[[ninit + k]] <-
        sim$epi[[base.vars[k + ncat*2]]] - (sim$epi[[base.vars[k]]] + sim$epi[[base.vars[k + ncat]]])
    }


  }
  ## return dataset
  return(sim)
}
