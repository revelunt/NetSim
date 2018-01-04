
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



## function to create compartment plot for SEI class model
comp_plot.SEI <- function (x, at = 1, digits = 3, ...) {

  nsteps <- x$control$nsteps
  dis.type <- "SEI"
  vital <- x$param$vital
  if (class(x) == "icm") {
    groups <- x$param$groups
  }
  if (class(x) == "netsim") {
    groups <- x$param$modes
  }
  if (groups != 1) {
    stop("Only 1-group/mode models currently supported",
         call. = FALSE)
  }
  if (at > nsteps | at < 1) {
    stop("Specify a timestep between 1 and ", nsteps, call. = FALSE)
  }

  df.mn <- as.data.frame(x, out = "mean")
  df.mn <- round(df.mn[at == df.mn$time, ], digits)
  df.sd <- as.data.frame(x, out = "sd")
  df.sd <- round(df.sd[at == df.sd$time, ], digits)
  ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)
  par(mar = c(0, 0, 2, 0))
  options(scipen = 10)
  plot(0:100, 0:100, type = "n", axes = FALSE)
  title(main = paste(dis.type, "Model Diagram"))
  mtext(paste0("Simulation means(sd) | time=", at), side = 3,
        cex = 0.8, line = -1)

  EpiModel:::mbox(5, 40, "Susceptible", paste0(df.mn$s.num, "(", df.sd$s.num, ")"))
  EpiModel:::mbox(40, 40, "Exposed", paste0(df.mn$e.num, "(", df.sd$e.num, ")"))
  EpiModel:::mbox(75, 40, "Infected", paste0(df.mn$i.num, "(", df.sd$i.num, ")"))
  EpiModel:::harrow(5, 40, "se.flow", df.mn$se.flow, dir = "right")
  EpiModel:::harrow(40, 40, "ei.flow", df.mn$ei.flow, dir = "right")
  if (vital == TRUE) {
    EpiModel:::varrow(5, 40, "ds.flow", df.mn$ds.flow, dir = "out")
    EpiModel:::varrow(40, 40, "de.flow", df.mn$de.flow, dir = "out")
    EpiModel:::varrow(75, 40, "di.flow", df.mn$di.flow, dir = "out")
    EpiModel:::varrow(5, 40, "b.flow", df.mn$b.flow, dir = "in")
  }

  on.exit(par(ops))
}



## function to create mixing matrix from the object type igraph
## calculate the mixing matrix of in igraph graph object 'mygraph', by some vertex attribute 'attrib'
## can change the default use.density=FALSE to return a matrix with raw number of edges rather than density

mixingmatrix <- function(nw, attr, use.density = TRUE) {

  require(igraph)

  # get unique list of characteristics of the attribute
  attlist <- sort(unique(get.vertex.attribute(nw,attr)))

  numatts <- length(attlist)

  # build an empty mixing matrix by attribute
  mm <- matrix(nrow=numatts,
               ncol=numatts,
               dimnames=list(attlist,attlist))

  # calculate edge density for each matrix entry by pairing type
  # lends itself to parallel if available
  el <- get.edgelist(nw,names=FALSE)
  for (i in 1:numatts) {
    for (j in 1:numatts) {
      mm[i,j] <- length(which(apply(el,1,function(x) {
        get.vertex.attribute(nw, attr, x[1] ) == attlist[i] &&
          get.vertex.attribute(nw, attr, x[2] ) == attlist[j]  } )))
    }
  }

  # convert to proportional mixing matrix if desired (by edge density)
  if (use.density) mm/ecount(nw) else mm

  # if undirected, set lower diagonal to NaN (upper diagonal is the heterophilous ties)
  if (!isTRUE(is_directed(nw))) mm[lower.tri(mm, diag = FALSE)] <- NaN
  mm
}
