
## helper functions
## parameters
nD <- 1258
nR <- 1142
n <- nD + nR
tau <- 0.5
max.time <- 1000 ## time unit is a day
nsims <- 100 ## no. of replicated simulations

## set no. of cores
ncores <- parallel::detectCores(logical = F)

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


## function to find the point of which s.num == i.num from object class netsim

find.point.to.plot <- function(netsim) {

  setDT(modelDT <- as.data.frame(netsim))
  overall <- modelDT[s.num - i.num < 0.05, time][1] - 1

    timeD <- modelDT[s.num.pidD - i.num.pidD < 0.05, time][1] - 1
    timeR <- modelDT[s.num.pidR - i.num.pidR < 0.05, time][1] - 1
    time <- c(overall = overall, timeD = timeD, timeR = timeR)

    return(time)
}

## function to print prevalence, network, and compartment plots from netsim into one pdf
print.plots.pdf <- function(netsim,
                            network.type = c("Bernoulli", "tree", "small-world", "homophilous"),
                            plot.title,
                            include.exposed = F) {

  cutoff <- find.point.to.plot(netsim)
  titles <- c(paste("Overall prevalence,", network.type, "network", sep = " "),
              paste("Prevalence among", c("Democrats,", "Republicans,"),  network.type, "network", sep = " "))
  ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)

  pdf(paste0(plot.title, ".pdf"), width = 10, height = 8)

  if (include.exposed == T) {
    plot(netsim, y = c("s.num", "e.num", "i.num"),
         mean.col = c('grey90', "grey50", "black"), qnts.col = c('grey90', "grey50", "black"),
         legend = F, popfrac = T, mean.smooth = T, qnts = 0.95)
    legend("right", legend = c("Infected: overall", "Suspected: overall", "Exposed: overall"), lwd = 3,
           col = c("black", "grey90", "grey50"), bg = "white", bty = "n"); title(titles[1])

    plot(netsim, y = c("s.num.pidD", "e.num.pidD", "i.num.pidD"),
         mean.col = c('grey90', "grey50", "black"), qnts.col = c('grey90', "grey50", "black"),
         legend = F, popfrac = T, mean.smooth = T, qnts = 0.95)
    legend("right", legend = c("Infected: Dem", "Suspected: Dem", "Exposed: Dem"), lwd = 3,
           col = c("black", "grey90", "grey50"), bg = "white", bty = "n"); title(titles[2])

    plot(netsim, y = c("s.num.pidR", "e.num.pidR", "i.num.pidR"),
         mean.col = c('grey90', "grey50", "black"), qnts.col = c('grey90', "grey50", "black"),
         legend = F, popfrac = T, mean.smooth = T, qnts = 0.95)
    legend("right", legend = c("Infected: Rep", "Suspected: Rep", "Exposed: Rep"), lwd = 3,
           col = c("black", "grey90", "grey50"), bg = "white", bty = "n"); title(titles[3])
  } else {
    ## only for first configuration ("SI model")
    plot(netsim, mean.col = c('grey', "black"), qnts.col = c('grey', "black"),
         legend = F, popfrac = T, qnts = 0.95)
    abline(v = cutoff[1], lty = 2); text(cutoff[1] - 40, 0.05, paste0("t = ", cutoff[1]), col = "black")
    legend("right", legend = c("Exposed: overall", "Suspected: overall"), lwd = 3,
           col = c("black", "grey"), bg = "white", bty = "n"); title(titles[1])

    plot(netsim, y = c("s.num.pidD", "i.num.pidD"), mean.col = c('grey', "black"), qnts.col = c('grey', "black"),
         legend = F, popfrac = T, qnts = 0.95)
    abline(v = cutoff[2], lty = 2);  text(cutoff[2] - 40, 0.05, paste0("t = ", cutoff[2]), col = "black")
    legend("right", legend = c("Exposed: Dem", "Suspected: Dem"), lwd = 3,
           col = c("black", "grey"), bg = "white", bty = "n"); title(titles[2])

    plot(netsim, y = c("s.num.pidR", "i.num.pidR"), mean.col = c('grey', "black"), qnts.col = c('grey', "black"),
         legend = F, popfrac = T, qnts = 0.95)
    abline(v = cutoff[3], lty = 2);  text(cutoff[3] - 40, 0.05, paste0("t = ", cutoff[3]), col = "black")
    legend("right", legend = c("Exposed: Rep", "Suspected: Rep"), lwd = 3,
           col = c("black", "grey"), bg = "white", bty = "n"); title(titles[3])
  }

  par(mfrow = c(1,2), mar = c(0,0,0,0))

  nw1 <- get_network(netsim, collapse = T, at = 1)
  cols <- ifelse(get.vertex.attribute.active(nw1, "testatus", at = 1) == "i", "grey20", "gray90")
  vertex.cex <- ifelse(get.vertex.attribute.active(nw1, "testatus", at = 1) == "i", 0.6, 0.4)
  vertex.cex[isolates(nw1)] <- 0.2; vertex.cex[cols == "grey20"] <- 0.6
  plot(nw1, mode = "fruchtermanreingold", displayisolates = T, vertex.lwd = 0.1, edge.lwd = 0.2,
       vertex.col = cols, vertex.border = "grey40", edge.col = "grey40", vertex.cex = vertex.cex)
  title("Prevalence at t1", line = -2)

  nw2 <- get_network(netsim, collapse = T, at = cutoff[1])
  cols <- ifelse(get.vertex.attribute.active(nw2, "testatus", at = cutoff[1]) == "i", "grey20", "gray90")
  vertex.cex <- ifelse(get.vertex.attribute.active(nw2, "testatus", at = cutoff[1]) == "i", 0.6, 0.4)
  vertex.cex[isolates(nw2)] <- 0.2; vertex.cex[cols == "grey20"] <- 0.6
  plot(nw2, mode = "fruchtermanreingold", displayisolates = T, vertex.lwd = 0.1, edge.lwd = 0.2,
       vertex.col = cols, vertex.border = "grey40", edge.col = "grey40", vertex.cex = vertex.cex)
  title(paste0("Prevalence at t", cutoff[1]), line = -2)

  dev.off()
  on.exit(par(ops))
}

## function to get mean and 95% CIs as a data.frame
## param: 'netsim' a object class netsim
## param: 'qvals' a vector of length 2, quantiles to probe level of CIs
get.summary.stats <- function(netsim, qvals = c(0.025, 0.5, 0.975)) {
  col.names <- colnames(as.data.frame(netsim))[-1]

  dat.mean <- as.data.frame(netsim, out = "qnt", qval = qvals[2]) ## median
  colnames(dat.mean) <- c("time", paste(col.names, "mean", sep = "."))

  dat.llci <- as.data.frame(netsim, out = "qnt", qval = qvals[1]) ## 95% lower
  colnames(dat.llci) <- c("time", paste(col.names, "llci", sep = "."))

  dat.ulci <- as.data.frame(netsim, out = "qnt", qval = qvals[3]) ## 95% upper
  colnames(dat.ulci) <- c("time", paste(col.names, "ulci", sep = "."))

  modelDT <- Reduce(function(x, y) merge(x, y, by = "time"), list(dat.mean,
                                                                  dat.llci,
                                                                  dat.ulci))
  data.table::setDT(modelDT)
  modelDT
}


## sourcing specific lines from the .R syntax
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}
