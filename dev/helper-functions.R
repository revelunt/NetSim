
## helper functions
## parameters
nD <- 1258
nR <- 1142
n <- nD + nR
tau <- 0.5
max.time <- 1000 ## time unit is a day
nsims <- 100 ## no. of replicated simulations

## initial status
status.vector <- c(rep(0, nD), rbinom(nR, 1, 0.15))
status.vector <- ifelse(status.vector == 1, "i", "s")

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

## function to normalize the range to 0 - 1
range01 <- function(x){
  normx <- (x-min(x)) / (max(x) - min(x))
  normx
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

## W Rank-Sum test
WRS.test <- function(dat1, dat2) {

  ## find min of distribution that do not different from each other
  test.min <- sapply(1:1000, function(i) {
    index <- seq(1, i)
    test <- suppressWarnings(wilcox.test(dat1[index, se.flow.mean/s.num.mean], dat2[index, se.flow.mean/s.num.mean]))
    test$p.value
  }, simplify = T)

  min <- min(which(test.min < 0.05))

  test.max <- sapply(min:1000, function(i) {
    index <- seq(min, i)
    test <- suppressWarnings(wilcox.test(dat1[index, se.flow.mean/s.num.mean], dat2[index, se.flow.mean/s.num.mean]))
    test$p.value
  }, simplify = T)

  max <- min + max(which(test.max < 0.05)) - 1

  return(c(min = min, max = max))
}

## ---------------------- FUNCTIONS FOR SIMULATIONS ----------------- ##

## infection module

infect <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  tea.status <- dat$control$tea.status
  nw <- dat$nw

  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)

  nElig <- length(idsInf)
  nInf <- 0

  if(!is.null(dat$param$pid.diff.rate)) {
    individual.diff.rate <- (dat$param$pid.diff.rate)/2
  } else {
    individual.diff.rate <- NULL
  }

  if (nElig > 0 && nElig < nActive) {
    del <- discord_edgelist(dat, at)
    if (!(is.null(del))) {
      ## if partisan differences in infection probabilities are assumed
      if (!is.null(individual.diff.rate)) {
        ## get individualized infection probabilities by adjusting from the baseline
        ## get pid vector of those who are infected
        Infpid <- (nw %v% "pid")[del$inf]
        del$transProb <- dat$param$inf.prob *
          ## if Republican, increase by the half of the size of pid.diff.rate,
          ## and if Democrat, decrease by the half of the szie of pid.diff.rate
          ifelse(Infpid == "R", (1 + individual.diff.rate), (1 - individual.diff.rate))
      } else {
        del$transProb <- dat$param$inf.prob
      }

      del$actRate <- dat$param$act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)

      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        dat$attr$infTime[idsNewInf] <- at

        if (tea.status == TRUE) {
          ## update network dynamic object
          nw <- activate.vertex.attribute(nw, prefix = "testatus", value = "e",
                                          onset = at, terminus = Inf, v = idsNewInf)
        }
      }
    }
  }

  ## update status for attr dataset
  dat$attr$status <- status

  if (at == 2) {
    dat$epi$se.flow <- c(0, nInf)
  }
  else {
    dat$epi$se.flow[at] <- nInf
  }
  dat$nw <- nw
  return(dat)
}

## Progression module
exposed.to.infectious.progress <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  nw <- dat$nw

  r.ei.rate <- dat$param$r.ei.rate
  d.ei.rate <- dat$param$d.ei.rate

  ## E to I progression
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "e")
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {

    Infpid <- (nw %v% "pid")[idsEligInf] ## get pid of those who are eligible for progression to infectious status
    ei.rate <- ifelse(Infpid == "D", d.ei.rate, r.ei.rate)

    ## discount the progression rates by accuracy (i.e., nfc)
    Infnfc <- (nw %v% "nfc")[idsEligInf]
    ## when nfc is 0.5, this gives identical ei.rate, and as nfc goes higher, ei.rate goes lower.
    ei.rate <- ei.rate^(Infnfc/0.5)
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }

  dat$attr$status <- status

  if (at == 2) {
    dat$epi$ei.flow <- c(0, nInf)
    dat$epi$e.num <- c(0, sum(active == 1 & status == "e"))

  }
  else {
    dat$epi$ei.flow[at] <- nInf
    dat$epi$e.num[at] <- sum(active == 1 & status == "e")
  }

  return(dat)
}

## Get_prev Module
get_prev.exposed.included <- function (dat, at) {

  active <- dat$attr$active
  l <- lapply(1:length(dat$attr), function(x) dat$attr[[x]][active ==
                                                              1])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL
  status <- l$status
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- dat$control$epi.by
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    assign(ebn, l[[ebn]])
  }

  if (at == 1) {
    dat$epi <- list()
    dat$epi$s.num <- sum(status == "s")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num", ebun[i])]] <- sum(status ==
                                                     "s" & get(ebn) == ebv[i])
      }
    }
    dat$epi$e.num <- sum(status == "e")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("e.num", ebun[i])]] <- sum(status ==
                                                     "e" & get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num <- sum(status == "i")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num", ebun[i])]] <- sum(status ==
                                                     "i" & get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num <- sum(status == "r")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num", ebun[i])]] <- sum(status ==
                                                       "r" & get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num <- length(status)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num", ebun[i])]] <- sum(get(ebn) ==
                                                   ebv[i])
      }
    }
  }
  else {
    dat$epi$s.num[at] <- sum(status == "s")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status ==
                                                         "s" & get(ebn) == ebv[i])
      }
    }
    dat$epi$e.num[at] <- sum(status == "e")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("e.num", ebun[i])]][at] <- sum(status ==
                                                         "e" & get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num[at] <- sum(status == "i")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num", ebun[i])]][at] <- sum(status ==
                                                         "i" & get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num[at] <- sum(status == "r")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num", ebun[i])]][at] <- sum(status ==
                                                           "r" & get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num[at] <- length(status)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num", ebun[i])]][at] <- sum(get(ebn) ==
                                                       ebv[i])
      }
    }
  }

  return(dat)
}


## revised module for recovery process
## this assumes homogenous recovery after average duration
recovery.delayed.random <- function (dat, at) {

  if (!(dat$control$type %in% c("SIR", "SIS"))) {
    return(dat)
  }

  ## control parameter...
  rec.start <- dat$param$rec.start

  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  tea.status <- dat$control$tea.status
  modes <- dat$param$modes
  mode <- idmode(dat$nw)
  type <- dat$control$type
  recovState <- ifelse(type == "SIR", "r", "s")
  rec.rand <- dat$control$rec.rand
  rec.rate <- dat$param$rec.rate
  rec.rate.m2 <- dat$param$rec.rate.m2
  nRecov <- nRecovM2 <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)
  infDur <- at - infTime[active == 1 & status == "i"]
  infDur[infDur == 0] <- 1
  lrec.rate <- length(rec.rate)
  if (lrec.rate == 1) {
    mElig <- mode[idsElig]
    rates <- c(rec.rate, rec.rate.m2)
    ratesElig <- rates[mElig]
  }
  else {
    mElig <- mode[idsElig]
    if (is.null(rec.rate.m2)) {
      rates <- ifelse(infDur <= lrec.rate, rec.rate[infDur],
                      rec.rate[lrec.rate])
    }
    else {
      rates <- ifelse(mElig == 1, ifelse(infDur <= lrec.rate,
                                         rec.rate[infDur], rec.rate[lrec.rate]), ifelse(infDur <=
                                                                                          lrec.rate, rec.rate.m2[infDur], rec.rate.m2[lrec.rate]))
    }
    ratesElig <- rates
  }

  ## we simply add this setup in order to bypass recovery when at is less than rec.start...
  if (at >= rec.start) {

    if (nElig > 0) {
      if (rec.rand == TRUE) {
        vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1) ## randomly recover from the infection status..
        if (length(vecRecov) > 0) {
          idsRecov <- idsElig[vecRecov]
          nRecov <- sum(mode[idsRecov] == 1)
          nRecovM2 <- sum(mode[idsRecov] == 2)
          status[idsRecov] <- recovState
          if (tea.status == TRUE) {
            dat$nw <- activate.vertex.attribute(dat$nw,
                                                prefix = "testatus", value = recovState,
                                                onset = at, terminus = Inf, v = idsRecov)
          }
        }
      }
      else {
        idsRecov <- idsRecovM2 <- NULL
        nRecov <- min(round(sum(ratesElig[mElig == 1])),
                      sum(mElig == 1))
        if (nRecov > 0) {
          idsRecov <- ssample(idsElig[mElig == 1], nRecov)
          status[idsRecov] <- recovState
        }
        if (modes == 2) {
          nRecovM2 <- min(round(sum(ratesElig[mElig ==
                                                2])), sum(mElig == 2))
          if (nRecovM2 > 0) {
            idsRecovM2 <- ssample(idsElig[mElig == 2],
                                  nRecovM2)
            status[idsRecovM2] <- recovState
          }
        }
        totRecov <- nRecov + nRecovM2
        if (tea.status == TRUE & totRecov > 0) {
          allids <- c(idsRecov, idsRecovM2)
          dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                              value = recovState, onset = at, terminus = Inf,
                                              v = allids)
        }
      }
    }
  }

  dat$attr$status <- status
  form <- get_nwparam(dat)$formation
  fterms <- get_formula_terms(form)
  if ("status" %in% fterms) {
    dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  }
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  outName[2] <- paste0(outName, ".m2")
  if (at == 2) {
    dat$epi[[outName[1]]] <- c(0, nRecov)
  }
  else {
    dat$epi[[outName[1]]][at] <- nRecov
  }
  if (modes == 2) {
    if (at == 2) {
      dat$epi[[outName[2]]] <- c(0, nRecovM2)
    }
    else {
      dat$epi[[outName[2]]][at] <- nRecovM2
    }
  }
  return(dat)
}

## threhold-based corrections
recovery.correction <- function(dat, at) {

  # nw <- dat$nw ## recover networks
  ## control parameter...
  rec.start <- dat$param$rec.start

  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  tea.status <- dat$control$tea.status
  modes <- dat$param$modes
  mode <- idmode(dat$nw)
  type <- dat$control$type
  recovState <- ifelse(type == "SIR", "r", "s")
  rec.rand <- dat$control$rec.rand
  rec.rate <- dat$param$rec.rate
  rec.rate.m2 <- dat$param$rec.rate.m2
  nRecov <- nRecovM2 <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)
  infDur <- at - infTime[active == 1 & status == "i"]
  infDur[infDur == 0] <- 1
  lrec.rate <- length(rec.rate)

  if (lrec.rate == 1) {
    mElig <- mode[idsElig]
    rates <- c(rec.rate, rec.rate.m2)
    ratesElig <- rates[mElig]
  } else {
    mElig <- mode[idsElig]
    if (is.null(rec.rate.m2)) {
      rates <- ifelse(infDur <= lrec.rate, rec.rate[infDur],
                      rec.rate[lrec.rate])
    }  else {
      rates <- ifelse(mElig == 1, ifelse(infDur <= lrec.rate,
                                         rec.rate[infDur], rec.rate[lrec.rate]), ifelse(infDur <=
                                                                                          lrec.rate, rec.rate.m2[infDur], rec.rate.m2[lrec.rate]))
    }
    ratesElig <- rates
  }

  idsInf <- which(active == 1 & status == "i")
  idsRec <- which(active == 1 & status == "r")

  correction.prob <- dat$param$correction.prob
  rho <- dat$param$rho

  ## I to R progression
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)

  ## revise rates based on social surroundings of an infected alters!
  if (nEligRec > 0) {

    ## see http://statnet.csde.washington.edu/workshops/SUNBELT/current/ndtv/ndtv_workshop.html#transmission-trees-and-constructed-animations
    ## loop through those who are infected ("idsEligRec")
    ## for every ego "id_EligRec"

    ratesElig.modified <- sapply(1:nEligRec, function(k) {

      ## get persistent ids
      id_EligRec <- idsEligRec[k]
      ## get alters of an ego who is not yet infected (based on active freindship ties of "infected" ego)
      active_alters <- get.neighborhood.active(dat$nw, v = id_EligRec, at = at)

      ## if active alters are not present
      if (length(active_alters) == 0) {
        rate.return <- ratesElig[k] ## copy the rates (nothing happens)

        ## if active alters are present
      } else {
        ## counts the number of alters who is not yet infected
        nCorrect <- length(active_alters[!(active_alters %in% idsInf)])
        ## RANDOM DRAW: n. of alters providing corrections, with probability of sending corrections = correction.prob
        vecCorrect <- which(rbinom(nCorrect, 1, correction.prob) == 1)
        nvecCorrect <- length(vecCorrect)

        ## if proportion of alters sending corrections relative to total neighbor is greater than (rho)
        if (nvecCorrect/length(active_alters) >= rho) {
          rate.return <- 1 - ratesElig[k] ## this person recovers from infection status (yet stochastically)

        } else if (nvecCorrect > 1) { ## if below rho and nvecCorrect > 0
          ## for each additional alters providing corrections (=nvecCorrect), it additionally increases the prob of recovery
          j <- (nvecCorrect-1)
          rate.return <- ratesElig[k]*exp(log(j+1)*1.1435)  ## this increases baseline rates. Numbers taken from Margolin et al. 2017 (mutual friendship coef, average of two values from Table 1/2, final models)
          ## since the recovery rate is the reciprocal of the disease duration,
          ## the increase in rate results to decrease in duration (-> fater recovery overall)
        } else { ## if nvecCorrect == 0
          rate.return <- ratesElig[k] ## copy the rates (nothing happens)
        }

      } ## end of presence of active alters
      ## return vector of revised rates
      return(rate.return)

    }, simplify = TRUE)

    ratesElig <- ratesElig.modified
  }

  ## we simply add this setup in order to bypass recovery when at is less than rec.start...
  if (at >= rec.start) {

    if (nElig > 0) {
      if (rec.rand == TRUE) {
        vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1) ## stochastically recover from the infection status..
        if (length(vecRecov) > 0) {
          idsRecov <- idsElig[vecRecov]
          nRecov <- sum(mode[idsRecov] == 1)
          nRecovM2 <- sum(mode[idsRecov] == 2)
          status[idsRecov] <- recovState
          if (tea.status == TRUE) {
            dat$nw <- activate.vertex.attribute(dat$nw,
                                                prefix = "testatus", value = recovState,
                                                onset = at, terminus = Inf, v = idsRecov)
          }
        } else {
          nRecov <- nRecovM2 <- 0
        }
      }
      else {
        idsRecov <- idsRecovM2 <- NULL
        nRecov <- min(round(sum(ratesElig[mElig == 1])),
                      sum(mElig == 1))
        if (nRecov > 0) {
          idsRecov <- ssample(idsElig[mElig == 1], nRecov)
          status[idsRecov] <- recovState
        }
        if (modes == 2) {
          nRecovM2 <- min(round(sum(ratesElig[mElig ==
                                                2])), sum(mElig == 2))
          if (nRecovM2 > 0) {
            idsRecovM2 <- ssample(idsElig[mElig == 2],
                                  nRecovM2)
            status[idsRecovM2] <- recovState
          }
        }
        totRecov <- nRecov + nRecovM2
        if (tea.status == TRUE & totRecov > 0) {
          allids <- c(idsRecov, idsRecovM2)
          dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                              value = recovState, onset = at, terminus = Inf,
                                              v = allids)
        }
      }
    }
  }

  dat$attr$status <- status
  form <- get_nwparam(dat)$formation
  fterms <- get_formula_terms(form)
  if ("status" %in% fterms) {
    dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  }
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  outName[2] <- paste0(outName, ".m2")
  if (at == 2) {
    dat$epi[[outName[1]]] <- c(0, nRecov)
  }
  else {
    dat$epi[[outName[1]]][at] <- nRecov
  }
  if (modes == 2) {
    if (at == 2) {
      dat$epi[[outName[2]]] <- c(0, nRecovM2)
    }
    else {
      dat$epi[[outName[2]]][at] <- nRecovM2
    }
  }
  return(dat)
}
