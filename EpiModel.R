
## simulation with epimodel
## install EpiModel and other required libraries if not already installed
list.of.packages <- c("EpiModel", "Rcpp","ergm","btergm","texreg","rstudioapi","data.table", "haven", "networkDynamic")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)

## load library
library("EpiModel")
library("sna")
library('magrittr')
source("helper-functions.R")

# In this case, we will focus on a "complex contagion", meaning something you
# have to hear about from multiple people before you will adopt it.
# We use a Susceptible-Infectious-Susceptible (SIS) model with no immunity
# (so individuals may transition back and forth between the susceptible and infected states)

# Setup
# 1. Construct a network (Bernoulli Random Graphs vs. small world graphs vs. more realistic social networks) with some key covariates
# (the overall prevalence of those covariates are taken from previous literature)
# 2. Choose one person at random, and set him/her and his/her neighbors as
#    initial adopters
# 3. Infect all the people who have more than tau neighbors infected
# 4. Set the people who were infected at the previous round as "recovered"
# 5. Repeat steps


### Set parameters in advance ###
nD <- 1258
nR <- 1142
n <- nD + nR
tau <- 2 / 3
max.time <- 1000 ## time unit is a day
nsims <- 100 ## no. of replicated simulations

## set no. of cores
ncores <- parallel::detectCores(logical = F)

### Step 1: Construct network ###
# Construct a random network
set.seed(12345)
sw.net <- network.initialize(n = n, directed = FALSE) ## or as.network(rgraph(n, 1, 0.5, mode = "graph"), directed = F)

# set some initial attributes
## partisanship -- one of the core mechanisms of directional motivation
pid <- c(rep("D", nD), rep("R", nR)) ## D == Democrats (50%), R == Republicans (50%)
#pid <- sample(pid, n, replace = F) ## randomly distribute pid values

## accuracy motivation
## accuracy_motive is ranged from 1 to 5
accuracy_motive <- sample(1:5, n, replace = T)

## Political sophistication
## important moderator: high sophistication tend to exhibit more direcitonal reasoning

## we set correlation of 0.5 for this sophistication vector
sophistication <- corr.vector(accuracy_motive, rho = 0.5)
## since there's some negative values, we round them to nearest integers
## and convert them to positive integers
sophistication <- round(sophistication, digits = 0)
sophistication <- sophistication - min(sophistication)
## check correlation (r = 0.49)
cor(accuracy_motive,sophistication)

## scale them to 0 - 1 range
library("scales")
accuracy_motive <- rescale(accuracy_motive, to = c(0,1))
sophistication <- rescale(sophistication, to = c(0,1))

## set network attributes
sw.net <- network::set.vertex.attribute(sw.net, "pid", pid)
sw.net <- network::set.vertex.attribute(sw.net, "accuracy_motive", accuracy_motive)
sw.net <- network::set.vertex.attribute(sw.net, "sophistication", sophistication)

## ---------------------------- ##
## Step 1: fit network dynamics ##
## ---------------------------- ##

formation <- list()
target.stats <- list()
coef.diss <- list()
## formation model for Bernoulli network (assumes	homogeneous	edge probability)
formation[[1]] <- ~ edges

## formation model for tree-like network (absence of triangles)
formation[[2]] <- ~ edges + concurrent + degree1.5

## formation model for scale-free network
formation[[3]] <- ~ edges + concurrent + degrange(2:8)

## formation model for more realistic network, conditioned by party identification
formation[[4]] <- ~ edges + nodematch("pid") + nodefactor("pid") + concurrent + degrange(2:8)

## for nodematch target statistics: from ANES 2008-2009 Panel Study
## see setups.R file for detail.

## target stats for formation statistics, each imply for
## edge (density 0.66 * n/2 = 1500),
## nodefactor (mean degree of 0.3 vs 0.3 in Dem vs. Rep group = 0.3*nR),
## nodematch (assuming 70% within the same pid group,
## the product of the number of edges and the probability of a same-group edge = 1500*0.7 = 1050),
## for concurrent (separately by pid), assume 35% of ties would have more than two edges,
## therefore 0.35*nD + 0.35*nR
## and degrange (no. of nodes with more than 5 edges, there's none, so 0)
## get.target.stats(nD, nR, 0.66, 0.73, 0.66, 0.26)
target.stats[[1]] <- 792
target.stats[[2]] <- c(792, 438, 2025)
target.stats[[3]] <- c(792, 285, 280, 72.675, 38.745, 27.605, 17.5, 11.38, 5.3)
target.stats[[4]] <- c(792, 578.16, 753.72, 285, 280, 72.675, 38.745, 27.605, 17.5, 11.38, 5.3)

## dissolution model
## see http://www.sciencedirect.com/science/article/pii/S0747563216309086#bib55
## almost 20% of the respondents reported that they have either hidden other people's comments
## or unfriended others because they do not share political views
coef.diss[[1]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[2]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[3]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[4]] <- dissolution_coefs(~offset(edges) + offset(nodematch("pid")), duration = c(80, 100))


## estimate dynamic network models
est.list <- lapply(1:4, function(i) {
  est <- netest(sw.net, formation[[i]], target.stats[[i]], coef.diss[[i]],
         set.control.ergm = control.ergm(MCMLE.maxit = 250), edapprox = T)
  est$fit <- logLik(est$fit, add = TRUE)
  est
})

## network dignostics
dx.list <- lapply(est.list, function(i) {
  dx <- netdx(i, nsims = 100, nsteps = max.time, ncores = ncores, verbose = T, keep.tedgelist = T)
  dx
})

## check the dignostic results
pdf("network target statistics.pdf", paper = 'a4r')
for (i in seq_len(length(dx.list))) {
  dx <- dx.list[[i]]
  plot(dx, type = "formation", legend = T)

  if (i == 4) {
    ## for heterogenous dissolution model, automatric plot is not available
    ## manually recover durations
    require(data.table)
    setDT(diss <- dx$tedgelist[[1]])
    dat1 <- diss[pid[diss$tail] != pid[diss$head], .(time = max.time - onset, duration)]
    dat2 <- diss[pid[diss$tail] == pid[diss$head], .(time = max.time - onset, duration)]

    require(ggplot2)
    ## heterogenous ties
    ggplot(dat1, aes(x = time, y = duration)) +
      geom_smooth(method = "auto", colour = "grey80", fill = "grey80", alpha = 0.2) +
      geom_hline(yintercept = 80, linetype = 2) + theme_bw() +
      ## and homogenous ties
      geom_smooth(data = dat2, method = "auto", colour = "grey20", fill = "grey20", alpha = 0.2) +
      geom_hline(yintercept = 100, linetype = 2) + theme_bw() +
      xlab("time") + ylab("Edge Age")
  } else {
    plot(dx, type = "duration")
    plot(dx, type = "dissolution")
  }
}
dev.off()
dev.off()

## plot networks and get some basic stats
require(igraph)
require(intergraph)

netstats <- vector(mode = "list", length = 4)
net.title <- list("Bernoulli network", "Tree network", "Realistic network", "Realistic homophilic network")
names(netstats) <- net.title

pdf("network.plots.pdf")
for (i in 1:4) {
  nw <- simulate.ergm(est.list[[i]]$fit, nsim = 1, seed = 12345)
  cols <- ifelse(network:::get.vertex.attribute(nw, "pid") == "D", "steelblue", "firebrick")
  netstats[[i]] <- list("degreedist" = table(network:::get.vertex.attribute(nw, "pid"),
                                    sna:::degree(nw, gmode = "graph", cmode = "freeman")),
                        "mixingmatrix" = mixingmatrix(asIgraph(nw), "pid"))
  vertex.cex <- rep(0.4, nw$gal$n)
  vertex.cex[sna::isolates(nw)] <- 0.1
  cols[sna::isolates(nw)] <- "gray80"
  nw <- asIgraph(nw)
  V(nw)$size <- vertex.cex*6
  plot(nw, edge.arrow.size = 0, vertex.label = NA, edge.curved = .1,
       layout = layout_with_kk, edge.width = 1.5,
       vertex.color = cols, vertex.frame.color = cols, edge.coloer = "grey30")
  title(net.title[[i]])
}
dev.off()

## save lists
save(est.list, dx.list, netstats, file = "net_and_dx_list.rda", compress = "bzip2")


## standard SI model
## set number of those who infected at the start
## we assume random 15% of Republicans are initially infected (n = 184)
status.vector <- c(rep(0, nD), rbinom(nR, 1, 0.15))
status.vector <- ifelse(status.vector == 1, "i", "s")
table(status.vector, pid)
init <- init.net(status.vector = status.vector)

## infection probability (risk of transmission), act rate (mean number of acts), recovery rates (prob. of recovery among infected)
## THIS NEEDS TO BE JUSTIFIED
## PEW data suggests approximately 16% of all partisans share FN news stories on their social networks,
## OSU political misperception data suggests "every day" or "almost every day" to "several times a week"
## is the modal value for the information sharing on social media
## we set arbitrary number of 2 here...
param <- param.net(inf.prob = 0.16, act.rate = 2)

## set control param. (model type and no of stpes, with n of replications)
control <- control.net(type = "SI", nsteps = max.time, nsims = nsims, epi.by = "pid",
                       ncores = ncores)

sim1 <- netsim(est1, param, init, control)
save(sim1, file = "sim1.SI.model.rda", compress = "bzip2")

## load the saved simulations
load("est1.rda")
load("sim1.SI.model.rda")

## summary of simulations
print(sim1)
summary(sim1, at = 500)
setDT(sim1DT <- as.data.frame(sim1))

## overall and party-specific prevalence
pdf("Prevalence_SI.model.pdf")
plot(sim1, mean.col = c('grey', "black"), qnts.col = c('grey', "black"),
     legend = F, popfrac = T, qnts = 0.95)
abline(v = 103, lty = 2);  text(130, 0.05, "t = 103", col = "black")
legend("right", legend = c("Infected: overall", "Suspected: overall"), lwd = 3,
       col = c("black", "grey"), bg = "white", bty = "n")
#par(mfrow = c(1,2), margin(0, 0, 0, 0)) ## plot side-by-side, with 95% CIs for 100 simulations
plot(sim1, y = c("s.num.pidD", "i.num.pidD"), mean.col = c('grey', "black"), qnts.col = c('grey', "black"),
     legend = F, popfrac = T, qnts = 0.95)
abline(v = 140, lty = 2);  text(170, 0.05, "t = 140", col = "black")
legend("right", legend = c("Infected: Dem", "Suspected: Dem"), lwd = 3,
       col = c("black", "grey"), bg = "white", bty = "n")
plot(sim1, y = c("s.num.pidR", "i.num.pidR"), mean.col = c('grey', "black"), qnts.col = c('grey', "black"),
     legend = F, popfrac = T, qnts = 0.95)
abline(v = 68, lty = 2); text(98, 0.05, "t = 68", col = "black")
legend("right", legend = c("Infected: Rep", "Suspected: Rep"), lwd = 3,
       col = c("black", "grey"), bg = "white", bty = "n")

par(mfrow = c(1,2), mar = c(0,0,0,0))

nw1 <- get_network(sim1, collapse = TRUE, at = 1)
cols <- ifelse(get.vertex.attribute.active(nw1, "testatus", at = 1) == "i", "grey20", "grey50")
vertex.cex <- ifelse(get.vertex.attribute.active(nw1, "testatus", at = 1) == "i", 0.6, 0.4)
vertex.cex[isolates(nw1)] <- 0.2
plot(nw1, mode = "fruchtermanreingold", displayisolates = T,
     vertex.col = cols, vertex.border = cols, edge.col = "grey40", vertex.cex = vertex.cex)
title("Prevalence at t1", line = -2)

nw103 <- get_network(sim1, collapse = TRUE, at = 103)
cols <- ifelse(get.vertex.attribute.active(nw103, "testatus", at = 103) == "i", "grey20", "grey50")
vertex.cex <- ifelse(get.vertex.attribute.active(nw103, "testatus", at = 103) == "i", 0.6, 0.4)
vertex.cex[isolates(nw103)] <- 0.2
plot(nw103, mode = "fruchtermanreingold", displayisolates = T,
     vertex.col = cols, vertex.border = cols, edge.col = "grey40", vertex.cex = vertex.cex)
title("Prevalence at t103", line = -2)

par(mfrow = c(1,1))
comp_plot(sim1, at = 2, digits = 2)
comp_plot(sim1, at = 103, digits = 2)
dev.off()

## extended SIS model with network structure as a moderator of epidemic process
## infection model, SIMPLE EXPOSURE FROM THOSE ALREADY INFECTED
## MODERATION BY INDIVIDUAL CHARACTERISTIC (SEE AGING MODULE FOR REF)
## param: pid.diff.rate (NULL or number between 0 and 1)
## this assumes Republicans have a higher infection prob than Democrats
## such that Inf_R.rate - Inf_D.rate = pid.diff.rate

infect <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
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
    del <- discord_edgelist(dat, idsInf, idsSus, at)
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
        dat$attr$status[idsNewInf] <- "e"
        dat$attr$infTime[idsNewInf] <- at
      }
    }
  }

  if (at == 2) {
    dat$epi$se.flow <- c(0, nInf)
  }
  else {
    dat$epi$se.flow[at] <- nInf
  }
  dat$nw <- nw
  return(dat)
}

## SET INITIAL INFECTION TO ONLY REP NODES?
## For initial conditions, one can use the i.num to set the initial number infected at the start,
## or pass in a vector with a disease status for each of the nodes in the network.
## EpiModel stores the individual-level disease status as a vector of lower-case letters:
## “s” for susceptible, “i” for infected, and “r” for recovered.
## Here, we specify that Republicans (R) has a baseline prevalence of 10% (randomly assigned),
## whereas there are no nodes in Democrats (D) infected.

status.vector <- c(rep(0, n/2), rbinom(n/2, 1, 0.1))
status.vector <- ifelse(status.vector == 1, "i", "s")

## PROGRESSION MODEL
## HOW TO MAKE THIS CONDITIONAL ON NETWORK?
## param: tau (proportion of already infected/susceptable alters as an infection threshold)

progress <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status

  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")

  ## recover simulated network at time t
  nw <- dat$nw

  ir.rate <- dat$param$ir.rate

  ## E to I progression  ## EXPOSED TO INFECTED
  nInf <- 0

  ## get ids of those who are active and eligible for infection ("exposed")
  ## as well as already infectious ("infected")

  idsEligInf <- which(active == 1 & status == "e")
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {

    ## see http://statnet.csde.washington.edu/workshops/SUNBELT/current/ndtv/ndtv_workshop.html#transmission-trees-and-constructed-animations
    ## loop through those who are exposed ("idsEligInf")
    ## for every ego "id_EligInf"
    for (id_EligInf in idsEligInf){
      ## get alters (based on active freindship ties of "exposed" ego
      active_alters <- get.neighborhood.active(nw, v = id_EligInf, at = at)
      ## counts the number of connected alters who are already infected
      ## and compare with the no. of total connected alters
      if (length(active_alters) > 0) {
      t <- length(active_alters[active_alters %in% idsInf]) / length(active_alters)
      ## if this is greater than certain proportion,
      if (t > tau) {
        nInf <- nInf + 1 ## increase the counter for nInf
        status[id_EligInf] <- "i"  ## change to "infected"
        }
      }
    }
  }

  ## I to R progression
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)

  if (nEligRec > 0) {

    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1) ## RANDOM DRAW

    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }

  dat$attr$status <- status

  if (at == 2) {
    dat$epi$ei.flow <- c(0, nInf)
    dat$epi$ir.flow <- c(0, nRec)
    dat$epi$e.num <- c(0, sum(active == 1 & status == "e"))
    dat$epi$r.num <- c(0, sum(active == 1 & status == "r"))
  }
  else {
    dat$epi$ei.flow[at] <- nInf
    dat$epi$ir.flow[at] <- nRec
    dat$epi$e.num[at] <- sum(active == 1 & status == "e")
    dat$epi$r.num[at] <- sum(active == 1 & status == "r")
  }

  return(dat)
}

## initial settings: infection prob = 0.16, difference betwwen R vs. D = 0.04 (0.18 vs. 0.14)
## assumes no recovery from infection (believing misperceptions),
## and progress to infection from exposure occurs when more than half of one's neighbors also
## believe the misperception that an ego is exposed to.
param <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 2, ir.rate = 0, tau = 0.5)
init <- init.net(status.vector = status.vector)

control <- control.net(type = "SI", nsteps = max.time, nsims = nsims, epi.by = "pid",
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL, skip.check = TRUE,
                       depend = F, verbose.int = 1)

sim2 <- netsim(est1, param, init, control)
save(sim2, file = "sim2.SEI.model.rda", compress = "bzip2")
print(sim2)

# load saved objects
load("sim2.SEI.model.rda")
print(sim2)

sim2 <- get.exposed.sim(sim2)
setDT(sim2DT <- as.data.frame(sim2))

pdf("Prevalence_SEI.model.pdf")
plot(sim2, y = c("s.num", "e.num", "i.num"),
     mean.col = c('grey90', "grey50", "black"), qnts.col = c('grey90', "grey50", "black"),
     legend = F, popfrac = T, mean.smooth = T, qnts = 0.95)
legend("topright", legend = c("Infected: overall", "Suspected: overall", "Exposed: overall"), lwd = 3,
       col = c("black", "grey90", "grey50"), bg = "white", bty = "n")

plot(sim2, y = c("s.num.pidD", "e.num.pidD", "i.num.pidD"),
     mean.col = c('grey90', "grey50", "black"), qnts.col = c('grey90', "grey50", "black"),
     legend = F, popfrac = T, mean.smooth = T, qnts = 0.95)
legend("topright", legend = c("Infected: Dem", "Suspected: Dem", "Exposed: Dem"), lwd = 3,
       col = c("black", "grey90", "grey50"), bg = "white", bty = "n")

plot(sim2, y = c("s.num.pidR", "e.num.pidR", "i.num.pidR"),
     mean.col = c('grey90', "grey50", "black"), qnts.col = c('grey90', "grey50", "black"),
     legend = F, popfrac = T, mean.smooth = T, qnts = 0.95)
legend("topright", legend = c("Infected: Rep", "Suspected: Rep", "Exposed: Rep"), lwd = 3,
       col = c("black", "grey90", "grey50"), bg = "white", bty = "n")

par(mfrow = c(1,2), mar = c(0,0,0,0))

nw1 <- get_network(sim2, collapse = TRUE, at = 1)
cols <- ifelse(get.vertex.attribute.active(nw1, "testatus", at = 1) == "i", "grey20", "grey50")
vertex.cex <- ifelse(get.vertex.attribute.active(nw1, "testatus", at = 1) == "i", 0.6, 0.4)
vertex.cex[isolates(nw1)] <- 0.2
plot(nw1, mode = "fruchtermanreingold", displayisolates = T,
     vertex.col = cols, vertex.border = cols, edge.col = "grey40", vertex.cex = vertex.cex)
title("Prevalence at t1", line = -2)

nw103 <- get_network(sim2, collapse = TRUE, at = 103)
cols <- ifelse(get.vertex.attribute.active(nw103, "testatus", at = 103) == "i", "grey20", "grey50")
vertex.cex <- ifelse(get.vertex.attribute.active(nw103, "testatus", at = 103) == "i", 0.6, 0.4)
vertex.cex[isolates(nw103)] <- 0.2
plot(nw103, mode = "fruchtermanreingold", displayisolates = T,
     vertex.col = cols, vertex.border = cols, edge.col = "grey40", vertex.cex = vertex.cex)
title("Prevalence at t103", line = -2)

par(mfrow = c(1,1))
comp_plot.SEI(sim2, at = 2, digits = 2)
comp_plot.SEI(sim2, at = 103, digits = 2)
dev.off()


## for some reason, summry function does not work with custom model
## see http://statnet.github.io/nme/d3-s2.html#data_extraction for manual data extraction
# summary(sim2, at = 5)
#
# summary.sim2 <- as.data.frame(sim2, out = "vals")
# setDT(summary.sim2)
# summary.sim2[, e.num.pidD := num.pidD - s.num.pidD - i.num.pidD]
# summary.sim2[, e.num.pidR := num.pidR - s.num.pidR - i.num.pidR]
#
# ### CHECK THIS CODE ##
# ggplot(summary.sim2, aes(x = time, y = e.num.pidD)) +
#   geom_smooth(method = "auto", colour = "steelblue") +
#   geom_smooth(aes(x = time, y = i.num.pidD), method = "auto", colour = "blue")
#
# ggplot(summary.sim2, aes(x = time, y = e.num.pidR)) +
#   geom_smooth(method = "auto", colour = "firebrick") +
#   geom_smooth(aes(x = time, y = i.num.pidD), method = "auto", colour = "red")


## examine simulated networks over time
require(ndtv)
nw10 <- get_network(sim2, sim = 10)
nw10 %n% "slice.par" <- list(start = 1, end = 120, interval = 10, aggregate.dur = 1, rule = 'latest')
compute.animation(nw10, animation.mode = 'MDSJ', chain.direction = 'reverse', verbose = FALSE)
render.d3movie(nw10, vertex.cex = 0.9, vertex.col = "pid",
               edge.col = "darkgrey",
               vertex.border = "lightgrey",
               displaylabels = FALSE,
               vertex.tooltip = function(slice){paste('name:',slice%v%'vertex.names','<br>',
                                                      'status:', slice%v%'testatus')})


## additional module for effects of correction (SEIR model)
## first module assumes single instances of "random" corrections
## this can be done by supply the parameter setups for "ir.rate" in param.net
param <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 2, ir.rate = 0.08, tau = 0.5)
control <- control.net(type = "SI", nsteps = max.time, nsims = nsims, epi.by = "pid",
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL, skip.check = TRUE,
                       depend = F, verbose.int = 1)

sim3 <- netsim(est1, param, init, control)
save(sim3, file = "sim3.SEIR.model.rda", compress = "bzip2")
print(sim3)






## network-based correction assumes that adoption of correcting information (therefore "recovering" from false beliefs)
## is more likely when mmultiple corrections are provided by one's immediate social contacts
## (socially contingent correction of false beliefs)
## param = correction.prob, recover.prob / in progress module, ir.rate should be set to zero unless random recovery is assumed

recovery.correction <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status

  idsSus <- which(active == 1 & status == "s")
  idsExp <- which(active == 1 & status == "e")
  idsInf <- which(active == 1 & status == "i")
  idsRec <- which(active == 1 & status == "r")

  ## recover simulated network at time t and relevant parameter
  nw <- dat$nw

  correction.prob <- dat$param$correction.prob
  recover.prob <- dat$param$recover.prob

  ## I to R progression
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)

  if (nEligRec > 0) {

    ## see http://statnet.csde.washington.edu/workshops/SUNBELT/current/ndtv/ndtv_workshop.html#transmission-trees-and-constructed-animations
    ## loop through those who are infected ("idsEligRec")
    ## for every ego "id_EligRec"
    for (id_EligRec in idsEligRec){

      ## get alters (based on active freindship ties of "infected" ego
      active_alters <- get.neighborhood.active(nw, v = id_EligRec, at = at)

      ## if active alters are present
      if (length(active_alters) > 0) {

        ## counts the number of alters who is exposed but not infected (idsExp) or already recovered (idsRec)
        nCorrect <- length(active_alters[(active_alters %in% idsExp) | active_alters %in% idsRec])

        ## RANDOM DRAW: n. of alters providing corrections, with probability of sending corrections = correction.prob
        vecCorrect <- which(rbinom(nCorrect, 1, correction.prob) == 1)
        nvecCorrect <- length(vecCorrect)

        if (nvecCorrect == 1) {

          ## if only one correction received, the final prob. of recovery is: recover.prob
          recovered <- (rbinom(1, 1, recover.prob) == 1) ## RANDOM DRAW

        } else if (nvecCorrect > 1) {

          ## for each additional alters providing corrections (=k), it additionally increases the prob of recovery by
          ## recover.prob + 0.06K - 0.018k^2 (THIS NEEDS TO BE JUSTIFIED...)
          recover.prob.temp <- recover.prob + 0.06*nvecCorrect - 0.018*(nvecCorrect^2)
          recovered <- (rbinom(1, 1, recover.prob.temp) == 1) ## RANDOM DRAW

        } else {
          ## receive no corrections, then recovery is not happening...
          recovered <- FALSE
        }

        ## if recovered, increase the counter for nRec and  change to "recovered"
        if (recovered == TRUE) {
          nRec <- nRec + 1
          status[id_EligRec] <- "r"
        }

      } ## if length(active_alters) == 0, there's no changes in the status
    }
  }


  if (at == 2) {
    dat$epi$ir.flow <- c(0, nRec)
    dat$epi$r.num <- c(0, sum(active == 1 & status == "r"))
  }
  else {
    dat$epi$ir.flow[at] <- nRec
    dat$epi$r.num[at] <- sum(active == 1 & status == "r")
  }

  return(dat)

  }









### TO DO LIST : FITTING EPIDEMIC MODELS ON OBSERVED NETWORK
## SEE http://statnet.github.io/gal/empnet.html
## and SEE
## http://statnet.csde.washington.edu/workshops/SUNBELT/current/ndtv/ndtv_workshop.html#importing-data-and-constructing-networkdynamic-objects

source("observednet.R")

# The first module that must be updated is the initialization module.
# The default module contained in the function initialize.net serves as our starting point for writing our own new module, net.init.mod.
# The default module does lots of work, such as simulating an initial network from the model fit,
# that we do not need to do with an observed network.
# Here, there are three key steps for the initialization module:
#   set up the master data list, dat, including its core data structures;
#   initialize infection among the nodes in the network;
#   and use the get_prev.net function to record summary epidemiological statistics.
# With this modules, whereas x would ordinarily be the fitted network model from netest,
# now it will just be the networkDynamic object detailed above.

net.init.mod <- function(x, param, init, control, s) {

  # Master Data List
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$temp <- list()

  # Network Parameters
  dat$nw <- x
  dat$param$modes <- 1

  # Initialization

  ## Infection Status and Time Modules
  n <- network.size(dat$nw)
  dat$attr$status <- init$status.vector

  dat$attr$active <- rep(1, n)
  dat$attr$entrTime <- rep(1, n)
  dat$attr$exitTime <- rep(NA, n)

  dat$attr$infTime <- rep(NA, n)
  dat$attr$infTime[dat$attr$status == "i"] <- 1

  ## Get initial prevalence
  dat <- get_prev.net(dat, at = 1)

  return(dat)
}


# The infection module must also be changed because of some features
# within the default function that depend on having a fitted network model.
# So the default module function, infection.net, serves as the basis of our new module, my.inf.mod.
# It is a stripped down version of the default that provides a much clearer picture of the processes within the module,
# but it is not general enough to handle all the epidemic modeling cases
# supported within EpiModel (e.g., time-vary infection probabilities or simulating epidemics over bipartite networks).
# The key element within both the default and this updated module is the discord_edgelist function
# that examines the current state of the network at 'at',
# and returns a matrix of disease discordant pairs (dyads over which disease may be transmitted).

my.inf.mod <- function(dat, at) {

  ## Variables ##
  active <- dat$attr$active
  status <- dat$attr$status

  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate

  nw <- dat$nw

  # Vector of infected and susceptible IDs
  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  nInf <- totInf <- 0

  ## Processes ##
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {

    # Get discordant edgelist
    del <- discord_edgelist(dat, idsInf, idsSus, at)

    # If some discordant edges, then proceed
    if (!(is.null(del))) {

      # Infection probabilities
      del$transProb <- inf.prob

      # Act rates
      del$actRate <- act.rate

      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - del$transProb) ^ del$actRate

      # Randomize transmissions and subset df
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Set new infections vector
      idsNewInf <- unique(del$sus)
      totInf <- length(idsNewInf)

      # Update attributes
      if (totInf > 0) {
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at
      }

    }
  }

  ## Summary statistics ##
  if (at == 2) {
    dat$epi$si.flow <- c(0, totInf)
  } else {
    dat$epi$si.flow[at] <- totInf
  }

  dat$nw <- nw
  return(dat)
}


# There are many other modules that could potentially run, and one in particular is the resim_nets.FUN module
# that handles the temporal network resimulation; we do not want network resimulation in this case.
# The  module.order argument, then, is a listing of the modules within the time loop that should be run at each time step.
# Since our model includes no demography (births or deaths), disease recovery, or network simulation,
# the only modules we need to run within the time loop are the infection module and the prevalence tracker module.
# Finally, the skip.check argument is set to false to bypass some standard error checking for built-in modules,
# and saving the network and network statistics are turned off
# (these are already contained in the networkDynamic object anyway).

control <- control.net(type = "SI", nsteps = 27, nsims = 100, ncores = 8,
                       initialize.FUN = net.init.mod, infection.FUN = my.inf.mod,
                       module.order = c("infection.FUN", "get_prev.FUN"),
                       skip.check = TRUE, save.nwstats = F, save.network = T)


# parameterization
param <- param.net(inf.prob = 0.9)

status.vector <- init.status.vector <- get.vertex.attribute.active(g, "partyid", at = 1)
status.vector[status.vector == 1] <- "s"
status.vector[status.vector == 0] <- sample(c("i", "s"), replace = T,
                                            size = length(status.vector[status.vector == 0]),
                                            prob = c(0.15, 0.85))
init <- init.net(status.vector = status.vector) ## total 15 infected, among Republican (coded as zero)

## simulation
sim3 <- netsim(g, param, init, control)
