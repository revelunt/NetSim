
## simulation with epimodel
## install EpiModel and other required libraries if not already installed
list.of.packages <- c("pbapply","EpiModel", "Rcpp","ergm","btergm","texreg","rstudioapi","data.table", "haven", "networkDynamic", "intergraph", "ape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)

## load library
library("EpiModel")
library("sna")
library('magrittr')
require(data.table)
source("dev/helper-functions.R")

# In this case, we will focus on a "complex contagion", meaning something you
# have to hear about from multiple people before you will adopt it.
# We use a Susceptible-Infectious-Susceptible (SIS) model with no immunity
# (so individuals may transition back and forth between the susceptible and infected states)


# ### Set parameters in advance ###
# nD <- 1258
# nR <- 1142
# n <- nD + nR
# tau <- 2 / 3
# max.time <- 1000 ## time unit is a day
# nsims <- 100 ## no. of replicated simulations

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
formation[[1]] <- ~ edges + nodefactor("pid")

## formation model for tree-like network (absence of triangles)
formation[[2]] <- ~ edges + nodefactor("pid") + concurrent + degree1.5

## formation model for scale-free network
formation[[3]] <- ~ edges + nodefactor("pid") + concurrent + degrange(2:8)

## formation model for more realistic network, conditioned by party identification
## crucial difference of this model is nodematch("pid"),
## which assumes partisan homophily in network formations
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
target.stats[[1]] <- c(792, 753.72)
target.stats[[2]] <- c(792, 753.72, 438, 2025)
target.stats[[3]] <- c(792, 753.72, 285, 280, 72.675, 38.745, 27.605, 17.5, 11.38, 5.3)
target.stats[[4]] <- c(792, 578.16, 799.4, 285, 280, 72.675, 38.745, 27.605, 17.5, 11.38, 5.3)

## dissolution model
## see http://www.sciencedirect.com/science/article/pii/S0747563216309086#bib55
## almost 20% of the respondents reported that they have either hidden other people's comments
## or unfriended others because they do not share political views
coef.diss[[1]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[2]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[3]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[4]] <- dissolution_coefs(~offset(edges) + offset(nodematch("pid")), duration = c(50, 100))


## estimate dynamic network models
est.list <- lapply(1:4, function(i) {
  est <- netest(sw.net, formation[[i]], target.stats[[i]], coef.diss[[i]],
                set.control.ergm = control.ergm(MCMLE.maxit = 250), edapprox = T)
  #est$fit <- logLik(est$fit, add = TRUE)
  est
})

## network dignostics
dx.list <- lapply(est.list, function(i) {
  dx <- netdx(i, nsims = 100, nsteps = max.time, ncores = ncores, verbose = T, keep.tedgelist = T)
  dx
})

## check the dignostic results
pdf("draft/network target statistics.pdf", paper = 'a4r', width = 12, height = 7)
for (i in 1:3) {
  dx <- dx.list[[i]]
  plot(dx, type = "formation", legend = T)
  plot(dx, type = "duration")
  plot(dx, type = "dissolution")
}

dx <- dx.list[[4]]
plot(dx, type = "formation", legend = T)
## for heterogenous dissolution model, automatric plot is not available
## manually recover durations
require(data.table)
dat1 <- dat2 <- data.frame()
for (j in seq_len(nsims)) {
  setDT(diss <- dx$tedgelist[[j]])
  dat1 <- rbind(dat1, diss[pid[diss$tail] != pid[diss$head], .(time = max.time - onset, duration)])
  dat2 <- rbind(dat2, diss[pid[diss$tail] == pid[diss$head], .(time = max.time - onset, duration)])
  rm(diss)
}
require(ggplot2)
## heterogenous ties
ggplot(dat1, aes(x = time, y = duration)) +
  geom_smooth(method = "auto", colour = "grey80", fill = "grey80", se = F) +
  geom_hline(yintercept = 50, linetype = 2) + theme_bw() +
  ## and homogenous ties
  geom_smooth(data = dat2, method = "auto", colour = "grey20", fill = "grey20", alpha = 0.2) +
  geom_hline(yintercept = 100, linetype = 2) + theme_bw() +
  xlab("time") + ylab("Edge Age")
dev.off()
dev.off()

## plot networks and get some basic stats
require(igraph)
require(intergraph)

netstats <- vector(mode = "list", length = 4)
net.title <- list("Bernoulli network", "Tree network", "Realistic network", "Realistic homophilic network")
names(netstats) <- net.title

pdf("draft/network.plots.pdf")
for (i in 1:4) {
  nw <- simulate.ergm(est.list[[i]]$fit, nsim = 1, seed = 345764)
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
save(est.list, dx.list, netstats, file = "results/net_and_dx_list.rda", compress = "bzip2")


## ------------------------- ##
## Step 2: Standard SI model ##
## ------------------------- ##

# load("results/net_and_dx_list.rda")
# rm(dx.list, netstats)

## set number of those who infected at the start
## For initial conditions, one can use the i.num to set the initial number infected at the start,
## or pass in a vector with a disease status for each of the nodes in the network.
## EpiModel stores the individual-level disease status as a vector of lower-case letters:
## “s” for susceptible, “i” for infected, and “r” for recovered.
## we assume random 15% of Republicans are initially infected (n = 184)
status.vector <- c(rep(0, nD), rbinom(nR, 1, 0.15))
status.vector <- ifelse(status.vector == 1, "i", "s")
table(status.vector, pid)
init <- init.net(status.vector = status.vector)

## infection probability (risk of transmission), act rate (mean number of acts)
## THIS NEEDS TO BE JUSTIFIED
## PEW data suggests approximately 16% of all partisans share FN news stories on their social networks,
## OSU political misperception data suggests "every day" or "almost every day" to "several times a week"
## is the modal value for the information sharing on social media
## we set arbitrary number of 2 here...
param1 <- param.net(inf.prob = 0.16, act.rate = 2)

## set control param. (model type and no of stpes, with n of replications)
control1 <- control.net(type = "SI", nsteps = max.time, nsims = nsims, epi.by = "pid",
                       ncores = ncores)

require(pbapply)
pblapply(seq_len(length(est.list)), function(i) {
  dat.name <- paste0("sim1.SI.model", i)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(542435)
  out <- netsim(est.list[[i]], param1, init, control1)
  assign(dat.name, out)
  save(dat.name, file = paste0('results/', dat.name, ".rda"))
})


## load the saved simulations
load("results/net_and_dx_list.rda")
load("results/sim1.SI.model1.rda")
load("results/sim1.SI.model2.rda")
load("results/sim1.SI.model3.rda")
load("results/sim1.SI.model4.rda")

## summary of simulations (using first simulation as an example)
print(sim1.SI.model1)
summary(sim1.SI.model1, at = 500)

## convert to data.frame for further processing
setDT(model1DT.SI <- get.summary.stats(sim1.SI.model1)); model1DT.SI
find.point.to.plot(sim1.SI.model1)

setDT(model2DT.SI <- get.summary.stats(sim1.SI.model2)); model2DT.SI
find.point.to.plot(sim1.SI.model2)

setDT(model3DT.SI <- get.summary.stats(sim1.SI.model3)); model3DT.SI
find.point.to.plot(sim1.SI.model3)

setDT(model4DT.SI <- get.summary.stats(sim1.SI.model4)); model4DT.SI
find.point.to.plot(sim1.SI.model4)

save(model1DT.SI, model2DT.SI, model3DT.SI, model4DT.SI, file = "results/sim1.SI.data.table.rda")
## overall and party-specific prevalence
print.plots.pdf(sim1.SI.model1, "Bernoulli", "Prevalence_SI.model1", F)
print.plots.pdf(sim1.SI.model2, "tree", "Prevalence_SI.model2", F)
print.plots.pdf(sim1.SI.model3, "small-world", "Prevalence_SI.model3", F)
print.plots.pdf(sim1.SI.model4, "homophilous", "Prevalence_SI.model4", F)


## -------------------------- ##
## Step 3: Extended SEI model ##
## -------------------------- ##

rm(list = ls())
source("dev/helper-functions.R")
load("results/net_and_dx_list.rda")
rm(dx.list, netstats)

## extended SE model with network structure as a moderator of epidemic process
## infection model, SIMPLE EXPOSURE FROM THOSE ALREADY INFECTED
## MODERATION BY INDIVIDUAL CHARACTERISTIC (SEE AGING MODULE FOR REF)
## param: pid.diff.rate (NULL or number between 0 and 1)
## this assumes Republicans have a higher infection prob than Democrats
## such that Inf_R.rate - Inf_D.rate = pid.diff.rate

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

## Here, we specify that Republicans (R) has a baseline prevalence of 15% (n = 184, randomly assigned),
## whereas there are no nodes in Democrats (D) infected.
status.vector <- c(rep(0, nD), rbinom(nR, 1, 0.15))
status.vector <- ifelse(status.vector == 1, "i", "s")

## PROGRESSION MODEL
## HOW TO MAKE THIS CONDITIONAL ON NETWORK?
## param: tau (proportion of already infected/susceptable alters as an infection threshold)

progress <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  tea.status <- dat$control$tea.status
  tau <- dat$param$tau

  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")
  idsRec <- which(active == 1 & status == "r")

  ## recover simulated network at time t
  nw <- dat$nw

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
          if (tea.status == TRUE) {
            ## update network dynamic object
            nw <- activate.vertex.attribute(nw,
                                            prefix = "testatus", value = "i",
                                            onset = at, terminus = Inf, v = id_EligInf)
          }
        }
      }
    }
  }

  ## update status for attr dataset
  dat$attr$status <- status

  if (at == 2) {
    dat$epi$ei.flow <- c(0, nInf)
    dat$epi$e.num <- c(0, sum(active == 1 & status == "e"))
    dat$epi$i.num <- c(0, sum(active == 1 & status == "i"))
  }
  else {
    dat$epi$ei.flow[at] <- nInf
    dat$epi$e.num[at] <- sum(active == 1 & status == "e")
    dat$epi$i.num[at] <- sum(active == 1 & status == "i")
  }
  dat$nw <- nw
  return(dat)
}

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

## initial settings: infection prob = 0.16, difference betwwen R vs. D = 0.04 (0.18 vs. 0.14)
## assumes no recovery from infection (believing misperceptions),
## and progress to infection from exposure occurs when more than half of one's neighbors also
## believe the misperception that an ego is exposed to.
param2 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 2, tau = 0.5)
init <- init.net(status.vector = status.vector)

control2 <- control.net(type = "SI", nsteps = max.time, nsims = nsims, epi.by = "pid",
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL,
                       get_prev.FUN = get_prev.exposed.included,
                       skip.check = TRUE,
                       depend = F, verbose.int = 1, save.other = "attr")

require(pbapply)
pblapply(seq_len(length(est.list)), function(i) {
  dat.name <- paste0("test.sim2.SEI.model", i)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(542435)
  out <- netsim(est.list[[i]], param2, init, control2)
  assign(dat.name, out)
  save(dat.name, file = paste0('results/', dat.name, ".rda"))
})

## load the saved simulations
load("results/net_and_dx_list.rda")
load("results/sim2.SEI.model1.rda")
load("results/sim2.SEI.model2.rda")
load("results/sim2.SEI.model3.rda")
load("results/sim2.SEI.model4.rda")

## summary of simulations (using first simulation as an example)
print(sim2.SEI.model1)
summary(sim2.SEI.model1, at = 500)

## convert to data.frame for further processing
setDT(model1DT.SEI <- get.summary.stats(sim2.SEI.model1)); model1DT.SEI
find.point.to.plot(sim2.SEI.model1)
# overall   timeD   timeR
# 466       494     430

setDT(model2DT.SEI <- get.summary.stats(sim2.SEI.model2)); model2DT.SEI
find.point.to.plot(sim2.SEI.model2)
# overall  timeD   timeR
# 458      487     421

setDT(model3DT.SEI <- get.summary.stats(sim2.SEI.model3)); model3DT.SEI
find.point.to.plot(sim2.SEI.model3)
# overall  timeD   timeR
# 506      537     466

setDT(model4DT.SEI <- get.summary.stats(sim2.SEI.model4)); model4DT.SEI
find.point.to.plot(sim2.SEI.model4)
# overall   timeD   timeR
# 399       447     334

save(model1DT.SEI, model2DT.SEI, model3DT.SEI, model4DT.SEI, file = "results/sim2.SEI.data.table.rda")
## overall and party-specific prevalence
print.plots.pdf(sim2.SEI.model1, "Bernoulli", "Prevalence_SEI.model1", include.exposed = T)
print.plots.pdf(sim2.SEI.model2, "tree", "Prevalence_SEI.model2", T)
print.plots.pdf(sim2.SEI.model3, "small-world", "Prevalence_SEI.model3", T)
print.plots.pdf(sim2.SEI.model4, "homophilous", "Prevalence_SEI.model4", T)


## cf. examine simulated networks over time
require(ndtv)
nw10 <- get_network(sim2.SEI.model1, sim = 10)
nw10 %n% "slice.par" <- list(start = 1, end = 1000, interval = 50, aggregate.dur = 1, rule = 'latest')
compute.animation(nw10, animation.mode = 'MDSJ', chain.direction = 'reverse', verbose = FALSE)
render.d3movie(nw10, vertex.cex = 0.9, vertex.col = "pid",
               edge.col = "darkgrey",
               vertex.border = "lightgrey",
               displaylabels = FALSE,
               vertex.tooltip = function(slice){paste('name:',slice%v%'vertex.names','<br>',
                                                      'status:', slice%v%'testatus')})


## ------------------ ##
## Step 4: SEIR model ##
## ------------------ ##

## module for effects of correction (SEIR model)
## this module assumes single instances of "random" corrections
## but introduced later in the time step
## this can be done by supply the parameter setups for "rec.rate" in param.net and set the model type to SIR
## but we need to change the default recovery module to accomodate delayed recovery.
## here we assume an arbitrary number of 168 for "when" corrections are introduced in the network (rec.start).
## The recovery rate implies that the average duration of disease is 336 hours (14 days),
## meaning on average, those who are infected are randomly receive corrections after 14 days,
## and they recover from the infected status
## (The recovery rate is the reciprocal of the disease duration: 1/336 = 0.00297619)
param3 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 2, tau = 0.5,
                    rec.rate = 0.00297619, rec.start = 168)
init <- init.net(status.vector = status.vector)
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

## control settings
control3 <- control.net(type = "SIR", nsteps = max.time, nsims = nsims, epi.by = "pid",
                        ncores = ncores,
                        infection.FUN = infect,
                        progress.FUN = progress,
                        recovery.FUN = recovery.delayed.random,
                        get_prev.FUN = get_prev.exposed.included,
                        skip.check = TRUE,
                        depend = F, verbose.int = 1)



RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEI.model1 <- netsim(est.list[[1]], param3, init, control3)
save(sim3.SEI.model1, file = "results/sim3.SEIR.model1.rda")

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEI.model2 <- netsim(est.list[[2]], param3, init, control3)
save(sim3.SEI.model2, file = "results/sim3.SEIR.model2.rda")

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEI.model3 <- netsim(est.list[[3]], param3, init, control3)
save(sim3.SEI.model3, file = "results/sim3.SEIR.model3.rda")

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEI.model4 <- netsim(est.list[[4]], param3, init, control3)
save(sim3.SEI.model4, file = "results/sim3.SEIR.model4.rda")



## network-based correction assumes that adoption of correcting information (therefore "recovering" from false beliefs)
## is more likely when mmultiple corrections are provided by one's immediate social contacts
## (socially contingent correction of false beliefs)
## set param = "correction.prob", "recover.prob" in param for progress module,
## and "rec.rate" should be set to NULL unless (additional) random recovery is assumed

recovery.correction <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  tea.status <- dat$control$tea.status
  nw <- dat$nw

  if (is.null(dat$param$rec.rate)) {
    rec.rate <- 0
  } else {
    rec.rate <- dat$param$rec.rate
  }

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

    vecRecov <- NULL

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
          ## following logistic-like function
          recover.prob.temp <- recover.prob*((1+exp(-1))/(1+exp(-nvecCorrect)))
          if (recover.prob.temp > 1) {
            recovered <- TRUE ## if greater than 1, change to recovered
          } else {
            recovered <- (rbinom(1, 1, recover.prob.temp) == 1) ## RANDOM DRAW
          }

        } else {
          ## receive no corrections, then recovery is not happening...
          recovered <- FALSE
        }

        ## if recovered, increase the counter for nRec and  change to "recovered"

        if (recovered == TRUE) {
          nRec <- nRec + 1
          status[id_EligRec] <- "r"
          vecRecov <- c(vecRecov, id_EligRec)
        }
      } ## if length(active_alters) == 0, there's no changes in the status
    }

    ## ir.progression, if rec.rate is supplied
    if (rec.rate != 0 ) {
      vecRecov <- unique(c(vecRecov, which(rbinom(nEligRec, 1, rec.rate) == 1)))
      idsRecov <- idsEligRec[vecRecov]
    }

    if (tea.status == TRUE) {
      ## update network dynamic object
      nw <- activate.vertex.attribute(nw,
                                      prefix = "testatus", value = "r",
                                      onset = at, terminus = Inf, v = idsRecov)
    }
  }

  ## update status for attr dataset
  dat$attr$status <- status

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
    del <- discord_edgelist(dat, at)

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
