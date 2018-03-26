
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
require(stargazer)

# In this case, we will focus on a "complex contagion", meaning something you
# have to hear about from multiple people before you will adopt it.
# We use a Susceptible-Infectious-Susceptible (SIS) model with no immunity
# (so individuals may transition back and forth between the susceptible and infected states)


# ### Set parameters in advance ###
# nD <- 1258
# nR <- 1142
# n <- nD + nR
# max.time <- 1000 ## time unit is a day
# nsims <- 100 ## no. of replicated simulations

## set no. of cores
ncores <- parallel::detectCores(logical = F)

### Step 1: Construct network ###
# Construct a random network
require(igraph)
set.seed(12345)
sw.net <- list()
sw.net[[1]] <- igraph::sample_gnm(n = n, m = 792) %>%  intergraph::asNetwork(.) ## n = node, m = edge
sw.net[[2]] <- igraph::sample_smallworld(1, n, p = 0, nei = 1, loops = FALSE, multiple = FALSE) %>%  intergraph::asNetwork(.) ## lattice network, since rewiring prob is zero
sw.net[[3]] <- igraph::sample_smallworld(1, n, p = .3, nei = 2, loops = FALSE, multiple = FALSE) %>%  intergraph::asNetwork(.) ## small-world
sw.net[[4]] <- igraph::sample_smallworld(1, n, p = .3, nei = 2, loops = FALSE, multiple = FALSE) %>%  intergraph::asNetwork(.) ## small-world

# set some initial attributes
## partisanship -- one of the core mechanisms of directional motivation
source_lines("dev/setups.R", 1:100)
pid <- anespanel0809[pid %in% c(1,2), ifelse(pid == 1, "D", "R")]
pidst <- as.vector(anespanel0809[pid %in% c(1,2), pidst])

## accuracy motivation: NFC
## NFC is ranged from -2 to 2
## with mean and sd by pid is:
## Dem: mean = 0.3526269, sd = 0.9299722
## Rep: mean = 0.3240918, sd = 0.8939632
nfc <- anespanel0809[pid %in% c(1,2), nfc]
nfc[is.na(nfc)] <- 0
## normalize the nfc values
nfc <- as.vector(range01(nfc))

attr <- data.frame(pid = pid, pidst = pidst, nfc = nfc)
attr <- attr[order(attr$pid), ]

## set network attributes
for (i in 1:4) {
  sw.net[[i]] <- network::set.vertex.attribute(sw.net[[i]], "pid", as.character(attr$pid))
  sw.net[[i]] <- network::set.vertex.attribute(sw.net[[i]], "pidst", as.vector(attr$pidst))
  sw.net[[i]] <- network::set.vertex.attribute(sw.net[[i]], "nfc", as.vector(attr$nfc))
}

## ---------------------------- ##
## Step 1: fit network dynamics ##
## ---------------------------- ##

formation <- list()
target.stats <- list()
coef.diss <- list()

## formation model for Bernoulli network (assumes	homogeneous	edge probability)
formation[[1]] <- ~ edges + nodefactor("pid")

## formation model for tree-like network (absence of triangles)
formation[[2]] <- ~ edges + nodefactor("pid") + concurrentties + degrange(3:12)

## formation model for scale-free network
formation[[3]] <- ~ edges + nodefactor("pid") + degrange(0) + concurrentties + triangle

## formation model for more realistic network, conditioned by party identification
## crucial difference of this model is nodematch("pid"),
## which assumes partisan homophily in network formations
formation[[4]] <- ~ edges + nodematch("pid") + nodefactor("pid") + degrange(0) + concurrentties + triangle

## for nodematch target statistics: from ANES 2008-2009 Panel Study
## see setups.R file for detail.

## target stats for formation statistics, each imply for
## edge (density 0.66 * n/2 = 1500),
## nodefactor (mean degree of 0.66 vs 0.66 in Dem vs. Rep group = 0.3*nR),
## nodematch (assuming 70% within the same pid group,
## the product of the number of edges and the probability of a same-group edge = 1500*0.7 = 1050),
## for concurrent (separately by pid), assume 35% of ties would have more than two edges,
## therefore 0.35*nD + 0.35*nR
## and degrange (no. of nodes with more than 5 edges, there's none, so 0)
## get.target.stats(nD, nR, 0.66, 0.73, 0.66, 0.26)
target.stats[[1]] <- c(792, 753.72)
target.stats[[2]] <- c(792, 753.72, 438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
target.stats[[3]] <- c(792, 753.72, 205.92, 500, 4)
target.stats[[4]] <- c(792, 578.16, 753.72, 205.92, 500, 3.5)

## dissolution model
coef.diss[[1]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[2]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[3]] <- dissolution_coefs(~offset(edges), duration = 100)
coef.diss[[4]] <- dissolution_coefs(~offset(edges) + offset(nodematch("pid")), duration = c(50, 100))


## estimate dynamic network models
est.list <- lapply(1:4, function(i) {
  est <- netest(sw.net[[i]], formation[[i]], target.stats[[i]], coef.diss[[i]],
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
require(intergraph)

netstats <- vector(mode = "list", length = 4)
net.title <- list("E-R", "Chain", "SW No-Homophily", "SW Homophily")
names(netstats) <- net.title

table1 <- matrix(NA, ncol = 4, nrow = 7)

layout <- list(layout_with_kk, layout_with_graphopt, layout_with_kk, layout_with_kk)


for (i in 1:4) {
  nw <- simulate.ergm(est.list[[i]]$fit, nsim = 1, seed = 12345)

  table1[1,i] <- format(sna::gden(nw), digits = 3, nsmall = 3)
  table1[2,i] <- format(summary(sna::degree(nw, gmode = "directed", cmode = "freeman"))[4], digits = 3, nsmall = 3)
  table1[3,i] <- format(summary(sna::degree(nw, gmode = "directed", cmode = "freeman"))[6], digits = 1, nsmall = 1)
  table1[4,i] <- format(igraph::transitivity(asIgraph(nw)), digits = 1, nsmall = 1)
  table1[5,i] <- format(igraph::mean_distance(asIgraph(nw), directed = T, unconnected = T), digits = 1, nsmall = 1)
  table1[6,i] <- sprintf("%.1f%%", 100*(sum(diag(mixingmatrix(asIgraph(nw), "pid")))/sum(mixingmatrix(asIgraph(nw), "pid"), na.rm = T)))
  table1[7,i] <- sprintf("%.1f%%", 100*(mixingmatrix(asIgraph(nw), "pid")[1,2]/sum(mixingmatrix(asIgraph(nw), "pid"), na.rm = T)))

  cols <- ifelse(network:::get.vertex.attribute(nw, "pid") == "D", "steelblue", "firebrick")
  netstats[[i]] <- list("degreedist" = table(network:::get.vertex.attribute(nw, "pid"),
                                             sna:::degree(nw, gmode = "graph", cmode = "freeman")),
                        "mixingmatrix" = mixingmatrix(asIgraph(nw), "pid"))
  vertex.cex <- rep(0.4, nw$gal$n)
  vertex.cex[sna::isolates(nw)] <- 0.1
  cols[sna::isolates(nw)] <- "white"
  nw <- asIgraph(nw)
  V(nw)$size <- vertex.cex*6

  #nw.com <- delete.vertices(nw, degree(nw)==0)
  #cols <- cols[cols != "gray80"]
  pdf(paste("draft/network_plots",i,".pdf",sep=""))
  plot(nw, edge.arrow.size = 0, vertex.label = NA, edge.curved = .1,
       layout = layout[[i]],
       edge.width = 1.5,
       vertex.color = cols, vertex.frame.color = cols, edge.coloer = "grey30")
  dev.off()
}

colnames(table1) <- net.title
rownames(table1) <- c("Graph density", "Mean degree", "Max degree", "Clustering coefficients",
                      "Mean distance", "Homophilous ties", "Heterophilous ties")
require(stargazer)
stargazer(table1, title = "Descriptive statistics of cross-sectional network from simulated networks")
#notes = "E-R: Erdos-Renyi random network. Chain: Chain network. SW No-Homophily: Small-world network absent of homophily. SW Homophily: Small-world network with partisan homophily"


## save lists
save(est.list, dx.list, netstats, file = "results/net_and_dx_list.rda", compress = "bzip2")


# ## ------------------------- ##
# ## Step 2: Standard SI model ##
# ## ------------------------- ##
#
# # load("results/net_and_dx_list.rda")
# # rm(dx.list, netstats)
#
# ## set number of those who infected at the start
# ## For initial conditions, one can use the i.num to set the initial number infected at the start,
# ## or pass in a vector with a disease status for each of the nodes in the network.
# ## EpiModel stores the individual-level disease status as a vector of lower-case letters:
# ## “s” for susceptible, “i” for infected, and “r” for recovered.
# ## we assume random 15% of Republicans are initially infected (n = 184)
# status.vector <- c(rep(0, nD), rbinom(nR, 1, 0.15))
# status.vector <- ifelse(status.vector == 1, "i", "s")
# table(status.vector, pid)
# init <- init.net(status.vector = status.vector)
#
# ## infection probability (risk of transmission), act rate (mean number of acts)
# ## PEW data suggests approximately 16% of all partisans share FN news stories on their social networks,
# ## while Guess et al. and Allcott & Gentzkow report 0.060 / 0.268 shares per unit time
# ## (0.060 + 0.268) / 2 = 0.164
# param1 <- param.net(inf.prob = 0.16, act.rate = 0.164)
#
# ## set control param. (model type and no of stpes, with n of replications)
# control1 <- control.net(type = "SI", nsteps = max.time, nsims = nsims, epi.by = "pid",
#                        ncores = ncores)
#
# require(pbapply)
# pblapply(seq_len(length(est.list)), function(i) {
#   dat.name <- paste0("sim1.SI.model", i)
#   RNGkind("L'Ecuyer-CMRG")
#   set.seed(542435)
#   out <- netsim(est.list[[i]], param1, init, control1)
#   assign(dat.name, out)
#   save(dat.name, file = paste0('results/', dat.name, ".rda"))
# })
#
#
# ## load the saved simulations
# load("results/net_and_dx_list.rda")
# load("results/sim1.SI.model1.rda")
# load("results/sim1.SI.model2.rda")
# load("results/sim1.SI.model3.rda")
# load("results/sim1.SI.model4.rda")
#
# ## summary of simulations (using first simulation as an example)
# print(sim1.SI.model1)
# summary(sim1.SI.model1, at = 500)
#
# ## convert to data.frame for further processing
# setDT(model1DT.SI <- get.summary.stats(sim1.SI.model1)); model1DT.SI
# find.point.to.plot(sim1.SI.model1)
#
# setDT(model2DT.SI <- get.summary.stats(sim1.SI.model2)); model2DT.SI
# find.point.to.plot(sim1.SI.model2)
#
# setDT(model3DT.SI <- get.summary.stats(sim1.SI.model3)); model3DT.SI
# find.point.to.plot(sim1.SI.model3)
#
# setDT(model4DT.SI <- get.summary.stats(sim1.SI.model4)); model4DT.SI
# find.point.to.plot(sim1.SI.model4)
#
# save(model1DT.SI, model2DT.SI, model3DT.SI, model4DT.SI, file = "results/sim1.SI.data.table.rda")
# ## overall and party-specific prevalence
# print.plots.pdf(sim1.SI.model1, "Bernoulli", "Prevalence_SI.model1", F)
# print.plots.pdf(sim1.SI.model2, "tree", "Prevalence_SI.model2", F)
# print.plots.pdf(sim1.SI.model3, "small-world", "Prevalence_SI.model3", F)
# print.plots.pdf(sim1.SI.model4, "homophilous", "Prevalence_SI.model4", F)


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
## such that Inf_R.rate - Inf_D.rate = pid.diff.rate (partisan differential infection probability)

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
## This is a constant hazard function (individual-level stochastic process)
## Republicans have a higher transition rates than Democrats
## and those who have higher accuracy motivations have a lower transition rates

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
param2 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                    r.ei.rate = 0.3, d.ei.rate = 0.07)
init <- init.net(status.vector = status.vector)

control2 <- control.net(type = "SI", nsteps = max.time, nsims = nsims, epi.by = "pid",
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = exposed.to.infectious.progress,
                       recovery.FUN = NULL, ## assumes no recovery process
                       get_prev.FUN = get_prev.exposed.included,
                       skip.check = TRUE,
                       depend = F, verbose.int = 1, save.other = "attr")

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim2.SEI.model1 <- netsim(est.list[[1]], param2, init, control2)
save(sim2.SEI.model1, file = "results/sim2.SEI.model1.rda")

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim2.SEI.model2 <- netsim(est.list[[2]], param2, init, control2)
save(sim2.SEI.model2, file = "results/sim2.SEI.model2.rda")

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim2.SEI.model3 <- netsim(est.list[[3]], param2, init, control2)
save(sim2.SEI.model3, file = "results/sim2.SEI.model3.rda")

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim2.SEI.model4 <- netsim(est.list[[4]], param2, init, control2)
save(sim2.SEI.model4, file = "results/sim2.SEI.model4.rda")



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
print.plots.pdf(sim2.SEI.model2, "tree", "Prevalence_SEI.model2", include.exposed = T)
print.plots.pdf(sim2.SEI.model3, "small-world", "Prevalence_SEI.model3", include.exposed = T)
print.plots.pdf(sim2.SEI.model4, "homophilous", "Prevalence_SEI.model4", include.exposed = T)

## get mean prevalence over entire time period
## those who are being exposed are: exposed + (infected|exposed)
dat.plot <- data.table(time = 1:1000,
                       model1 = model1DT.SEI[, (e.num.mean + i.num.mean)/num.mean],
                       model2 = model2DT.SEI[, (e.num.mean + i.num.mean)/num.mean],
                       model3 = model3DT.SEI[, (e.num.mean + i.num.mean)/num.mean],
                       model4 = model4DT.SEI[, (e.num.mean + i.num.mean)/num.mean],

                       model1.llci = model1DT.SEI[, (e.num.llci + i.num.llci)/num.mean],
                       model2.llci = model2DT.SEI[, (e.num.llci + i.num.llci)/num.mean],
                       model3.llci = model3DT.SEI[, (e.num.llci + i.num.llci)/num.mean],
                       model4.llci = model4DT.SEI[, (e.num.llci + i.num.llci)/num.mean],

                       model1.ulci = model1DT.SEI[, (e.num.ulci + i.num.ulci)/num.mean],
                       model2.ulci = model2DT.SEI[, (e.num.ulci + i.num.ulci)/num.mean],
                       model3.ulci = model3DT.SEI[, (e.num.ulci + i.num.ulci)/num.mean],
                       model4.ulci = model4DT.SEI[, (e.num.ulci + i.num.ulci)/num.mean],

                       model1.suspect = model1DT.SEI[, s.num.mean/num.mean],
                       model2.suspect = model2DT.SEI[, s.num.mean/num.mean],
                       model3.suspect = model3DT.SEI[, s.num.mean/num.mean],
                       model4.suspect = model4DT.SEI[, s.num.mean/num.mean])

dat.plot2 <- data.table(time = 1:1000,
                        model1.D = model1DT.SEI[, (e.num.pidD.mean + i.num.pidD.mean)/num.mean],
                        model1.R = model1DT.SEI[, (e.num.pidR.mean + i.num.pidR.mean)/num.mean],
                        model2.D = model2DT.SEI[, (e.num.pidD.mean + i.num.pidD.mean)/num.mean],
                        model2.R = model2DT.SEI[, (e.num.pidR.mean + i.num.pidR.mean)/num.mean],
                        model3.D = model3DT.SEI[, (e.num.pidD.mean + i.num.pidD.mean)/num.mean],
                        model3.R = model3DT.SEI[, (e.num.pidR.mean + i.num.pidR.mean)/num.mean],
                        model4.D = model4DT.SEI[, (e.num.pidD.mean + i.num.pidD.mean)/num.mean],
                        model4.R = model4DT.SEI[, (e.num.pidR.mean + i.num.pidR.mean)/num.mean],

                        model1.D.llci = model1DT.SEI[, (e.num.pidD.llci + i.num.pidD.llci)/num.mean],
                        model1.R.llci = model1DT.SEI[, (e.num.pidR.llci + i.num.pidR.llci)/num.mean],
                        model2.D.llci = model2DT.SEI[, (e.num.pidD.llci + i.num.pidD.llci)/num.mean],
                        model2.R.llci = model2DT.SEI[, (e.num.pidR.llci + i.num.pidR.llci)/num.mean],
                        model3.D.llci = model3DT.SEI[, (e.num.pidD.llci + i.num.pidD.llci)/num.mean],
                        model3.R.llci = model3DT.SEI[, (e.num.pidR.llci + i.num.pidR.llci)/num.mean],
                        model4.D.llci = model4DT.SEI[, (e.num.pidD.llci + i.num.pidD.llci)/num.mean],
                        model4.R.llci = model4DT.SEI[, (e.num.pidR.llci + i.num.pidR.llci)/num.mean],

                        model1.D.ulci = model1DT.SEI[, (e.num.pidD.ulci + i.num.pidD.ulci)/num.mean],
                        model1.R.ulci = model1DT.SEI[, (e.num.pidR.ulci + i.num.pidR.ulci)/num.mean],
                        model2.D.ulci = model2DT.SEI[, (e.num.pidD.ulci + i.num.pidD.ulci)/num.mean],
                        model2.R.ulci = model2DT.SEI[, (e.num.pidR.ulci + i.num.pidR.ulci)/num.mean],
                        model3.D.ulci = model3DT.SEI[, (e.num.pidD.ulci + i.num.pidD.ulci)/num.mean],
                        model3.R.ulci = model3DT.SEI[, (e.num.pidR.ulci + i.num.pidR.ulci)/num.mean],
                        model4.D.ulci = model4DT.SEI[, (e.num.pidD.ulci + i.num.pidD.ulci)/num.mean],
                        model4.R.ulci = model4DT.SEI[, (e.num.pidR.ulci + i.num.pidR.ulci)/num.mean])

# p2_1 <- ggplot(dat.plot2, aes(time, model1.D)) + geom_line(color = "steelblue") + ## model1 no. of exposed
#   geom_ribbon(aes(ymin = model1.D.llci, ymax = model1.D.ulci), fill = "steelblue", alpha = 0.1) +
#   geom_line(aes(time, model1.R), linetype = "dotted", color = "firebrick") +
#   geom_ribbon(aes(ymin = model1.R.llci, ymax = model1.R.ulci), fill = "firebrick", alpha = 0.1) +
#   theme_bw() + ylab("prevalence over time (proportions)") + xlab("time") +
#   ggtitle("E-R random network")
#
# p2_2 <- ggplot(dat.plot2, aes(time, model2.D)) + geom_line(color = "steelblue") + ## model1 no. of exposed
#   geom_ribbon(aes(ymin = model2.D.llci, ymax = model2.D.ulci), fill = "steelblue", alpha = 0.1) +
#   geom_line(aes(time, model2.R), linetype = "dotted", color = "firebrick") +
#   geom_ribbon(aes(ymin = model2.R.llci, ymax = model2.R.ulci), fill = "firebrick", alpha = 0.1) +
#   theme_bw() + ylab("prevalence over time (proportions)") + xlab("time") +
#   ggtitle("Chain network")
#
# p2_3 <- ggplot(dat.plot2, aes(time, model3.D)) + geom_line(color = "steelblue") + ## model1 no. of exposed
#   geom_ribbon(aes(ymin = model3.D.llci, ymax = model3.D.ulci), fill = "steelblue", alpha = 0.1) +
#   geom_line(aes(time, model3.R), linetype = "dotted", color = "firebrick") +
#   geom_ribbon(aes(ymin = model3.R.llci, ymax = model3.R.ulci), fill = "firebrick", alpha = 0.1) +
#   theme_bw() + ylab("prevalence over time (proportions)") + xlab("time") +
#   ggtitle("SW absent homophily")
#
# p2_4 <- ggplot(dat.plot2, aes(time, model4.D)) + geom_line(color = "steelblue") + ## model1 no. of exposed
#   geom_ribbon(aes(ymin = model4.D.llci, ymax = model4.D.ulci), fill = "steelblue", alpha = 0.1) +
#   geom_line(aes(time, model4.R), linetype = "dotted", color = "firebrick") +
#   geom_ribbon(aes(ymin = model4.R.llci, ymax = model4.R.ulci), fill = "firebrick", alpha = 0.1) +
#   theme_bw() + ylab("prevalence over time (proportions)") + xlab("time") +
#   ggtitle("SW homophily")

pdf("draft/Figure 2. Prevalence.pdf", paper = 'a4r', width = 12, height = 7)
ggplot(dat.plot, aes(time, model1)) + geom_line() + ## model1 no. of exposed
  #geom_ribbon(aes(ymin = model1.llci, ymax = model1.ulci), alpha = 0.1) +
  geom_line(aes(time, model1.suspect), linetype = "dotted") +
  geom_line(aes(time, model2), color = "grey") +
  #geom_ribbon(aes(x = time, ymin = model2.llci, ymax = model2.ulci), fill = "grey", alpha = 0.1) +
  geom_line(aes(time, model2.suspect), color = "grey", linetype = "dotted") +
  geom_line(aes(time, model3), color = "steelblue") +
  #geom_ribbon(aes(ymin = model3.llci, ymax = model3.ulci), fill = "blue", alpha = 0.1) +
  geom_line(aes(time, model3.suspect), color = "steelblue", linetype = "dotted") +
  geom_line(aes(time, model4), color = "firebrick") +
  #geom_ribbon(aes(ymin = model4.llci, ymax = model4.ulci), fill = "red", alpha = 0.1) +
  geom_line(aes(time, model4.suspect), color = "firebrick", linetype = "dotted") +
  theme_bw() + ylab("Cumulative Fraction Exposed to Misinformation") + xlab("Time (Hours)")

## adoption rates:
ggplot(model1DT.SEI, aes(time, se.flow.mean)) +
  geom_smooth(fill = "black", colour = "black", method = "loess", alpha = 0.1) +
  geom_smooth(data = model2DT.SEI, aes(time, se.flow.mean),
              fill = "grey", colour = "grey", method = "loess", alpha = 0.1) +
  geom_smooth(data = model3DT.SEI, aes(time, se.flow.mean),
              fill = "steelblue", colour = "steelblue", method = "loess", alpha = 0.1) +
  geom_smooth(data = model4DT.SEI, aes(time, se.flow.mean),
              fill = "firebrick", colour = "firebrick", method = "loess", alpha = 0.1) +
  theme_bw() + ylab("Suspected to Exposed Flow") + xlab("Time (Hours)")
dev.off()



## prevalence table
tb3 <- matrix(NA, nrow = 6, ncol = 4)
rownames(tb3) <- c("overall", "95% CIs", "Among Dem", "95% CIs", "Among Rep", "95% CIs")
colnames(tb3) <- c("E-R", "Chain", "SW No-Homophily", "SW Homophily")

datlist <- list(model1DT.SEI, model2DT.SEI, model3DT.SEI, model4DT.SEI)
for (i in 1:4) {
dat <- datlist[[i]]
tb3[,i] <- c(format(dat[, mean((e.num.mean + i.num.mean)/num.mean)],
                    trim = T, digits = 3, nsmall = 3, justify = "centre"),
             paste0("[", format(dat[, mean((e.num.llci + i.num.llci)/num.mean)],
                                trim = T, digits = 3, nsmall = 3),
                 ", ", format(dat[, mean((e.num.ulci + i.num.ulci)/num.mean)],
                              trim = T, digits = 3, nsmall = 3), "]"),
             format(dat[, mean((e.num.pidD.mean + i.num.pidD.mean)/num.mean, na.rm = T)],
                    trim = T, digits = 3, nsmall = 3, justify = "centre"),
             paste0("[", format(dat[, mean((e.num.pidD.llci + i.num.pidD.llci)/num.mean, na.rm = T)],
                                trim = T, digits = 3, nsmall = 3),
                    ", ", format(dat[, mean((e.num.pidD.ulci + i.num.pidD.ulci)/num.mean, na.rm = T)],
                                 trim = T, digits = 3, nsmall = 3), "]"),
             format(dat[, mean((e.num.pidR.mean + i.num.pidR.mean)/num.mean, na.rm = T)],
                    trim = T, digits = 3, nsmall = 3, justify = "centre"),
             paste0("[", format(dat[, mean((e.num.pidR.llci + i.num.pidR.llci)/num.mean, na.rm = T)],
                                trim = T, digits = 3, nsmall = 3),
                    ", ", format(dat[, mean((e.num.pidR.ulci + i.num.pidR.ulci)/num.mean, na.rm = T)],
                                 trim = T, digits = 3, nsmall = 3), "]"))
}

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

tb4 <- matrix(NA, nrow = 3, ncol = 4)
rownames(tb4) <- c("Chain", "SW No-Homophily", "SW Homophily")
colnames(tb4) <- c("E-R", "Chain", "SW No-Homophily", "SW Homophily")
tb4[,1] <- c(paste(WRS.test(model1DT.SEI, model2DT.SEI), collapse = ":"),
             paste(WRS.test(model1DT.SEI, model3DT.SEI), collapse = ":"),
             paste(WRS.test(model1DT.SEI, model4DT.SEI), collapse = ":"))
tb4[,2] <- c("-",
             paste(WRS.test(model2DT.SEI, model3DT.SEI), collapse = ":"),
             paste(WRS.test(model2DT.SEI, model4DT.SEI), collapse = ":"))
tb4[,3] <- c(NA,
             "-",
             paste(WRS.test(model3DT.SEI, model4DT.SEI), collapse = ":"))
tb4[,4] <- c(NA, NA, "-")
stargazer(tb4)


## cf. examine simulated networks over time
require(ndtv)
nw10 <- get_network(sim2.SEI.model1, sim = 10)
nw10 %n% "slice.par" <- list(start = 1, end = 1000, interval = 100, aggregate.dur = 1, rule = 'latest')
compute.animation(nw10, animation.mode = 'kamadakawai', chain.direction = 'reverse', verbose = FALSE)
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
## meaning on average, those who are infected are randomly receive corrections after 14 days (two weeks),
## and they recover from the infected status
## (The recovery rate is the reciprocal of the disease duration: 1/336 = 0.00297619)
## this is multiplied by hypothetical correction.prob (0.5)
## since recovery rate of 0.00297619 here means correction.prob is always 1, so we need to discount this number by
## the probability of recieving correction (correction.prob) that is not explicitly accounted for.
param3 <-  param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                     r.ei.rate = 0.3, d.ei.rate = 0.07,
                     rec.rate = 0.001488095, rec.start = 168)
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
                        progress.FUN = exposed.to.infectious.progress,
                        recovery.FUN = recovery.delayed.random,
                        get_prev.FUN = get_prev.exposed.included,
                        skip.check = TRUE,
                        depend = F, verbose.int = 1, save.other = "attr")

## See Table 2 of the main manuscript for details regarding numbering of the scenarios.
## This batch represents asocial correction (random recovery) in Table 2
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model1 <- netsim(est.list[[1]], param3, init, control3) ## model 1 in Table 2
save(sim3.SEIR.model1, file = "results/sim3.SEIR.model1.rda")
rm(sim3.SEIR.model1)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model5 <- netsim(est.list[[2]], param3, init, control3) ## model 5 in Table 2
save(sim3.SEIR.model5, file = "results/sim3.SEIR.model5.rda")
rm(sim3.SEIR.model5)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model9 <- netsim(est.list[[3]], param3, init, control3) ## model 9 in Table 2
save(sim3.SEIR.model9, file = "results/sim3.SEIR.model9.rda")
rm(sim3.SEIR.model9)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model13 <- netsim(est.list[[4]], param3, init, control3) ## model 13 in Table 2
save(sim3.SEIR.model13, file = "results/sim3.SEIR.model13.rda")
rm(sim3.SEIR.model13)


## network-based correction assumes that adoption of correcting information (therefore "recovering" from false beliefs)
## is more likely when mmultiple corrections are provided by one's immediate social contacts
## (socially contingent correction of false beliefs)
## set param = "correction.prob", "rec.rate" in param for progress module,
## and "rec.rate" should be set to baseline, which then modified by the number of alters who sends correction

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

## control settings
init <- init.net(status.vector = status.vector)
control4 <- control.net(type = "SIR", nsteps = max.time, nsims = nsims, epi.by = "pid",
                        ncores = ncores,
                        infection.FUN = infect,
                        progress.FUN = exposed.to.infectious.progress,
                        recovery.FUN = recovery.correction,
                        get_prev.FUN = get_prev.exposed.included,
                        skip.check = TRUE,
                        depend = F, verbose.int = 1, save.other = "attr")

## thresholds (rho) of 0.25 batch -- model 2, 6, 10, 14
param4_1 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                    r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                    rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.5, rho = 0.25)

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model2 <- netsim(est.list[[1]], param4_1, init, control4)
save(sim3.SEIR.model2, file = "results/sim3.SEIR.model2.rda")
rm(sim3.SEI.model2)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model6 <- netsim(est.list[[2]], param4_1, init, control4)
save(sim3.SEIR.model6, file = "results/sim3.SEIR.model6.rda")
rm(sim3.SEIR.model6)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model10 <- netsim(est.list[[3]], param4_1, init, control4)
save(sim3.SEIR.model10, file = "results/sim3.SEIR.model10.rda")
rm(sim3.SEIR.model10)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model14 <- netsim(est.list[[4]], param4_1, init, control4)
save(sim3.SEIR.model14, file = "results/sim3.SEIR.model14.rda")
rm(sim3.SEIR.model14)

## thresholds (rho) of 0.5 batch -- model 3, 7, 11, 15
param4_2 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                      r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                      rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.5, rho = 0.5)

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model3 <- netsim(est.list[[1]], param4_2, init, control4)
save(sim3.SEIR.model3, file = "results/sim3.SEIR.model3.rda")
rm(sim3.SEIR.model3)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model7 <- netsim(est.list[[2]], param4_2, init, control4)
save(sim3.SEIR.model7, file = "results/sim3.SEIR.model7.rda")
rm(sim3.SEIR.model7)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model11 <- netsim(est.list[[3]], param4_2, init, control4)
save(sim3.SEIR.model11, file = "results/sim3.SEIR.model11.rda")
rm(sim3.SEIR.model11)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model15 <- netsim(est.list[[4]], param4_2, init, control4)
save(sim3.SEIR.model15, file = "results/sim3.SEIR.model15.rda")
rm(sim3.SEIR.model15)

## thresholds (rho) of 0.75 batch -- model 4, 8, 12, 16
param4_3 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                      r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                      rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.5, rho = 0.75)

RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model4 <- netsim(est.list[[1]], param4_3, init, control4)
save(sim3.SEIR.model4, file = "results/sim3.SEIR.model4.rda")
rm(sim3.SEIR.model4)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model8 <- netsim(est.list[[2]], param4_3, init, control4)
save(sim3.SEIR.model8, file = "results/sim3.SEIR.model8.rda")
rm(sim3.SEIR.model8)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model12 <- netsim(est.list[[3]], param4_3, init, control4)
save(sim3.SEIR.model12, file = "results/sim3.SEIR.model12.rda")
rm(sim3.SEIR.model12)
RNGkind("L'Ecuyer-CMRG")
set.seed(542435)
sim3.SEIR.model16 <- netsim(est.list[[4]], param4_3, init, control4)
save(sim3.SEIR.model16, file = "results/sim3.SEIR.model16.rda")
rm(sim3.SEIR.model16)


## trends in overall efficacy
## In E-R network
files.ER <- paste0("results/", "sim3.SEIR.model", 1:4, ".rda")
for (i in files.ER) {
  load(i, .GlobalEnv)
}

library(gridGraphics)
library(grid)

plot(sim3.SEIR.model1, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "firebrick", qnts.col = "firebrick", qnts.alpha = 0.2, ylim = c(0, 0.4), ylab = "Number of Recovered (Proportions)")
plot(sim3.SEIR.model2, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "steelblue", qnts.col = "steelblue", qnts.alpha = 0.2, add = T, ylim = c(0, 0.4))
plot(sim3.SEIR.model3, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#41AB5D", qnts.col = "#41AB5D", qnts.alpha = 0.2, add = T, ylim = c(0, 0.4))
plot(sim3.SEIR.model4, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#AA4488", qnts.col = "#AA4488", qnts.alpha = 0.2, add = T, ylim = c(0, 0.4))
abline(v = 168, lty = 2, col = "red") ## start of the recovery
grid.echo()
p1 <- grid.grab() ## grid.draw(p1)
rm(list = ls(pattern = "^sim3.SEIR.model[0-9]"))


files.chain <- paste0("results/", "sim3.SEIR.model", 5:8, ".rda")
for (i in files.chain) {
  load(i, .GlobalEnv)
}

plot(sim3.SEIR.model5, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "firebrick", qnts.col = "firebrick", qnts.alpha = 0.2, ylim = c(0, 0.4), ylab = "Number of Recovered (Proportions)")
plot(sim3.SEIR.model6, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "steelblue", qnts.col = "steelblue", qnts.alpha = 0.2, add = T, ylim = c(0, 0.4))
plot(sim3.SEIR.model7, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#41AB5D", qnts.col = "#41AB5D", qnts.alpha = 0.2, add = T, ylim = c(0, 0.4))
plot(sim3.SEIR.model8, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#AA4488", qnts.col = "#AA4488", qnts.alpha = 0.2, add = T, ylim = c(0, 0.4))
abline(v = 168, lty = 2, col = "red") ## start of the recovery
grid.echo()
p2 <- grid.grab() ## grid.draw(p1)
rm(list = ls(pattern = "^sim3.SEIR.model[0-9]"))


## load all saves simulations
files <- list.files("results/", "sim3.*")
#files <- c(files, list.files("results/", "sim2.SEI.model*"))

require(gtools)
files <- mixedsort(files)

counterfactuals <- list.files("results/", "sim2.SEI.model*")

## get percent infections averted (compared to no-recovery cases)
require(parallel)
PIA.stats <- lapply(seq_len(length(files)), function(i) {
  fl <- files[i]
  load(paste0(getwd(), "/results/", fl))
  fl <- get(gsub(".rda", "", fl)) ## pointer to a loaded file

  if (i %in% c(1:4)) {
    cf <- counterfactuals[1]
  } else if (i %in% c(5:8)) {
    cf <- counterfactuals[2]
  } else if (i %in% c(9:12)) {
    cf <- counterfactuals[3]
  } else {
    cf <- counterfactuals[4]
  }

  ## get counterfactual model (no recovery model)
  load(paste0(getwd(), "/results/", cf))
  cf <- get(gsub(".rda", "", cf))

  ## without intervensions minus with intervensions, divided by baseline (so proportion), over time by models
  out.overall <- sapply(1:1000, function(j) fl$epi$i.num[j,]/cf$epi$i.num[j,], simplify = T)
  out.dem <- sapply(1:1000, function(j) fl$epi$i.num.pidD[j,]/cf$epi$i.num.pidD[j,], simplify = T)
  out.rep <- sapply(1:1000, function(j) fl$epi$i.num.pidR[j,]/cf$epi$i.num.pidR[j,], simplify = T)

  rm(list = ls(pattern = "^sim[0-9]")); rm(cf, fl)
  return(list(out.overall = as.vector(out.overall),
           out.dem = as.vector(out.dem),
           out.rep = as.vector(out.rep)))
})
names(PIA.stats) <- paste0("Model", 1:16)
save(PIA.stats, file = "results/PIA.stats.Rdata")

out.overall <- out.dem <- out.rep <- vector()
for (i in 1:16) {
  out.overall <- cbind(out.overall, apply(t(PIA.stats[[i]]$out.overall)[168:1000,],
                                          2, function(j) median(unlist(j)))) ## median across all time points
  out.dem <- cbind(out.dem, apply(t(PIA.stats[[i]]$out.dem)[168:1000,],
                                          2, function(j) median(unlist(j))))
  out.rep <- cbind(out.rep, apply(t(PIA.stats[[i]]$out.rep)[168:1000,],
                                          2, function(j) median(unlist(j))))
}

colnames(out.overall) <- paste0("Model", 1:16)
ggplot(data = melt(out.overall), aes(x = Var2, y = value)) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted, Overall") +
  theme_bw()

colnames(out.dem) <- paste0("Model", 1:16)
ggplot(data = melt(out.dem), aes(x = Var2, y = value)) + geom_boxplot() +
  xlab("") + ylab("Percent Infections Averted, Among Democrats") +
  theme_bw()

colnames(out.rep) <- paste0("Model", 1:16)
ggplot(data = melt(out.rep), aes(x = Var2, y = value)) + geom_boxplot() +
  xlab("") + ylab("Percent Infections Averted, Among Republicans") +
  theme_bw()






























sim4.SEI.model4 <- netsim(est.list[[4]], param4, init, control4)

model4DT.SEIR2 <- get.summary.stats(sim4.SEI.model4)


p1 <- ggplot(model4DT.SEIR, aes(x = time, y = r.num.mean)) +
  geom_line(color = "grey") +
  geom_line(aes(x = time, y = r.num.pidD.mean), color = "blue") +
  geom_ribbon(aes(x = time, ymin = r.num.pidD.llci, ymax = r.num.pidD.ulci), alpha = 0.2, fill = "blue") +
  geom_line(aes(x = time, y = r.num.pidR.mean), color = "red") +
  geom_ribbon(aes(x = time, ymin = r.num.pidR.llci, ymax = r.num.pidR.ulci), alpha = 0.2, fill = "red")

p2 <- ggplot(model4DT.SEIR2, aes(x = time, y = r.num.mean)) +
  geom_line(color = "grey") +
  geom_line(aes(x = time, y = r.num.pidD.mean), color = "blue") +
  geom_ribbon(aes(x = time, ymin = r.num.pidD.llci, ymax = r.num.pidD.ulci), alpha = 0.2, fill = "blue") +
  geom_line(aes(x = time, y = r.num.pidR.mean), color = "red") +
  geom_ribbon(aes(x = time, ymin = r.num.pidR.llci, ymax = r.num.pidR.ulci), alpha = 0.2, fill = "red")

p1 + p2 + plot_layout(nrow = 1)

