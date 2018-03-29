
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
infect

## Here, we specify that Republicans (R) has a baseline prevalence of 15% (n = 184, randomly assigned),
## whereas there are no nodes in Democrats (D) infected.
status.vector <- c(rep(0, nD), rbinom(nR, 1, 0.15))
status.vector <- ifelse(status.vector == 1, "i", "s")

## PROGRESSION MODEL
## This is a constant hazard function (individual-level stochastic process)
## Republicans have a higher transition rates than Democrats
## and those who have higher accuracy motivations have a lower transition rates
exposed.to.infectious.progress
get_prev.exposed.included

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
load("results/net_and_dx_list.rda", verbose = T)
load("results/sim2.SEI.model1.rda", verbose = T)
load("results/sim2.SEI.model2.rda", verbose = T)
load("results/sim2.SEI.model3.rda", verbose = T)
load("results/sim2.SEI.model4.rda", verbose = T)

## summary of simulations (using first simulation as an example)
# print(sim2.SEI.model1)
# summary(sim2.SEI.model1, at = 500)

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


pdf("draft/Figure 2. Prevalence.pdf", paper = 'a4r', width = 12, height = 7)
ggplot(dat.plot, aes(time, model1)) + geom_line(color = "firebrick") + ## model1 no. of exposed
  #geom_ribbon(aes(ymin = model1.llci, ymax = model1.ulci), alpha = 0.1) +
  geom_line(aes(time, model1.suspect), color = "firebrick", linetype = "dotted") +

  geom_line(aes(time, model2), color = "steelblue") +
  #geom_ribbon(aes(x = time, ymin = model2.llci, ymax = model2.ulci), fill = "grey", alpha = 0.1) +
  geom_line(aes(time, model2.suspect), color = "steelblue", linetype = "dotted") +

  geom_line(aes(time, model3), color = "#41AB5D") +
  #geom_ribbon(aes(ymin = model3.llci, ymax = model3.ulci), fill = "blue", alpha = 0.1) +
  geom_line(aes(time, model3.suspect), color = "#41AB5D", linetype = "dotted") +

  geom_line(aes(time, model4), color = "#AA4488") +
  #geom_ribbon(aes(ymin = model4.llci, ymax = model4.ulci), fill = "red", alpha = 0.1) +
  geom_line(aes(time, model4.suspect), color = "#AA4488", linetype = "dotted") +
  theme_bw() + ylab("Cumulative Fraction Exposed to Misinformation") + xlab("Time (Hours)")

## adoption rates:
ggplot(model1DT.SEI, aes(time, se.flow.mean)) +
  geom_smooth(fill = "firebrick", colour = "firebrick", method = "loess", alpha = 0.1) +
  geom_smooth(data = model2DT.SEI, aes(time, se.flow.mean),
              fill = "steelblue", colour = "steelblue", method = "loess", alpha = 0.1) +
  geom_smooth(data = model3DT.SEI, aes(time, se.flow.mean),
              fill = "#41AB5D", colour = "#41AB5D", method = "loess", alpha = 0.1) +
  geom_smooth(data = model4DT.SEI, aes(time, se.flow.mean),
              fill = "#AA4488", colour = "#AA4488", method = "loess", alpha = 0.1) +
  theme_bw() + ylab("Suspected to Exposed Flow") + xlab("Time (Hours)")
dev.off()

## prevalence table
tb3 <- matrix(NA, nrow = 6, ncol = 4)
rownames(tb3) <- c("overall", "95% CIs", "Among Dem", "95% CIs", "Among Rep", "95% CIs")
colnames(tb3) <- c("E-R", "Chain", "SW No-Homophily", "SW Homophily")

datlist <- list(model1DT.SEI, model2DT.SEI, model3DT.SEI, model4DT.SEI)
for (i in 1:4) {
dat <- datlist[[i]]
tb3[,i] <- c(format(dat[, mean((e.num.mean + i.num.mean)/num.mean*100)],
                    trim = T, digits = 3, nsmall = 3, justify = "centre"),
             paste0("[", format(dat[, mean((e.num.llci + i.num.llci)/num.mean*100)],
                                trim = T, digits = 3, nsmall = 3),
                 ", ", format(dat[, mean((e.num.ulci + i.num.ulci)/num.mean*100)],
                              trim = T, digits = 3, nsmall = 3), "]"),
             format(dat[, mean((e.num.pidD.mean + i.num.pidD.mean)/num.mean*100, na.rm = T)],
                    trim = T, digits = 3, nsmall = 3, justify = "centre"),
             paste0("[", format(dat[, mean((e.num.pidD.llci + i.num.pidD.llci)/num.mean*100, na.rm = T)],
                                trim = T, digits = 3, nsmall = 3),
                    ", ", format(dat[, mean((e.num.pidD.ulci + i.num.pidD.ulci)/num.mean*100, na.rm = T)],
                                 trim = T, digits = 3, nsmall = 3), "]"),
             format(dat[, mean((e.num.pidR.mean + i.num.pidR.mean)/num.mean*100, na.rm = T)],
                    trim = T, digits = 3, nsmall = 3, justify = "centre"),
             paste0("[", format(dat[, mean((e.num.pidR.llci + i.num.pidR.llci)/num.mean*100, na.rm = T)],
                                trim = T, digits = 3, nsmall = 3),
                    ", ", format(dat[, mean((e.num.pidR.ulci + i.num.pidR.ulci)/num.mean*100, na.rm = T)],
                                 trim = T, digits = 3, nsmall = 3), "]"))
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
# require(ndtv)
# nw10 <- get_network(sim2.SEI.model1, sim = 10)
# nw10 %n% "slice.par" <- list(start = 1, end = 1000, interval = 100, aggregate.dur = 1, rule = 'latest')
# compute.animation(nw10, animation.mode = 'kamadakawai', chain.direction = 'reverse', verbose = FALSE)
# render.d3movie(nw10, vertex.cex = 0.9, vertex.col = "pid",
#                edge.col = "darkgrey",
#                vertex.border = "lightgrey",
#                displaylabels = FALSE,
#                vertex.tooltip = function(slice){paste('name:',slice%v%'vertex.names','<br>',
#                                                       'status:', slice%v%'testatus')})


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

param3 <-  param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                     r.ei.rate = 0.3, d.ei.rate = 0.07,
                     rec.rate = 0.00297619, rec.start = 168)
init <- init.net(status.vector = status.vector)

## revised module for recovery process
## this assumes homogenous recovery after average duration
recovery.delayed.random

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
recovery.correction

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
  load(i, .GlobalEnv, verbose = T)
}

setDT(model1DT.sim3 <- get.summary.stats(sim3.SEIR.model1))
setDT(model2DT.sim3 <- get.summary.stats(sim3.SEIR.model2))
setDT(model3DT.sim3 <- get.summary.stats(sim3.SEIR.model3))
setDT(model4DT.sim3 <- get.summary.stats(sim3.SEIR.model4))
datlist <- list(model1DT.sim3, model2DT.sim3, model3DT.sim3, model4DT.sim3)

tb5 <- matrix(NA, nrow = 8, ncol = 4)
rownames(tb5) <- c("E-R", "",  "Chain", "", "SW No-Homophily", "", "SW Homophily", "")
colnames(tb5) <- c("Asocial", "Threshods (0.25)", "Threshods (0.50)", "Threshods (0.75)")
for (i in 1:4) {
  dat <- datlist[[i]]
  tb5[1,i] <- paste0(format(dat[, mean((e.num.mean + i.num.mean)/num.mean*100)],
                      trim = T, digits = 3, nsmall = 3, justify = "centre"),
               paste0(" [", format(dat[, mean((e.num.llci + i.num.llci)/num.mean*100)],
                                  trim = T, digits = 3, nsmall = 3),
                      ", ", format(dat[, mean((e.num.ulci + i.num.ulci)/num.mean*100)],
                                   trim = T, digits = 3, nsmall = 3), "]"))
  tb5[2,i] <- paste0(format(dat[, mean(ir.flow.mean)],
                              trim = T, digits = 3, nsmall = 3, justify = "centre"),
                       paste0(" [", format(dat[, mean(ir.flow.llci)],
                                          trim = T, digits = 3, nsmall = 3),
                              ", ", format(dat[, mean(ir.flow.ulci)],
                                           trim = T, digits = 3, nsmall = 3), "]"))
}

pdf("draft/Fig3a.pdf", paper = 'a4r', width = 12, height = 7)
plot(sim3.SEIR.model1, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "firebrick", qnts.col = "firebrick", qnts.alpha = 0.2, ylab = "Number of Recovered (Proportions)")
plot(sim3.SEIR.model2, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "steelblue", qnts.col = "steelblue", qnts.alpha = 0.2, add = T)
plot(sim3.SEIR.model3, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#41AB5D", qnts.col = "#41AB5D", qnts.alpha = 0.2, add = T)
plot(sim3.SEIR.model4, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#AA4488", qnts.col = "#AA4488", qnts.alpha = 0.2, add = T)
abline(v = 168, lty = 2, col = "red") ## start of the recovery
dev.off()
rm(list = ls(pattern = "^sim3.SEIR.model[0-9]"))


files.chain <- paste0("results/", "sim3.SEIR.model", 5:8, ".rda")
for (i in files.chain) {
  load(i, .GlobalEnv, verbose = T)
}
setDT(model5DT.sim3 <- get.summary.stats(sim3.SEIR.model5))
setDT(model6DT.sim3 <- get.summary.stats(sim3.SEIR.model6))
setDT(model7DT.sim3 <- get.summary.stats(sim3.SEIR.model7))
setDT(model8DT.sim3 <- get.summary.stats(sim3.SEIR.model8))
datlist <- list(model5DT.sim3, model6DT.sim3, model7DT.sim3, model8DT.sim3)
for (i in 1:4) {
  dat <- datlist[[i]]
  tb5[3,i] <- paste0(format(dat[, mean((e.num.mean + i.num.mean)/num.mean*100)],
                            trim = T, digits = 3, nsmall = 3, justify = "centre"),
                     paste0(" [", format(dat[, mean((e.num.llci + i.num.llci)/num.mean*100)],
                                         trim = T, digits = 3, nsmall = 3),
                            ", ", format(dat[, mean((e.num.ulci + i.num.ulci)/num.mean*100)],
                                         trim = T, digits = 3, nsmall = 3), "]"))
  tb5[4,i] <- paste0(format(dat[, mean(ir.flow.mean)],
                            trim = T, digits = 3, nsmall = 3, justify = "centre"),
                     paste0(" [", format(dat[, mean(ir.flow.llci)],
                                         trim = T, digits = 3, nsmall = 3),
                            ", ", format(dat[, mean(ir.flow.ulci)],
                                         trim = T, digits = 3, nsmall = 3), "]"))
}

pdf("draft/Fig3b.pdf", paper = 'a4r', width = 12, height = 7)
plot(sim3.SEIR.model5, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "firebrick", qnts.col = "firebrick", qnts.alpha = 0.2, ylab = "Number of Recovered (Proportions)")
plot(sim3.SEIR.model6, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "steelblue", qnts.col = "steelblue", qnts.alpha = 0.2, add = T)
plot(sim3.SEIR.model7, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#41AB5D", qnts.col = "#41AB5D", qnts.alpha = 0.2, add = T)
plot(sim3.SEIR.model8, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#AA4488", qnts.col = "#AA4488", qnts.alpha = 0.2, add = T)
abline(v = 168, lty = 2, col = "red") ## start of the recovery
dev.off()
rm(list = ls(pattern = "^sim3.SEIR.model[0-9]"))


files.SWHNO <- paste0("results/", "sim3.SEIR.model", 9:12, ".rda")
for (i in files.SWHNO) {
  load(i, .GlobalEnv, verbose = T)
}
setDT(model9DT.sim3 <- get.summary.stats(sim3.SEIR.model9))
setDT(model10DT.sim3 <- get.summary.stats(sim3.SEIR.model10))
setDT(model11DT.sim3 <- get.summary.stats(sim3.SEIR.model11))
setDT(model12DT.sim3 <- get.summary.stats(sim3.SEIR.model12))
datlist <- list(model9DT.sim3, model10DT.sim3, model11DT.sim3, model12DT.sim3)
for (i in 1:4) {
  dat <- datlist[[i]]
  tb5[5,i] <- paste0(format(dat[, mean((e.num.mean + i.num.mean)/num.mean*100)],
                            trim = T, digits = 3, nsmall = 3, justify = "centre"),
                     paste0(" [", format(dat[, mean((e.num.llci + i.num.llci)/num.mean*100)],
                                         trim = T, digits = 3, nsmall = 3),
                            ", ", format(dat[, mean((e.num.ulci + i.num.ulci)/num.mean*100)],
                                         trim = T, digits = 3, nsmall = 3), "]"))
  tb5[6,i] <- paste0(format(dat[, mean(ir.flow.mean)],
                            trim = T, digits = 3, nsmall = 3, justify = "centre"),
                     paste0(" [", format(dat[, mean(ir.flow.llci)],
                                         trim = T, digits = 3, nsmall = 3),
                            ", ", format(dat[, mean(ir.flow.ulci)],
                                         trim = T, digits = 3, nsmall = 3), "]"))
}
pdf("draft/Fig3c.pdf", paper = 'a4r', width = 12, height = 7)
plot(sim3.SEIR.model9, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "firebrick", qnts.col = "firebrick", qnts.alpha = 0.2, ylab = "Number of Recovered (Proportions)")
plot(sim3.SEIR.model10, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "steelblue", qnts.col = "steelblue", qnts.alpha = 0.2, add = T)
plot(sim3.SEIR.model11, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#41AB5D", qnts.col = "#41AB5D", qnts.alpha = 0.2, add = T)
plot(sim3.SEIR.model12, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#AA4488", qnts.col = "#AA4488", qnts.alpha = 0.2, add = T)
abline(v = 168, lty = 2, col = "red") ## start of the recovery
dev.off()
rm(list = ls(pattern = "^sim3.SEIR.model[0-9]"))


files.SWHYES <- paste0("results/", "sim3.SEIR.model", 13:16, ".rda")
for (i in files.SWHYES) {
  load(i, .GlobalEnv, verbose = T)
}
setDT(model13DT.sim3 <- get.summary.stats(sim3.SEIR.model13))
setDT(model14DT.sim3 <- get.summary.stats(sim3.SEIR.model14))
setDT(model15DT.sim3 <- get.summary.stats(sim3.SEIR.model15))
setDT(model16DT.sim3 <- get.summary.stats(sim3.SEIR.model16))
datlist <- list(model13DT.sim3, model14DT.sim3, model15DT.sim3, model16DT.sim3)
for (i in 1:4) {
  dat <- datlist[[i]]
  tb5[7,i] <- paste0(format(dat[, mean((e.num.mean + i.num.mean)/num.mean*100)],
                            trim = T, digits = 3, nsmall = 3, justify = "centre"),
                     paste0(" [", format(dat[, mean((e.num.llci + i.num.llci)/num.mean*100)],
                                         trim = T, digits = 3, nsmall = 3),
                            ", ", format(dat[, mean((e.num.ulci + i.num.ulci)/num.mean*100)],
                                         trim = T, digits = 3, nsmall = 3), "]"))
  tb5[8,i] <- paste0(format(dat[, mean(ir.flow.mean)],
                            trim = T, digits = 3, nsmall = 3, justify = "centre"),
                     paste0(" [", format(dat[, mean(ir.flow.llci)],
                                         trim = T, digits = 3, nsmall = 3),
                            ", ", format(dat[, mean(ir.flow.ulci)],
                                         trim = T, digits = 3, nsmall = 3), "]"))
}
pdf("draft/Fig3d.pdf", paper = 'a4r', width = 12, height = 7)
plot(sim3.SEIR.model13, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "firebrick", qnts.col = "firebrick", qnts.alpha = 0.2, ylab = "Number of Recovered (Proportions)")
plot(sim3.SEIR.model14, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "steelblue", qnts.col = "steelblue", qnts.alpha = 0.2, add = T)
plot(sim3.SEIR.model15, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#41AB5D", qnts.col = "#41AB5D", qnts.alpha = 0.2, add = T)
plot(sim3.SEIR.model16, popfrac = T, y = "r.num", qnts = 0.9, mean.smooth = F, qnts.smooth = F,
     mean.col = "#AA4488", qnts.col = "#AA4488", qnts.alpha = 0.2, add = T)
abline(v = 168, lty = 2, col = "red") ## start of the recovery
dev.off()
rm(list = ls(pattern = "^sim3.SEIR.model[0-9]"))

correction.models <- lapply(ls(pattern = "^model[0-9]*DT"), function(i) get(i))
save(correction.models, file = "results/correction.models.Rdata")

tb5 <- cbind(rep(c("I: Prevalence", "R: Incidence"), 4), tb5)
require(xtable); print.xtable(tb5)

## load all saves simulations
files <- list.files("results/", "sim3.*")

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
  out.overall <- t(sapply(1:1000, function(j) 1 - (fl$epi$i.num[j,]/cf$epi$i.num[j,]), simplify = T))
  out.dem <- t(sapply(1:1000, function(j) 1 - (fl$epi$i.num.pidD[j,]/cf$epi$i.num.pidD[j,]), simplify = T))
  out.rep <- t(sapply(1:1000, function(j) 1 - (fl$epi$i.num.pidR[j,]/cf$epi$i.num.pidR[j,]), simplify = T))

  rm(list = ls(pattern = "^sim[0-9]")); rm(cf, fl)
  return(list(out.overall = as.vector(out.overall),
           out.dem = as.vector(out.dem),
           out.rep = as.vector(out.rep)))
})
names(PIA.stats) <- paste0("Model", 1:16)
save(PIA.stats, file = "results/PIA.stats.Rdata")

out.overall <- out.dem <- out.rep <- vector()
for (i in 1:16) {
  out.overall <- cbind(out.overall, apply(PIA.stats[[i]]$out.overall,
                                          2, function(j) mean(unlist(j)))) ## median across all time points
  out.dem <- cbind(out.dem, apply(PIA.stats[[i]]$out.dem,
                                          2, function(j) mean(unlist(j), na.rm = T)))
  out.rep <- cbind(out.rep, apply(PIA.stats[[i]]$out.rep,
                                          2, function(j) mean(unlist(j), na.rm = T)))
}

fig4.dat <- melt(out.overall)
fig4.dat$type <- "overall"
colnames(fig4.dat) <- c("sim", "Model", "PIA", "type")

fig4.dat <- rbind(fig4.dat, data.frame(
  sim = melt(out.dem)$Var1,
  Model = melt(out.dem)$Var2,
  PIA = melt(out.dem)$value,
  type = "Democrats"
))

fig4.dat <- rbind(fig4.dat, data.frame(
  sim = melt(out.rep)$Var1,
  Model = melt(out.rep)$Var2,
  PIA = melt(out.rep)$value,
  type = "Republicans"
))

setDT(fig4.dat)
## In E-R network
p4_1 <- ggplot(fig4.dat[Model %in% c("Model1", "Model2", "Model3", "Model4"), ],
               aes(x = Model, y = PIA, fill = factor(type))) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted") +
  theme_bw() + scale_x_discrete(labels = c("Model1" = "Model 1 \n(Asocial)",
                                           "Model2" = "Model 2 \n(Thresholds = 0.25)",
                                           "Model3" = "Model 3 \n(Thresholds = 0.50)",
                                           "Model4" = "Model 4 \n(Thresholds = 0.75)")) +
  guides(fill = FALSE) + ggtitle("Panel A: E-R Network") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))

## in Chain network
p4_2 <- ggplot(fig4.dat[Model %in% c("Model5", "Model6", "Model7", "Model8"), ],
               aes(x = Model, y = PIA, fill = factor(type))) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted") +
  theme_bw() + scale_x_discrete(labels = c("Model5" = "Model 5 \n(Asocial)",
                                           "Model6" = "Model 6 \n(Thresholds = 0.25)",
                                           "Model7" = "Model 7 \n(Thresholds = 0.50)",
                                           "Model8" = "Model 8 \n(Thresholds = 0.75)")) +
  guides(fill = FALSE) + ggtitle("Panel B: Chain Network") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))


## in SW Network without homophily
p4_3 <- ggplot(fig4.dat[Model %in% c("Model9", "Model10", "Model11", "Model12"), ],
               aes(x = Model, y = PIA, fill = factor(type))) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted") +
  theme_bw() + scale_x_discrete(labels = c("Model9" = "Model 9 \n(Asocial)",
                                           "Model10" = "Model 10 \n(Thresholds = 0.25)",
                                           "Model11" = "Model 11 \n(Thresholds = 0.50)",
                                           "Model12" = "Model 12 \n(Thresholds = 0.75)")) +
  guides(fill = FALSE) + ggtitle("Panel C: SW No-Homophily Network") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))

## in SW Network with homophily (baseline)
p4_4 <- ggplot(fig4.dat[Model %in% c("Model13", "Model14", "Model15", "Model16"), ],
               aes(x = Model, y = PIA, fill = factor(type))) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted") +
  theme_bw() + scale_x_discrete(labels = c("Model13" = "Model 13 \n(Asocial)",
                                           "Model14" = "Model 14 \n(Thresholds = 0.25)",
                                           "Model15" = "Model 15 \n(Thresholds = 0.50)",
                                           "Model16" = "Model 16 \n(Thresholds = 0.75)")) +
  guides(fill = guide_legend(title = "")) + ggtitle("Panel D: SW Homophily Network") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))
p4_4 <- p4_4 + theme(legend.position = c(0.8, 0.2), legend.text=element_text(size = 11))

## all in one plot
require(patchwork)
p4_1 + p4_2 + p4_3 + p4_4 + plot_layout(nrow = 2, byrow = T)



## ------------------------------------------------------- ##
## SENSITIVITY ANALYSIS, VARYING DEGREE OF CORRECTION PROB ##
## ------------------------------------------------------- ##

rm(list = ls())
source("dev/helper-functions.R")
load("results/net_and_dx_list.rda")
rm(dx.list, netstats)

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

param5_1 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                      r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                      rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.25, rho = 0.25)
## simulations for correction.prob = .25, rho = .25
names <- paste0("sensitivity.model", c(2,6,10,14), ".correction.prob.0.25")
require(pbapply)
pblapply(1:4, function(i) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(542435)
  sim.out <- netsim(est.list[[i]], param5_1, init, control4)
  assign(names[i], sim.out)
  save(sim.out, list = names[i], file = paste0("results/", names[i], ".rda"))
  rm(sim.out); print(paste("finished processing files,", i, "of 4"))
})

param5_2 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                      r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                      rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.25, rho = 0.5)
## simulations for correction.prob = .25, rho = .5
names <- paste0("sensitivity.model", c(3,7,11,15), ".correction.prob.0.25")
pblapply(1:4, function(i) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(542435)
  sim.out <- netsim(est.list[[i]], param5_2, init, control4)
  assign(names[i], sim.out)
  save(sim.out, list = names[i], file = paste0("results/", names[i], ".rda"))
  rm(sim.out); print(paste("finished processing files,", i, "of 4"))
})

param5_3 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                      r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                      rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.25, rho = 0.75)
## simulations for correction.prob = .25, rho = .75
names <- paste0("sensitivity.model", c(4,8,12,16), ".correction.prob.0.25")
pblapply(1:4, function(i) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(542435)
  sim.out <- netsim(est.list[[i]], param5_3, init, control4)
  assign(names[i], sim.out)
  save(sim.out, list = names[i], file = paste0("results/", names[i], ".rda"))
  rm(sim.out); print(paste("finished processing files,", i, "of 4"))
})


## correction prob = 0.75
param5_4 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                      r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                      rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.75, rho = 0.25)
## simulations for correction.prob = .75, rho = .25
names <- paste0("sensitivity.model", c(2,6,10,14), ".correction.prob.0.75")
require(pbapply)
pblapply(1:4, function(i) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(542435)
  sim.out <- netsim(est.list[[i]], param5_4, init, control4)
  assign(names[i], sim.out)
  save(sim.out, list = names[i], file = paste0("results/", names[i], ".rda"))
  rm(sim.out); print(paste("finished processing files,", i, "of 4"))
})

param5_5 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                      r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                      rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.75, rho = 0.5)
## simulations for correction.prob = .25, rho = .5
names <- paste0("sensitivity.model", c(3,7,11,15), ".correction.prob.0.75")
pblapply(1:4, function(i) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(542435)
  sim.out <- netsim(est.list[[i]], param5_5, init, control4)
  assign(names[i], sim.out)
  save(sim.out, list = names[i], file = paste0("results/", names[i], ".rda"))
  rm(sim.out); print(paste("finished processing files,", i, "of 4"))
})

param5_6 <- param.net(inf.prob = 0.16, pid.diff.rate = 0.04, act.rate = 0.164,
                      r.ei.rate = 0.3, d.ei.rate = 0.07, rec.rand = TRUE,
                      rec.rate = 0.00297619, rec.start = 168, correction.prob = 0.75, rho = 0.75)
## simulations for correction.prob = .25, rho = .75
names <- paste0("sensitivity.model", c(4,8,12,16), ".correction.prob.0.75")
pblapply(1:4, function(i) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(542435)
  sim.out <- netsim(est.list[[i]], param5_6, init, control4)
  assign(names[i], sim.out)
  save(sim.out, list = names[i], file = paste0("results/", names[i], ".rda"))
  rm(sim.out); print(paste("finished processing files,", i, "of 4"))
})


## load saved simulations
## correction.prob = 0.25

files <- list.files("results/", "sensitivity.model[0-9]*.correction.prob.0.[0-9]*.rda")
require(gtools)
files <- gtools::mixedsort(files)

files2 <- gtools::mixedsort(list.files("results/", "sim3.SEIR.model[0-9]*"))
files2 <- files2[-c(1,5,9,13)] ## remove single-shot, asocial correction with correction.prob = 0.5

vec <- matrix(c(files[grepl("sensitivity.model[0-9]*.correction.prob.0.25.rda", files)], files2,
  files[grepl("sensitivity.model[0-9]*.correction.prob.0.75.rda", files)]), nrow = 12, byrow = F)

vec <- vector()
for (i in 1:length(files2)) {
  vec <- c(vec, files[(2*i-1)], files2[i], files[2*i])
}

files <- vec; rm(files2)
counterfactuals <- list.files("results/", "sim2.SEI.model*")

PIA.stats2 <- lapply(seq_len(length(files)), function(i) {
  fl <- files[i]
  load(paste0(getwd(), "/results/", fl))
  fl <- get(gsub(".rda", "", fl)) ## pointer to a loaded file

  if (i %in% c(1:9)) {
    cf <- counterfactuals[1]
  } else if (i %in% c(10:18)) {
    cf <- counterfactuals[2]
  } else if (i %in% c(19:27)) {
    cf <- counterfactuals[3]
  } else {
    cf <- counterfactuals[4]
  }

  ## get counterfactual model (no recovery model)
  load(paste0(getwd(), "/results/", cf))
  cf <- get(gsub(".rda", "", cf))

  ## without intervensions minus with intervensions, divided by baseline (so proportion), over time by models
  out.overall <- t(sapply(1:1000, function(j) 1 - (fl$epi$i.num[j,]/cf$epi$i.num[j,]), simplify = T))
  out.dem <- t(sapply(1:1000, function(j) 1 - (fl$epi$i.num.pidD[j,]/cf$epi$i.num.pidD[j,]), simplify = T))
  out.rep <- t(sapply(1:1000, function(j) 1 - (fl$epi$i.num.pidR[j,]/cf$epi$i.num.pidR[j,]), simplify = T))

  rm(list = ls(pattern = "^sim[0-9]")); rm(cf, fl)
  return(list(out.overall = as.vector(out.overall),
              out.dem = as.vector(out.dem),
              out.rep = as.vector(out.rep)))
})
out.overall2 <- out.dem2 <- out.rep2 <- vector()
for (i in 1:36) {
  out.overall2 <- cbind(out.overall2, apply(PIA.stats2[[i]]$out.overall,
                                          2, function(j) mean(unlist(j)))) ## median across all time points
  out.dem2 <- cbind(out.dem2, apply(PIA.stats2[[i]]$out.dem,
                                  2, function(j) mean(unlist(j), na.rm = T)))
  out.rep2 <- cbind(out.rep2, apply(PIA.stats2[[i]]$out.rep,
                                  2, function(j) mean(unlist(j), na.rm = T)))
}

files <- gsub(".rda", "", files); files <- gsub("sensitivity.", "", files); files <- gsub("sim3.SEIR.", "", files)
colnames(out.overall2) <- files

fig5.dat <- melt(out.overall2)
setDT(fig5.dat)
colnames(fig5.dat) <- c("sim", "Model", "PIA")

fig5.dat[grep("model2", Model, value = F), threshods := 0.25]
fig5.dat[grep("model6", Model, value = F), threshods := 0.25]
fig5.dat[grep("model10", Model, value = F), threshods := 0.25]
fig5.dat[grep("model14", Model, value = F), threshods := 0.25]

fig5.dat[grep("model3", Model, value = F), threshods := 0.5]
fig5.dat[grep("model7", Model, value = F), threshods := 0.5]
fig5.dat[grep("model11", Model, value = F), threshods := 0.5]
fig5.dat[grep("model15", Model, value = F), threshods := 0.5]

fig5.dat[grep("model4", Model, value = F), threshods := 0.75]
fig5.dat[grep("model8", Model, value = F), threshods := 0.75]
fig5.dat[grep("model12", Model, value = F), threshods := 0.75]
fig5.dat[grep("model16", Model, value = F), threshods := 0.75]

fig5.dat[grep("correction.prob.0.25", Model, value = F), correction.prob := 0.25]
fig5.dat[grep("correction.prob.0.75", Model, value = F), correction.prob := 0.75]
fig5.dat[is.na(correction.prob), correction.prob := 0.5]

fig5.dat[, Model := gsub(".correction.prob.0.[0-9]*", "", Model)]

## E-R network
p5_1 <- ggplot(fig5.dat[Model %in% c("model2", "model3", "model4"), ],
       aes(x = Model, y = PIA, fill = factor(correction.prob))) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted") +
  theme_bw() + scale_x_discrete(labels = c("model2" = "Model 2 \n(Thresholds = 0.25)",
                                           "model3" = "Model 3 \n(Thresholds = 0.50)",
                                           "model4" = "Model 4 \n(Thresholds = 0.75)")) +
  guides(fill = FALSE) + ggtitle("Panel A: E-R Network") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))

## Chain network
p5_2 <- ggplot(fig5.dat[Model %in% c("model6", "model7", "model8"), ],
       aes(x = Model, y = PIA, fill = factor(correction.prob))) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted") +
  theme_bw() + scale_x_discrete(labels = c("model6" = "Model 6 \n(Thresholds = 0.25)",
                                           "model7" = "Model 7 \n(Thresholds = 0.50)",
                                           "model8" = "Model 8 \n(Thresholds = 0.75)")) +
  guides(fill = FALSE) + ggtitle("Panel B: Chain Network") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))

## SW No-homophily
p5_3 <- ggplot(fig5.dat[Model %in% c("model10", "model11", "model12"), ],
       aes(x = Model, y = PIA, fill = factor(correction.prob))) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted") +
  theme_bw() + scale_x_discrete(labels = c("model10" = "Model 10 \n(Thresholds = 0.25)",
                                           "model11" = "Model 11 \n(Thresholds = 0.50)",
                                           "model12" = "Model 12 \n(Thresholds = 0.75)")) +
  guides(fill = FALSE) + ggtitle("Panel C: SW No-Homophily Network") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))

## SW Homophily
p5_4 <- ggplot(fig5.dat[Model %in% c("model14", "model15", "model16"), ],
       aes(x = Model, y = PIA, fill = factor(correction.prob))) + geom_boxplot() + ## across all simulations, plot medians
  xlab("") + ylab("Proportion Infections Averted") +
  theme_bw() + scale_x_discrete(labels = c("model14" = "Model 14 \n(Thresholds = 0.25)",
                                           "model15" = "Model 15 \n(Thresholds = 0.50)",
                                           "model16" = "Model 16 \n(Thresholds = 0.75)")) +
  guides(fill = guide_legend(title = "Probability of alters \nsending corrections",
                             title.position = "top")) + ggtitle("Panel D: SW Homophily Network") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        legend.position = c(0.2, 0.2), legend.text=element_text(size = 12),
        legend.direction = "horizontal", legend.title.align = 0.5)

## all in one plot
p5_1 + p5_2 + p5_3 + p5_4 + plot_layout(nrow = 2, byrow = T)


