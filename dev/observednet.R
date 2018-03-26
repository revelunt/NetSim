
##
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
