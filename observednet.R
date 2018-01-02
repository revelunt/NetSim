
## EpiModel with observed social network

## we first need a panel of observations, and convert them to a networkDynamic object.
## after converting them to a list of static network objects,

library(foreign)
library(plyr)

if(!("haven" %in% installed.packages()[,"Package"])) install.packages("haven")
if(!("data.table" %in% installed.packages()[,"Package"])) install.packages("data.table")
library(haven)
library(data.table)
require(igraph)
require(sna)
require(car)
require(parallel)

file_path <- "/Users/songh95/Dropbox/GitHub/Korean2012ElectionProject"

## load the dataset from SPSS and set as data.table object
dat <- haven::read_spss(paste(file_path, "/Dat/DiscussionForumThreeWavePanel(N=341).sav", sep = ""))
setDT(dat)

## load the node dataset for subsetting the network
dat2 <- haven::read_spss(paste(file_path, "/Dat/Survey Data for Network Matrix Data.sav", sep = ""))
setDT(dat2)
dat2 <- na.omit(dat2[, c('r_id', 'p_image', 'm_image', 'ide_self', 'evalcrit1', 'evalcrit2',
                         'policy_c', 'policy_l', 'knowtotal', 'talk', 'interest', "female", "edu", "age", "opleader", "follower", "income")])

## this yeilds a total of 312 cases
lengthn <- dat2[, .N]
lengthn
vids <- dat2$r_id

## remove "NaN" in data -- W1
dat[is.na(pv311), pv311 := 0 ]
dat[is.na(pv313), pv313 := 0 ]
dat[is.na(pv317), pv317 := 0 ]
## remove "NaN" in data -- W2
dat[is.na(kv194), kv194 := 0 ]
dat[is.na(kv196), kv196 := 0 ]
dat[is.na(kv200), kv200 := 0 ]
## remove "NaN" in data -- W3
dat[is.na(hv253), hv253 := 0 ]
dat[is.na(hv255), hv255 := 0 ]
dat[is.na(hv259), hv259 := 0 ]

## motivation for using online forum
consistency.motivation <- dat[vids, .(as.numeric(pv18),
                                      as.numeric(pv19),
                                      as.numeric(pv20),
                                      as.numeric(pv21),
                                      as.numeric(pv23),
                                      as.numeric(pv24))]
understanding.motivation <- dat[vids, pv13:pv16]
hedomic.motivation <- dat[vids, pv27:pv29]

##---------------------------------##
## create a dependent network list ##
##---------------------------------##

## pre-wave: Nov 13 to Nov 26; Wave 1 survey: Nov 27 to Nov 29,
## Wave 2 survey: Dec 11th to 13th, and Wave 3 survey: Dec 21th to 23th

net <- read.csv(paste(file_path, "/Dat/Reading_1113-1126_Participants(N=341)_Count(N=160836).csv", sep = ""))
net2 <- read.csv(paste(file_path, "/Dat/Reading_1127-1219_Participants(N=341)_Count(N=160836).csv", sep = ""))
net2 <- data.frame(reading.time = net2$Reading.Time, reader.id = net2$Reader.Id, poster.id = net2$Poster.Id)
net <- rbind(net, net2)
setDT(net)
net[, reading.date := as.Date(reading.time, format = "%Y-%m-%d %H:%M:%S")]

date.range <- unique(sort(net[, reading.date]))

thresholds <- sapply(1:length(date.range), function(i) {
  g[[i]] <- net[reading.date %in% date.range[i],]
  g[[i]] <- data.frame(g[[i]][,2], g[[i]][,3])
  setDT(g[[i]]); g[[i]][, count :=1]
  g[[i]][, sum(count), by = c("poster.id", "reader.id")][poster.id %in% vids & reader.id %in% vids, mean(V1)]
})

g <- list()
for (i in 1:length(date.range)) {

  g[[i]] <- net[reading.date %in% date.range[i],]
  g[[i]] <- data.frame(g[[i]][,2], g[[i]][,3])
  g[[i]] <- graph.data.frame(g[[i]], directed = TRUE, vertices = 1:341)
  g[[i]] <- induced_subgraph(g[[i]], vids = vids)
  g[[i]] <- as.matrix(as_adj(g[[i]]))
  rownames(g[[i]]) <- colnames(g[[i]]) <- vids
  g[[i]] <- sna::event2dichot(g[[i]], method = "absolute", thresh = thresholds[i])
  g[[i]] <- as.network(g[[i]])

}

## convert the list to the networkDynamic object
## this takes some time....
require(networkDynamic)
g <- networkDynamic(network.list = g)

## outputs ##
# Neither start or onsets specified, assuming start = 0
# Onsets and termini not specified, assuming each network in network.list should have a discrete spell of length 1
# Argument base.net not specified, using first element of network.list instead
# Created net.obs.period to describe network
# Network observation period info:
# Number of observation spells: 1
# Maximal time range observed: 0 until 27
# Temporal mode: discrete
# Time unit: step
# Suggested time increment: 1

## set some network attributes
g %n% 'net.obs.period' <- list(
  observations = list(c(0,27)),
  mode = "discrete",
  time.increment = 1,
  time.unit = "day")

## time-invariant covariates
activate.vertex.attribute(g, "interest", value = dat[vids, pv165:pv166][, rowMeans(.SD), by = vids][,V1])
activate.vertex.attribute(g, "internal.efficacy", value = dat[vids, pv126:pv129][, rowMeans(.SD), by = vids])
activate.vertex.attribute(g, "consistency.motivation", value = apply(consistency.motivation, 1, mean))
activate.vertex.attribute(g, "understanding.motivation", value = apply(understanding.motivation, 1, mean))
activate.vertex.attribute(g, "age", value = dat[vids, as.numeric(age)/10])
activate.vertex.attribute(g, "gender", value = dat[vids, as.numeric(sex) - 1])
activate.vertex.attribute(g, "edu", value = dat[vids, as.numeric(edu)])

## time-variying covariates
activate.vertex.attribute(g, "partyid", value = split(
  partyid <- cbind(
    rep.col(dat[vids, as.numeric(canpref1)], 7),
    rep.col(dat[vids, as.numeric(canpref2)], 14),
    rep.col(dat[vids, as.numeric(canpref3)], 7)),  seq(nrow(partyid))))



