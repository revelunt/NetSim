
library(magrittr)
library(statnet)
library(ggnetwork)
library(tidyverse)

# First, in our simple example, let's set our parameters ahead of time
n.people <- 100
p.infection <- 0.5
pct.starting.infected <- 0.01
max.time <- 5000
contact.rate <- 0.05 / 2  # prob. of contact between any 2 people in the network

### Step 1: Set up world ###
# Create a random graph, where edges are placed randomly.  This is called a
# Bernoulli or an Erdos-Renyi random graph
set.seed(919)
net <- rgraph(n = n.people, tprob = contact.rate) %>%
  symmetrize %>%  # Make symmetric -- doesn't matter who is connected to who
  network(matrix.type = "adjacency")  # convert to statnet

# Chose a subset of people who are infected initially
infected <- sample(
  x = c(T, F),      # people can be infected (T) or susceptible (F)
  size = n.people,  # create a vector that is n.people long
  replace = T,      # with replacement, so that >1 person can be infected
  prob = c(pct.starting.infected, 1 - pct.starting.infected)
)

### Step 2: Choose an edge ###
# For each step, we're going to choose an edge at random, and then, if the edge
# is discordant, flip a coin to determine whether the susceptible person gets
# infected.

# First, create an edgelist...
el <- as.edgelist(net) %>% as.data.frame
colnames(el) <- c("from", "to")

  # ... attach the values of infected...
el <- el %>% mutate(from.infected = infected[from], to.infected = infected[to],
         # ... and create a discordant variable
         discordant = (from.infected != to.infected))

# Next, choose an edge at random
random.edge <- sample(nrow(el), size = 1)

# Check if the edge is discordant
el[random.edge, "discordant"]  # it's not, so we do nothing.

# For the example, I'm going to speed this up by choosing a discordant edge
discordant.edge <- sample(which(el$discordant), size = 1)



### Step 3: Flip a coin to see if the person gets infected ###

# Now, flip a coin to see if the uninfected person gets infected
el[discordant.edge, ]  # Person 62 is the suceptible, but we will want to be
# able to determine that without looking manually

# A little tricky indexing to pull out the ID of the susceptible person...
who.susceptible <- with(
  el[discordant.edge, ],
  c(from, to)[!c(from.infected, to.infected)]
)

# Flip the coin to determine if infection spreads
infected[who.susceptible] <- sample(c(T, F), size = 1,
                                     prob = c(p.infection, 1 - p.infection))

### Step 4: Repeat ###
# To repeat this process, we actually embed steps 1 and 2 in a loop.

# Set up a list with the output
infections <- vector(length = max.time, mode = "list")

# Save what we already did as the first iteration
infections[[1]] <- infected

# Quick aside -- what did we create?
head(infections)

# Drop the "from.infected", "to.infected", and "discordant" columns from el,
# because they'll actually change with every iteration
el %<>% select(-from.infected, -to.infected, -discordant)

# Next, run the loop
set.seed(27708)
for (t in 2:max.time) {
  infections[[t]] <- infections[[t - 1]]

  # Pick an edge at random
  random.edge <- sample(nrow(el), size = 1)

  # If the edge is discordant, flip a coin to decide if the infection spreads
  if (with(el[random.edge, ],
           infections[[t]][from] != infections[[t]][to])) {

    who.susceptible <- with(
      el[random.edge, ],
      c(from, to)[!c(infections[[t]][from], infections[[t]][to])]
    )

    infections[[t]][who.susceptible] <- sample(
      c(T, F),
      size = 1,
      prob = c(p.infection, 1 - p.infection)
    )
  }
}


