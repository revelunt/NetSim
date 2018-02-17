
## setup file for target statistics

## Derivation for:
## proportion of cross-party discussion ties
## (data from 2008-2009 ANES panel)

require(haven)
require(data.table)
anespanel0809 <- read_spss("dev/anes2008_2009panel.sav")
setDT(anespanel0809)

## "w9l1" and "w9l3" = ego partisanship
## "w9zd12_1" = alter partisanship

anespanel0809[w9l1 == "-1", pid := w9l3]
anespanel0809[w9l3 == "-1", pid := w9l1]
## 1 = Democrates (n = 980), 2 = Republican (n = 866), 3 = independents (n = 768), & 4 = something else (n=119)
anespanel0809[, table(pid)]
## set strengths of pid (3 = strong D/R, 2 = not so strong D/R) for party identifiers
anespanel0809[pid %in% c(1,2) & !(w9l5 %in% c(-7,-5)), pidst := 4 - w9l5]
anespanel0809[pid %in% c(1,2) & (w9l5 %in% c(-7,-5)), pidst := 2]

anespanel0809[(pid %in% c(3,4)), w9l6] ## 1. Closer to the Republican Party, 2. Closer to the Democratic Party, & 3. Neither
anespanel0809[(pid %in% c(3,4)) & (w9l6 %in% c(1)), c("pid", "pidst") := list(2, 1)]
anespanel0809[(pid %in% c(3,4)) & (w9l6 %in% c(2)), c("pid", "pidst") := list(1, 1)]
anespanel0809[(pid %in% c(3,4)) & (w9l6 %in% c(3)), c("pid", "pidst") := list(0, 0)]
anespanel0809[(pid %in% c(0,3,-7)), c("pid", "pidst") := list(0, 0)]

## distribution of ego partyids
anespanel0809[, table(pid, pidst, exclude = NULL)]

## process alter partyid
anespanel0809[w9zd12_1 == 1, alter1pid := 1]
anespanel0809[w9zd12_1 == 2, alter1pid := 2]
anespanel0809[w9zd13_1 == 1, alter1pid := 2]
anespanel0809[w9zd13_1 == 2, alter1pid := 1]
anespanel0809[w9zd16_1 == 1, alter1pid := 1]
anespanel0809[w9zd16_1 == 2, alter1pid := 2]

anespanel0809[w9zd12_2 == 1, alter2pid := 1]
anespanel0809[w9zd12_2 == 2, alter2pid := 2]
anespanel0809[w9zd13_2 == 1, alter2pid := 2]
anespanel0809[w9zd13_2 == 2, alter2pid := 1]
anespanel0809[w9zd16_2 == 1, alter2pid := 1]
anespanel0809[w9zd16_2 == 2, alter2pid := 2]

anespanel0809[w9zd12_3 == 1, alter3pid := 1]
anespanel0809[w9zd12_3 == 2, alter3pid := 2]
anespanel0809[w9zd13_3 == 1, alter3pid := 2]
anespanel0809[w9zd13_3 == 2, alter3pid := 1]
anespanel0809[w9zd16_3 == 1, alter3pid := 1]
anespanel0809[w9zd16_3 == 2, alter3pid := 2]


pid_mixing <- anespanel0809[, table(pid, alter1pid, exclude = c(0, NA))]
pid_mixing <- pid_mixing + anespanel0809[, table(pid, alter2pid, exclude = c(0, NA))]
pid_mixing <- pid_mixing + anespanel0809[, table(pid, alter3pid, exclude = c(0, NA))]

colnames(pid_mixing) <- c("alter_D", "alter_R")
rownames(pid_mixing) <- c("ego_D", "ego_R")
pid_mixing.prob <- pid_mixing/sum(pid_mixing)

## proportion of cross-partisan ties ( = .23)
sum(pid_mixing.prob[lower.tri(pid_mixing.prob)|upper.tri(pid_mixing.prob)])

## this number is similar to Eveland et al (2017) where
## they report 19.5% of respondents reporting exposure to difference in standard name generator
## yet their data suggests that this number might underestimate the true exposure
## see http://www.sciencedirect.com/science/article/pii/S0378873317302113#bib0215
## if we assume we just add one more alter for those who report no cross-partisan discussion ties,
## and additional 39% of ties are cross-partisan ties (by their calculation),
## this would increase the observed proportion a bit.
## there are a total of 2400 individuals either R or D.
anespanel0809[pid %in% c(1,2), .N]
## among them, 1258 are democrats and 1142 are republicans
anespanel0809[pid %in% c(1:2), table(pid)]

## among 2400, 626 individuals do not have any R or D discussants as their primary discussants (first alter)
## (that does not facter into pid_mixing matrix above)
anespanel0809[pid %in% c(1,2) & is.na(alter1pid), .N] ## 626
## among 1774 individuals who have at least one of either D or R alters,
## 148 have only one discussants with the same preference
anespanel0809[pid == alter1pid & is.na(alter2pid), .N] ## 148
## 110 have two discussants with same preferences
anespanel0809[pid == alter1pid & pid == alter2pid & is.na(alter3pid), .N] # 110
## and 775 has three discussants with same preferences
anespanel0809[pid == alter1pid & pid == alter2pid & pid == alter3pid, .N] # 775
## therefore 34.4% of respondents reporting at least one exposure to difference

## including those who have no R/D discussants at all (n = 626), there are 1659 potential dyadic cases to be added.
## out of that, 647 dyadic interactions would be cross-party ties (39%, by Eveland et al's calculation),
## therefore ((585 + 571) + 647) / (5001 + 1659) = 27.07% are cross-party heterogenous dyads, and
## ((2079 + 1766) + 1012) / (5001 + 1659) = 72.92% are in-party homogenous dyads




## for nodefactor target statistics
## we look at the mean number of alters by pid group
anespanel0809[pid %in% c(1,2) & is.na(alter1pid), n_partisan_alter := 0]
anespanel0809[pid %in% c(1,2) & !is.na(alter1pid), n_partisan_alter := 1]
anespanel0809[pid %in% c(1,2) & !is.na(alter1pid) & !is.na(alter2pid), n_partisan_alter := 2]
anespanel0809[pid %in% c(1,2) & !is.na(alter1pid) & !is.na(alter2pid) & !is.na(alter3pid), n_partisan_alter := 3]

## degree distributions appears to be no difference by pid, therefore set to both mean degree
## also, approximately 26% of alters has more than two degrees, therefore set concurrent to be 0.26*n
anespanel0809[pid %in% 1:2, table(pid, n_partisan_alter, exclude = NULL)]
chisq.test(anespanel0809[pid %in% 1:2, table(pid, n_partisan_alter, exclude = NULL)], correct = FALSE)

## mean degree = average degree / maximun possible degree
anespanel0809[pid %in% 1:2, mean(n_partisan_alter)] / 3 ## = 0.66



## Fake news spread and consumption
require(haven)
require(data.table)
FN.Pew <- read_spss("dev/FNdata_Pew.sav")
setDT(FN.Pew)

## get PID of respondents (including leaners)
FN.Pew[, table(party, partyln, exclude = NULL)]
FN.Pew[party %in% 1:2, partyid := 3 - party] ## partyid 1 = Democrats, 2 = Republicans
FN.Pew[is.na(party), partyid := 3 - partyln]
FN.Pew[partyid == -6, partyid := car::recode(partyid, "-6 = NA")]
FN.Pew[, table(partyid, exclude = NULL)]

FN.Pew[, table(partyid, pew3) + table(partyid, pew4)] /
  rowSums(FN.Pew[, table(partyid, pew3) + table(partyid, pew4)]) ## 1 == Yes, 2 == No

## Among democrats, 14% have ever shared a poliical story online, either completely made up or later they found to be not true
## Democrats are bit higher for this tendency, approximately 18% of them shared such a story online
