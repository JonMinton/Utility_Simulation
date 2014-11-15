rm(list=ls())

require(plyr)
require(ggplot2)
require(reshape2)
require(MCMCpack)


source("scripts/functions.r")

##############################################################################################
##############################################################################################
#### DATA
##############################################################################################
# using data from table 1 of Rivero-Arias 
# complete cases
# Dead, MRS 0, MRS1, MRS2, MRS3, MRS4, MRS5
# Dead/nondead following stroke:
# from Table 1: 319 dead out of 1283
ra_tab01_following_stroke <- data.frame(
  dead=319,
  total=1283
  )

# mRS following stroke, if not dead
# from Table 1, 24 months

ra_tab01_24_months <- c(
  mrs0=61,
  mrs1=143,
  mrs2=111,
  mrs3=82,
  mrs4=24,
  mrs5=4
  )


# Using table 3 (24 months column) from Rivero-Arias

ra_tab03_24months <- data.frame(
  var=c(
    "s0",
    "s1",
    "s2",
    "s3",
    "s4",
    "s5"
    ),
  mean=c(
    .959,
    .812,
    .656,
    .545,
    .248,
    .020
    ),
  sd=c(
    .074,
    .181,
    .218,
    .277,
    .281,
    .046
    )
  )

#####################################
######################################
n_psa <- 10000

attach(ra_tab01_following_stroke)
prop_dead <- rbinom(
    n=n_psa,
    size= total,
    prob= dead/total
  )/total
detach(ra_tab01_following_stroke)


mrs_following_stroke <- rdirichlet(
  n_psa, 
  ra_tab01_24_months
)

colnames(mrs_following_stroke) <- names(ra_tab01_24_months)


# three state reduction:
dep_indep_following_nonfatal_stroke <- cbind(
  apply(mrs_following_stroke[,1:3], 1, sum), 
  apply(mrs_following_stroke[,4:6], 1, sum)
  )
colnames(dep_indep_following_nonfatal_stroke) <- c("independent", "dependent")

prop_alive <- 1-prop_dead

dead_dep_indep_following_stroke <- cbind(
  prop_dead, 
  prop_alive * dep_indep_following_nonfatal_stroke[,1], 
  prop_alive * dep_indep_following_nonfatal_stroke[,2]
  )
colnames(dead_dep_indep_following_stroke) <- c("dead", "independent", "dependent")
dead_dep_indep_following_stroke <- data.frame(dead_dep_indep_following_stroke)
dead_dep_indep_following_stroke$psa <- 1:n_psa 


tmp_long <- melt(dead_dep_indep_following_stroke, id.var="psa")

g1 <- ggplot(tmp_long, aes(x=value, group=variable))
g1 + geom_density(aes(fill=variable), alpha=0.5) 


write.csv(DeadDepInd_followingStroke, "C:/temp/tmp5.csv")


######################################################
# Using table 3 (24 months column) from Rivero-Arias
attach(ra_tab03_24months)
util_ests <- data.frame(
  s0 = rnorm(n_psa, mean[var=="s0"], sd[var=="s0"]),
  s1 = rnorm(n_psa, mean[var=="s1"], sd[var=="s1"]),
  s2 = rnorm(n_psa, mean[var=="s2"], sd[var=="s2"]),
  s3 = rnorm(n_psa, mean[var=="s3"], sd[var=="s3"]),
  s4 = rnorm(n_psa, mean[var=="s4"], sd[var=="s4"]),
  s5 = rnorm(n_psa, mean[var=="s5"], sd[var=="s5"])
  )
detach(ra_tab03_24months)

attach(util_ests)
util_mult_ests <- data.frame(
  mult_s1 = s1/s0,
  mult_s2 = s2/s0,
  mult_s3 = s3/s0,
  mult_s4 = s4/s0,
  mult_s5 = s5/s0
  )
detach(util_ests)

tmp <- melt(util_mult_ests)

g1 <- ggplot(tmp, aes(x=value, group=variable, fill=variable, colour=variable))
g1 + geom_density(alpha=0.5)

#######################################################

stroke_ind <- mrs_following_stroke[,1:3]
stroke_dep <- mrs_following_stroke[,4:6]

stroke_dep_sums <- apply(stroke_dep, 1, sum)
stroke_ind_sums <- apply(stroke_ind, 1, sum)

stroke_dep <- apply(stroke_dep, 2, function (x) x / stroke_dep_sums)
stroke_ind <- apply(stroke_ind, 2, function (x) x / stroke_ind_sums)


attach(util_mult_ests)
stroke_ind_utils <-   stroke_ind[,"mrs0"] * 1        + 
                      stroke_ind[,"mrs1"] * mult_s1  + 
                      stroke_ind[,"mrs2"] * mult_s2

stroke_dep_utils <-   stroke_dep[,"mrs3"] * mult_s3  + 
                      stroke_dep[,"mrs4"] * mult_s4  +
                      stroke_dep[,"mrs5"] * mult_s5

detach(util_mult_ests)


###

n_boots <- 10000

stroke_ind_utils_mean <- bootstrap(stroke_ind_utils)
stroke_dep_utils_mean <- bootstrap(stroke_dep_utils)

plot(density(stroke_ind_utils), ylim=c(0, max(density(stroke_ind_utils)$y) * 1.5), xlim=c(0,1), main="Utility multipliers following stroke")

lines(density(stroke_dep_utils), lty="dashed")
lines(density(stroke_ind_utils_mean), lwd=2)
lines(density(stroke_dep_utils_mean), lwd=2, lty="dashed")

write.csv(Stroke.Ind.utils.mean[1:1000], "clipboard")
write.csv(Stroke.Dep.utils.mean[1:1000], "clipboard")

########################################################################
#####
# Now to do something similar for utilities following a bleed
# Using following mappying between Glasgow Outcome Scale and modified Rankin Scale

# GOS 1: Dead                 -> mRS 6 (Dead)
# GOS 2: Vegetative State     -> mRS 6 (Dead)
# GoS 3: Severely Disabled    -> mRS 4-5
# GoS 4: Moderately Disabled  -> mRS 2-3
# GOS 5: Good Recovery        -> mRS 0-1

# Independent State : mRS 0-2
# Dependent State : mRS 3-5

GOS_5 <- mRS_followingStroke[,1:2]
GOS_4 <- mRS_followingStroke[,3:4]
GOS_3 <- mRS_followingStroke[,5:6]


GOS_5.sums <- apply(GOS_5, 1, sum)
GOS_4.sums <- apply(GOS_4, 1, sum)
GOS_3.sums <- apply(GOS_3, 1, sum)

GOS_5 <- apply(GOS_5, 2, function (x) x / GOS_5.sums)
GOS_4 <- apply(GOS_4, 2, function (x) x / GOS_4.sums)
GOS_3 <- apply(GOS_3, 2, function (x) x / GOS_3.sums)


GOS_5.utils <- GOS_5[,1] * 1        + GOS_5[,2] * mult.s1
GOS_4.utils <- GOS_4[,1] * mult.s2  + GOS_4[,2] * mult.s3
GOS_3.utils <- GOS_3[,1] * mult.s4  + GOS_3[,2] * mult.s5


plot(density(GOS_5.utils), ylim=c(0, max(density(GOS_5.utils)$y) * 1.5), xlim=c(0,1), main="Utility multipliers following IC Bleed")

lines(density(GOS_4.utils), col="red")
lines(density(GOS_3.utils), col="green")

GOS_5.mean <- vector("numeric", n.bootstraps)
GOS_4.mean <- vector("numeric", n.bootstraps)
GOS_3.mean <- vector("numeric", n.bootstraps)

for (i in 1:n.bootstraps){GOS_5.mean[i] <- mean(GOS_5.utils[sample(1:N.PSA, n.bootstraps, replace=T)])}
for (i in 1:n.bootstraps){GOS_4.mean[i] <- mean(GOS_4.utils[sample(1:N.PSA, n.bootstraps, replace=T)])}
for (i in 1:n.bootstraps){GOS_3.mean[i] <- mean(GOS_3.utils[sample(1:N.PSA, n.bootstraps, replace=T)])}

lines(density(GOS_5.mean), lwd=2)
lines(density(GOS_4.mean), lwd=2, col="red")
lines(density(GOS_3.mean), lwd=2, col="green")

quantile(GOS_5.mean, c(0.025, 0.5, 0.975))
quantile(GOS_4.mean, c(0.025, 0.5, 0.975))
quantile(GOS_3.mean, c(0.025, 0.5, 0.975))

#write.csv(Stroke.Ind.utils.mean[1:1000], "clipboard")
#write.csv(Stroke.Dep.utils.mean[1:1000], "clipboard")

write.csv(cbind(gos3=GOS_3.mean, gos4=GOS_4.mean, gos5=GOS_5.mean), "C:/temp/tmp8.csv")
mean(GOS_5.mean)

# mean util following ICH
# Frequencies: GOS 2, 3, 4,  5
gos_followingICH <- rdirichlet(N.PSA, c(115.5, 140, 79.3, 665.1))
util_followingICH <- gos_followingICH[,1] * 0  + gos_followingICH[,2] * GOS_3.mean + gos_followingICH[,3] * GOS_4.mean + gos_followingICH[,4] * GOS_5.mean
# Not doing this this way for now. Will do for each GOS


