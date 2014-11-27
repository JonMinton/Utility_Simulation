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
n_boots <- 10000

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
  independent=apply(mrs_following_stroke[,1:3], 1, sum), 
  dependent=apply(mrs_following_stroke[,4:6], 1, sum)
  )


prop_alive <- 1-prop_dead

dead_dep_indep_following_stroke <- cbind(
  dead=prop_dead, 
  independent=prop_alive * dep_indep_following_nonfatal_stroke[,"independent"], 
  dependent=prop_alive * dep_indep_following_nonfatal_stroke[,"dependent"]
  )

dead_dep_indep_following_stroke <- data.frame(dead_dep_indep_following_stroke)
dead_dep_indep_following_stroke$psa <- 1:n_psa 
dead_dep_indep_following_stroke <- dead_dep_indep_following_stroke[
  c("psa", "dead", "dependent", "independent")
  ]

tmp_long <- melt(dead_dep_indep_following_stroke, 
                 variable.name="state", 
                 id.var="psa"
                 )

# g1 <- ggplot(tmp_long, aes(x=value, group=state))
# g1 + geom_density(aes(fill=state, linetype=state), alpha=0.5) + coord_cartesian(xlim=c(0, 1))

# geom_area

g1 <- qplot(x=psa, y=value, group=state, colour=state, fill=state, data=tmp_long, geom="area")
g2 <- g1 + scale_fill_grey() + scale_colour_grey() + labs(x="PSA number", y="Cumulative proportion")

tiff("figures/fig_02.tiff", 1000, 1000)
print(g2)
dev.off()

# write.csv(dead_dep_indep_following_stroke, "data/generated/dead_dep_indep_following_stroke.csv")


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
util_ests$psa <- 1:n_psa

util_ests <- util_ests[
  c("psa", "s0", "s1", "s2", "s3", "s4", "s5")
  ]

util_mult_ests <- transform(
  util_ests,
  s1=s1/s0,
  s2=s2/s0,
  s3=s3/s0,
  s4=s4/s0,
  s5=s5/s0,
  s0=NULL
  )

util_mult_ests <- util_mult_ests[
  c("psa", "s1", "s2","s3", "s4", "s5")
  ]


tmp <- melt(util_mult_ests, id.var="psa", variable.name="state")

g1 <- ggplot(tmp, aes(x=value, group=state, fill=state, linetype=state, colour=state))
g2 <- g1 + geom_density(alpha=0.5) + scale_fill_grey() + scale_colour_grey() + coord_flip()
g3 <- g2 + labs(x="HRQL estimates", y="Density of estimates") + geom_vline(xintercept=c(0,1), linetype="dashed")
print(g3)

tiff("figures/fig_03.tiff", 1000, 1000)
print(g3)
dev.off()

#######################################################

stroke_ind <- mrs_following_stroke[,1:3]
stroke_dep <- mrs_following_stroke[,4:6]

stroke_dep_sums <- apply(stroke_dep, 1, sum)
stroke_ind_sums <- apply(stroke_ind, 1, sum)

stroke_dep <- apply(stroke_dep, 2, function (x) x / stroke_dep_sums)
stroke_ind <- apply(stroke_ind, 2, function (x) x / stroke_ind_sums)

stroke_dep <- as.data.frame(stroke_dep)

stroke_dep <- data.frame(
  psa=1:n_psa,
  stroke_dep
  )

stroke_ind <- as.data.frame(stroke_ind)
stroke_ind <- data.frame(
  psa=1:n_psa,
  stroke_ind
  )

tmp_dep <- melt(stroke_dep,
                id.var="psa",
                variable.name="state"
                )

tmp_ind <- melt(stroke_ind,
                id.var="psa",
                variable.name="state"
                )

g1 <- qplot(x=psa, y=value, group=state, colour=state, fill=state, data=tmp_dep, geom="area")
g2 <- g1 + scale_fill_grey() + scale_colour_grey() + labs(x="PSA number", y="Cumulative proportion")

tiff("figures/fig_04.tiff", 1000, 1000)
print(g2)
dev.off()


g1 <- qplot(x=psa, y=value, group=state, colour=state, fill=state, data=tmp_ind, geom="area")
g2 <- g1 + scale_fill_grey() + scale_colour_grey() + labs(x="PSA number", y="Cumulative proportion")

tiff("figures/fig_05.tiff", 1000, 1000)
print(g2)
dev.off()

rm(tmp_dep, tmp_ind)


stroke_util_mult_ind_dep <- data.frame(
  psa=1:n_psa,
  independent = stroke_ind[,"mrs0"] * 1        + 
                stroke_ind[,"mrs1"] * util_mult_ests$s1  + 
                stroke_ind[,"mrs2"] * util_mult_ests$s2,
  dependent   = stroke_dep[,"mrs3"] * util_mult_ests$s3  + 
                stroke_dep[,"mrs4"] * util_mult_ests$s4  +
                stroke_dep[,"mrs5"] * util_mult_ests$s5
  )


tmp <- melt(stroke_util_mult_ind_dep, id.var="psa", variable.name="state")

g1 <- ggplot(data=tmp) + geom_density(
  aes(x=value, group=state, colour=state, fill=state),
  alpha=0.5
  )
g2 <- g1 + scale_fill_grey() + scale_colour_grey() + labs(x="HRQL multiplier estimate") 
g3 <- g2 + geom_vline(xintercept=c(0,1), linetype="dashed")

tiff("figures/fig_06.tiff", 1000, 1000)
print(g3)
dev.off()

rm(tmp)
###


utils_bootstrapped_means <- data.frame(
  boot= 1:10000,                                    
  independent= bootstrap(stroke_util_mult_ind_dep$independent),
  dependent = bootstrap(stroke_util_mult_ind_dep$dependent)
  )

tmp <- melt(utils_bootstrapped_means, id.var="boot", variable.name="state")

g1 <- ggplot(data=tmp) + geom_density(
  aes(x=value, group=state, colour=state, fill=state),
  alpha=1
)
g2 <- g1 + scale_fill_grey() + scale_colour_grey() + labs(x="Bootstrapped mean of HRQL multiplier estimate") 
g3 <- g2 + geom_vline(xintercept=c(0,1), linetype="dashed")
tiff("figures/fig_07.tiff", 1000, 1000)
print(g3)
dev.off()


rm(tmp)


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

gos_5 <- mrs_following_stroke[,c("mrs0", "mrs1")]
gos_4 <- mrs_following_stroke[,c("mrs2", "mrs3")]
gos_3 <- mrs_following_stroke[,c("mrs4", "mrs5")]


gos_5_sums <- apply(gos_5, 1, sum)
gos_4_sums <- apply(gos_4, 1, sum)
gos_3_sums <- apply(gos_3, 1, sum)

gos_5 <- apply(gos_5, 2, function (x) x / gos_5_sums)
gos_4 <- apply(gos_4, 2, function (x) x / gos_4_sums)
gos_3 <- apply(gos_3, 2, function (x) x / gos_3_sums)

gos_utils <- data.frame(
  gos5 = gos_5[,"mrs0"] * 1                   + gos_5[,"mrs1"]  * util_mult_ests$s1,
  gos4 = gos_4[,"mrs2"] * util_mult_ests$s2   + gos_4[,"mrs3"]  * util_mult_ests$s3,
  gos3 = gos_3[,"mrs4"] * util_mult_ests$s4   + gos_3[,"mrs5"]  * util_mult_ests$s5
  )

gos_utils_bootstrapped_means <- data.frame(
  gos5 = bootstrap(gos_utils$gos5),
  gos4 = bootstrap(gos_utils$gos4),
  gos3 = bootstrap(gos_utils$gos3)
  )


