# Bespoke functions: Load these first

#####
 # convert from RR with CIs to simulations

# This function takes the central estimate plus lower and upper 95% CIs of a log normal distribution and 
# produces 10000 simulated values from this distribution.
RR.simulate <- function(central, lower, upper, simulates=10000){

  mu <- log(central)
  sigma.l <- (1/ 1.96) * (mu - log(lower))
  sigma.h <- (1/ 1.96) * (log(upper) - mu )
  sigma <- (sigma.l + sigma.h) / 2
  
  sims <- rnorm(n=simulates, mean=mu, sd=sigma)
  return(exp(sims))
}

# This function produces bootstrapped CIs of means of a vector
Bootstrapper <- function(inputs, simulates = 10000){
  X.mean <- vector("numeric", simulates)
  N.inputs <- length(inputs)
  for (i in 1:simulates) {X.mean[i] <- mean(inputs[sample(1:N.inputs, replace=T)])}
  return(X.mean)
}

# This function provides an estimated mu and sigma for a normal distribution 
# given summary estimates for the log normal distribution associated with it
RR2.simulate <- function(central, lower, upper){
  mu <- log(central)
  sigma.l <- (1/ 1.96) * (mu - log(lower))
  sigma.h <- (1/ 1.96) * (log(upper) - mu )
  sigma <- (sigma.l + sigma.h) / 2
  
  return(list(mu=mu, sigma=sigma))
}

#####################################
######################################


require(MCMCpack)

# using data from table 1 of Rivero-Arias 
# complete cases
# Dead, MRS 0, MRS1, MRS2, MRS3, MRS4, MRS5
N.PSA <- 10000
# Dead/nondead following stroke:
# from Table 1: 319 dead out of 1283
Dead_nonDead <- rbinom(N.PSA, 1283, (319/1283)) / 1283

# mRS following stroke, if not dead
# from Table 1, 24 months
mRS_followingStroke <- rdirichlet(N.PSA, c(61, 143, 111, 82, 24, 4))

# three state reduction:
DepInd_followingStroke <- cbind(apply(mRS_followingStroke[,1:3], 1, sum), apply(mRS_followingStroke[,4:6], 1, sum))

DeadDepInd_followingStroke <- cbind(Dead_nonDead, (1 - Dead_nonDead) * DepInd_followingStroke[,1], (1-Dead_nonDead) * DepInd_followingStroke[,2])

apply(DeadDepInd_followingStroke,2, mean)
colnames(DeadDepInd_followingStroke) <- c("Dead", "Independent", "Dependent")

# for producing figure
d.dead <- density(DeadDepInd_followingStroke[,1])
d.ind <- density(DeadDepInd_followingStroke[,2])
d.dep <- density(DeadDepInd_followingStroke[,3])

maxval <- max(c(d.dead$y, d.dep$y, d.ind$y))
d.dead$y <- d.dead$y / maxval
d.dep$y <- d.dep$y / maxval
d.ind$y <- d.ind$y / maxval

d.dead <- cbind(x = d.dead$x, y=d.dead$y)
d.dep <- cbind(x = d.dep$x, y=d.dep$y)
d.ind <- cbind(x = d.ind$x, y=d.ind$y)

plot(y ~ x, data=d.dead, type="l", ylab="", ylim=c(0,1), xlim=c(0,0.8), xlab="Probability", main="Outcome probabilities following a stroke")
lines(y ~ x, data=d.dep, lty="dashed")
lines(y ~ x, data=d.ind, lty="dashed", lwd=2)
legend("topright", legend=c("dead", "dependent state", "independent state"), lwd=c(1,1,2), lty=c("solid", "dashed", "dashed"))

write.csv(DeadDepInd_followingStroke, "C:/temp/tmp5.csv")

# Using table 3 (24 months column) from Rivero-Arias

s0 <- rnorm(N.PSA, .959, .074)
s1 <- rnorm(N.PSA, .812 , .181)
s2 <- rnorm(N.PSA, .656, .218)
s3 <- rnorm(N.PSA, .545, .277)
s4 <- rnorm(N.PSA, .248, .281)
s5 <- rnorm(N.PSA, .020, .046)

mult.s1 <- s1/s0
mult.s2 <- s2/s0
mult.s3 <- s3/s0
mult.s4 <- s4/s0
mult.s5 <- s5/s0

mean(mult.s1)
mean(mult.s2)
mean(mult.s3)
mean(mult.s4)
mean(mult.s5)


plot(density(mult.s1), xlim=c(0, 1))
lines(density(mult.s2), col="red")
lines(density(mult.s3), col="green")
lines(density(mult.s4), col="blue")
lines(density(mult.s5), col="black", lty="dashed")

# Independent State : mRS 0-2
# Dependent State : mRS 3-5

Stroke.Ind <- mRS_followingStroke[,1:3]
Stroke.Dep <- mRS_followingStroke[,4:6]

Stroke.Dep.sums <- apply(Stroke.Dep, 1, sum)
Stroke.Ind.sums <- apply(Stroke.Ind, 1, sum)

Stroke.Dep <- apply(Stroke.Dep, 2, function (x) x / Stroke.Dep.sums)
Stroke.Ind <- apply(Stroke.Ind, 2, function (x) x / Stroke.Ind.sums)


Stroke.Ind.utils <-   Stroke.Ind[,1] * 1          + Stroke.Ind[,2] * mult.s1   + Stroke.Ind[,3] * mult.s2

Stroke.Dep.utils <-   Stroke.Dep[,1] * mult.s3    + Stroke.Dep[,2] * mult.s4   + Stroke.Dep[,3] * mult.s5



###

n.bootstraps <- 10000
Stroke.Ind.utils.mean <- vector("numeric", n.bootstraps)
Stroke.Dep.utils.mean <- vector("numeric", n.bootstraps)

for (i in 1:n.bootstraps){Stroke.Ind.utils.mean[i] <- mean(Stroke.Ind.utils[sample(1:N.PSA, n.bootstraps, replace=T)])}
for (i in 1:n.bootstraps){Stroke.Dep.utils.mean[i] <- mean(Stroke.Dep.utils[sample(1:N.PSA, n.bootstraps, replace=T)])}
Stroke.Ind.utils.mean <- Bootstrapper(Stroke.Ind.utils)
Stroke.Dep.utils.mean <- Bootstrapper(Stroke.Dep.utils)

plot(density(Stroke.Ind.utils), ylim=c(0, max(density(Stroke.Ind.utils)$y) * 1.5), xlim=c(0,1), main="Utility multipliers following stroke")

lines(density(Stroke.Dep.utils), lty="dashed")
lines(density(Stroke.Ind.utils.mean), lwd=2)
lines(density(Stroke.Dep.utils.mean), lwd=2, lty="dashed")

write.csv(Stroke.Ind.utils.mean[1:1000], "clipboard")
write.csv(Stroke.Dep.utils.mean[1:1000], "clipboard")

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

####### PSAs#######

# Risk (Stroke  | LABN)  

n.PSA <- 10000
# Using Bernoulli simulation approach...
X <- rbinom(n.PSA, 1, 4/50) # 4/50 events per year from original study
n.bootstraps <- 10000

X.means <- Bootstrapper(X)

write.csv(X.means[1:1000], "clipboard")


# NOTE : The following have been updated June 2012 using Friberg et al data

# Use the updated values in the future. 
# From table 3 of Friberg et al

# CHADS SCORE /  Rate (events/100 patient years)  / N
# 0  0.6  13258            
# 1  3.0  23041
# 2  4.2  25813
# 3  7.1  15527
# 4 11.1   8767
# 5 12.5   3315
# 6 13.0    769

Friberg.Data <- data.frame(
  Cscore=0:6,
  rate=c(
    0.6,
    3.0,
    4.2,
    7.1,
    11.1,
    12.5,
    13.0
    ),
  N=c(
    13258,
    23041,
    25813,
    15527,
     8767,
     3315,
      769
    )
  )

# parameters for beta distributions
shape1 <- (Friberg.Data$rate/100) * Friberg.Data$N
shape2 <- Friberg.Data$N - shape1

Friberg.Data <- data.frame(Friberg.Data, shape1=shape1, shape2=shape2)

# Risk (Stroke | CHADs O)	
# Risk (Stroke | CHADS = 1)	
# Risk (Stroke | CHADS = 2)	
# Risk (Stroke | CHADS = 3 )	
# Risk (Stroke | CHADS = 4)	

C.0 <- rbeta(1000, shape1[1], shape2[1])
C.1 <- rbeta(1000, shape1[2], shape2[2])
C.2 <- rbeta(1000, shape1[3], shape2[3])
C.3 <- rbeta(1000, shape1[4], shape2[4])
C.4 <- rbeta(1000, shape1[5], shape2[5])

CHADS <- data.frame(c0 = C.0, c1 = C.1, c2 = C.2, c3 = C.3, c4 = C.4)

#jpeg(filename="graphs/updatedStrokeRisk.jpeg", height=800, width=800)
plot(density(C.0), xlim=c(0,0.15), main="Estimated distribution of annual stroke risk by score", xlab="Annual risk (Proportion)")
lines(density(C.1), lty="dashed")
lines(density(C.2), lwd=2)
lines(density(C.3), lwd=2, lty="dashed")
lines(density(C.4), lty="dotdash")
legend("topright", cex=1.3, legend=c("Score of 0", "Score of 1", "Score of 2", "Score of 3", "Score of 4"), lwd=c(1,1,2,2,1), lty=c(1,2,1,2,4))
#dev.off()



write.csv(CHADS, file="J:/Echo AF/tmp8a.csv")

#####################################################################################
#####################################################################################

# HR (Bleed | Age <75 & Treat = Dab)	
# HR (Bleed | Age =>75 & Treat = Dab)	

# Slightly stumped with this... CIs not given around the percentages, and denominators 
# not provided by age.
# Will do the following but uncertain whether it's the right approach:

N.under75 <- round(10855/3) # assume one third (Dab 150, cf Dab 110 or Warfarin)
N.75orover <- round(7258/3) # assume one third (Dab 150, cf Dab 110 or Warfarin)

probBleed.under75 <- 0.0212
probBleed.75orover <- 0.0510

X.under75 <- rbinom(10000, N.under75, probBleed.under75)
X.75orover <- rbinom(10000, N.75orover, probBleed.75orover)

X.under75 <- X.under75 / N.under75
X.75orover <- X.75orover / N.75orover

quantile(X.under75, c(0.025, 0.5, 0.975))
quantile(X.75orover, c(0.025, 0.5, 0.975))

# producing something pretty for the report...

D.under75 <- density(X.under75)
D.75orover <- density(X.75orover)
maxval <- max(c(D.under75$y, D.75orover$y))
D.under75 <- cbind(x=100 * D.under75$x, y=D.under75$y / maxval)
D.75orover <- cbind(x=100 * D.75orover$x, y=D.75orover$y / maxval)

plot(y ~ x, data=D.under75, ylim=c(0,1), xlab="Risk per year (%)", 
     ylab="", xlim=c(0, 8), main="Annual Risk of Major Bleed", type="l", lwd=2)

lines(y ~ x, data=D.75orover, lwd=2, lty="dashed")
legend("topright", legend=c("Under 75", "75 or older"), lwd=2, lty=c("solid", "dashed"))

write.csv(cbind(X.under75, X.75orover), "C:/temp/tmp4.csv")

# HR (Bleed | Age <75 & Treat = War)	
# HR (Bleed | Age =>75 & Treat =War)	

# Slightly stumped with this... CIs not given around the percentages, and denominators 
# not provided by age.
# Will do the following but uncertain whether it's the right approach:

N.under75 <- round(10855/3) # assume one third (Dab 150, cf Dab 110 or Warfarin)
N.75orover <- round(7258/3) # assume one third (Dab 150, cf Dab 110 or Warfarin)

probBleed.under75 <- 0.0304 # Now using Warfarin numbers instead of Dab numbers
probBleed.75orover <- 0.0437 # Now using Warfarin numbers instead of Dab numbers

X.under75 <- rbinom(10000, N.under75, probBleed.under75)
X.75orover <- rbinom(10000, N.75orover, probBleed.75orover)

X.under75 <- X.under75 / N.under75
X.75orover <- X.75orover / N.75orover

quantile(X.under75, c(0.025, 0.5, 0.975))
quantile(X.75orover, c(0.025, 0.5, 0.975))

write.csv(cbind(X.under75, X.75orover), "C:/temp/tmp4a.csv")


# HR (Bleed | Age < 75 & Treat = Riv)	
# HR (Bleed | Age => 75 & Treat = Riv)	

# Using Patel data, table 3
# HR of major bleed, 1.04 (0.90 to 1.20)
# not differentiated by age

# So
Risk.Riv.under75 <- RR.simulate(1.04, 0.90, 1.20) * X.under75
Risk.Riv.75orover <- RR.simulate(1.04, 0.90, 1.20) * X.75orover

mean(Risk.Riv.under75)
mean(Risk.Riv.75orover)

quantile(Risk.Riv.under75, c(0.025, 0.5, 0.975))
quantile(Risk.Riv.75orover, c(0.025, 0.5, 0.975))

# HR (Stroke | Dab)	

RR.sims.Warf_vs_placebo <- RR.simulate(0.33, 0.24, 0.45) # Lip 2006
RR.sims.D150_vs_Warf <- RR.simulate(0.66, 0.53, 0.82) # Connolly 2009
RR.sims.D150_vs_placebo <- RR.sims.D150_vs_Warf * RR.sims.Warf_vs_placebo

mean(RR.sims.D150_vs_placebo)

quantile(RR.sims.D150_vs_placebo, c(0.025, 0.5, 0.975))

plot(density(RR.sims.D150_vs_placebo), main="", xlab="utility", lwd=2)

# A simplish plotting function...

Make.plot <- function(A=density(RR.sims.Warf_vs_placebo), B=density(RR.sims.D150_vs_Warf), C=density(RR.sims.D150_vs_placebo)){
  Standardise <- function(X) {
    X$y <- X$y / max(X$y)
    tmp <- cbind(y=X$y, x=X$x)
    return(tmp)
    }
  
  sA <- Standardise(A)
  sB <- Standardise(B)
  sC <- Standardise(C)
  
  plot(y ~ x, data=sA, type="l", main="Distribution of RRs", xlab="RR", ylab="", xlim=c(0,1), ylim=c(0,1))
  lines(y ~ x, data=sB, lty="dashed")
  lines(y ~ x, data=sC, lwd=2, lty="dashed")
  legend("bottomright", legend=c("Warfarin vs Placebo", "Dabigatran vs Warfarin", "Dabigatran vs Placebo"), lwd=c(1,1,2), lty=c("solid", "dashed", "dashed"))
}


plot(Standardise(density(RR.sims.Warf_vs_placebo)), main="Distribution of RRs", xlab="utility", xlim=c(0,1))
lines(Standardise(density(RR.sims.D150_vs_Warf)), lty="dashed")
lines(Standardise(density(RR.sims.D150_vs_placebo)), lwd=2, lty="dashed")


# HR (Stroke | Warfarin)	

RR.sims.Warf_vs_placebo <- RR.simulate(0.33, 0.24, 0.45) # Lip 2006
# n.b. these aren't quite *annual* risk ratios as not all studies are one year. Please replace with better 
# estimate if you can find one.

write.csv(RR.sims.Warf_vs_placebo, file="C:/temp/tmpWarStr.csv")


# HR (Stroke | Rivaroxoban)	

#Patel et al 2011

RR.sims.Warf_vs_placebo <- RR.simulate(0.33, 0.24, 0.45) # Lip 2006
RR.sims.Riv_vs_Warf <- RR.simulate(0.88, 0.74, 1.03) # Patel 2011
RR.sims.Riv_vs_placebo <- RR.sims.Riv_vs_Warf * RR.sims.Warf_vs_placebo

mean(RR.sims.Riv_vs_placebo)


# Prob (Death | Stroke)	
# Prob (Dep state | Stroke)	
# Prob (Indep State | Stroke)	 

# See above... Rivero-Arias estimates

# Prob (Death | Bleed & Treat = Dab)	

# Calculated based on Simpson 2009 values (from Holmes' table)

# mean 

#Outcomes: 
#nonfatal GI  : 0.795
#nonfatal IC  : 0.091
#fatal        : 0.114   

# Values for Dirichlet distribution: 
# A =  22.7
# B =  28.4
# C = 198.9
require(MCMCpack)
IC_Dead_GI <- rdirichlet(10000, c(22.7, 28.4, 198.9))
colnames(IC_Dead_GI) <- c("ic", "dead", "gi")

write.csv(IC_Dead_GI, "C:/temp/tmp6.csv")



# Prob (GOS 2 | ICH)	
# Prob (GOS 3 | ICH)	
# Prob (GOS 4 | ICH)	
# Prob (GOS 5 | ICH)	

# GOS 2 0.116
# GOS 3 0.140
# GOS 4 0.079
# GOS 5 0.665

# Dirichlet 
# A = 115.5 # GOS 2
# B = 140 # GOS 3
# C = 79.3 # GOS 4
# D = 665.1 # GOS 5

GOS_states <- rdirichlet (10000, c(115.5, 140, 79.3, 665.1))
colnames(GOS_states) <- c("gos2", "gos3", "gos4", "gos5")

write.csv(GOS_states, "C:/temp/tmp7.csv")


# Util (Dep Stroke)	
# Util (Indep Stroke)	

# Above

# Util (GOS 2)	
# Util (GOS 3)	
# Util (GOS 4)	
# Util (GOS 5)	

# Util (NICH)	
# uniform distribution ranging from 0.997 to 1.000 ?
# slightly different: still want 0.997 as mean value
# so 0.997 +/- 0.003

X <- runif(1000, 0.994, 1.000)

# Cost (Death | Stroke)	
# Sandercock, Table 6
# Mean length of stay: 33-34 days
# cost per night: 200 (150-500) # assuming log normal distribution

# using RR.simulate function below
los <- runif(10000, 33,34)

# Inflation adjustment from 2000 to 2011 levels
# from Curtis PSSRU, RPI figures
inf.adjust2000_2010 <- 222.7/167.7

cpn <- RR.simulate(200, 150, 500)
cpn <- cpn * inf.adjust2000_2010

cost.death.stroke <- los*cpn

mean(cost.death.stroke)

cost.stroke.death.mean <- Bootstrapper(cost.death.stroke)


plot(density(cost.stroke.death.mean))
# Cost (Dep State, Inst)	
# Cost (Dep State, Cont)	
Dep.cost.init <- rnorm(1000, 2830, 62)
Dep.cost.cont <- rnorm(1000, 6386, 325)

# Cost (Indep State, Inst)	
# Cost (Indep State, Cont)	
Ind.cost.init <- rnorm(1000, 542, 15)
Ind.cost.cont <- rnorm(1000, 3195, 165)

strokeCost <- data.frame(dep.init=Dep.cost.init, dep.cont=Dep.cost.cont, ind.init=Ind.cost.init, ind.cont=Ind.cost.cont)

write.csv(strokeCost, "C:/temp/temp11.csv")

# From Holmes' table 4

# all GOS incorporate an initial GI haemorrhage costs 
# GI haemorrhage initial cost Y ~ N( 1261, 25)

# Cost( GOS 2, Inst)	
# GOS 2 intensive care costs 
# mean 15469
# gamma a = 165, b = 94
# GOS 2 rehabilitation cost
# mean 27960
# gamma a = 250, b = 120


GOS_2.cost.init <- rnorm(10000, 1261, 25) + rgamma(10000, scale=165, shape=94) + rgamma(10000, scale=250, shape=120) 
# GI equivalent costs + GOS 2 intensitve care costs + GOS 2 rehabilitation costs

# all GOS incorporate an initial GI haemorrhage costs 
# GI haemorrhage initial cost Y ~ N( 1261, 25)
# Cost( GOS 2, Cont)	
# GOS 2 weekly nursing home cost
# mean 893
# gamma a = 159, b = 6

GOS_2.cost.cont <- 52 * rgamma(10000, scale=159, shape=6)
GOS_2.cost.cont <- Bootstrapper(GOS_2.cost.cont) 
# GOS 2 nursing home costs



# Cost( GOS 3, Inst)	

GOS_3.cost.init <- rnorm(10000, 1261, 25) + rnorm(10000, 8829, 633)
# GI initial cost + intracranial procedures except trauma with haemorrhagic cerebrovascular


# Cost( GOS 3, Cont)	
GOS_3.cost.cont <- rgamma(10000, scale=326, shape=104)
GOS_3.cost.cont <- Bootstrapper(GOS_3.cost.cont)
# GOS_3 annual care costs


# Cost( GOS 4, Inst)	
GOS_4.cost.init <- rnorm(10000, 1261, 25) + rnorm(10000, 8829, 633) + rgamma(10000, scale=385, shape=45)
# GI initial costs + intracranial procedures except trauma with haemorrhagic cerebrovascular + GOS 4 rehabilitation costs

# Cost( GOS 4, Cont)	
GOS_4.cost.cont <- rep(0,10000)
# none

# Cost( GOS 5, Inst)	
GOS_5.cost.init <- rnorm(10000, 1261, 25)
# GI inital cost 

# Cost( GOS 5, Cont)	
GOS_5.cost.cont <- rep(0,10000)
# None


# Bringing the GOSes together

GOS <- data.frame(g2.init=GOS_2.cost.init, g2.cont=GOS_2.cost.cont, 
                  g3.init=GOS_3.cost.init, g3.cont=GOS_3.cost.cont,
                  g4.init=GOS_4.cost.init, g4.cont=GOS_4.cost.cont,
                  g5.init=GOS_5.cost.init, g5.cont=GOS_5.cost.cont)

write.csv(GOS, "C:/temp/temp10.csv")

# Cost (Death | Bleed)	
# Assume equal to death from stroke

# Cost (NICH)	
NICH.cost <- rnorm(1000, 1261, 25)
# NOTE: This should probably be subcategorised into instantaneous and continuous costs
# These are instantaneous costs
# continuous costs are 0
write.csv(NICH.cost, "clipboard")

# Cost (Dabigatran)	
# fixed

# Cost (Warfarin)	

# Cost (Rivaroxoban)



#### Confidence intervals around 



# RR for Warfarin vs Control
#  From Connolly 1994
# central = 0.32
# lower = 0.21
# upper = 0.50

# RR for Warfarin vs Control
#  From Ezekowitz 1992
# central = 0.79
# lower = 0.52
# upper = 0.90

# RR For Warfarin vs Placebo # USING THIS ONE
# from Lip 2006
# Central = 0.33
# lower = 0.24
# upper = 0.45


# Connolly 2009
# RR Dab 150 cf Warfarin, Stroke
# central = 0.66 
# lower = 0.53
# higher = 0.82

# Recalculated RR from Cochrane Review Numbers (Warfarin vs antiplatelet)
## central = 0.6877
# lower = 0.5534
# upper = 0.8548

#RR.sims.Warf_vs_C <- RR.simulate(0.6877, 0.5534, 0.8548) # Cochrane Re-analysis




# > quantile(RR.sims.D110_vs_C.under75, c(0.025, 0.5, 0.975))
#      2.5%       50%     97.5% 
# 0.4483713 0.6394948 0.9157832 
# > quantile(RR.sims.D110_vs_C.75plus, c(0.025, 0.5, 0.975))
#      2.5%       50%     97.5% 
# 0.4216202 0.6068090 0.8678763 
# > 
# > quantile(RR.sims.D150_vs_C.under75, c(0.025, 0.5, 0.975))
#      2.5%       50%     97.5% 
# 0.2963574 0.4325622 0.6346056 
#> quantile(RR.sims.D150_vs_C.75plus, c(0.025, 0.5, 0.975))
#     2.5%       50%     97.5% 
#0.3144194 0.4607246 0.6660244 

###
# Numbers for meta-analysis, from Cochrane Review
# http://onlinelibrary.wiley.com/doi/10.1002/14651858.CD006186.pub2/pdf

#Trial      nT	  NT	  nC	  NC
#ACTIVE W	  64	  3371	106	  3335
#ADASAK I	  9	    335	  16	  336
#AFASAK II	10	  170	  9	    169
#ATHENS	    1	    16	  2	    15
#NASFPEAK	  6	    237	  11	  242
#PATAF	    3	    131	  4	    141
#SPAF Iia	  19	  358	  21	  357
#SPAF IIb	  20	  197	  21	  188

trial.name  <- c("ACTIVE", "ADASAK I" ,"ADASAK II", "ATHENS", "NASFPEAK", "PATAF", "SPAF IIa", "SPAF IIb")
trial.nT <- c(64, 9, 10, 1, 6, 3, 19, 20)
trial.NT <- c(3371, 335, 170, 16, 237, 131, 358, 197)
trial.nC <- c(106, 16, 9, 2, 11, 4, 21, 21)
trial.NC <- c(3335, 336, 169, 15, 242, 141, 357, 188)

Data.Cochrane <- data.frame(study=trial.name, nT=trial.nT, NT=trial.NT, nC=trial.nC, NC=trial.NC)
Data.Cochrane <- data.frame(study=trial.name, nT=trial.nT, NT=trial.NT, nC=trial.nC, NC=trial.NC)

require(metafor)



results.OR <- escalc(measure="OR", ai=nT, n1i=NT, ci=nC, n2i=NC, data=Data.Cochrane)
rma.peto(ai=nT, n1i=NT, ci=nC, n2i=NC, data=Data.Cochrane)
# reproduces results in Cochrane
rma.mh(measure="RR", ai=nT, n1i=NT, ci=nC, n2i=NC, data=Data.Cochrane)
# central = 0.6877
# lower = 0.5534
# upper = 0.8548


X <- read.delim("clipboard")
X

require(metafor)
rma.mh(measure="RR", ai=Wn, n1i=WN, ci=Pn, n2i=PN, data=X)

###############################################################################################
###############################################################################################
# Alternative approach idea

# According to Gage 2001:
#The stroke rate per 100 patient-years without antithrombotic therapy
#increased by a factor of 1.5 (95% CI, 1.3-1.7) for each 1-point increase in the CHADS2
#score:

plus1 <- rnorm(1000, 1.5, (1.5-1.3)/1.96)
plus2 <- Bootstrapper(plus1) * rnorm(1000, 1.5, (1.5-1.3)/1.96)
plus3 <- Bootstrapper(plus2) * rnorm(1000, 1.5, (1.5-1.3)/1.96)
plus4 <- Bootstrapper(plus3) * rnorm(1000, 1.5, (1.5-1.3)/1.96)
plus5 <- Bootstrapper(plus4) * rnorm(1000, 1.5, (1.5-1.3)/1.96)


plot(density(plus1), xlim=c(0,5))
lines(density(plus2))
lines(density(plus3))
lines(density(plus4))
lines(density(plus5))

Additions <- data.frame(plus1, plus2, plus3, plus4, plus5)
write.csv(Additions, "C:/temp/tmp9.csv")


#####################################################################
# Original calculation for CHADS scores (From Gage et al)


# Risk (Stroke | CHADs O)  
# Risk (Stroke | CHADS = 1)	
# Risk (Stroke | CHADS = 2)	
# Risk (Stroke | CHADS = 3 )	
# Risk (Stroke | CHADS = 4)	

# Only information to go on : CHADS2
#Score   No. of Patients   No. of Strokes  NRAF Crude            NRAF Adjusted
#        (n = 1733)        (n = 94)        Stroke Rate per       Stroke Rate,
#                                          100 Patient-Years     (95% CI)
#
#0       120               2               1.2                   1.9 (1.2-3.0)
#1       463               17              2.8                   2.8 (2.0-3.8)
#2       523               23              3.6                   4.0 (3.1-5.1)
#3       337               25              6.4                   5.9 (4.6-7.3)
#4       220               19              8.0                   8.5 (6.3-11.1)
#5       65                6               7.7                   12.5 (8.2-17.5)
#6       5                 2               44.0                  18.2 (10.5-27.4)


#"Confidence intervals calculated from using the binomial method"

# using RR.simulates function produced below
C.0 <- RR.simulate(0.019, 0.012, 0.030, 1000)
C.1 <- RR.simulate(0.028, 0.020, 0.038, 1000)
C.2 <- RR.simulate(0.040, 0.031, 0.051, 1000)
C.3 <- RR.simulate(0.059, 0.046, 0.073, 1000)
C.4 <- RR.simulate(0.085, 0.063, 0.111, 1000)

unsorted.CHADS <- data.frame(c0 = C.0, c1 = C.1, c2 = C.2, c3 = C.3, c4 = C.4)
sorted.CHADS <- data.frame(c0 = sort(C.0), c1 = sort(C.1), c2 = sort(C.2), c3 = sort(C.3), c4 = sort(C.4))

plot(unsorted.CHADS$c1 - unsorted.CHADS$c0, ylim=c(-0.05,0.05))
points(unsorted.CHADS$c2 - unsorted.CHADS$c1, col="red")
points(unsorted.CHADS$c3 - unsorted.CHADS$c2, col="green")
points(unsorted.CHADS$c4 - unsorted.CHADS$c3, col="blue")

plot(sorted.CHADS$c1 - sorted.CHADS$c0, ylim=c(-0.05,0.05))
points(sorted.CHADS$c2 - sorted.CHADS$c1, col="red")
points(sorted.CHADS$c3 - sorted.CHADS$c2, col="green")
points(sorted.CHADS$c4 - sorted.CHADS$c3, col="blue")

plot(density(C.0), xlim=c(0,0.15), main="Estimated distribution of annual stroke risk by score", xlab="Annual risk (Proportion)")
lines(density(C.1), lty="dashed")
lines(density(C.2), lwd=2)
lines(density(C.3), lwd=2, lty="dashed")
lines(density(C.4), lty="dotdash")
legend("topright", legend=c("Score of 0", "Score of 1", "Score of 2", "Score of 3", "Score of 4"), lwd=c(1,1,2,2,1), lty=c(1,2,1,2,4))

# Summary statistics of effect of sorting:
quantile(sorted.CHADS[,5] - sorted.CHADS[,4], c(0, 0.025, 0.975, 1))
quantile(sorted.CHADS[,4] - sorted.CHADS[,3], c(0, 0.025, 0.975, 1))
quantile(sorted.CHADS[,3] - sorted.CHADS[,2], c(0, 0.025, 0.975, 1))
quantile(sorted.CHADS[,2] - sorted.CHADS[,1], c(0, 0.025, 0.975, 1))
quantile(unsorted.CHADS[,5] - unsorted.CHADS[,4], c(0, 0.025, 0.975, 1))
quantile(unsorted.CHADS[,4] - unsorted.CHADS[,3], c(0, 0.025, 0.975, 1))
quantile(unsorted.CHADS[,3] - unsorted.CHADS[,2], c(0, 0.025, 0.975, 1))
quantile(unsorted.CHADS[,2] - unsorted.CHADS[,1], c(0, 0.025, 0.975, 1))
length(which(unsorted.CHADS[,5] < unsorted.CHADS[,4])) / dim(unsorted.CHADS)[1]
length(which(unsorted.CHADS[,4] < unsorted.CHADS[,3])) / dim(unsorted.CHADS)[1]
length(which(unsorted.CHADS[,3] < unsorted.CHADS[,2])) / dim(unsorted.CHADS)[1]
length(which(unsorted.CHADS[,2] < unsorted.CHADS[,1])) / dim(unsorted.CHADS)[1]

length(which(sorted.CHADS[,5] < sorted.CHADS[,4])) / dim(sorted.CHADS)[1]
length(which(sorted.CHADS[,4] < sorted.CHADS[,3])) / dim(sorted.CHADS)[1]
length(which(sorted.CHADS[,3] < sorted.CHADS[,2])) / dim(sorted.CHADS)[1]
length(which(sorted.CHADS[,2] < sorted.CHADS[,1])) / dim(sorted.CHADS)[1]

write.csv(sorted.CHADS, file="C:/temp/tmp8a.csv")



