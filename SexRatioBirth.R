#############################################################################
#                                                                           #  
# THE DYNAMICS OF SON PREFERENCE, TECHNOLOGY DIFFUSION, AND FERTILITY       #
# DECLINE UNDERLYING DISTORTED SEX RATIOS AT BIRTH: A SIMULATION APPROACH   #
#                                                                           #
#                                                                           #
# Ridhi KASHYAP - kashyap@demogr.mpg.de/ridhi.kashyap@nuffield.ox.ac.uk     #
# Francisco VILLAVICENCIO - villavicencio@imada.sdu.dk                      #                  
#                                                                April 2016 #
#############################################################################


############################################################
### 1)  PACKAGES, LIBRARIES AND WORKING DIRECTORY        ###  
############################################################

# Clear workspace
rm(list = ls())

# Packages
# install.packages("snowfall")
# install.packages("parallel")
# install.packages(c("snowfall","parallel","RColorBrewer"), Ncpus=2)

library(snowfall)
library(parallel)

# Set the working directory for the code and the subfolders 'Data' and 'Sensitivity'
setwd(".")



#### Try a testrun #######################################################################
# Temporarily disable everything below this point
if (FALSE) {
## Define globals used by FunctionsSRB.R BEFORE sourcing it
iniYear   <- 1980
niter     <- 10              # for smoke test
minAge    <- 15
maxAge    <- 50
cohortLen <- 5
periodStudy <- seq(iniYear, iniYear + niter - 1, 1)
cohortNames <- seq(iniYear - maxAge,
                   iniYear + trunc(niter / cohortLen, 0) * cohortLen,
                   cohortLen)

sensitivity <- FALSE
model <- "Empirical"
country <- "Korea"

iniPop    <- 50000
startYear <- 1945
nsim      <- 2
srb       <- 0.5122

alphaPar <- 0.075; gammaPar <- 0.20
rhoPar   <- 0.50;  phiPar   <- 7
sigmaPar <- 1.70;  betaPar  <- 0.20

inputData     <- get.data(country)
spData        <- funcSonPref(country)
sonPrefMx     <- spData$sonPrefMx
sonPrefTheory <- spData$sonPrefTheory
startPop      <- initialization(iniPop, startYear)

# Quick error fixes for a test run
fertMod <- 1   # use the standard UN fertility schedule (“fertility5”)
sonTheoPar <- 1L      # dummy value; only used in Theory mode

#
ncpus <- min(detectCores(), 2)
sfInit(parallel = TRUE, cpus = ncpus)
sfExport("iniYear","nsim","niter","srb","minAge","maxAge",
         "alphaPar","gammaPar","rhoPar","phiPar","sigmaPar","betaPar",
         "cohortLen","periodStudy","cohortNames","startPop","inputData",
         "sonPrefMx","sonPrefTheory","model","country",
         "sensitivity","fertMod","sonTheoPar")
#sfSource("Code/FunctionsSRB.R")
sfSource("FunctionsSRB.R")

Start <- Sys.time()
out <- sfLapply(1:ncpus, run.simulations)
End <- Sys.time(); print(End-Start)
sfStop()

# dir.create("SimulationResults", showWarnings = FALSE)
fileName <- "Korea_SMOKETEST"
save(out, file = file.path("SimulationResults", paste0(fileName,".RData")))
cat("Saved:", file.path("SimulationResults", paste0(fileName,".RData")), "\n")


## Visual Check for Testrun ##
# 1) Load your smoke-test results
load("SimulationResults/Korea_SMOKETEST.RData")

# 2) Build SRB and TFR matrices from 'out' (works for any niter)
#    SRB in boys per 100 girls = p/(1-p) * 100, where p = propBoys
srb_runs <- lapply(out, function(z) t((z$propBoys / (1 - z$propBoys)) * 100))
srb_mat  <- do.call(rbind, srb_runs)

# Years are the rownames of propBoys
years <- as.integer(rownames(out[[1]]$propBoys))

# TFR is either in $TFR (non-sensitivity) or embedded in fertilityData
tfr_runs <- lapply(out, function(z) {
  if (!is.null(z$TFR)) {
    t(z$TFR)
  } else {
    # fallback if needed
    t(z$fertilityData[, ncol(z$fertilityData), , drop = FALSE])
  }
})
tfr_mat <- do.call(rbind, tfr_runs)

# 3) Smooth and compute means/quantiles
ma <- function(x, n=5) stats::filter(x, rep(1/n, n), sides = 2)

srb_ma <- t(apply(srb_mat, 1, ma))
tfr_ma <- t(apply(tfr_mat, 1, ma))

srb_mean <- apply(srb_ma, 2, mean, na.rm = TRUE)
srb_q10  <- apply(srb_ma, 2, quantile, 0.10, na.rm = TRUE)
srb_q90  <- apply(srb_ma, 2, quantile, 0.90, na.rm = TRUE)

tfr_mean <- apply(tfr_ma, 2, mean, na.rm = TRUE)
tfr_q10  <- apply(tfr_ma, 2, quantile, 0.10, na.rm = TRUE)
tfr_q90  <- apply(tfr_ma, 2, quantile, 0.90, na.rm = TRUE)

# 4) Plots
library(RColorBrewer)
cols <- brewer.pal(8, "Greys")
ni <- c(1, 2, length(years)-1, length(years))  # trim MA edge NAs if present
ni <- ni[ni >= 1 & ni <= length(years)]        # guard for very short series

# SRB
plot(years, srb_mean, type = "n", ylim = c(104, 118),
     xlab = "Year", ylab = "SRB (boys per 100 girls)", frame.plot = FALSE)
grid(lwd = 2)
if (length(years) > 4) {
  polygon(c(years[-ni], rev(years[-ni])),
          c(srb_q10[-ni], rev(srb_q90[-ni])),
          col = cols[4], border = NA)
}
lines(years, srb_mean, lwd = 3, col = cols[7])
title("South Korea — Simulated SRB (smoke test)")

# TFR
plot(years, tfr_mean, type = "n",
     ylim = range(tfr_q10, tfr_q90, na.rm = TRUE),
     xlab = "Year", ylab = "TFR", frame.plot = FALSE)
grid(lwd = 2)
if (length(years) > 4) {
  polygon(c(years[-ni], rev(years[-ni])),
          c(tfr_q10[-ni], rev(tfr_q90[-ni])),
          col = cols[4], border = NA)
}
lines(years, tfr_mean, lwd = 3, col = cols[7])
title("South Korea — Simulated TFR (smoke test)")

### End of test run ############################################################ 
}



###########################################################
### 2)  INITIAL PARAMETERS                              ###  
###########################################################

# SENSITIVITY: TRUE to carry out sensitivity analysis
sensitivity <- FALSE
# sensitivity <- TRUE

# Parameters for the sensitivity analysis
if (sensitivity) {
  
  # In order to carry out a sensitivity analysis, the user may choose one of the 4 tables below with sets of parameters.
  # Alternatively, the user may generate its own table with different parameter combinations.
  
  lhsPar <- read.table("Sensitivity/beta_sigma_sensitivity.txt", header = T)
  # lhsPar <- read.table("Sensitivity/rho_phi_sensitivity.txt", header = T)
  # lhsPar <- read.table("Sensitivity/lhs_44x5sims_244set.txt", header = T)
  # lhsPar <- read.table("Sensitivity/lhs_180x1sim_244set.txt", header = T)
 
  lhsPar <- lhsPar[, order(names(lhsPar))]
}

# CHOOSE A MODEL (Theory or empirical)
model <- "Empirical"
# model <- "Theory"

# GLOBAL INITIAL PARAMETERS
# Number of iterations
if (model == "Empirical") niter <- 70 else niter <- 50
# Initial year of the simulation
iniYear <- 1980
# Initial population
iniPop <- 100000
# Starting year for population initialization
startYear <- 1945
# Number of simulations
nsim <- 2
# Sex Ratio at Birth 
srb <- 0.5122
# Min age reproduction
minAge <- 15
# Max age reproduction
maxAge <- 50

# MODEL PARAMETERS: 
# Length of the vector parameter !!! WHEN MODEL = "Theory" it only considers the first value of each vector
n <- 1    
# parity scaling parameter
alphaPar <- rep(0.075, n)
# Birth risk proportional expansion parameter without meeting son preference
gammaPar <- rep(0.2, n)
# Techonolgy difussion parameters
rhoPar <- rep(0.5, n)
phiPar <- rep(7, n)
# Abortion probabilities scaling parameter
sigmaPar <- rep(1.7, n)
betaPar <- rep(.2, n)

# Son preferences for the THEORETICAL scenario: 1 = 0%, 2 = 10%, 3 = 20%,... 11 = 100%)
sonTheoPar <- c(1, 6, 11)
# sonTheoPar <- 3

# Fertility declines to be tested on the THEORETICAL model: 
#   1 = normal decline, 
#   2 = 25% faster, 3 = 66.7% faster
#   4 = 28.5% slower, 5 = 50% slower
# fertMod <- 1:3 (The first value is used by DEFAULT)
fertMod <- 1

# Period of study
periodStudy <- seq(iniYear, iniYear + niter - 1, 1)
# Cohorts considered (5-year cohort)
cohortLen <- 5
cohortNames <- seq(iniYear-maxAge, iniYear+trunc(niter/cohortLen,0)*cohortLen, cohortLen)

# Source code 
source("Code/FunctionsSRB.R")

###########################################################
### 3)  DATA AND INITIAILIZATION OF THE MODEL           ###  
###########################################################

# COUNTRY: Data from 2 countries are available to run the model: South Korea and India
country <- "Korea"
# country <- "India"

# INPUT DATA: Female Population, Fertility rates, and death rates
inputData <- get.data(country)

# SON PREFERENCE DATA
dataSonPref <- funcSonPref(country)
sonPrefMx <- dataSonPref$sonPrefMx
sonPrefTheory <- dataSonPref$sonPrefTheory

# INITIALIZE THE MODEL: STARTING POPULATION
startPop <- initialization(iniPop, startYear)


###########################################################
### 4)  RUN THE MODEL WITH PARALLEL COMPUTING           ###  
###########################################################

# Number of CPUS for parallel computing depending on the test
# 1) Empirical model 
if (model == "Empirical") {
  fertMod <- fertMod[1]
  # When only ONE value per parameter, split the simulaitons among cores
  if (length(alphaPar) == 1 & !sensitivity) {
    ncpus <- min(detectCores(), 24)
    nsim <- max(1, trunc(nsim / ncpus))
    alphaPar <- rep(alphaPar, ncpus)
    gammaPar <- rep(gammaPar, ncpus)
    rhoPar <- rep(rhoPar, ncpus)
    phiPar <- rep(phiPar, ncpus)
    sigmaPar <- rep(sigmaPar, ncpus)
    betaPar <- rep(betaPar, ncpus)
  } else ncpus <- min(detectCores(), length(alphaPar)) 
} 

# 2) Theoretical model with endogenous abortion probabilities
if (model == "Theory") {
  # Vary across different son preference values
  if (length(sonTheoPar) > 1) {
    ncpus <- min(detectCores(), length(sonTheoPar))
    fertMod <- fertMod[1]
  } else ncpus <- min(detectCores(), length(fertMod))
}

# SENSITIVITY ANALYSIS
if (sensitivity) {

  # Number of CPUS
  ncpus <- min(detectCores(), nrow(lhsPar))
  
  # Model parameters from the input table 'lhsPar'
  alphaPar <- lhsPar[, 1][1:ncpus]
  betaPar <- lhsPar[, 2][1:ncpus]
  gammaPar <- lhsPar[, 3][1:ncpus]
  phiPar <- lhsPar[, 4][1:ncpus]
  rhoPar <- lhsPar[, 5][1:ncpus]
  sigmaPar <- lhsPar[, 6][1:ncpus]
  if (model == "Theory") sonTheoPar <- lhsPar[, 7][1:ncpus]
}
  
# Start cpus
sfInit(parallel = TRUE, cpus = ncpus)

# Export variables to cpus
sfExport("iniYear", "nsim", "niter", "srb", "minAge", "maxAge", "alphaPar", 
         "gammaPar", "rhoPar", "phiPar", "sigmaPar", "betaPar", "cohortLen", 
         "periodStudy", "cohortNames", "startPop", "inputData", "sonPrefMx",
         "sonPrefTheory", "model", "sonTheoPar", "fertMod", "sensitivity")

# Run common prep functions on all cpus:
sfSource("Code/FunctionsSRB.R")

# Run process in cpus for "niter" iterations:
Start <- Sys.time()
out <- sfLapply(c(1:ncpus), run.simulations)
End <- Sys.time()
print(End - Start)

# stop parallel processes
sfStop()


###########################################################
### 5)  RESULTS - SAVE OUTPUT                           ###  
###########################################################

# This code gives the option to save the output of the simulations in a RData file

# Name of the RData file (extension no needed)
fileName <- "NameOfTheOuputFile"

# Save Output of the Simulations in an RData file
#save(out, file = paste(fileName, ".RData", sep = ""))

save(out, file = file.path("SimulationResults", paste0(fileName,".RData")))
cat("Saved:", file.path("SimulationResults", paste0(fileName,".RData")), "\n")


###########################################################
### 6) PLOTS RESULTS                                    ###  
###########################################################

# Here the code to replicate Figures 4, 6 and 7 from the Demography paper is provided.
# The RData files with the results to plot the figures are stored in the SimulationResults folder

# Packages
# install.packages("RColorBrewer")
library(RColorBrewer)


#----------#
# FIFGURE 4: Empirical, with sex-selective abortion simulation data
#----------#

# DATA
korea_emp_ssa <- extractdata_all("Korea_Empirical_newPar_240sim_1M", "empirical")

# TFR
korea_tfr_ssa <- korea_emp_ssa[,c(77:146)]
korea_tfr_ssa <- apply(korea_tfr_ssa, 1, ma)
korea_tfr_ssa_mean <- apply(korea_tfr_ssa, 1, mean, na.rm= TRUE)
korea_tfr_ssa_q10 <- apply(korea_tfr_ssa, 1, quantile, 0.025, na.rm= TRUE)
korea_tfr_ssa_q90 <- apply(korea_tfr_ssa, 1, quantile, 0.975, na.rm= TRUE)

# SRB
korea_srb_ssa <- korea_emp_ssa[,c(7:76)]
korea_srb_ssa <- apply(korea_srb_ssa, 1, ma)
korea_srb_ssa_mean <- apply(korea_srb_ssa, 1, mean, na.rm= TRUE)
korea_srb_ssa_q10 <- apply(korea_srb_ssa, 1, quantile, 0.025, na.rm= TRUE)
korea_srb_ssa_q90 <- apply(korea_srb_ssa, 1, quantile, 0.975, na.rm= TRUE)


# PLOT

niter <- 70
x <- seq(1980, 1980+niter-1, 1)

srb.colours <- brewer.pal(n = 8, name = "Greys")

def.par <- par(no.readonly = TRUE)
par(def.par)
#layout.show()
layout(matrix(c(1,2), nrow = 2, ncol =1, byrow = TRUE), 
       widths = c(2, 5.5), heights = c(2, 7))
par(mar=c(0, 0, 0, 0)+.1)
plot.new()

par(mar=c(0, 0, 0, 0)+.1)
legend(x = .1, y = 0.7, title = expression(bold("Total Fertility Rate")), col = c(srb.colours[6], 
                                                                                  srb.colours[3]), 
       lty = c(1, 1, 3), lwd = c(3,3, 3), y.intersp = 1,
       legend = c(expression(paste("Simulated with sex-selective abortion")),  
                  
                  "UN estimates"), cex = 0.6, bty = "n")
legend(x =0.6, y = 0.7, title = expression(bold("Sex Ratio at Birth")), col = c(srb.colours[6], srb.colours[3]), 
       pch=c(1,1,1), lwd = c(3, 3,3), y.intersp = 1,
       legend = c(expression(paste("Simulated with sex-selective abortion")),
                  "UN estimates"), cex = 0.6, bty = "n")

par(mar=c(4,4,1.2,4)+.1)
plot(x, korea_srb_ssa_mean, typ = "n", col=srb.colours[5], lwd = 2,
     xlab="Years", xlim = c(1980, 2050), yaxt = 'n', ylab = '', ylim = c(104, 116),
     font.lab = 2, main = "South Korea", frame.plot = FALSE)
grid(lwd = 2)
noInclude <- c(1, 2, length(x)-1, length(x))
polygon(x = c(x[-noInclude], rev(x[-noInclude])), 
        y = c(korea_srb_ssa_q10[-noInclude], rev(korea_srb_ssa_q90[-noInclude])), col=srb.colours[4], border=NA)
lines(x=x, y=korea_srb_ssa_mean, typ = "b", lwd =4, col = srb.colours[6])
lines(x[1:33], y = ma(srb.approx.plot.sk)[1:33], col = srb.colours[3], typ = 'b', lwd = 3)

axis(4)
mtext("Sex Ratio at Birth", side=4, line=3, font = 2)

par(new=TRUE)
plot(x, korea_tfr_ssa_mean, type="n", col="red", lwd = 2,
     xlab="Years", ylab="Total Fertility Rate", xlim = c(1980, 2050),
     font.lab = 2, main = "South Korea", frame.plot = FALSE)
lines(x=x, y=korea_tfr_ssa_mean, lty = 1, lwd =3, col = srb.colours[6])

lines(x[1:33], ma(tfr.approx.plot.sk)[1:33], col = srb.colours[4], lty = 3, lwd = 4)
lines(x[33:length(tfr.approx.plot.sk)], ma(tfr.approx.plot.sk)[33:length(tfr.approx.plot.sk)], col = srb.colours[3], lty = 3, lwd = 3)


#----------#
# FIGURE 6 plot: SRB at different theoretical son preference values
#----------#

# load data 
dat1 <- extractdata_all("Korea_Theory_sonPref0-10_500K", "theoretical")
y_ssa <- korea_srb_ssa_mean[1:50]
y_ssa_lower <- korea_srb_ssa_q10[1:49]
y_ssa_upper <- korea_srb_ssa_q90[1:49]

# Son preference percentage: 1) 0%, 2) 10%, ... , 11) 100%
# We are interested in 0, 50, and 100% here.
sonPref <- c(1, 6, 11)
nsim <- nrow(dat1) / 11

# SRB of each son preference scenario
SRB <- matrix(NA, nrow = length(sonPref), ncol = 50)
SRB_upper <- SRB_lower <- matrix(NA, nrow = length(sonPref), ncol = 50)

#Son Preference at 0% 
SRB_b <- dat1[1:nsim+nsim*(sonPref[1]-1), 8:57]

SRB_b <- t(apply(SRB_b, 1, ma))
SRB[1, ] <- apply(SRB_b, 2, mean, na.rm = TRUE)
SRB_upper[1, ] <- apply(SRB_b, 2, quantile, 0.1, na.rm = TRUE)
SRB_lower[1, ] <- apply(SRB_b, 2, quantile, 0.9, na.rm = TRUE)

#Son preference at 50% (1980 levels)
SRB_b <- dat1[1:nsim+nsim*(sonPref[2]-1), 8:57]
SRB_b <- t(apply(SRB_b, 1, ma))
SRB[2, ] <- apply(SRB_b, 2, mean, na.rm = TRUE)
SRB_upper[2, ] <- apply(SRB_b, 2, quantile, 0.1, na.rm = TRUE)
SRB_lower[2, ] <- apply(SRB_b, 2, quantile, 0.9, na.rm = TRUE)

#Son preference at 100% 
SRB_b <- dat1[1:nsim+nsim*(sonPref[3]-1), 8:57]
SRB_b <- t(apply(SRB_b, 1, ma))
SRB[3, ] <- apply(SRB_b, 2, mean, na.rm = TRUE)
SRB_upper[3, ] <- apply(SRB_b, 2, quantile, 0.1, na.rm = TRUE)
SRB_lower[3, ] <- apply(SRB_b, 2, quantile, 0.9, na.rm = TRUE)

# PLOT
years <- seq(1980, 2029, 1) 
# setEPS()
# postscript("Fig6.eps", height=8, width=10)
par(mfrow = c(1,1))
plot(years, SRB[1, ], xlab = 'time (years)', ylab = 'SRB',
     ylim = c(100, max(SRB, na.rm = TRUE)),
     type = 'n', frame.plot = FALSE)
grid(lwd = 2)

noInclude <- c(1, 2, length(years)-1, length(years))
polygon(x = c(years[-noInclude], rev(years[-noInclude])), 
        y = c(SRB_lower[2, -noInclude], rev(SRB_upper[2, -noInclude])), col="grey", border=NA)

lines(years, SRB[2, ], lty = 3, col = "black", lwd = 2)

noInclude <- c(1, 2, length(years)-1, length(years))
polygon(x = c(years[-noInclude], rev(years[-noInclude])), 
        y = c(SRB_lower[3, -noInclude], rev(SRB_upper[3, -noInclude])), col="grey", border=NA)

lines(years, SRB[3, ], lty = 3, col = "black", lwd = 2)


noInclude <- c(1, 2, length(years)-1, length(years))
polygon(x = c(years[-noInclude], rev(years[-noInclude])), 
        y = c(y_ssa_lower[-noInclude], rev(y_ssa_upper[-noInclude])), col="grey", border=NA)

lines(years[1:48], y_ssa[1:48], lty = 4, col = "black", lwd = 2)

legend(x =2000, y = 140, col = "black", 
       lty = c(3,2, 4), y.intersp = 1, lwd = 2,
       legend = c("Son Preference Constant at 100%",
                  "Son Preference Constant at 1980 levels",
                  "Declining Son Preference"), cex = 0.7, bty = "n")
#dev.off()

#----------#
# FIGURE 7: Different ability and readiness assumptions
#----------#

#Load data
dat1 <- extractdata_all("Korea_Empirical_technology05_500K", "empirical")
dat2 <- extractdata_all("Korea_Empirical_beta0_500K", "empirical")

korea_tech05 <- dat1[,7:57]
korea_tech05 <- t(apply(korea_tech05, 1, ma))
korea_tech05_mean <- apply(korea_tech05, 2, mean, na.rm= TRUE)
korea_tech05_lower <- apply(korea_tech05, 2, quantile, 0.1, na.rm= TRUE)
korea_tech05_upper <- apply(korea_tech05, 2, quantile, 0.9, na.rm= TRUE)


korea_beta0 <- dat2[,7:57]
korea_beta0 <- t(apply(korea_beta0, 1, ma))
korea_beta0_mean <- apply(korea_beta0, 2, mean, na.rm = TRUE)
korea_beta0_lower <- apply(korea_beta0, 2, quantile, 0.1, na.rm = TRUE)
korea_beta0_upper <- apply(korea_beta0, 2, quantile, 0.9, na.rm = TRUE)

# PLOT
years <- seq(1980, 2030, 1) 
#setEPS()
#postscript("Fig7.eps", height=8, width=10)
plot(years, korea_beta0_mean, xlab = 'time (years)', ylab = 'SRB',
     ylim = c(105, 120),
     #ylim = c(100, max(SRB, na.rm = TRUE)),
     type = 'n', frame.plot = FALSE)
grid(lwd = 2)
noInclude <- c(1, 2, length(years)-1, length(years))
polygon(x = c(years[-noInclude], rev(years[-noInclude])), 
        y = c(korea_beta0_lower[-noInclude], rev(korea_beta0_upper[-noInclude])), col="grey", border=NA)

lines(years, korea_beta0_mean, lty = 1, col = "black", lwd = 2)

noInclude <- c(1, 2, length(years)-1, length(years))
polygon(x = c(years[-noInclude], rev(years[-noInclude])), 
        y = c(korea_tech05_lower[-noInclude], rev(korea_tech05_upper[-noInclude])), col="grey", border=NA)
lines(years, korea_tech05_mean, lty = 3, col = "black", lwd = 2)

legend(x =2000, y = 120, col = "black", 
       lty = c(3,1), y.intersp = 1, lwd = 2,
       legend = c(expression(paste("Technology Constant at 50% (", rho, "=0, ", phi, "=0)")),
                  expression(paste("Lower readiness (", sigma, "=1, ", beta, "=0)"))),
       cex = 0.7, bty = "n")
#dev.off()

#####################################################################
############### the case of India 1 ###################################
#####################################################################
# Clear workspace
rm(list = ls())


library(snowfall)
library(parallel)



# --- core model toggles ---
sensitivity <- FALSE
model <- "Empirical"
country <- "India"

# --- simulation settings (small–moderate load first) ---
iniYear   <- 1980
niter     <- 60       # 1980–2039
iniPop    <- 100000   # starting synthetic population
startYear <- 1945
nsim      <- 10       # number of stochastic replications per CPU
srb       <- 0.5122   # baseline SRB (≈105 boys per 100 girls)

# --- global constants ---
minAge <- 15
maxAge <- 50
cohortLen <- 5
periodStudy <- seq(iniYear, iniYear + niter - 1, 1)
cohortNames <- seq(iniYear - maxAge,
                   iniYear + trunc(niter / cohortLen, 0) * cohortLen,
                   cohortLen)

# --- fertility schedule selector ---
fertMod <- 1    # use “fertility5” (UN standard decline)

# --- theoretical dummy (not used here but required) ---
sonTheoPar <- 1L

# --- source function definitions ---
source("Code/FunctionsSRB.R")


# 2) Parameters -------------------------------------------

n <- 1
alphaPar <- rep(0.075, n)   # parity scaling
gammaPar <- rep(0.35,  n)   # extra fertility if no son
rhoPar   <- rep(0.30,  n)   # technology diffusion slope
phiPar   <- rep(11,    n)   # technology diffusion midpoint
sigmaPar <- rep(1.0,   n)   # abortion probability scaling
betaPar  <- rep(0.10,  n)   # readiness scaling

# 3) Load data & initialize -------------------------------
inputData     <- get.data(country)
spData        <- funcSonPref(country)
sonPrefMx     <- spData$sonPrefMx
sonPrefTheory <- spData$sonPrefTheory
startPop      <- initialization(iniPop, startYear)

# 4) Parallel computing setup -----------------------------
ncpus <- min(detectCores(), 4)  # adjust if you want more
sfInit(parallel = TRUE, cpus = ncpus)
sfExport("iniYear","nsim","niter","srb","minAge","maxAge",
         "alphaPar","gammaPar","rhoPar","phiPar","sigmaPar","betaPar",
         "cohortLen","periodStudy","cohortNames","startPop","inputData",
         "sonPrefMx","sonPrefTheory","model","country",
         "sensitivity","fertMod","sonTheoPar")
sfSource("Code/FunctionsSRB.R")

# 5) Run simulations --------------------------------------
Start <- Sys.time()
out <- sfLapply(1:ncpus, run.simulations)
End <- Sys.time(); print(End - Start)
sfStop()

# 6) Save results -----------------------------------------
dir.create("SimulationResults", showWarnings = FALSE)
fileName <- "India_Empirical_newPar_nsim10"
save_path <- file.path("SimulationResults", paste0(fileName, ".RData"))
save(out, file = save_path)
cat("Saved simulation output to:", save_path, "\n")

# 7) Quick check ------------------------------------------
load(save_path)
str(out, max.level = 1)

# 8) Plots 
library(RColorBrewer)
dat <- extractdata_all("India_Empirical_newPar_nsim10", "empirical")

years <- 1980:(1980 + 59)
srb_mat <- as.matrix(dat[, paste0("srb", years)])
srb_ma  <- t(apply(srb_mat, 1, ma))
srb_mean <- apply(srb_ma, 2, mean, na.rm = TRUE)
srb_q10  <- apply(srb_ma, 2, quantile, 0.10, na.rm = TRUE)
srb_q90  <- apply(srb_ma, 2, quantile, 0.90, na.rm = TRUE)
cols <- brewer.pal(8, "Greys")
ni <- c(1,2,length(years)-1,length(years))

plot(years, srb_mean, type = "n", ylim = c(104,118),
     xlab = "Year", ylab = "SRB (boys per 100 girls)", frame.plot = FALSE)
grid(lwd = 2)
polygon(c(years[-ni], rev(years[-ni])),
        c(srb_q10[-ni], rev(srb_q90[-ni])), col = cols[4], border = NA)
lines(years, srb_mean, lwd = 3, col = cols[7])
title("India — Simulated SRB (Empirical model)")

#####################################################################
############### the case of India 2 ###################################
#####################################################################

#Multiple set of scenarios 

model   <- "Empirical"
country <- "India"
niter   <- 60
iniPop  <- 100000
nsim    <- 8
# India-appropriate scenarios (start with 3–4)
alphaPar <- c(0.075,0.075,0.075,0.075)
gammaPar <- c(0.30,0.35,0.40,0.30)
rhoPar   <- c(0.30,0.30,0.25,0.35)
phiPar   <- c(11,12,12,10)
sigmaPar <- c(1.0,0.9,1.2,1.1)
betaPar  <- c(0.10,0.10,0.08,0.12)

ncpus <- min(detectCores(), length(alphaPar))
nsim  <- max(1, floor(nsim / ncpus))  # split simulations across cores

# re-run the same pipeline and save:
fileName <- "India_Scenarios_K4"


#####################################################################
############### the case of Vietnam 1 ###################################
#####################################################################
# Clear workspace
rm(list = ls())

library(snowfall)
library(parallel)

# --- core model toggles ---
model   <- "Empirical"
country <- "Vietnam"
niter   <- 60
iniPop  <- 100000
nsim    <- 8

alphaPar <- 0.075
gammaPar <- 0.27
rhoPar   <- 0.45
phiPar   <- 17
sigmaPar <- 1.4
betaPar  <- 0.17

ncpus <- min(detectCores(), 4)
# re-run pipeline, save:
fileName <- "Vietnam_SeedRun"

# --- simulation settings (small–moderate load first) ---
iniYear   <- 1980
niter     <- 60       # 1980–2039
iniPop    <- 100000   # starting synthetic population
startYear <- 1945
nsim      <- 10       # number of stochastic replications per CPU
srb       <- 0.5122   # baseline SRB (≈105 boys per 100 girls)

# --- global constants ---
minAge <- 15
maxAge <- 50
cohortLen <- 5
periodStudy <- seq(iniYear, iniYear + niter - 1, 1)
cohortNames <- seq(iniYear - maxAge,
                   iniYear + trunc(niter / cohortLen, 0) * cohortLen,
                   cohortLen)

# --- fertility schedule selector ---
fertMod <- 1    # use “fertility5” (UN standard decline)

# --- theoretical dummy (not used here but required) ---
sonTheoPar <- 1L

# --- source function definitions ---
source("Code/FunctionsSRB.R")

# 3) Load data & initialize -------------------------------
inputData     <- get.data(country)
spData        <- funcSonPref(country)
sonPrefMx     <- spData$sonPrefMx
sonPrefTheory <- spData$sonPrefTheory
startPop      <- initialization(iniPop, startYear)

# 4) Parallel computing setup -----------------------------
ncpus <- min(detectCores(), 4)  # adjust if you want more
sfInit(parallel = TRUE, cpus = ncpus)
sfExport("iniYear","nsim","niter","srb","minAge","maxAge",
         "alphaPar","gammaPar","rhoPar","phiPar","sigmaPar","betaPar",
         "cohortLen","periodStudy","cohortNames","startPop","inputData",
         "sonPrefMx","sonPrefTheory","model","country",
         "sensitivity","fertMod","sonTheoPar")
sfSource("Code/FunctionsSRB.R")

# 5) Run simulations --------------------------------------
Start <- Sys.time()
out <- sfLapply(1:ncpus, run.simulations)
End <- Sys.time(); print(End - Start)
sfStop()

# 6) Save results -----------------------------------------
dir.create("SimulationResults", showWarnings = FALSE)
fileName <- "India_Empirical_newPar_nsim10"
save_path <- file.path("SimulationResults", paste0(fileName, ".RData"))
save(out, file = save_path)
cat("Saved simulation output to:", save_path, "\n")

# 7) Quick check ------------------------------------------
load(save_path)
str(out, max.level = 1)

# 8) Plots 
library(RColorBrewer)
dat <- extractdata_all("Vietnam_Empirical_newPar_nsim10", "empirical")

years <- 1980:(1980 + 59)
srb_mat <- as.matrix(dat[, paste0("srb", years)])
srb_ma  <- t(apply(srb_mat, 1, ma))
srb_mean <- apply(srb_ma, 2, mean, na.rm = TRUE)
srb_q10  <- apply(srb_ma, 2, quantile, 0.10, na.rm = TRUE)
srb_q90  <- apply(srb_ma, 2, quantile, 0.90, na.rm = TRUE)
cols <- brewer.pal(8, "Greys")
ni <- c(1,2,length(years)-1,length(years))

plot(years, srb_mean, type = "n", ylim = c(104,118),
     xlab = "Year", ylab = "SRB (boys per 100 girls)", frame.plot = FALSE)
grid(lwd = 2)
polygon(c(years[-ni], rev(years[-ni])),
        c(srb_q10[-ni], rev(srb_q90[-ni])), col = cols[4], border = NA)
lines(years, srb_mean, lwd = 3, col = cols[7])
title("Vietnam — Simulated SRB (Empirical model)")
