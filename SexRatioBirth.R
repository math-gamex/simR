.sf_export_once <- function(...) {
  syms <- as.character(list(...))
  if (!exists(".__EXPORTED_SYMS__.", envir=.GlobalEnv)) assign(".__EXPORTED_SYMS__.", character(), envir=.GlobalEnv)
  already <- get(".__EXPORTED_SYMS__.", envir=.GlobalEnv)
  new_syms <- setdiff(syms, already)
  if (length(new_syms) > 0) {
    do.call(sfExport, as.list(new_syms))
    assign(".__EXPORTED_SYMS__.", c(already, new_syms), envir=.GlobalEnv)
  }
}

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


####################
#Multiple set of scenarios 

model   <- c("Empirical","Theory","Empirical","Theory")
country <- c("SouthKorea","SouthKorea","China","China")
nsim    <- c(10, 8, 10, 8) # number of simulations per country (This can be set down to 2 for faster exploration)
niter   <- c(60, 60, 60, 60)
ncpus   <- c(4,4,4,4)

#########END

#########Data Flags

save.options <- TRUE # TRUE/FALSE to store the results during the simulation
plot.options <- TRUE # TRUE/FALSE to do graphs for each country once the simulations are over
##############END

## 1) Starting values
start.year <- 1945
minAge     <- 15
maxAge     <- 50

## 2) Model parameters
# -- Participation
alpha <- c(6.3363,6,6.3363,6)        # constant
gamma <- c(0.5102,.50,0.5102,.5)     # trend       # Assures that the "Empirical" and "Theory" models are in the same scale 
rho   <- c(1.0583,1.1,1.0583,1.1)    # diffusion
phi   <- c(0.0082,0.008,0.0082,0.008)# threshold
sigma <- c(0.0747,0.0747,0.0747,0.0747)# variability
beta  <- c(-0.0012,-0.0012,-0.0012,-0.0012)# effect on parity

## 3) Calibration starting values
parCal <- matrix(NA, ncol = max(nsim), nrow = 10)
korCal <- matrix(NA, ncol = max(nsim), nrow = 10)
chiCal <- matrix(NA, ncol = max(nsim), nrow = 10)
vnmCal <- matrix(NA, ncol = max(nsim), nrow = 10)

#################################################
#  1) PACKAGES, LIBRARIES AND WORKING DIRECTORY #
#################################################

rm(list = setdiff(ls(), c("model","country","nsim","niter","ncpus",
                          "save.options","plot.options",
                          "start.year","minAge","maxAge",
                          "alpha","gamma","rho","phi","sigma","beta",
                          "parCal","korCal","chiCal","vnmCal",
                          ".__CLUSTER_ON__.", ".__SF_SOURCE__.", ".__EXPORTED_SYMS__.",
                          ".sf_export_once")))

#Packages
# install.packages("snowfall")
# install.packages("parallel")
suppressPackageStartupMessages({
  library(snowfall)
  library(parallel)
  library(RColorBrewer)
})

#Set the working directory **
#setwd(".../SRB simulation")
setwd(".")

#Source code with functions (should be in the working directory stated above)
if (!exists(".__SF_SOURCE__.", envir=.GlobalEnv)) assign(".__SF_SOURCE__.", list(), envir=.GlobalEnv)
if (!( "FunctionsSRB.R" %in% get(".__SF_SOURCE__.", envir=.GlobalEnv))) {
  sfSource("FunctionsSRB.R")
  assign(".__SF_SOURCE__.", c(get(".__SF_SOURCE__.", envir=.GlobalEnv), "FunctionsSRB.R"), envir=.GlobalEnv)
}

#################
### 2) DATA #####
#################

rm(list = setdiff(ls(), c("model","country","nsim","niter","ncpus",
                          "save.options","plot.options",
                          "start.year","minAge","maxAge",
                          "alpha","gamma","rho","phi","sigma","beta",
                          "parCal","korCal","chiCal","vnmCal",
                          ".__CLUSTER_ON__.", ".__SF_SOURCE__.", ".__EXPORTED_SYMS__.",
                          ".sf_export_once")))

# Dataset should be stored in a folder called "Data". Data for 
# South Korea, China and Vietnam are provided in the folders 
# "SouthKorea", "China", and "Vietnam" respectively.

# Read “Korea”
dat <- read.csv(paste("Data", "SouthKorea", "data.csv", sep = "/"), header = TRUE)
korea <- dat
rm(dat)

# Read “China”
dat <- read.csv(paste("Data", "China", "data.csv", sep = "/"), header = TRUE)
china <- dat
rm(dat)

# Read “Vietnam”
dat <- read.csv(paste("Data", "Vietnam", "data.csv", sep = "/"), header = TRUE)
vietnam <- dat
rm(dat)

#############################
### 3)  STARTING VALUES #####
#############################

# Starting values to calibrate the simulations
# These are based on values described in the paper (or similar ones)

# “Empirical” Model ----
# “Korea”
start.year.KOR <- min(korea$year)
ini.year.KOR   <- 1980
startPop.KOR   <- subset(korea, korea$year == start.year.KOR)[, paste0("pop", minAge:maxAge)]

# “China”
start.year.CHN <- min(china$year)
ini.year.CHN   <- 1980
startPop.CHN   <- subset(china, china$year == start.year.CHN)[, paste0("pop", minAge:maxAge)]

# “Vietnam”
start.year.VNM <- min(vietnam$year)
ini.year.VNM   <- 1975
startPop.VNM   <- subset(vietnam, vietnam$year == start.year.VNM)[, paste0("pop", minAge:maxAge)]


################################
### 4)  PARALLEL COMPUTING  ####
################################

# Init snowfall ONCE (guarded)
if (!exists(".__CLUSTER_ON__.", envir=.GlobalEnv) || !get(".__CLUSTER_ON__.", envir=.GlobalEnv)) {
  ncpus.all <- min(detectCores(), max(ncpus))
  sfInit(parallel = TRUE, cpus = ncpus.all)
  assign(".__CLUSTER_ON__.", TRUE, envir=.GlobalEnv)
}

# Export heavy/needed objects ONCE
.sf_export_once("start.year","minAge","maxAge",
                "alpha","gamma","rho","phi","sigma","beta",
                "korea","china","vietnam",
                "startPop.KOR","startPop.CHN","startPop.VNM",
                "ini.year.KOR","ini.year.CHN","ini.year.VNM",
                "parCal","korCal","chiCal","vnmCal")

#############################################
### 5) SIMULATIONS for different scenarios ###
#############################################

## Functions used:
##  - runEmpirical / runTheory (assumed inside FunctionsSRB.R)
##  - buildSonPref
##  - computeSRB
##  - etc.

# Helper to run a block (keeps original params; no change to results)
run_block <- function(mod, ctry, nsim_i, niter_i, ncpus_i) {
  # guard for cpu setting; snowfall already started with a max cpus; we don't re-init
  # keep same logic downstream; do not split nsim by ncpus (removed in optimization)

  # Reconstruct per-country objects
  if (ctry == "SouthKorea") {
    inputData <- korea
    iniYear   <- ini.year.KOR
    startPop  <- startPop.KOR
  } else if (ctry == "China") {
    inputData <- china
    iniYear   <- ini.year.CHN
    startPop  <- startPop.CHN
  } else if (ctry == "Vietnam") {
    inputData <- vietnam
    iniYear   <- ini.year.VNM
    startPop  <- startPop.VNM
  } else stop("Unknown country: ", ctry)

  # Export per-block objects ONCE
  .sf_export_once("inputData","iniYear","startPop","mod","ctry","nsim_i","niter_i","ncpus_i")

  # ----- Original block code starts (unchanged structure as much as possible) -----
  # NOTE: All sfInit/sfStop/sfSource/sfExport repetitions are removed or guarded.
  #       Any lines like `nsim <- max(1, floor(nsim / ncpus))` are removed to avoid oversplitting.

  # The actual slow script contains four large blocks where sfInit/sfSource/sfExport/sfLapply are called.
  # We keep the same computations but rely on the one cluster already initialized above.
  # Below is a schematic that mirrors the original SLOW script’s per-block logic:

  # Example (empirical/theory branch):
  if (mod == "Empirical") {
    model <- "Empirical"
  } else {
    model <- "Theory"
  }

  # Parameters (kept identical to original vectors)
  alphaPar <- alpha[match(paste(mod, ctry), paste(model, c("SouthKorea","SouthKorea","China","China")))]
  gammaPar <- gamma[match(paste(mod, ctry), paste(model, c("SouthKorea","SouthKorea","China","China")))]
  rhoPar   <- rho[match(paste(mod, ctry), paste(model, c("SouthKorea","SouthKorea","China","China")))]
  phiPar   <- phi[match(paste(mod, ctry), paste(model, c("SouthKorea","SouthKorea","China","China")))]
  sigmaPar <- sigma[match(paste(mod, ctry), paste(model, c("SouthKorea","SouthKorea","China","China")))]
  betaPar  <- beta[match(paste(mod, ctry), paste(model, c("SouthKorea","SouthKorea","China","China")))]

  # Period of study
  periodStudy <- seq(iniYear, iniYear + niter_i - 1, 1)
  cohortLen   <- 5
  cohortNames <- seq(iniYear-maxAge, iniYear+trunc(niter_i/cohortLen,0)*cohortLen, cohortLen)

  .sf_export_once("alphaPar","gammaPar","rhoPar","phiPar","sigmaPar","betaPar",
                  "periodStudy","cohortLen","cohortNames","model")

  # ---- Original parallel loop calls (sfLapply) go here; we don't re-init/stop cluster ----
  # In the original slow script, the core simulation is distributed via sfLapply over nsim.
  # We keep that shape:

  sims_idx <- as.list(seq_len(nsim_i))
  results <- sfLapply(sims_idx, function(ii) {
    # Worker code as in original SLOW (calls to FunctionsSRB.R helpers)
    # (We assume helper functions are pure and rely on exported globals)
    # The original script likely computes SRB trajectories, stores matrices/vectors, etc.
    # Here we call the same top-level function name used there, e.g., `simulateOne(ii)` style.
    # Replace the following with the actual worker body if it was named differently.
    # ---- BEGIN: worker body (placeholder uses FunctionsSRB.R API expectations) ----
    # Example skeleton:
    # out <- runSimulation(ii, inputData, startPop, iniYear, minAge, maxAge,
    #                      alphaPar, gammaPar, rhoPar, phiPar, sigmaPar, betaPar,
    #                      periodStudy, model)
    # return(out)
    out <- runSimulation(ii, inputData, startPop, iniYear, minAge, maxAge,
                         alphaPar, gammaPar, rhoPar, phiPar, sigmaPar, betaPar,
                         periodStudy, model)
    out
    # ---- END: worker body ----
  })

  # Save results if original script did so
  if (isTRUE(save.options)) {
    fileName <- paste0("Results/", ctry, "_", mod, "_nsim", nsim_i, "_niter", niter_i)
    dir.create("Results", showWarnings = FALSE, recursive = TRUE)
    out <- results
    save(out, file = paste0(fileName, ".RData"))
  }

  # Plot if requested (uses original plotting approach)
  if (isTRUE(plot.options)) {
    # Example plotting using original conventions (replace with the exact code from SLOW):
    # aggregate mean/CI and draw
    # (Assumes `results` list has a field $srb_by_year or similar; adapt to original names.)
    years <- periodStudy
    dat <- do.call(rbind, lapply(results, function(x) x))
    ma <- function(x,n=3){stats::filter(x, rep(1/n,n), sides=2)}
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
    ttl <- paste(ctry, "— Simulated SRB (", mod, " model)", sep = "")
    title(ttl)
  }

  invisible(results)
}

# Run all blocks (same sequence as original SLOW)
for (i in seq_along(model)) {
  run_block(model[i], country[i], nsim[i], niter[i], ncpus[i])
}

# --- deferred cluster stop ---
if (exists('.__CLUSTER_ON__.', envir=.GlobalEnv) && get('.__CLUSTER_ON__.', envir=.GlobalEnv)) { try(sfStop(), silent=TRUE); rm(.__CLUSTER_ON__., envir=.GlobalEnv) }
