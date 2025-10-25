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
source("FunctionsSRB.R")


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

sfSource("FunctionsSRB.R")

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