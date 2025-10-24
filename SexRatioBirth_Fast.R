#############################################################################
# Fast driver for SexRatioBirth_slow.R (same simulation, faster execution)
# - Initializes snowfall cluster ONCE
# - Exports heavy globals ONCE
# - Runs simulations in parallel by assigning a batch of nsim per worker
# - Concatenates worker outputs along the simulation dimension
#############################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(snowfall)
  library(parallel)
  library(RColorBrewer)
})

# -------------------------- User parameters (same defaults as slow) --------------------------
model   <- c("Empirical","Theory","Empirical","Theory")
country <- c("Korea","Korea","China","China")
nsim    <- c(10, 8, 10, 8)    # total simulations per scenario (keep as in slow script)
niter   <- c(60, 60, 60, 60)
ncpus   <- c(4,4,4,4)         # choose modest cpus to reduce export overhead

# Save / Plot options
save.options <- TRUE
plot.options <- TRUE

# Global constants (as in slow)
iniPop    <- 50000
startYear <- 1945
srb       <- 0.5122
minAge    <- 15
maxAge    <- 50

# Parameters (examples carried over from slow script header)
alphaPar <- 0.075; gammaPar <- 0.20
rhoPar   <- 0.50;  phiPar   <- 7
sigmaPar <- 1.70;  betaPar  <- 0.20

# Model choices
sensitivity <- FALSE
fertMod <- 1
sonTheoPar <- 1

# --------------------------------------------------------------------------------------------
# Utility: data getters (mirror slow expectations; adapt if you already have get.data)
# --------------------------------------------------------------------------------------------
get.data <- function(country) {
  dat <- read.csv(file.path("Data", country, "data.csv"), header = TRUE)
  # NOTE: If your original slow script builds a richer 'inputData' list,
  # replace this stub with that constructor to match FunctionsSRB expectations.
  return(list(raw = dat))
}

# Load core simulation functions (uses initialization(), funcSonPref(), run.simulations(), etc.)
source("FunctionsSRB.R")

# Moving average helper to match slow visuals
ma <- function(x, n=5) stats::filter(x, rep(1/n, n), sides = 2)

# -------------------------- Parallel bootstrap --------------------------
# Initialize snowfall ONCE
ncpus_all <- min(detectCores(), max(ncpus))
sfInit(parallel = TRUE, cpus = ncpus_all)
on.exit({ try(sfStop(), silent = TRUE) }, add = TRUE)

# Export invariant globals once
sfExport("iniPop","startYear","srb","minAge","maxAge",
         "alphaPar","gammaPar","rhoPar","phiPar","sigmaPar","betaPar",
         "sensitivity","fertMod","sonTheoPar","ma")

sfSource("FunctionsSRB.R")

# -------------------------- Worker wrapper --------------------------
# Keep FunctionsSRB::run.simulations() intact; set batch nsim in worker and run.
worker_run <- function(nsim_batch,
                       iniYear, periodStudy, cohortLen, cohortNames,
                       startPop, inputData, sonPrefMx, sonPrefTheory,
                       model, country) {
  assign("nsim", nsim_batch, envir = .GlobalEnv)
  assign("niter", length(periodStudy), envir = .GlobalEnv)
  assign("iniYear", iniYear, envir = .GlobalEnv)
  assign("cohortLen", cohortLen, envir = .GlobalEnv)
  assign("periodStudy", periodStudy, envir = .GlobalEnv)
  assign("cohortNames", cohortNames, envir = .GlobalEnv)
  assign("startPop", startPop, envir = .GlobalEnv)
  assign("inputData", inputData, envir = .GlobalEnv)
  assign("sonPrefMx", sonPrefMx, envir = .GlobalEnv)
  assign("sonPrefTheory", sonPrefTheory, envir = .GlobalEnv)
  assign("model", model, envir = .GlobalEnv)
  assign("country", country, envir = .GlobalEnv)
  run.simulations(1)
}
sfExport("worker_run")

# -------------------------- Combine utility --------------------------
concat_along_sim <- function(res_list, total_nsim) {
  take <- function(field) {
    xs <- lapply(res_list, function(x) x[[field]])
    if (is.null(xs[[1]])) return(NULL)
    if (length(dim(xs[[1]])) == 2) {
      out <- do.call(cbind, xs)                   # (niter x nsim)
      return(out[, seq_len(total_nsim), drop = FALSE])
    } else if (length(dim(xs[[1]])) == 3) {
      if (!requireNamespace("abind", quietly = TRUE)) stop("Please install.packages('abind')")
      out <- abind::abind(xs, along = 3)          # (niter x K x nsim)
      return(out[, , seq_len(total_nsim), drop = FALSE])
    } else {
      return(xs[[length(xs)]])
    }
  }
  last <- res_list[[length(res_list)]]
  list(
    propBoys    = take("propBoys"),
    TFR         = take("TFR"),
    abortData   = take("abortData"),
    alpha       = last$alpha, gamma = last$gamma, rho = last$rho, phi = last$phi,
    sigma       = last$sigma, beta = last$beta,
    sonPrefTheo = if (!is.null(last$sonPrefTheo)) last$sonPrefTheo else NA,
    nsim        = total_nsim,
    niter       = if (!is.null(last$niter)) last$niter else nrow(last$TFR),
    iniYear     = if (!is.null(last$iniYear)) last$iniYear else NA,
    endYear     = if (!is.null(last$endYear)) last$endYear else NA,
    iniPop      = last$iniPop, endPop = last$endPop
  )
}
sfExport("concat_along_sim")

# -------------------------- Scenario runner --------------------------
run_scenario <- function(model_i, country_i, nsim_i, niter_i, ncpus_i) {
  # 1) Inputs
  inputData     <- get.data(country_i)
  spData        <- funcSonPref(country_i)
  sonPrefMx     <- spData$sonPrefMx
  sonPrefTheory <- spData$sonPrefTheory
  startPop      <- initialization(iniPop, startYear)
  iniYear       <- if (country_i == "Vietnam") 1975 else 1980
  periodStudy   <- seq(iniYear, iniYear + niter_i - 1, 1)
  cohortLen     <- 5
  cohortNames   <- seq(iniYear - maxAge,
                       iniYear + trunc(niter_i / cohortLen, 0) * cohortLen,
                       cohortLen)

  # 2) Export per-scenario globals once
  sfExport("inputData","sonPrefMx","sonPrefTheory","startPop",
           "iniYear","periodStudy","cohortLen","cohortNames","model_i","country_i")

  # 3) Build batches (exactly nsim_i in total)
  cpus_use <- min(ncpus_all, ncpus_i)
  batch_sizes <- rep(floor(nsim_i / cpus_use), cpus_use)
  remainder <- nsim_i - sum(batch_sizes)
  if (remainder > 0) batch_sizes[seq_len(remainder)] <- batch_sizes[seq_len(remainder)] + 1

  # 4) Parallel run
  Start <- Sys.time()
  res_parts <- sfLapply(batch_sizes, function(nsim_batch) {
    worker_run(nsim_batch,
               iniYear, periodStudy, cohortLen, cohortNames,
               startPop, inputData, sonPrefMx, sonPrefTheory,
               model_i, country_i)
  })
  End <- Sys.time(); print(End - Start)

  # 5) Concatenate
  out <- concat_along_sim(res_parts, nsim_i)

  # 6) Save / Plot (same outputs/conventions as slow)
  if (isTRUE(save.options)) {
    dir.create("SimulationResults", showWarnings = FALSE, recursive = TRUE)
    fileName <- paste(country_i, model_i, "nsim", nsim_i, "niter", niter_i, sep="_")
    save(out, file = file.path("SimulationResults", paste0(fileName, ".RData")))
  }

  if (isTRUE(plot.options)) {
    years <- periodStudy
    dat <- t(out$propBoys) * 100  # assuming propBoys in [0,1]
    ma_dat <- t(apply(dat, 1, ma))
    srb_mean <- apply(ma_dat, 2, mean, na.rm = TRUE)
    srb_q10  <- apply(ma_dat, 2, quantile, 0.10, na.rm = TRUE)
    srb_q90  <- apply(ma_dat, 2, quantile, 0.90, na.rm = TRUE)
    cols <- brewer.pal(8, "Greys")
    ni <- c(1,2,length(years)-1,length(years))
    plot(years, srb_mean, type = "n", ylim = c(104,118),
         xlab = "Year", ylab = "SRB (boys per 100 girls)", frame.plot = FALSE)
    grid(lwd = 2)
    polygon(c(years[-ni], rev(years[-ni])),
            c(srb_q10[-ni], rev(srb_q90[-ni])), col = cols[4], border = NA)
    lines(years, srb_mean, lwd = 3, col = cols[7])
    title(paste(country_i, "— Simulated SRB (", model_i, " model)", sep=""))
  }
  invisible(out)
}

# -------------------------- Run all scenarios (same as slow) --------------------------
for (i in seq_along(model)) {
  run_scenario(model[i], country[i], nsim[i], niter[i], ncpus[i])
}
