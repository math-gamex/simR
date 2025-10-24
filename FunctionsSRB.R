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


###########################################################
### FUNCTIONS OF THE MODEL                              ###  
###########################################################

# 1) UPLOAD EMPIRICAL INPUT DATA
get.data <- function(country) {
  
  # 5-year interval
  n <- 5
  
  # INDIA
  if (country == "India") {
    pop1950 <- read.table("Data/UN_IndiaFemPop1950.txt", header=TRUE) 
    pop1980 <- read.table("Data/IndiaFemPop1980.txt", header=FALSE)
    fertRates <- read.table("Data/UN_IndiaASFR.txt", header = FALSE)
    lifeTable <- read.table("Data/UN_India_LifeTable.txt", header=TRUE)
    
    # DEATH RATES
    lifeTable <- lifeTable[!(lifeTable$Sex == "Both"),]
    lifeTable <- lifeTable[lifeTable$AgeGrpStart < 55,]
    # Age-specific mortality rates
    DRall <- lifeTable[, c("Time", "MidPeriod", "Sex", "AgeGrpStart", "AgeGrpSpan", "mx")]
    yearVec <- seq(1950, 2050, 5)
    year <- 1953
    DRmales <- DRall[DRall$Sex == "Male",]
    DRmales <- DRmales[1:(12*length(yearVec)),]
    DRfemales <- DRall[DRall$Sex == "Female",]
    DRfemales <- DRfemales[1:(12*length(yearVec)),]
    
    for (i in 1:length(yearVec)) {
      
      # Females
      fdr <- c(DRfemales[DRfemales$MidPeriod == year, "mx"][1],
               rep(DRfemales[DRfemales$MidPeriod == year, "mx"][2], n-1),
               rep(DRfemales[DRfemales$MidPeriod == year, "mx"][-c(1,2)], each=n))
      # Males
      mdr <- c(DRmales[DRmales$MidPeriod == year, "mx"][1],
               rep(DRmales[DRmales$MidPeriod == year, "mx"][2], n-1),
               rep(DRmales[DRmales$MidPeriod == year, "mx"][-c(1,2)], each=n))
      dr <- cbind(c(rep(year - 2, maxAge + 2), rep(year - 1, maxAge + 2), rep(year, maxAge + 2), 
                    rep(year + 1, maxAge + 2), rep(year + 2, maxAge + 2)),
                  rep(0:(maxAge + 1), n), rep(fdr[1:(maxAge + 2)], n), 
                  rep(mdr[1:(maxAge + 2)], n))
      if (i == 1) {
        deathRates <- dr
      } else {
        deathRates <- rbind(deathRates, dr) 
      }
      year <- year + 5
    }
    deathRates <- as.data.frame(deathRates)
    colnames(deathRates) <- c("year", "age", "female", "male")
  } 
  
  # SOUTH KOREA
  if (country == "Korea") {
    pop1950 <- read.table("Data/UN_KoreaFemPop1950.txt", header=TRUE) 
    pop1980 <- read.table("Data/UN_KoreaFemPop1980.txt", header=TRUE)
    fertRates <- read.table("Data/UN_KoreaASFR.txt", header = FALSE)
    deathRates <- read.table("Data/UN_KoreaDeathRates.txt", header = TRUE)
    deathRates <- deathRates[deathRates$age < 55,]
    for (i in 1:length(unique(deathRates$period))) {
      drmale <- deathRates[deathRates$period == unique(deathRates$period)[i], "male"]
      drfemale <- deathRates[deathRates$period == unique(deathRates$period)[i], "female"]
      if (i == 1) {
        mdr <- rep(c(drmale[1],
                     rep(drmale[2], n-1),
                     rep(drmale[-c(1,2)], each=n)), n)
        fdr <- rep(c(drfemale[1],
                     rep(drfemale[2], n-1),
                     rep(drfemale[-c(1,2)], each=n)), n)
      } else {
        mdr <- c(mdr, rep(c(drmale[1],
                            rep(drmale[2], n-1),
                            rep(drmale[-c(1,2)], each=n)), n))
        fdr <- c(fdr, rep(c(drfemale[1],
                            rep(drfemale[2], n-1),
                            rep(drfemale[-c(1,2)], each=n)), n))
      }
    }
    agevec <- 0:(max(deathRates$age)+4)
    yearvec <- (min(deathRates$period)-4):max(deathRates$period)
    deathRates <- cbind(rep(yearvec, each=length(agevec)),
                        rep(agevec, length(yearvec)),
                        fdr, mdr)
    deathRates <- as.data.frame(deathRates)
    colnames(deathRates) <- c("year", "age", "female", "male")  
  }
  
  # Female population in 1950
  pop1950 <- rep(pop1950$female/n, each=n)
  pop1950 <- cumsum(pop1950)
  femPop1950 <- pop1950 / pop1950[length(pop1950)]
  
  # Female population in 1980
  pop1980 <- rep(pop1980[,2]/n, each=n)
  pop1980 <- cumsum(pop1980[1:50])
  femPop1980 <- pop1980 / pop1980[length(pop1980)]
  
  # FERTILITY RATES
  # Tables for different theoretical scenarios: 'fertility5' is the standard (United Nations)
  # 'fertility3' means that the decline observed in 'fertility5' is compressed in 3-year periods
  # 'fertility7' is the opposite, the 5-year decline is observed in 7 years, etc.
  for (i in 1:nrow(fertRates)) {
    if (i == 1) {
      fertility5 <- rep(unlist(rep(fertRates[i,-1], each=n)), n)
    } else {
      if (i == 6) {
        fertility3 <- rep(unlist(rep(fertRates[i,-1], each=n)), 1)
        fertility4 <- rep(unlist(rep(fertRates[i,-1], each=n)), 1)
        fertility7 <- rep(unlist(rep(fertRates[i,-1], each=n)), 1)
        fertility10 <- rep(unlist(rep(fertRates[i,-1], each=n)), 1)
      } else if (i > 6) {
        fertility3 <- c(fertility3, rep(unlist(rep(fertRates[i,-1], each=n)), n-2))
        fertility4 <- c(fertility4, rep(unlist(rep(fertRates[i,-1], each=n)), n-1))
        fertility7 <- c(fertility7, rep(unlist(rep(fertRates[i,-1], each=n)), n+2))
        fertility10 <- c(fertility10, rep(unlist(rep(fertRates[i,-1], each=n)), n+5))
      }
      fertility5 <- c(fertility5, rep(unlist(rep(fertRates[i,-1], each=n)), n))
    }
  }
  fertility3 <- unname(fertility3 / 1000)
  fertility4 <- unname(fertility4 / 1000)
  fertility5 <- unname(fertility5 / 1000)
  fertility7 <- unname(fertility7 / 1000)
  fertility10 <- unname(fertility10 / 1000)  
  maxYear <- length(fertility3)/(maxAge - minAge)-1
  fertility3 <- cbind(rep(0:maxYear + iniYear,  each = maxAge - minAge),
                      rep(minAge:(maxAge - 1), 1+(n-2)*(nrow(fertRates)-6)), 
                      fertility3)
  fertility3 <- rbind(fertility3, 
                      cbind(rep(maxYear + iniYear + 1, maxAge - minAge),
                            minAge:(maxAge - 1),
                            fertility3[(nrow(fertility3) - 34):nrow(fertility3),3]))
  maxYear <- length(fertility4)/(maxAge - minAge) - 1
  fertility4 <- cbind(rep(0:maxYear + iniYear,  each = maxAge - minAge),
                      rep(minAge:(maxAge - 1), 1+(n-1)*(nrow(fertRates)-6)), 
                      fertility4)
  fertility5 <- cbind(rep((min(fertRates[,1])+1):(max(fertRates[,1])+n), 
                          each = maxAge - minAge),
                      rep(minAge:(maxAge - 1), n*nrow(fertRates)), fertility5)
  maxYear <- length(fertility7)/(maxAge - minAge) - 1
  fertility7 <- cbind(rep(0:maxYear + iniYear,  each = maxAge - minAge),
                      rep(minAge:(maxAge - 1), 1+(n+2)*(nrow(fertRates)-6)), 
                      fertility7)  
  maxYear <- length(fertility10)/(maxAge - minAge) - 1
  fertility10 <- cbind(rep(0:maxYear + iniYear,  each = maxAge - minAge),
                       rep(minAge:(maxAge - 1), 1+(n+5)*(nrow(fertRates)-6)), 
                       fertility10)
  fertility3 <- as.data.frame(fertility3)
  colnames(fertility3) <- c("year", "age", "ASFR")
  fertility4 <- as.data.frame(fertility4)
  colnames(fertility4) <- c("year", "age", "ASFR")
  fertility5 <- as.data.frame(fertility5)
  colnames(fertility5) <- c("year", "age", "ASFR")
  fertility7 <- as.data.frame(fertility7) 
  colnames(fertility7) <- c("year", "age", "ASFR")
  fertility10 <- as.data.frame(fertility10)
  colnames(fertility10) <- c("year", "age", "ASFR")

  # Output
  return(list(femPop1950 = femPop1950, femPop1980 = femPop1980, 
              fertRates = list(fertility5 = fertility5, fertility4 = fertility4, fertility3 = fertility3,
                               fertility7 = fertility7, fertility10 = fertility10), 
              deathRates = deathRates))
}


# 2) OBTAIN SON PREFERENCE
funcSonPref <- function(country) {
  
  # Son preference matrix for the theoretical scenarios
  sonPrefTheory <- rbind(seq(1, 0, -.1), rep(1, 11))
  
  # INDIA
  if (country == "India") {  
    sonPrefData <- read.table("Data/India_SonPrefData_NFHS.txt", header=TRUE, quote="\"")
    dichotomoussonreg_logit <- glm(cbind((sonPrefData[,1]*sonPrefData[,4]), 
                                         ((1 - sonPrefData[,1])* sonPrefData[,4])) ~ (sonPrefData[,3]) + (sonPrefData[,2]), 
                                   family = "binomial")
    
    dichotomoussonreg_logit_onlytime <- glm(cbind((sonPrefData[,1]*sonPrefData[,4]), 
                                                  ((1 - sonPrefData[,1])* sonPrefData[,4])) ~ (sonPrefData[,3]), 
                                            family = "binomial")
    
    years <- seq(1980, 2050, 1)
    cohorts <- c(1:length(cohortNames)) - 6
  
    # Predicted values across time and period effects from logistic regression
    predCohorts <- matrix(NA, nrow = length(years), ncol = length(cohorts))
    for(i in 1: length(years)) {
      for (j in 1:length(cohorts)) {
        predCohorts[i, j] <- exp(dichotomoussonreg_logit$coeff[2]*years[i] + dichotomoussonreg_logit$coeff[3]*cohorts[j] + dichotomoussonreg_logit$coeff[1])/
          (1+exp(dichotomoussonreg_logit$coeff[2]*years[i] + dichotomoussonreg_logit$coeff[1] + dichotomoussonreg_logit$coeff[3]*cohorts[j]))
        
      }
    }
    colnames(predCohorts) <- cohortNames
    rownames(predCohorts) <- years
    spMeanPeriod <- matrix(NA, nrow = 1, ncol = length(years))
    colnames(spMeanPeriod) <- years
    
    for (i in 1: length(years)) { 
      cohortstoconsider <- intersect(which(years[i] - as.numeric(colnames(predCohorts)) >= minAge),
                                     which(years[i] - as.numeric(colnames(predCohorts)) <= maxAge))
      spMeanPeriod[1,i] <- mean(predCohorts[i, cohortstoconsider])
    }
    
    observed <- cbind(c(1992.5,   1998.5,   2005.5),
                      c(0.9002342, 0.8585739, 0.7901512))  
    
    return(list(sonPrefMx = predCohorts, sonPrefTheory = sonPrefTheory,
                meanSonPref = spMeanPeriod, obsData = observed))  
  }
  
  # KOREA
  if (country == "Korea") {
    # Proportion stating "must have son"
    sp <- c(.477, .418, 0.263, 0.248, 0.192, 0.168)
    # Time points for "must have son"
    time <- c(1985, 1991, 1994, 1997, 2000, 2003)
    # Logistic regression with one covariate (time)
    dichotomoussonreg_logit <- glm(cbind(sp*1000, (1000 - (sp*1000))) ~ time, 
                                   family = "binomial")
    
    # Years to predict for = periodStudy
    years <- seq(1980, 2050, 1) 
    
    # As we have only one covariate (time), we create prediction by time only
    predTime <- rep(NA, length(years))
    
    for(i in 1: length(years)) {
        predTime[i] <- exp(dichotomoussonreg_logit$coeff[2]*years[i] + dichotomoussonreg_logit$coeff[1])/
          (1+exp(dichotomoussonreg_logit$coeff[2]*years[i] + dichotomoussonreg_logit$coeff[1]))
        
    }
    # predTime generates predictions for each year of the simulation (periodStudy)
    predCohorts <- matrix(rep(predTime, each = length(cohortNames)), 
                          nrow=length(predTime), ncol = length(cohortNames))
    colnames(predCohorts) <- cohortNames
    rownames(predCohorts) <- years
    
    return(list(sonPrefMx = predCohorts, sonPrefTheory = sonPrefTheory) ) 
  }
}

# 3) INITIALIZATION OF THE MODEL
initialization <- function(pop, startYear) {
  
  # INPUT DATA
  # Initial distribution of the population
  femPop <- inputData$femPop1950
  # Death rates
  deathRates <- inputData$deathRates
  # Fertility rates
  fertRates <- inputData$fertRates[[1]]
  # Starting year to initialize the simulations
  year <- startYear
  
  # INITIAL VARIABLES
  id <- 1:pop
  age <- rep(NA, pop)
  age <- findInterval(runif(pop), femPop, rightmost.closed = TRUE)
  cohort <- trunc((startYear - age)/cohortLen,0)*cohortLen
  sex <- rep(1, pop)
  children <- rep(0, pop)
  sons <- rep(0, pop)
  mother <- rep(0, pop)
  dr <- rep(0, pop)
  fr <- rep(0, pop)
  idrep <- which(age >= minAge & age < maxAge)
  if (year < min(deathRates$year)) {
    dr <- deathRates[deathRates$year==min(deathRates$year), "female"][age+1]
    fr[idrep] <- fertRates[fertRates$year==min(deathRates$year), "ASFR"][age[idrep] - minAge + 1]
  } else {
    dr <- deathRates[deathRates$year==year, "female"][age+1]
    fr[idrep] <- fertRates[fertRates$year==year, "ASFR"][age[idrep] - minAge + 1]
  }
  
  # Population evolution from 1945 to 1979
  for (i in 1:(iniYear-startYear)) {
    # AGEING
    idold <- which(age >= maxAge)
    if (length(idold) != 0) {
      id <- id[-idold]
      age <- age[-idold]
      sex <- sex[-idold]
      cohort <- cohort[-idold]
      children <- children[-idold]
      sons <- sons[-idold]
      dr <- dr[-idold]
      fr <- fr[-idold]
      mother <- mother[-idold]
      pop <- length(id)
    }
    # Identify those who will die
    iddie <- which(runif(pop) < dr)
    if (length(iddie) != 0) {
      # Mothers of individuals to die
      idmother <- which(id %in% mother[iddie])
      children[idmother] <- children[idmother] - 1
      # Males to die
      malesToDie <- which(sex[iddie] == 2)
      # Which mothers of males to die are alive
      if (length(malesToDie) != 0) {
        idmother <- which(id %in% mother[iddie][malesToDie])
        sons[idmother] <- sons[idmother] - 1
      }
      id <- id[-iddie]
      age <- age[-iddie]
      sex <- sex[-iddie]
      cohort <- cohort[-iddie]
      children <- children[-iddie]
      sons <- sons[-iddie]
      dr <- dr[-iddie]
      fr <- fr[-iddie]
      mother <- mother[-iddie]
      pop <- length(id)
    }
    
    # REPRODUCING    
    # Women in reproductive ages
    idrep <- which(sex == 1 & age >= minAge & age < maxAge)
    # Women having a child
    idborn <- which(runif(length(idrep)) < fr[idrep])
    if (length(idborn) != 0) {
      # Select sex of child
      sexNew <- sample(c(1, 2), length(idborn), replace = T, prob = c(1 - srb, srb))
      # Update number of children and sons
      children[idrep][idborn] <- children[idrep][idborn] + 1
      sons[idrep][idborn][which(sexNew == 2)] <- sons[idrep][idborn][which(sexNew == 2)] + 1
      # Add new individuals to the data set
      id <- c(id, 1:length(idborn) + max(id))
      sex <- c(sex, sexNew)
      age <- c(age + 1, rep(0, length(idborn)))
      cohort <- c(cohort, rep(trunc(year/cohortLen, 0)*cohortLen, length(idborn)))
      children <- c(children, rep(0, length(idborn)))
      sons <- c(sons, rep(0, length(idborn)))
      mother <- c(mother, id[idrep][idborn])
    }
    
    # Update population size and year
    pop <- length(id)
    year <- year + 1
    
    # Update fertility rates
    fr <- rep(0, pop)
    if (year < min(deathRates$year)) {
      fr[idrep] <- fertRates[fertRates$year==min(deathRates$year), "ASFR"][age[idrep] - minAge + 1]
    } else {
      fr[idrep] <- fertRates[fertRates$year==year, "ASFR"][age[idrep] - minAge + 1]
    }
    
    # Update death rates vector
    idfemale <- which(sex == 1)
    idmale <- which(sex == 2)
    dr <- rep(0, pop)
    if (year < min(deathRates$year)) {
      dr[idfemale] <- deathRates[deathRates$year==min(deathRates$year), "female"][age[idfemale]+1]
      dr[idmale] <- deathRates[deathRates$year==min(deathRates$year), "male"][age[idmale]+1]    
    } else {
      dr[idfemale] <- deathRates[deathRates$year==year, "female"][age[idfemale]+1]
      dr[idmale] <- deathRates[deathRates$year==year, "male"][age[idmale]+1]    
    }
  }
  
  # Son preference
  sonPrefYear <- rbind(1 - sonPrefMx[rownames(sonPrefMx) == iniYear,],
                       rep(1, ncol(sonPrefMx)))
  sonPref <- runif(pop)
  sonPref <- mapply(findSonPref, sonPref, cohort, MoreArgs = list(M = sonPrefYear))
  
  # Output
  return(list(id = id, age = age, sex = sex, cohort = cohort, children = children,
              sons = sons, dr = dr, fr = fr, mother = mother, 
              sonPref = sonPref))  
}


# 4) RUN THE MODEL
run.simulations <- function(k) {
  
  # Parameters
  alpha <- alphaPar[1]
  gamma <- gammaPar[1]
  rho <- rhoPar[1]
  phi <- phiPar[1]
  sigma <- sigmaPar[1]
  beta <- betaPar[1]
  # Fertility rates
  fertRates <- inputData$fertRates[[fertMod[1]]]
  # Son preference theoretical scenarios
  sonTheo <- sonTheoPar[1]
  
  # MODIFY INPUT DEPENDING ON THE MODEL
  # EMPIRICAL
  if (model == "Empirical") {
    
    # Initial Son Preference
    sonPrefIni <- startPop$sonPref
    
    # Parameters: k-th element of each of the parameter's vectors
    alpha <- alphaPar[k]
    gamma <- gammaPar[k]
    rho <- rhoPar[k]
    phi <- phiPar[k]
    sigma <- sigmaPar[k]
    beta <- betaPar[k]
  }
  # THEORETICAL
  if (model == "Theory") {
    
    # Initial Son Preference
    if (length(sonTheoPar) > 1) sonTheo <- sonTheoPar[k]
    if (!sensitivity) {
      sonPrefIni <- findInterval(runif(length(startPop$id)), sonPrefTheory[, sonTheo],
                                 rightmost.close=TRUE)
    }
    
    # Fertility rates: Different scenarios of fertility decline
    if (length(sonTheoPar) > 1) {
      fertRates <- inputData$fertRates[[fertMod]]      
    } else {
      fertRates <- inputData$fertRates[[k]]
      fertMod <- fertMod[k]
    }    
  }
      
  # SENSITIVITY ANALYSIS
  if (sensitivity) {
    alpha <- alphaPar[k]
    gamma <- gammaPar[k]
    rho <- rhoPar[k]
    phi <- phiPar[k]
    sigma <- sigmaPar[k]
    beta <- betaPar[k]
    if (model == "Theory") {
      sonPrefIni <- findInterval(runif(length(startPop$id)), c(1 - sonTheo, 1),
                                 rightmost.close=TRUE)
    }
  }
    
  # DEATH RATES
  deathRates <- inputData$deathRates
    
  # OUTPUT TABLES
  # Number of boys born per cohort and year
  boysCohort <- array(0, dim=c(niter, length(cohortNames), nsim))
  colnames(boysCohort) <- cohortNames
  rownames(boysCohort) <- periodStudy
  
  # Parity per year (not per cohort)
  parityMale <- array(0, dim=c(niter, 6, nsim))
  rownames(parityMale) <- periodStudy
  
  # Proportion boys born per year (not per cohort)
  propBoys <- array(0, dim=c(niter, nsim))
  rownames(propBoys) <- periodStudy
  
  # Average child per women alife in age of reproduction, per cohort and total per year
  avgChildCohort <- array(0, dim=c(niter, length(cohortNames)+1, nsim))
  colnames(avgChildCohort) <- c(cohortNames, "AvgChildYear")
  rownames(avgChildCohort) <- periodStudy
  
  # Abortion data
  abortData <- array(0, dim=c(niter, nsim))
  rownames(abortData) <- periodStudy
  
  # Fertility data
  fertilityData <- array(0, dim=c(niter, 10, nsim))
  colnames(fertilityData) <- c("15-19", "20-24", "25-29", "30-34", "35-39",
                               "40-44", "45-49", "totalBirths","birthRate", "TFR")
  rownames(fertilityData) <- periodStudy
  
  # RUN SIMULATIONS
  for (j in 1:nsim) {
    
    # INITIAL POPULATION: Data stored in vectors; each cell corresponds to 1 individual
    # Individual ID
    id <- startPop$id
    # Population size
    pop <- length(id)
    # Age
    age <- startPop$age
    # Cohort
    cohort <- startPop$cohort
    # Sex
    sex <- startPop$sex
    # Number of children
    children <- startPop$children
    # Number of sons
    sons <- startPop$sons
    # Fertility rate
    fr <- startPop$fr
    # Death rate
    dr <- startPop$dr
    # Mother ID
    mother <- startPop$mother
    # Son preference (1 or 0)
    sonPref <- sonPrefIni
    # Abortion rate
    ar <- rep(0, pop)
    # Number of abortions
    abortion <- rep(0, pop)
    # Acces to hosopital / technology
    hospital <- rep(0, pop)
    
    # Initial year for the simulations
    year <- iniYear
    
    # START SIMULATIONS
    for (i in 1:niter) {
          
      # DYING
      idold <- which(age >= maxAge)
      if (length(idold) != 0) {
        id <- id[-idold]
        age <- age[-idold]
        sex <- sex[-idold]
        cohort <- cohort[-idold]
        children <- children[-idold]
        sons <- sons[-idold]
        sonPref <- sonPref[-idold]
        dr <- dr[-idold]
        fr <- fr[-idold]
        ar <- ar[-idold]
        mother <- mother[-idold]
        abortion <- abortion[-idold]
        hospital <- hospital[-idold]
        pop <- length(age)
      }
      # a) Identify those who will die
      iddie <- which(runif(pop) < dr)
      if (length(iddie) != 0) {
        # Mothers of individuals to die
        idmother <- which(id %in% mother[iddie])
        children[idmother] <- children[idmother] - 1
        # Males to die
        malesToDie <- which(sex[iddie] == 2)
        # Which mothers of males to die are alive
        if (length(malesToDie) != 0) {
          idmother <- which(id %in% mother[iddie][malesToDie])
          sons[idmother] <- sons[idmother] - 1
        }
        id <- id[-iddie]
        age <- age[-iddie]
        sex <- sex[-iddie]
        cohort <- cohort[-iddie]
        children <- children[-iddie]
        sons <- sons[-iddie]
        sonPref <- sonPref[-iddie]
        dr <- dr[-iddie]
        fr <- fr[-iddie]
        ar <- ar[-iddie]
        mother <- mother[-iddie]
        abortion <- abortion[-iddie]
        hospital <- hospital[-iddie]
        pop <- length(age)
      }
      # b) Update age of individuals alive
      sonYear <- rep(0, pop)
      childYear <- rep(0, pop)

      #### ACCESS TO TECHNOLOGY ####
      # Diffusion parameter
      diffPar <- exp(rho * (i - phi)) / (1 + exp(rho*(i - phi)))
      # Update access to technology
      idnoTech <- which(age >= minAge & sex == 1 & hospital == 0)
      idNewTech <- which(runif(length(idnoTech)) < diffPar)
      hospital[idnoTech][idNewTech] <- 1

      #### REPRODUCTION ####
      # Women in reproductive ages
      idrep <- which(sex == 1 & age >= minAge & age < maxAge)
      # Women having a child
      idborn <- which(runif(length(idrep)) < fr[idrep])
      if (length(idborn) != 0) {
        # Select sex of child
        sexNew <- sample(c(1, 2), length(idborn), replace = T, prob = c(1 - srb, srb))
        idgirls <- which(sexNew == 1) 
        # Abortion procedure
        idabort <- which(runif(length(idgirls)) < ar[idrep][idborn][idgirls] &
                           sons[idrep][idborn][idgirls] < sonPref[idrep][idborn][idgirls] &
                           hospital[idrep][idborn][idgirls] == 1)
        if (length(idabort) != 0) {
          abortion[idrep][idgirls][idabort] <- abortion[idrep][idgirls][idabort] + 1
          sexNew[idgirls][idabort] <- 9
          idborn <- idborn[-which(sexNew == 9)]
          sexNew <- sexNew[-which(sexNew == 9)]
          abortYear <- length(idabort)
        }
        # Update number of children and sons
        children[idrep][idborn] <- children[idrep][idborn] + 1
        childYear[idrep][idborn] <- 1
        sons[idrep][idborn][which(sexNew == 2)] <- sons[idrep][idborn][which(sexNew == 2)] + 1
        sonYear[idrep][idborn][which(sexNew == 2)] <- 1
        # Add new individuals to the data set
        id <- c(id, 1:length(idborn) + max(id))
        sex <- c(sex, sexNew)
        age <- c(age, rep(0, length(idborn)))
        cohort <- c(cohort, 
                    rep(trunc((year)/cohortLen, 0)*cohortLen, length(idborn)))
        children <- c(children, rep(0, length(idborn)))
        childYear <- c(childYear, rep(0, length(idborn)))
        sons <- c(sons, rep(0, length(idborn)))
        sonYear <- c(sonYear, rep(0, length(idborn)))
        hospital <- c(hospital, rep(0, length(idborn)))
        mother <- c(mother, id[idrep][idborn])
        abortion <- c(abortion, rep(0, length(idborn)))
      } else idabort <- NULL
      
      #### STORE RESULTS ####
      # Abortion data and proportion of boys born per year
      if (sum(childYear) != 0) {
        abortData[i, j] <- round(100*length(idabort) / sum(childYear), 6)
        propBoys[i, j] <- round(sum(sonYear) / sum(childYear), 6)
      }
      
      # Period Fertility Rate
      id15 <- which(sex == 1 & age >= 15 & age < 20)
      id20 <- which(sex == 1 & age >= 20 & age < 25)
      id25 <- which(sex == 1 & age >= 25 & age < 30)
      id30 <- which(sex == 1 & age >= 30 & age < 35)
      id35 <- which(sex == 1 & age >= 35 & age < 40)
      id40 <- which(sex == 1 & age >= 40 & age < 45)
      id45 <- which(sex == 1 & age >= 45 & age < 50)
      fertilityData[i, 1, j] <- sum(childYear[id15]) / length(id15)      
      fertilityData[i, 2, j] <- sum(childYear[id20]) / length(id20)      
      fertilityData[i, 3, j] <- sum(childYear[id25]) / length(id25)      
      fertilityData[i, 4, j] <- sum(childYear[id30]) / length(id30)      
      fertilityData[i, 5, j] <- sum(childYear[id35]) / length(id35)      
      fertilityData[i, 6, j] <- sum(childYear[id40]) / length(id40)      
      fertilityData[i, 7, j] <- sum(childYear[id45]) / length(id45)      
      fertilityData[i, 8, j] <- sum(childYear)      
      fertilityData[i, 9, j] <- sum(childYear)/length(idrep)
      fertilityData[i, 10, j] <- round(sum(fertilityData[i, 1:7, j])*cohortLen, 6)

      # Update year and population size
      if (year < max(deathRates$year)) year <- year + 1
      pop <- length(id)
      
      # Son preferences
      if (model == "Empirical") {
        sonPrefYear <- rbind(1 - sonPrefMx[rownames(sonPrefMx) == year,],
                             rep(1, ncol(sonPrefMx)))
        sonPref <- mapply(findSonPref, runif(pop), cohort, 
                          MoreArgs = list(M = sonPrefYear))
      }
      if (model == "Theory") {
        if (!sensitivity) {
          sonPref <- c(sonPref,
                       findInterval(runif(length(idborn)), sonPrefTheory[, sonTheo],
                                    rightmost.close=TRUE))
        } else sonPref <- c(sonPref,
                            findInterval(runif(length(idborn)), c(1 - sonTheo, 1),
                                         rightmost.close=TRUE))
      }
            
      # UPDATE AGE, except for new borns
      age[1:(pop-length(idborn))] <- age[1:(pop-length(idborn))] + 1
            
    
      #### UPDATE VARIABLES ####
      
      # Indicator: sex
      idfemale <- which(sex == 1)
      idmale <- which(sex == 2)
      # Indicator: women in reproductive ages
      idrep <- which(sex == 1 & age >= minAge & age < maxAge)
      # Indicator: women that have reached son and parity preferences
      idsucc <- which(sex == 1 & sons >= sonPref & sonPref != 0)
      # Indicator: women that have NOT reached son preferences
      idunsu <- which(sex == 1 & sons < sonPref)
            
      # Update death rates vector
      dr <- rep(0, pop)
      dr[idfemale] <- deathRates[deathRates$year==year, "female"][age[idfemale]+1]
      dr[idmale] <- deathRates[deathRates$year==year, "male"][age[idmale]+1]
      
      # Update fertility rates
      fr <- rep(0, pop)
      fr[idrep] <- fertRates[fertRates$year==year, "ASFR"][age[idrep] - minAge + 1]
      fr[idsucc] <- fr[idsucc]*(1 - alpha)
      fr[idunsu] <- fr[idunsu]*(1 + gamma)
      
      # Update abortion rate
      ar <- rep(0, pop)
      # Identify candidates for abortion
      idabort <- which(sex == 1 & sons < sonPref)
      TFR <- sum(fertilityData[i, 1:7, j]) * cohortLen
      ar[idabort] <- children[idabort] / TFR
      ar <- ar * sigma
      ar[which(sex == 1 & sons < sonPref & children == 0)] <- beta / TFR 
      ar[which(ar > 1)] <- 1
    }  
  }
  
  # Return results
  if (!(sensitivity)) {
    TFR <- fertilityData[, 10, 1:nsim]
    return(list(propBoys = propBoys, TFR = TFR, abortData = abortData,
                alpha = alpha, gamma = gamma, rho = rho, sigma = sigma, phi = phi, beta = beta, 
                sonPrefTheo = 1 - sonPrefTheory[1, sonTheo], fertMod = fertMod,
                nsim = nsim, niter = niter, iniYear = iniYear, endYear = year,
                iniPop = length(startPop$id), endPop = pop))
  } else {
    # Compute means accross simulations
    TFR <- fertilityData[, 10, 1:nsim]
    return(list(propBoys = propBoys, TFR = TFR, abortData = abortData,
                alpha = alpha, gamma = gamma, rho = rho, phi = phi, 
                sigma = sigma, beta = beta, sonPrefTheo = sonTheo,
                nsim = nsim, iniPop = length(startPop$id), endPop = pop))
  }
}


# 5) FIND SON PREFERENCE
findSonPref <- function(x, y, M) {
  M <- unlist(M)
  findInterval(x, M[, paste(y)], rightmost.closed = TRUE)
}

#FUNCTIONS to extract data for plotting results and study model behavior
#There are four functions: 
#(1) extractdata_all -- this recovers TFR and SRB values for all years for different simulation types (empirical, lhs (sensitivity analysis), or theoretical)
#(2) extractdata_1990 -- SRB value in 1990 in simulation run
#(3) extractdata_maxsrb -- extracts maximum SRB value attained in a simulation run
#(4) extractdata_rmse -- calculates RMSE for each simulation run (see eq. 5 in paper)

# Moving average function, n = number of years for plotting
ma <- function(x,n=5){
  filter(x,rep(1/n,n), sides=2)
}


extractdata_all <- function(file, type) {
  load(file = paste("SimulationResults/", file, ".RData", sep = ""))
  x <- out
  outputall <- c()
  for (i in 1:length(x)) {
    if (type == "empirical") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$rho)),  length(x[[i]]$propBoys['1990', ])), nrow = length(x[[i]]$propBoys['1990', ]), 
                      ncol = 6, byrow = TRUE) }
    if (type == "lhs") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$rho)),  length(x[[i]]$propBoys['1990', ])), nrow = length(x[[i]]$propBoys['1990', ]), 
                      ncol = 6, byrow = TRUE) }
    
    if (type == "theoretical") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$sonPrefTheo, x[[i]]$rho)),  length(x[[i]]$propBoys['1990', ])), nrow = length(x[[i]]$propBoys['1990', ]), 
                      ncol = 7, byrow = TRUE) }
    
    
    propboys <- t((x[[i]]$propBoys / (1 - x[[i]]$propBoys)) * 100)
    
    
    if (length(x[[i]]$fertilityData) == 0) {
      tfr <- t(x[[i]]$TFR)
    } else  tfr <- t(x[[i]]$fertilityData[,ncol(x[[i]]$fertilityData),])
    
    
    output <- cbind(pars, propboys, tfr)
    outputall <- rbind(outputall, output) 
  }
  
  if (type == "empirical") { 
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'rho', paste('srb', seq(1980, 2049,1), sep = ''), paste('tfr', seq(1980, 2049,1), sep = '')) }
  if (type == "lhs") { 
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'rho', paste('srb', seq(1980, 2049,1), sep = ''), paste('tfr', seq(1980, 2049,1), sep = '')) }
  if (type == "theoretical") {
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'sonpref', 'rho', paste('srb', seq(1980, 2029,1), sep = ''), paste('tfr', seq(1980, 2029,1), sep = ''))}
  
  return(outputall)
}


#Function to extract SRB in 1990 all values (rather than mean)
extractdata_1990_v2 <- function(file, type) {
  load(file = paste("SimulationResults/", file, ".RData", sep = ""))
  x <- out
  outputall <- c()
  for (i in 1:length(x)) {
    if (type == "empirical") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$rho)),  length(x[[i]]$propBoys['1990', ])), nrow = length(x[[i]]$propBoys['1990', ]), 
                      ncol = 6, byrow = TRUE) }
    if (type == "theoretical") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$sonPrefTheo, x[[i]]$rho)),  length(x[[i]]$propBoys['1990', ])), nrow = length(x[[i]]$propBoys['1990', ]), 
                      ncol = 7, byrow = TRUE) }
    propboys <- apply(x[[i]]$propBoys, 2, ma)[11,] / (1 - apply(x[[i]]$propBoys, 2, ma)[11,]) * 100
    #propboys <- (x[[i]]$propBoys['1990',] / (1 - x[[i]]$propBoys['1990',])) * 100
    output <- cbind(pars, propboys)
    
    outputall <- rbind(outputall, output) }
  if (type == "empirical") { 
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'rho', 'srb') }
  if (type == "theoretical") {
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'sonpref', 'rho', 'srb')}
  return(outputall)
  
}

#Function to extract max SRB
extractdata_maxsrb <- function(file, type) {
  load(file = paste("SimulationResults/", file, ".RData", sep = ""))
  x <- out
  outputall <- c()
  for (i in 1:length(x)) {
    if (type == "empirical") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$rho)),  length(x[[i]]$propBoys['1990', ])), nrow = length(x[[i]]$propBoys['1990', ]), 
                      ncol = 6, byrow = TRUE)}
    
    if (type == "theoretical") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$sonPrefTheo, x[[i]]$rho)),  length(x[[i]]$propBoys['1990', ])), nrow = length(x[[i]]$propBoys['1990', ]), 
                      ncol = 7, byrow = TRUE)}
    propboysint <- apply(x[[i]]$propBoys, 2, max)
    #propboysint <- apply(apply(x[[1]]$propBoys, 2, ma), 2, max, na.rm = TRUE)
    propboys <- (propboysint / (1 - propboysint)) * 100 
    output <- cbind(pars, propboys)
    
    outputall <- rbind(outputall, output) }
  if (type == "empirical") {
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'rho', 'srb') }
  if (type == "theoretical") {
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'sonpref', 'rho', 'srb') }
  return(outputall)
  
}



#extractdata RMSE

extractdata_rmse <- function(file, type) {
  load(file = paste("/SimulationResults/", file, ".RData", sep = ""))
  x <- out
  outputall <- c()
  for (i in 1:length(x)) {
    if (type == "empirical") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$rho)),  1), nrow = 1, 
                      ncol = 6, byrow = TRUE) }
    if (type == "lhs") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$rho)),  1), nrow = 1, 
                      ncol = 6, byrow = TRUE) }
    
    if (type == "theoretical") {
      pars <-  matrix(rep(rbind(c(x[[i]]$alpha, x[[i]]$gamma, x[[i]]$sigma, 
                                  x[[i]]$beta, x[[i]]$phi, x[[i]]$sonPrefTheo, x[[i]]$rho)),  1), nrow = 1, 
                      ncol = 7, byrow = TRUE) }
    
    
    propboys <- as.matrix(t((x[[i]]$propBoys / (1 - x[[i]]$propBoys)) * 100))
    propboysts <- apply(propboys[,1:30], 1, ma)
    
    
    srb.skorea.wpp2012 <- c(1.07,   1.07,	 1.14,	 1.14,	 1.10,	 1.10,	 1.07,
                            1.07,	 1.07,	 1.07,	 1.07,	 1.07,	 1.07,	 1.07,	 1.07,	 1.07)
    
    time <- seq(1978, 2053, 5)
    
    srb.approx.sk <- approx(x = time, y = srb.skorea.wpp2012, n = 76)$y * 100
    srb.approx.plot.sk <- srb.approx.sk[3:(length(srb.approx.sk)-4)]
    
    srb.output <- ma(srb.approx.plot.sk[1:30])
    rmse <- apply(propboysts, 2, FUN = function(x) { sqrt(sum((srb.output - x)^2, na.rm = TRUE) / 30)})
    rmse_mean <- mean(rmse)
    
    
    output <- cbind(pars, rmse_mean)
    
    outputall <- rbind(outputall, output) }
  if (type == "empirical") { 
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'rho', 'rmse') }
  if (type == "lhs") { 
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'rho', 'rmse') }
  if (type == "theoretical") {
    colnames(outputall) <- c('alpha', 'gamma', 'sigma', 'beta', 'phi', 'sonpref', 'rho', 'rmse')}
  return(outputall)
  
}

#DATA FOR CALIBRATION for South Korea
#UN estimates of TFR and SRB
time <- seq(1978, 2053, 5) #retrieved from wppexplorer. 
#data from 1980 - 2055
tfr.skorea.wpp2012 <- c(2.92, 2.23,	 1.60,	 1.70,	 1.51,	 1.22,	 1.23,	 1.32,	
                        1.39,	 1.46,	 1.52,	 1.57,	 1.61,	 1.65,	 1.68,	 1.71)

tfr.approx.sk <- approx(x = time, y = tfr.skorea.wpp2012, n = 76)$y
tfr.approx.plot.sk <- tfr.approx.sk[3:(length(tfr.approx.sk)-4)]

srb.skorea.wpp2012 <- c(1.07,   1.07,	 1.14,	 1.14,	 1.10,	 1.10,	 1.07,
                        1.07,	 1.07,	 1.07,	 1.07,	 1.07,	 1.07,	 1.07,	 1.07,	 1.07)

srb.approx.sk <- approx(x = time, y = srb.skorea.wpp2012, n = 76)$y * 100
srb.approx.plot.sk <- srb.approx.sk[3:(length(srb.approx.sk)-4)]