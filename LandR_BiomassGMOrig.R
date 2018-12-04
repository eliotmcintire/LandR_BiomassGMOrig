
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "LandR_BiomassGMOrig",
  description = "insert module description here",
  keywords = c("insert key words here"),
  authors = person("First", "Last", email = "first.last@example.com", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("1.3.1.9035"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "LandR_BiomassGMOrig"),
  reqdPkgs = list("data.table", "raster", "PredictiveEcology/pemisc@development"), ## TODO: update package list
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("calibrate", "logical", TRUE, NA, NA, "should the model have detailed outputs?"),
    defineParameter("growthInitialTime", "numeric", default = 0, min = NA_real_, max = NA_real_,
                    desc = "Initial time for the growth event to occur"),
    defineParameter("successionTimestep", "numeric", 10, NA, NA,
                    desc = "defines the simulation time step, default is 10 years"),
    defineParameter(".plotInitialTime", "numeric", default = 0, min = NA, max = NA,
                    desc = "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".saveInitialTime", "numeric", default = 0, min = NA, max = NA,
                    desc = paste("This describes the simulation time at which the first save event should occur.",
                                 "Set to NA if no saving is desired.")),
    defineParameter(".useParallel", "ANY", default = parallel::detectCores(),
                    desc = paste("Used only in seed dispersal. If numeric, it will be passed to data.table::setDTthreads,",
                                 "if logical and TRUE, it will be passed to parallel::makeCluster,",
                                 "and if cluster object it will be passed to parallel::parClusterApplyLB"))
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput("cohortData", "data.table",
                 desc = "age cohort-biomass table hooked to pixel group map by pixelGroupIndex at
                 succession time step",
                 sourceURL = NA),
    expectsInput("lastReg", "numeric",
                 desc = "time at last regeneration", sourceURL = NA),
    expectsInput("species", "data.table",
                 desc = "a table that has species traits such as longevity...",
                 sourceURL = "https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession-Archive/master/biomass-succession-archive/trunk/tests/v6.0-2.0/species.txt"),
    expectsInput("speciesEcoregion", "data.table",
                 desc = "table defining the maxANPP, maxB and SEP, which can change with both ecoregion and simulation time",
                 sourceURL = "https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession-dynamic-inputs_test.txt")
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput("calculateAgeMortality", "function",
                  desc = "function to calculate aging and mortality"),
    createsOutput("calculateANPP", "function",
                  desc = "function to calculate ANPP"),
    createsOutput("calculateCompetition", "function",
                  desc = "function to calculate competition for light"),
    createsOutput("calculateGrowthMortality", "data.table",
                  desc = "function to calculate growth and mortality"),
    createsOutput("calculateSumB", "function",
                  desc = "function to sum biomass"),
    createsOutput("cohortData", "function",
                  desc = "tree-level data by pixel group"),
    createsOutput("simulationTreeOutput", "data.table",
                  desc = "Summary of several characteristics about the stands, derived from cohortData"),
    createsOutput("updateSpeciesAttributes", "function",
                  desc = "function to add/update species attributes in species cohort table"),
    createsOutput("updateSpeciesEcoregionAttributes", "function",
                  desc = "function to add/update species ecoregion attributes in species cohort table")
   )
))

## event types
#   - type `init` is required for initialiazation
doEvent.LandR_BiomassGMOrig = function(sim, eventTime, eventType, debug = FALSE) {
  if (is.numeric(P(sim)$useParallel)) {
    a <- data.table::setDTthreads(P(sim)$useParallel)
    message("Mortality and Growth should be using >100% CPU")
    on.exit(setDTthreads(a))
  }
  switch(eventType,
         init = {
           ## do stuff for this event
           sim <- Init(sim)

           sim <- scheduleEvent(sim, start(sim) + P(sim)$growthInitialTime,
                                "LandR_BiomassGMOrig", "mortalityAndGrowth", eventPriority = 5)
         },

         mortalityAndGrowth = {
           sim <- MortalityAndGrowth(sim)
           sim <- scheduleEvent(sim, time(sim) + 1, "LandR_BiomassGMOrig", "mortalityAndGrowth",
                                eventPriority = 5)
         },
         warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                       "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
Init <- function(sim) {
  ## export local functions to simList
  sim$updateSpeciesEcoregionAttributes <- updateSpeciesEcoregionAttributes
  sim$updateSpeciesAttributes <- updateSpeciesAttributes
  sim$calculateSumB <- calculateSumB
  sim$calculateAgeMortality <- calculateAgeMortality
  sim$calculateANPP <- calculateANPP
  sim$calculateGrowthMortality <- calculateGrowthMortality
  sim$calculateCompetition <- calculateCompetition

  return(invisible(sim))
}

MortalityAndGrowth <- function(sim) {
  if (is.numeric(P(sim)$useParallel)) {
    data.table::setDTthreads(P(sim)$useParallel)
    message("Mortality and Growth should be using >100% CPU")
  }

  sim$cohortData <- sim$cohortData[, .(pixelGroup, ecoregionGroup,
                                       speciesCode, age, B, mortality, aNPPAct)]
  cohortData <- sim$cohortData
  sim$cohortData <- cohortData[0, ]
  pixelGroups <- data.table(pixelGroupIndex = unique(cohortData$pixelGroup),
                            temID = 1:length(unique(cohortData$pixelGroup)))
  cutpoints <- sort(unique(c(seq(1, max(pixelGroups$temID), by = 10^4), max(pixelGroups$temID))))
  #cutpoints <- c(1,max(pixelGroups$temID))
  if (length(cutpoints) == 1) cutpoints <- c(cutpoints, cutpoints + 1)
  pixelGroups[, groups := cut(temID, breaks = cutpoints,
                              labels = paste("Group", 1:(length(cutpoints) - 1), sep = ""),
                              include.lowest = T)]
  for (subgroup in paste("Group", 1:(length(cutpoints) - 1), sep = "")) {
    subCohortData <- cohortData[pixelGroup %in% pixelGroups[groups == subgroup, ]$pixelGroupIndex, ]
    #   cohortData <- sim$cohortData
    set(subCohortData, NULL, "age", subCohortData$age + 1)
    subCohortData <- updateSpeciesEcoregionAttributes(speciesEcoregion = sim$speciesEcoregion,
                                                      time = round(time(sim)),
                                                      cohortData = subCohortData)
    subCohortData <- updateSpeciesAttributes(species = sim$species, cohortData = subCohortData)
    subCohortData <- calculateSumB(cohortData = subCohortData,
                                       lastReg = sim$lastReg,
                                       simuTime = time(sim),
                                       successionTimestep = P(sim)$successionTimestep)
    subCohortData <- subCohortData[age <= longevity,]
    subCohortData <- calculateAgeMortality(cohortData = subCohortData)
    set(subCohortData, NULL, c("longevity", "mortalityshape"), NULL)
    subCohortData <- calculateCompetition(cohortData = subCohortData)
    if (!P(sim)$calibrate) {
      set(subCohortData, NULL, "sumB", NULL)
    }
    #### the below two lines of codes are to calculate actual ANPP
    subCohortData <- calculateANPP(cohortData = subCohortData)
    set(subCohortData, NULL, "growthcurve", NULL)
    set(subCohortData, NULL, "aNPPAct", pmax(1, subCohortData$aNPPAct - subCohortData$mAge))
    subCohortData <- calculateGrowthMortality(cohortData = subCohortData)
    set(subCohortData, NULL, "mBio", pmax(0, subCohortData$mBio - subCohortData$mAge))
    set(subCohortData, NULL, "mBio", pmin(subCohortData$mBio, subCohortData$aNPPAct))
    set(subCohortData, NULL, "mortality", subCohortData$mBio + subCohortData$mAge)
    set(subCohortData, NULL, c("mBio", "mAge", "maxANPP", "maxB", "maxB_eco", "bAP", "bPM"), NULL)
    if (P(sim)$calibrate) {
      set(subCohortData, NULL, "deltaB", as.integer(subCohortData$aNPPAct - subCohortData$mortality))
      set(subCohortData, NULL, "B", subCohortData$B + subCohortData$deltaB)
      tempcohortdata <- subCohortData[,.(pixelGroup, Year = time(sim), siteBiomass = sumB, speciesCode,
                                         Age = age, iniBiomass = B - deltaB, ANPP = round(aNPPAct, 1),
                                         Mortality = round(mortality,1), deltaB, finBiomass = B)]

      tempcohortdata <- setkey(tempcohortdata, speciesCode)[
        setkey(sim$species[, .(species, speciesCode)], speciesCode),
        nomatch = 0][, ':='(speciesCode = species, species = NULL, pixelGroup = NULL)]
      setnames(tempcohortdata, "speciesCode", "Species")
      sim$simulationTreeOutput <- rbind(sim$simulationTreeOutput, tempcohortdata)
      set(subCohortData, NULL, c("deltaB", "sumB"), NULL)
    } else {
      set(subCohortData, NULL, "B",
          subCohortData$B + as.integer(subCohortData$aNPPAct - subCohortData$mortality))
    }
    sim$cohortData <- rbindlist(list(sim$cohortData, subCohortData))
    rm(subCohortData)
    gc() # TODO: use .gc()
  }
  rm(cohortData, cutpoints, pixelGroups)
  return(invisible(sim))
}

## other functions
updateSpeciesEcoregionAttributes <- function(speciesEcoregion, time, cohortData) {
  # the following codes were for updating cohortdata using speciesecoregion data at current simulation year
  # to assign maxB, maxANPP and maxB_eco to cohortData
  specieseco_current <- speciesEcoregion[year <= time]
  specieseco_current <- setkey(specieseco_current[year == max(year),
                                                  .(speciesCode, maxANPP,
                                                    maxB, ecoregionGroup)],
                               speciesCode, ecoregionGroup)
  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]

  cohortData <- setkey(cohortData, speciesCode, ecoregionGroup)[specieseco_current, nomatch = 0]
  return(cohortData)
}

updateSpeciesAttributes <- function(species, cohortData) {
  # to assign longevity, mortalityshape, growthcurve to cohortData
  species_temp <- setkey(species[, .(speciesCode, longevity, mortalityshape,
                                     growthcurve)], speciesCode)
  setkey(cohortData, speciesCode)
  cohortData <- cohortData[species_temp, nomatch = 0]
  return(cohortData)
}

calculateSumB <- function(cohortData, lastReg, simuTime, successionTimestep) {
  # this function is used to calculate total stand biomass that does not include the new cohorts
  # the new cohorts are defined as the age younger than simulation time step
  # reset sumB
  pixelGroups <- data.table(pixelGroupIndex = unique(cohortData$pixelGroup),
                            temID = 1:length(unique(cohortData$pixelGroup)))
  cutpoints <- sort(unique(c(seq(1, max(pixelGroups$temID), by = 10^4), max(pixelGroups$temID))))
  if (length(cutpoints) == 1) {cutpoints <- c(cutpoints, cutpoints + 1)}
  pixelGroups[, groups := cut(temID, breaks = cutpoints,
                              labels = paste("Group", 1:(length(cutpoints) - 1), sep = ""),
                              include.lowest = TRUE)]
  for (subgroup in paste("Group",  1:(length(cutpoints) - 1), sep = "")) {
    subCohortData <- cohortData[pixelGroup %in% pixelGroups[groups == subgroup, ]$pixelGroupIndex, ]
    set(subCohortData, NULL, "sumB", 0L)
    if (simuTime == lastReg + successionTimestep - 2) {
      sumBtable <- subCohortData[age > successionTimestep,
                                 .(tempsumB = as.integer(sum(B, na.rm=TRUE))), by = pixelGroup]
    } else {
      sumBtable <- subCohortData[age >= successionTimestep,
                                 .(tempsumB = as.integer(sum(B, na.rm=TRUE))), by = pixelGroup]
    }
    subCohortData <- merge(subCohortData, sumBtable, by = "pixelGroup", all.x = TRUE)
    subCohortData[is.na(tempsumB), tempsumB := as.integer(0L)][, ':='(sumB = tempsumB, tempsumB = NULL)]
    if (subgroup == "Group1") {
      newcohortData <- subCohortData
    } else {
      newcohortData <- rbindlist(list(newcohortData, subCohortData))
    }
    rm(subCohortData, sumBtable)
  }
  rm(cohortData, pixelGroups, cutpoints)
  return(newcohortData)
}

calculateAgeMortality <- function(cohortData, stage = "nonSpinup", spinupMortalityfraction) {
  # for age-related mortality calculation
  if (stage == "spinup") {
    cohortData[age > 0, mAge := B*(exp((age)/longevity*mortalityshape)/exp(mortalityshape))]
    cohortData[age > 0, mAge := mAge+B*spinupMortalityfraction]
    cohortData[age > 0, mAge := pmin(B, mAge)]
  } else {
    set(cohortData, NULL, "mAge",
        cohortData$B*(exp((cohortData$age)/cohortData$longevity*cohortData$mortalityshape)/exp(cohortData$mortalityshape)))
    set(cohortData, NULL, "mAge",
        pmin(cohortData$B,cohortData$mAge))
  }
  return(cohortData)
}

calculateANPP <- function(cohortData, stage = "nonSpinup") {
  if (stage == "spinup") {
    cohortData[age > 0, aNPPAct := maxANPP*exp(1)*(bAP^growthcurve)*exp(-(bAP^growthcurve))*bPM]
    cohortData[age > 0, aNPPAct := pmin(maxANPP*bPM,aNPPAct)]
  } else {
    set(cohortData, NULL, "aNPPAct",
        cohortData$maxANPP*exp(1)*(cohortData$bAP^cohortData$growthcurve)*exp(-(cohortData$bAP^cohortData$growthcurve))*cohortData$bPM)
    set(cohortData, NULL, "aNPPAct",
        pmin(cohortData$maxANPP*cohortData$bPM,cohortData$aNPPAct))
  }
  return(cohortData)
}

calculateGrowthMortality <- function(cohortData, stage = "nonSpinup") {
  if (stage == "spinup") {
    cohortData[age > 0 & bAP %>>% 1.0, mBio := maxANPP*bPM]
    cohortData[age > 0 & bAP %<=% 1.0, mBio := maxANPP*(2*bAP) / (1 + bAP)*bPM]
    cohortData[age > 0, mBio := pmin(B, mBio)]
    cohortData[age > 0, mBio := pmin(maxANPP*bPM, mBio)]
  } else {
    cohortData[bAP %>>% 1.0, mBio := maxANPP*bPM]
    cohortData[bAP %<=% 1.0, mBio := maxANPP*(2*bAP)/(1 + bAP)*bPM]
    set(cohortData, NULL, "mBio",
        pmin(cohortData$B, cohortData$mBio))
    set(cohortData, NULL, "mBio",
        pmin(cohortData$maxANPP*cohortData$bPM, cohortData$mBio))
  }
  return(cohortData)
}

calculateCompetition <- function(cohortData, stage = "nonSpinup") {
  # two competition indics are calculated bAP and bPM
  if (stage == "spinup") {
    cohortData[age > 0, bPot := pmax(1, maxB - sumB + B)]
    cohortData[age > 0, bAP := B/bPot]
    set(cohortData, NULL, "bPot", NULL)
    cohortData[, cMultiplier := pmax(as.numeric(B^0.95), 1)]
    cohortData[age > 0, cMultTotal := sum(cMultiplier), by = pixelGroup]
    cohortData[age > 0, bPM := cMultiplier / cMultTotal]
    set(cohortData, NULL, c("cMultiplier", "cMultTotal"), NULL)
  } else {
    set(cohortData, NULL, "bPot", pmax(1, cohortData$maxB - cohortData$sumB + cohortData$B))
    set(cohortData, NULL, "bAP", cohortData$B/cohortData$bPot)
    set(cohortData, NULL, "bPot", NULL)
    set(cohortData, NULL, "cMultiplier", pmax(as.numeric(cohortData$B^0.95), 1))
    cohortData[, cMultTotal := sum(cMultiplier), by = pixelGroup]
    set(cohortData, NULL, "bPM", cohortData$cMultiplier/cohortData$cMultTotal)
    set(cohortData, NULL, c("cMultiplier", "cMultTotal"), NULL)
  }
  return(cohortData)
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects", "function:spades")
  dPath <- asPath(dataPath(sim))

  # read species txt and convert it to data table
  if (!suppliedElsewhere("species", sim)) {
    mainInput <- prepInputsMainInput(url = NULL, dPath, cacheTags) ## uses default URL

    sim$species <- prepInputsSpecies(url = extractURL("species"), dPath, cacheTags)
  }

  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    sim$speciesEcoregion <- prepInputsSpeciesEcoregion(url = extractURL("speciesEcoregion"),
                                                       dPath = dPath, cacheTags = cacheTags)
  }

  return(invisible(sim))
}
