
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "LBMR_GMOrig",
  description = "insert module description here",
  keywords = c("insert key words here"),
  authors = person("First", "Last", email = "first.last@example.com", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("1.3.1.9035"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "LBMR_GMOrig.Rmd"),
  reqdPkgs = list(),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter(name = "growthInitialTime", class = "numeric", default = 0,
                    min = NA_real_, max = NA_real_,
                    desc = "Initial time for the growth event to occur"),
    defineParameter(name = ".plotInitialTime", class = "numeric", default = 0,
                    min = NA, max = NA,
                    desc = "This describes the simulation time at which the
                    first plot event should occur"),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = 0,
                    min = NA, max = NA,
                    desc = "This describes the simulation time at which the first save event should occur.
                    Set to NA if no saving is desired."),
    defineParameter("useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated?"),
    defineParameter("successionTimestep", "numeric", 10, NA, NA, "defines the simulation time step, default is 10 years"),
    defineParameter("calibrate", "logical", TRUE, NA, NA, "should the model have detailed outputs?")
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "cohortData", objectClass = "data.table",
                 desc = "tree-level data by pixel group",
                 sourceURL = NA),
    expectsInput(objectName = "lastReg", objectClass = "numeric",
                 desc = "time at last regeneration", sourceURL = NA),
    expectsInput(objectName = "species", objectClass = "data.table", 
                 desc = "species attribute table", sourceURL = NA),
    expectsInput(objectName = "speciesEcoregion", objectClass = "data.table",
                 desc = "species ecoregion data", sourceURL = NA)
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = "cohortData", objectClass = "data.table", 
                  desc = "tree-level data by pixel group")
  )
))

## event types
#   - type `init` is required for initialiazation
doEvent.LBMR_GMOrig = function(sim, eventTime, eventType, debug = FALSE) {
  if (is.numeric(sim$useParallel)) {
    a <- data.table::setDTthreads(P(sim)$useParallel)
    message("Mortality and Growth should be using >100% CPU")
    on.exit(setDTthreads(a))
  }
  switch(eventType,
         init = {
           ## do stuff for this event
           sim <- Init(sim)
           
           sim <- scheduleEvent(sim, start(sim) + P(sim)$growthInitialTime,
                                "LBMR_GMOrig", "mortalityAndGrowth", eventPriority = 5)
           },
         
         mortalityAndGrowth = {
           sim <- mortalityAndGrowth(sim)
           sim <- scheduleEvent(sim, time(sim) + 1, "LBMR_GMOrig", "mortalityAndGrowth",
                                eventPriority = 5)
           },
          warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                        "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
  return(invisible(sim))
}

### template for your event1
mortalityAndGrowth <- function(sim) {
  cohortData <- sim$cohortData
  sim$cohortData <- cohortData[0,]
  pixelGroups <- data.table(pixelGroupIndex = unique(cohortData$pixelGroup), 
                            temID = 1:length(unique(cohortData$pixelGroup)))
  cutpoints <- sort(unique(c(seq(1, max(pixelGroups$temID), by = 10^4), max(pixelGroups$temID))))
  if(length(cutpoints) == 1){cutpoints <- c(cutpoints, cutpoints+1)}
  pixelGroups[, groups:=cut(temID, breaks = cutpoints,
                            labels = paste("Group", 1:(length(cutpoints)-1),
                                           sep = ""),
                            include.lowest = T)]
  for(subgroup in paste("Group",  1:(length(cutpoints)-1), sep = "")){
    subCohortData <- cohortData[pixelGroup %in% pixelGroups[groups == subgroup, ]$pixelGroupIndex, ]
    #   cohortData <- sim$cohortData
    set(subCohortData, ,"age", subCohortData$age + 1)
    subCohortData <- updateSpeciesEcoregionAttributes_GMM(speciesEcoregion = sim$speciesEcoregion,
                                                          time = round(time(sim)), cohortData = subCohortData)
    subCohortData <- updateSpeciesAttributes_GMM(species = sim$species, cohortData = subCohortData)
    subCohortData <- calculateSumB_GMM(cohortData = subCohortData, 
                                       lastReg = sim$lastReg, 
                                       simuTime = time(sim),
                                       successionTimestep = P(sim)$successionTimestep)
    subCohortData <- subCohortData[age <= longevity,]
    subCohortData <- calculateAgeMortality_GMM(cohortData = subCohortData)
    set(subCohortData, , c("longevity", "mortalityshape"), NULL)
    subCohortData <- calculateCompetition_GMM(cohortData = subCohortData)
    if(!P(sim)$calibrate){
      set(subCohortData, , "sumB", NULL)
    }
    #### the below two lines of codes are to calculate actual ANPP
    subCohortData <- calculateANPP_GMM(cohortData = subCohortData)
    set(subCohortData, , "growthcurve", NULL)
    set(subCohortData, ,"aNPPAct",
        pmax(1, subCohortData$aNPPAct - subCohortData$mAge))
    subCohortData <- calculateGrowthMortality_GMM(cohortData = subCohortData)
    set(subCohortData, ,"mBio",
        pmax(0, subCohortData$mBio - subCohortData$mAge))
    set(subCohortData, ,"mBio",
        pmin(subCohortData$mBio, subCohortData$aNPPAct))
    set(subCohortData, ,"mortality",
        subCohortData$mBio + subCohortData$mAge)
    set(subCohortData, ,c("mBio", "mAge", "maxANPP",
                          "maxB", "maxB_eco", "bAP", "bPM"),
        NULL)
    if(P(sim)$calibrate){
      set(subCohortData, ,"deltaB",
          as.integer(subCohortData$aNPPAct - subCohortData$mortality))
      set(subCohortData, ,"B",
          subCohortData$B + subCohortData$deltaB)
      tempcohortdata <- subCohortData[,.(pixelGroup, Year = time(sim), siteBiomass = sumB, speciesCode,
                                         Age = age, iniBiomass = B - deltaB, ANPP = round(aNPPAct, 1),
                                         Mortality = round(mortality,1), deltaB, finBiomass = B)]
      
      tempcohortdata <- setkey(tempcohortdata, speciesCode)[setkey(sim$species[,.(species, speciesCode)],
                                                                   speciesCode),
                                                            nomatch = 0][, ':='(speciesCode = species,
                                                                                species = NULL,
                                                                                pixelGroup = NULL)]
      setnames(tempcohortdata, "speciesCode", "Species")
      sim$simulationTreeOutput <- rbind(sim$simulationTreeOutput, tempcohortdata)
      set(subCohortData, ,c("deltaB", "sumB"), NULL)
    } else {
      set(subCohortData, ,"B",
          subCohortData$B + as.integer(subCohortData$aNPPAct - subCohortData$mortality))
    }
    sim$cohortData <- rbindlist(list(sim$cohortData, subCohortData))
    rm(subCohortData)
    gc()
  }
  rm(cohortData, cutpoints, pixelGroups)
  return(invisible(sim))
  
}


.inputObjects = function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create an named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can use 'sim$.userSuppliedObjNames' in their function below to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call. e.g.,
  # if (!('defaultColor' %in% sim$.userSuppliedObjNames)) {
  #  sim$defaultColor <- 'red'
  # }
  # ! ----- EDIT BELOW ----- ! #
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}


updateSpeciesEcoregionAttributes_GMM <- function(speciesEcoregion, time, cohortData){
  # the following codes were for updating cohortdata using speciesecoregion data at current simulation year
  # to assign maxB, maxANPP and maxB_eco to cohortData
  specieseco_current <- speciesEcoregion[year <= time]
  specieseco_current <- setkey(specieseco_current[year == max(specieseco_current$year),
                                                  .(speciesCode, maxANPP,
                                                    maxB, ecoregionGroup)],
                               speciesCode, ecoregionGroup)
  specieseco_current[, maxB_eco:=max(maxB), by = ecoregionGroup]
  
  cohortData <- setkey(cohortData, speciesCode, ecoregionGroup)[specieseco_current, nomatch=0]
  return(cohortData)
}

updateSpeciesAttributes_GMM <- function(species, cohortData){
  # to assign longevity, mortalityshape, growthcurve to cohortData
  species_temp <- setkey(species[,.(speciesCode, longevity, mortalityshape,
                                    growthcurve)], speciesCode)
  setkey(cohortData, speciesCode)
  cohortData <- cohortData[species_temp, nomatch=0]
  return(cohortData)
}

calculateSumB_GMM <- function(cohortData, lastReg, simuTime, successionTimestep){
  # this function is used to calculate total stand biomass that does not include the new cohorts
  # the new cohorts are defined as the age younger than simulation time step
  # reset sumB
  pixelGroups <- data.table(pixelGroupIndex = unique(cohortData$pixelGroup), 
                            temID = 1:length(unique(cohortData$pixelGroup)))
  cutpoints <- sort(unique(c(seq(1, max(pixelGroups$temID), by = 10^4), max(pixelGroups$temID))))
  pixelGroups[, groups:=cut(temID, breaks = cutpoints,
                            labels = paste("Group", 1:(length(cutpoints)-1),
                                           sep = ""),
                            include.lowest = T)]
  for(subgroup in paste("Group",  1:(length(cutpoints)-1), sep = "")){
    subCohortData <- cohortData[pixelGroup %in% pixelGroups[groups == subgroup, ]$pixelGroupIndex, ]
    set(subCohortData, ,"sumB", 0L)
    if(simuTime == lastReg + successionTimestep - 2){
      sumBtable <- subCohortData[age > successionTimestep,
                                 .(tempsumB = as.integer(sum(B, na.rm=TRUE))), by = pixelGroup]
    } else {
      sumBtable <- subCohortData[age >= successionTimestep,
                                 .(tempsumB = as.integer(sum(B, na.rm=TRUE))), by = pixelGroup]
    }
    subCohortData <- merge(subCohortData, sumBtable, by = "pixelGroup", all.x = TRUE)
    subCohortData[is.na(tempsumB), tempsumB:=as.integer(0L)][,':='(sumB = tempsumB, tempsumB = NULL)]
    if(subgroup == "Group1"){
      newcohortData <- subCohortData
    } else {
      newcohortData <- rbindlist(list(newcohortData, subCohortData))
    }
    rm(subCohortData, sumBtable)
  }
  rm(cohortData, pixelGroups, cutpoints)
  gc()
  return(newcohortData)
}


calculateAgeMortality_GMM <- function(cohortData){
  set(cohortData, ,"mAge",
      cohortData$B*(exp((cohortData$age)/cohortData$longevity*cohortData$mortalityshape)/exp(cohortData$mortalityshape)))
  set(cohortData, ,"mAge",
      pmin(cohortData$B,cohortData$mAge))
  return(cohortData)
}

calculateANPP_GMM <- function(cohortData){
  set(cohortData, ,"aNPPAct",
      cohortData$maxANPP*exp(1)*(cohortData$bAP^cohortData$growthcurve)*exp(-(cohortData$bAP^cohortData$growthcurve))*cohortData$bPM)
  set(cohortData, ,"aNPPAct",
      pmin(cohortData$maxANPP*cohortData$bPM,cohortData$aNPPAct))
  return(cohortData)
}

calculateGrowthMortality_GMM <- function(cohortData){
  cohortData[bAP %>>% 1.0, mBio := maxANPP*bPM]
  cohortData[bAP %<=% 1.0, mBio := maxANPP*(2*bAP)/(1 + bAP)*bPM]
  set(cohortData, , "mBio",
      pmin(cohortData$B, cohortData$mBio))
  set(cohortData, , "mBio",
      pmin(cohortData$maxANPP*cohortData$bPM, cohortData$mBio))
  return(cohortData)
}

calculateCompetition_GMM <- function(cohortData){
  set(cohortData, , "bPot", pmax(1, cohortData$maxB - cohortData$sumB + cohortData$B))
  set(cohortData, , "bAP", cohortData$B/cohortData$bPot)
  set(cohortData, , "bPot", NULL)
  set(cohortData, , "cMultiplier", pmax(as.numeric(cohortData$B^0.95), 1))
  cohortData[, cMultTotal := sum(cMultiplier), by = pixelGroup]
  set(cohortData, , "bPM", cohortData$cMultiplier/cohortData$cMultTotal)
  set(cohortData, , c("cMultiplier", "cMultTotal"), NULL)
  return(cohortData)
}

