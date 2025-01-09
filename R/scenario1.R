## This function runs a single instance of SiMPore Scenario 1
scenario1 <- function(evCount = 1e5,
                      fovDist = rep(1, 31)/31,
                      prevRates = rep(.5, 10),
                      bindRates = rep(1, length(prevRates)),
                      nonSpecBindRate = .005,
                      maxProteins = 10,
                      maxBinds = 10,
                      typeIRate = .01,
                      bloodDrawEV = 1e8,
                      assumedNonSpecBindRate = nonSpecBindRate){
  #Get arguments, generate FOVs from those arguments
  argg <- as.list(match.call.defaults())
  fovs <- do.call('genFOVs', argg)

  #For scenario 1 we don't care about FOVs specifically, so collapse into one matrix
  evs <- do.call('rbind', fovs)

  #Get percentage estimates of prevalence per protein
  protEsts <- colMeans(evs)

  #Get Bonferroni adjusted alpha (can change this if would help)
  alphBon <- typeIRate/length(prevRates)

  #Get test statistics for H_0:\hat{pi} = assumedNonSpecBindRate
  p <- assumedNonSpecBindRate
  zs <- (protEsts - p)/sqrt(p*(1-p)/nrow(evs))

  #Calculate test result, return results
  res <- data.frame('TruePrev' = prevRates, 'ObservedPrev' = protEsts, 'ProteinExists' = prevRates != 0,
                    'TestExists' = zs > -qnorm(alphBon))
  res$Type1Error <- !res$ProteinExists & res$TestExists
  res$Type2Error <- res$ProteinExists & !res$TestExists
  ret <- list('Results' = res, 'Params' = argg)

  return(ret)
}
