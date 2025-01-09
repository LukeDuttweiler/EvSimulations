#' Run a Single Instance of SiMPore Scenario 1
#'
#' This function runs a simulation for SiMPore Scenario 1, which analyzes field of view (FOV) data, and tests for protein prevalence using a Bonferroni-adjusted alpha threshold for hypothesis testing.
#'
#' @inheritParams genFOVs
#' @param typeIRate Numeric value indicating the family-wise type one error rate desired. Defaults to .01.
#' @param assumedNonSpecBindRate Numeric value indicating the assumed non-specific binding rate used in the hypothesis test. Defaults to `nonSpecBindRate`.
#'
#' @return A list with two components:
#' \describe{
#'   \item{Results}{A data frame containing the true prevalence, observed prevalence, protein existence flag, test result flag, type 1 error flag, and type 2 error flag for each protein.}
#'   \item{Params}{A list of the parameters used for generating the simulation results.}
#' }
#'
#' @seealso \link{genFOVs}
#'
#' @export
scenario1 <- function(evCount = 1e5,
                      fovDist = rep(1, 31)/31,
                      prevRates = rep(.5, 10),
                      bindRates = rep(1, length(prevRates)),
                      nonSpecBindRate = .005,
                      maxProteins = 10,
                      maxBinds = 10,
                      bloodDrawEV = 1e8,
                      typeIRate = .01,
                      assumedNonSpecBindRate = nonSpecBindRate){
  #Get arguments, generate FOVs from those arguments
  argg <- as.list(match.call.defaults())
  #fovs <- do.call('genFOVs', argg)
  fovs <- genFOVs(evCount = evCount, fovDist = fovDist, prevRates = prevRates,
                  bindRates = bindRates, nonSpecBindRate = nonSpecBindRate,
                  maxProteins = maxProteins, maxBinds = maxBinds, bloodDrawEV = bloodDrawEV)

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
