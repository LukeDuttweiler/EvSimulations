##### This file contains helper functions to support the 4 Simulation scenarios #####
##### Each function is documented #####

#' Generate Field of Vision (FOV) Data
#'
#' This function generates information available for each of the 31 Fields of Vision (FOVs)
#' by simulating extracellular vesicles (EVs) with certain protein binding profiles. The function
#' utilizes parameters such as the number of EVs, prevalence rates, binding rates, and others to
#' generate a list of binary matrices indicating the binding profiles for each FOV.
#'
#' @param particleCount A numeric value representing the total number of particles to be simulated across all FOVs. EVs represent about .15 percent of these particles. Default is \code{1e5}.
#' @param fovDist A numeric vector representing the distribution of EVs across the FOVs. The values should sum to 1. The number of FOVs desired is length(fovDist). Default is \code{rep(1, 31)/31}.
#' @param prevRates A numeric vector representing the prevalence rates of proteins. Number of different proteins is length(prevRates) Default is \code{rep(.5, 10)}.
#' @param bindRates A numeric vector representing the binding rates of antibodies. Number of antibodies is length(bindRates). length(bindRates) and length(prevRates) must be equal. Default is \code{rep(1, 10)}.
#' @param nonSpecBindRate A numeric value representing the non-specific binding rate of antibodies. Default is \code{.005}.
#' @param maxProteins A numeric value specifying the maximum number of proteins that can exist in a single EV. Default is \code{15}.
#' @param maxBinds A numeric value specifying the maximum number of bindings allowed in a single EV. Default is \code{15}
#' @param bloodDrawEV A numeric value representing the total number of EVs in the blood draw, from which evCount is a small portion that is experimented on. Defaults to 1e8. Really shouldn't try to set higher than this.
#' @param ... Catch unecessary arguments from scenario functions.
#'
#' @return A list of binary matrices, each corresponding to a FOV. In each matrix, rows represent EVs and columns represent the binding of antibodies.
#'
#' @export
genFOVs <- function(particleCount = 1e5,
                    fovDist = rep(1, 31)/31,
                    prevRates = rep(.5, 10),
                    bindRates = rep(1, length(prevRates)),
                    nonSpecBindRate = .005,
                    maxProteins = 10,
                    maxBinds = 10,
                    bloodDrawEV = 1e8,
                    ...){
  #Make sure fovDist sums to 1
  if(!all.equal(sum(fovDist), 1)){
    stop('fovDist must sum to 1')
  }

  #Make sure prevRates and bindRates are of the same length
  if(length(prevRates) != length(bindRates)){
    stop('prevRates and bindRates must have the same number of entries')
  }

  #Make sure enough particles are loaded
  if(particleCount <= 1e3){
    stop('not enough particles loaded to detect anything.')
  }

  #Get protein values for full blood draw
  bloodDrawMat <- sapply(prevRates, function(r){
    rbinom(bloodDrawEV, 1, r)
  })

  #Get evCount as a random count (from a scaled beta distribution) of the particles input into the system
  lowLim <- particleCount/1e3
  uppLim <- particleCount/1e2
  evCount <- (uppLim - lowLim)*rbeta(1, 2, 5) + lowLim

  #Get number of EVs per FOV
  evFOV <- ceiling(evCount*fovDist)

  #Get which EVs are in each FOV (sampling from bloodDrawMat rows without replacement)
  evWhich <- vector('list', length(evFOV))
  rowsUsed <- vector('integer')
  for(i in 1:length(evWhich)){
    rowSamp <- sample(bloodDrawEV, evFOV[i])
    stp <- 0 #Emergency stopper for while loop
    while(any(rowSamp) %in% rowsUsed & stp < 100){
      badRows <- which(rowSamp %in% rowsUsed)
      rowSamp[badRows] <- rowSamp[badRows] + 1
      stp <- stp + 1
    }
    rowsUsed <- c(rowsUsed, rowSamp)
    evWhich[[i]] <- rowSamp
  }

  #Per FOV create EVs with protein profiles
  fovs <- lapply(1:length(fovDist), function(i){
    #Number of EVs in this FOV
    n <- evFOV[i]

    #Binary matrix sampled from bloodDraw, rows are evs i, columns represent existence of protein j in ev i
    protMat <- bloodDrawMat[evWhich[[i]],]
    #Make sure number of proteins does not exceed max
    protMat <- t(apply(protMat, 1, rowShrink, maxSum = maxProteins))

    #Binary matrix, rows are evs i, columns represent possibility of binding antibody j to ev i
    possMat <- sapply(bindRates, function(r){
      rbinom(n, 1, r)
    })

    #Binary matrix, rows are evs i, columns represent specific binding of antibody j to ev i
    specMat <- protMat * possMat

    #Binary matrix, rows are evs i, columns represent possibility of non-specific binding of antibody j to i
    nonSpecMat <- sapply(rep(nonSpecBindRate, length(bindRates)), function(r){
      rbinom(n, 1, r)
    })

    #Binary matrix, rows are evs i, columns are binding of antibody j to ev i (regardless of truth)
    bindMat <- matrix(as.integer(as.logical(specMat + nonSpecMat)), nrow = n)
    #Ensure binding maximum is not surpassed anywhere
    bindMat <- t(apply(bindMat, 1, rowShrink, maxSum = maxBinds))

    return(bindMat)
  })

  rm(bloodDrawMat)
  gc()
  return(fovs)
}

#' Row Shrink Function
#'
#' This function takes a binary vector as input and reduces the number of `1`s in the vector at random positions
#' such that the sum of the vector does not exceed a specified maximum sum.
#'
#' @param binVec A binary vector (numeric vector with values 0 and 1).
#' @param maxSum A numeric value specifying the upper bound for the sum of the elements in `binVec`.
#'
#' @return A binary vector where the number of `1`s is reduced to meet the specified maximum sum, if necessary.
#' If the sum of `binVec` is already less than or equal to `maxSum`, the original vector is returned.
#'
rowShrink <- function(binVec, maxSum){
  currSum <- sum(binVec)
  sumDiff <- currSum - maxSum

  if(sumDiff <= 0){
    return(binVec)
  }else{
    whichOnes <- sample(which(binVec==1), sumDiff, replace = F)
    binVec[whichOnes] <- 0
    return(binVec)
  }
}

#Helper function from stackOverflow, catches all calls including defaults
match.call.defaults <- function(...) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))

  for(i in setdiff(names(formals), names(call)))
    call[i] <- list( formals[[i]] )


  match.call(sys.function(sys.parent()), call)
}
