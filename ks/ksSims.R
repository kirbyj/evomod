require(plyr)
require(foreach)

#source("ksPlots.R")
 
##
## run simulations for journal ms
##
##
## arguments
## _________
## (defaults are shown in the function definition)
## 
## nGens: Number of generations
##
## nExamples: Number of examples received by each learner
##
## nLearners: Number of learners per generation
## 
## muA, sigA: parameters of /a/ normal distribution
## 
## muI, sigI: parameters of /i/ normal distribution
##
## teachers: how many teachers from previous generation? (character string; options: 'single', 'some', 'all')
##
## nTeachers: exact num of teachers from previous generation (only works with teachers = 'some', must be 1< __ < nLearners)
##
## initialDist, distParams: distribution of c in the initial generation (string, list)
##              - initialDist='gaussian' : normal dist with params pMuStart and pSigStart
##              - initialDist='custom' : initial c values supplied as a vector of length nLearners
##
## algorithm, algParams: learner's algorithm applied to data (A in the paper; string, list)
##              - simple : Anaive (no prior)
##              - gaussianPrior : MAP/EV using gaussian prior with SD param tau
##              - quadraticPrior : MAP using quadratic prior with
##                shape param a, lenition param lambda, SD omega,
##                lenition applied fraction biasProp (0-1) of
##                productions
##              
## storeIvl: store estimated PDF of c every storeIvl generations
##
## use.optimize: use optimize() to maximize posterior in quadPrior
## algorithm (TRUE), vs. using grid search (FALSE). optimize() faster,
## but only \approx guaranteed to find optimum of a *unimodal* function.
##
## stopTest : If TRUE, for a given run, after window or more
## generations, check every generation whether range of mean, 5th
## percentile, and 95th percentile are (all) less than 'thresh'. If
## so, abort this run.
##
## window : see 'stopTest'
##
## thresh : see 'stopTest'
##
## writeFile: write the data out as a file
##
## writeFinal: write just the params and the final results
##
## NOTE: can't both be true, since two .params files with different timestamps then written out.
##
## fPrefix: directory to write the file to (NB: need the trailing slash)
runSim <- function(nGens=500, nExamples = 100, nLearners = 1000,
                   muA = 730, sigA = 50, muI=530, sigI=50,
                   teachers = 'single', nTeachers = NA, initialDist = 'gaussian',
                   distParams=list(cMuStart = 10, cSigStart = 10),
                   algorithm = 'simple', algParams = list(lambda=0.0, omega=0.0, biasProp=0),
                   storeIvl = 5, use.optimize=TRUE,
                   stopTest = FALSE, window = 100, thresh = 2, 
                   writeFile = FALSE, writeFinal=FALSE, fPrefix = '/Users/jkirby/Desktop/runs/'){
  
  stopifnot(is.numeric(nLearners) && round(nLearners)==nLearners && nLearners>0,
            is.numeric(nExamples) && round(nExamples)==nExamples && nExamples>0,
            is.numeric(nGens) && round(nGens)==nGens && nGens>0,
            is.numeric(muA) && muA>0, is.numeric(sigA) && sigA>0,
            is.numeric(muI) && muI>0, is.numeric(sigI) && sigI>0,
            muA > muI,
            initialDist %in% c('gaussian', 'custom'),
            algorithm %in% c('simple', 'gaussianPrior', 'quadraticPrior'),
            teachers %in% c('single', 'some', 'all'),
            is.numeric(storeIvl) && round(storeIvl)==storeIvl && storeIvl>0,
            is.logical(writeFile),
            is.character(fPrefix)
            )

  if(writeFile && writeFinal){stop("writeFile and writeFinal can't both be TRUE")}

  if(teachers=='single'){
    nTeachers = 1
  }
  else if(teachers == 'all'){
    nTeachers = nLearners
  }
  else if(teachers == 'some'){
    stopifnot(is.numeric(nTeachers) && round(nTeachers)==nTeachers && nTeachers > 1 && nTeachers < nLearners)
  }
  else{
    stop("shouldn't get here")
  }

  if(algorithm == 'simple'){stopifnot('lambda' %in% names(algParams))}
  if(algorithm == 'gaussianPrior'){stopifnot('tau' %in% names(algParams), 'lambda' %in% names(algParams))}
  if(algorithm == 'quadraticPrior'){stopifnot('a' %in% names(algParams), 'lambda' %in% names(algParams), 'biasProp' %in% names(algParams), 'omega' %in% names(algParams))}
  
  ## generation 0 distribution over p

  ## gaussian distribution
  if(initialDist == 'gaussian'){
    stopifnot('cMuStart' %in% names(distParams), 'cSigStart' %in% names(distParams))

    cMuStart <- distParams$cMuStart
    cSigStart <- distParams$cSigStart
    teacherC <- rnorm(nLearners, cMuStart, cSigStart)
  }
  ## other 'distribution' (actually just list of c values for each
  ## member of G0, which you pre-generated from some dist)
  else if(initialDist == 'custom'){
    stopifnot('teacherC' %in% names(distParams))
    
    teacherC <- distParams$teacherC
  }

  ## estimated PDF over c in G0
  dEst <- density(teacherC, from=0, to=muA - muI)
  
  ##  which timesteps to storecp density at
  storeIts <- seq(from=storeIvl, to=nGens, by=storeIvl)
  if(!(nGens %in% storeIts)){ storeIts <- c(storeIts, n) }
  
  ## dataframe for same
  densityDf <- data.frame(t=0, p=dEst$x, probDist=dEst$y)

  ## dataframe to store summary stats about c distribution in each gen
  summaryStatsDf <- data.frame(t=0, mean=mean(teacherC), sd=sd(teacherC),
                               q05 = as.numeric(quantile(teacherC, 0.05)),
                               q25 = as.numeric(quantile(teacherC, 0.25)),
                               q50 = as.numeric(quantile(teacherC, 0.50)),
                               q75 = as.numeric(quantile(teacherC, 0.75)),
                               q95 = as.numeric(quantile(teacherC, 0.95)))
  
  ## dataframe to store simulation run arguments
  #if(length(algParams) > 0) {
      simParamsDf <- data.frame(nGens = nGens, nExamples = nExamples, nLearners = nLearners,
                                muA = muA, sigA = sigA, muI = muI, sigI = sigI,
                                teachers = teachers, nTeachers = nTeachers, initialDist = initialDist,
                                distParams = distParams, algorithm = algorithm,
                                algParams = algParams, storeIvl = storeIvl
                               )
  #} else {
  #      simParamsDf <- data.frame(nGens = nGens, nExamples = nExamples, nLearners = nLearners,
  #                              muA = muA, sigA = sigA, muI = muI, sigI = sigI,
  #                              teachers = teachers, nTeachers = nTeachers, initialDist = initialDist,
  #                              distParams = distParams, algorithm = algorithm,
  #                              storeIvl = storeIvl
  #                             )  
  #}
  print(summary(simParamsDf))
  ## create progress bar
  progBar <- txtProgressBar(min = 0, max = nGens, style = 3)


  ## run for nGens generations
  for(i in 1:nGens) {
    learners <- rep(0, nLearners)

    #learnerC <- foreach(j=1:nLearners, .combine=cbind) %do%{

    learnerC <- sapply(1:nLearners, function(x){
      ## sample F1 values from previous gen
      y <- sampleData(teachers, nTeachers, teacherC,  nExamples, muA, sigA)
     
      ## add bias (simple and gaussian cases; should make this consistent eventually)
      if (algorithm == 'simple' | algorithm == 'gaussianPrior') {
          ##  y <- y - algParams$lambda  WAS THIS: CHANGED FOR CONSISTENCY WITH PAPER 9/9/14
          y <- y- c(rnorm(as.integer(nExamples*algParams$biasProp), algParams$lambda, algParams$omega), rep(0,nExamples - as.integer(nExamples*algParams$biasProp)))
      }
 
      ## calculate p_hat estimate
      if(algorithm == 'simple'){
        muA - mean(y)
      }
      else if(algorithm == 'gaussianPrior'){
        (muA - mean(y))/(1 + (sigA^2)/(nExamples*(algParams$tau)^2))
      }
      else if(algorithm == 'quadraticPrior'){
        quadPriorAlg(y, nExamples,  muA, sigA, muI, algParams, use.optimize=use.optimize)
      }
    }
                       )
                       
    learnerC <- as.numeric(learnerC)
    dEst <- density(learnerC, from=0, to=(muA - muI))

    if(i %in% storeIts) {densityDf <- rbind(densityDf, data.frame(t=i, p=dEst$x, probDist=dEst$y))}
    
    summaryStatsDf <- rbind(summaryStatsDf, data.frame(t=i,mean=mean(learnerC), sd=sd(learnerC),
                                                       q05 = as.numeric(quantile(learnerC, 0.05)),
                                                       q25 = as.numeric(quantile(learnerC, 0.25)),
                                                       q50 = as.numeric(quantile(learnerC, 0.50)),
                                                       q75 = as.numeric(quantile(learnerC, 0.75)),
                                                       q95 = as.numeric(quantile(learnerC, 0.95))
                                                       )
                            )

    doStop = FALSE
    if(stopTest && i > window){
      ## check if mean, q05, q95 have remained stable
      tempDf <- summaryStatsDf[(i-window+1):i, ]

      ## what's the maximum of max - min over past 'window' runs of:
      ## mean, 5th percentile, 95th percentile.
      maxDiff <- with(tempDf, max(max(mean) - min(mean), max(q05) - min(q05), max(q95) - min(q95)))

      cat("maxDiff = ", maxDiff, "\n")

      ## if there's been little enough change, abort the simulation
      if(maxDiff < thresh){
        doStop = TRUE
      }
    }
    
    ## these learners are the teachers for the next gen
    teacherC <- learnerC
    
    ## update progress bar
    setTxtProgressBar(progBar, i)

    ## if we're aborting this simulation
    if(doStop){
      ## store the PDF, because i probably wasn't in storeIts
      densityDf <- rbind(densityDf, data.frame(t=i, p=dEst$x, probDist=dEst$y))
      cat("\nStopping run at t=", i, " : maxDiff = ", maxDiff, "\n", sep="")
      break
    }
  }
  

  ## calculate CDF for each stored PDF for c
  mySum <- function(x){y <- rep(0,length(x)); for(i in seq(1,length(x))){y[i] <- sum(x[1:i])}; return(y/max(y))}
  densityDf <- ddply(densityDf, .(t), function(x){data.frame(x, cumDist = mySum(x$probDist))})

  ## store sim parameters and also final values
  simParamsDf <- data.frame(nGens = nGens, nExamples = nExamples, nLearners = nLearners,
                            muA = muA, sigA = sigA, muI = muI, sigI = sigI,
                            teachers = teachers, nTeachers = nTeachers, initialDist = initialDist,
                            distParams = distParams, algorithm = algorithm,
                            algParams = algParams, storeIvl = storeIvl, 
                            finalMean = summaryStatsDf[nGens,'mean'], finalSd = summaryStatsDf[nGens,'sd'])
  
  ## store results
  resultList <- list(summaryStatsDf=summaryStatsDf, densityDf=densityDf, simParamsDf=simParamsDf)

  ## write a list containing summary stats and density dataframes to a
  ## compressed (.xz) file, which can be reloaded later
  if(writeFile || writeFinal){
    ## note that filename is uninformative: need to load .params file
    ## (this is different to how it was done in the original ksSims.R)
    body <- paste0(format(Sys.time(), "%Y-%m-%d-%H%M%S"), '-', round(runif(1,100,999),digits=0))
    paramsFile <- paste0(fPrefix, algorithm, 'Run-', body, '.params')
    
    ## write .param file
    write.table(resultList[[3]], file = paramsFile)
    
    ## write .xz file
    if(writeFile){
      resultsFile <- paste0(fPrefix, algorithm, 'Run-', body, '.xz')
      save(resultList, file=resultsFile, compress="xz")
    }
    
  }
  
  return(resultList)
}


## sample nExamples F1 values from 'teachers' agents' values of teacherC
sampleData <- function(teachers, nTeachers, teacherC, nExamples, muA, sigA){
  
  if(teachers == 'all'){
    exampleCs <- sample(teacherC, nExamples, replace=TRUE)
    y <- rnorm(nExamples, (muA - exampleCs), rep(sigA, nExamples))
  }
  else if(teachers == 'some'){
    teacherSubset <- sample(teacherC, nTeachers, replace=FALSE)
    #cat(teacherSubset, "\n")
    exampleCs <- sample(teacherSubset, nExamples, replace=TRUE)
    #cat(exampleCs, "\n")
    y <- rnorm(nExamples, (muA - exampleCs), rep(sigA, nExamples))
  }
  else{
    y <- rnorm(nExamples, muA - sample(teacherC,nTeachers), sigA)
  }
  
  return(y)
}

## take F1 values (y), infer p_hat for quadratic prior case. bias
## applied *here* (for convenience), even though technically it's part
## of production by teachers
##
quadPriorAlg <- function(y, nExamples,  muA, sigA, muI, algParams, use.optimize=TRUE){

  ## add bias (should eventually wrap this and complex prior case into one block of code)
    y <- y- c(rnorm(as.integer(nExamples*algParams$biasProp), algParams$lambda, algParams$omega), rep(0,nExamples - as.integer(nExamples*algParams$biasProp)))
    
  tempFun <- function(p){
      temp <- sum((y - muA + p)^2)  
    return(-temp/(2*sigA^2) + log(algParams$a + ((p/(muA - muI)) - 0.5)^2))
  }

  if(use.optimize){
    ## note: default tol on parmesan is 0.0001220703, so tol below is
    ## much more accurate.  nonetheless, optimize() only guaranteed to
    ## find the global max (+- epislon) if tempFun is unimodal, which
    ## I *think* is true for reasonable values of nExamples, sigA, and
    ## (muA - muI) (i.e., nExamples in hundreds, sigA/(muA - muI) < 1)

    optimize(tempFun, c(0, muA-muI), maximum=TRUE, tol=0.0000001)$maximum
  }
  else{
    ## calculate posterior
    pVec <- seq(0, (muA - muI), 1)
    lpVec <- vapply(pVec, FUN=function(p){sum((y - muA + p)^2)}, FUN.VALUE=0)
    lpVec <- -lpVec/(2*sigA^2) + log(algParams$a + ((pVec/(muA - muI)) - 0.5)^2)
    
    ## find MAP estimate of p_hat
    temp <- pVec[which(lpVec==max(lpVec))]
    if(length(temp)>1){warning("Non-unique MAP estimate")}
    
    pVec[which(lpVec==max(lpVec))][1]
  }   
}


## 10 gens:
## after change:
## 1.263   0.460   1.727

## before change:
## 1.388   0.185   1.582

## 100 gens:
## after: 13.037   4.706  17.963
##

## before: 13.664   1.580  15.360
##

## vapply instead of sapply: similar
