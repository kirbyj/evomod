##
## Functions to simulate evolution of a single-category Pierrehumbert (2001) exemplar model
##
## LI 308
## 7/2015
##

##
## EXAMPLES
## --------
##
## 0. First, run this file in R (either using source(), or Source File in the File window
## > source("pierrehumbertOne.R")
##
## 1. run a simulation with tau=2000, window width=0.05, amount of
## noise=0.1 (as in Pierrehumbert single-category simulations), no
## entrenchment, no lenition, for 10000 steps. Store in outputDf.
##
## > 
##


library(plyr)
library(ggplot2)

##
## do one-category pierrehumbert simulation
##
## arguments
## -------
## tau: memory decay time (positive integer; P: 2000)
##
## windowThresh: window threshold (positive number; P: 0.05)
##
## lambda: lenition bias (number)
##
## E: production noise from uniform dist between -E and E (positive number; P: 0.1)
##
## n: number of iterations (positive integer)
##
## entrenchment: do entrenchment? (TRUE or FALSE)
##
## ntrench: num of closest exemplars to consider (positive integer; P: 500)
##
## cutoff: num of most recent exemplars to keep. set to -1 to keep all
## exemplars, otherwise a positive integer (P: -1)
##
## initialValue: initial value of the 'phonetic character' for the category (number; P: 1)
##
## storeIvl: store category structure every storeIvl steps, for plotting distributions later (positive integer)
##
pierrehumbertOneCat <- function(tau=2000, windowThresh=0.05, lambda=0,
                                E=0.1, n=10000, entrenchment=FALSE, ntrench=500,
                                cutoff=10000, initialValue=0,
                                storeIvl=2000){
  ## check arguments
  try(stopifnot(is.numeric(tau) && tau>0))
  try(stopifnot(is.numeric(windowThresh) && windowThresh>0))
  try(stopifnot(is.numeric(lambda)))
  try(stopifnot(is.numeric(E) && E>0))
  try(stopifnot(is.numeric(n) && n>0))
  try(stopifnot(is.logical(entrenchment)))
  try(stopifnot(is.numeric(ntrench) && ntrench>0))
  try(stopifnot(is.numeric(cutoff) && (cutoff==-1 || cutoff>0)))
  try(stopifnot(is.numeric(initialValue)))
  try(stopifnot(is.numeric(storeIvl) && storeIvl>0))
  

  ## figure out which timesteps to store category structure at
  storeIts <- seq(from=storeIvl, to=n, by=storeIvl)

  if(!(n %in% storeIts)){ storeIts <- c(storeIts, n) }
  
  ## set up data frame where simulation results will be kept
  summaryDf <- data.frame(category = character(0), cue = numeric(0), strength = numeric(0), weight = numeric(0),  t = numeric(0))

  
  ## initialize category 
  category <- data.frame(cue = initialValue, T_i = 1, strength = 1)
  initialCategory <- category

  ## mean/var df
  meanVarDf <- data.frame(t = numeric(0), mean = numeric(0), var = numeric(0))
  
  ## agent talks to herself n times
  for(i in 2:n){
    if(i %% 100 == 0) cat(i, "\n")
    
    ## update strength of all exemplars currently in list
    category$strength <-  category$strength*exp(-1/tau)

    ## store mean and variance
    meanVarDf <- rbind(meanVarDf, data.frame(t=i,mean=weighted.mean(category$cue, category$strength), var=weighted.var(category$cue, category$strength)))
    
    ## if we're using a cutoff, check if there are >cutoff exemplars
    ## in each category, and discard oldest exemplar if so
    if(cutoff != -1){
      if(nrow(category) > cutoff){
        category <- category[2:nrow(category),]
        stopifnot(nrow(category)==cutoff)
        rownames(category) <- seq(1:cutoff)
      }
    }

    ## produce an exemplar
    x <- produceToken1Cat(category, p, entrenchment,ntrench, E, lambda)

    ## classify the exemplar
    score <- classifyToken1Cat(x, category,windowThresh)

    
    ## if score=NA, don't store utterance
    if(!is.na(score)){
      category <- rbind(category, data.frame(cue = x, T_i = i, strength=1))
    }
    
    ## store current category (all exemplars) in summaryDf, every so often
    if(i %in% storeIts){
      df <- category

      ## catStrength: relative activation among exemplars in *same* category
      df$catStrength <- df$strength/sum(df$strength)

      ## overallStrength: relative activation among *all* exemplars in
      ## all categories
      df$overallStrength <- df$strength / sum(df$strength) 

      df$T_i <- NULL
      df$t <- i
      
      summaryDf <- rbind(summaryDf, df)
    }
    
  }
  
  
  return(list(summaryDf=summaryDf, meanVarDf=meanVarDf))
  
}




##
## produce a token from the exemplar list 'category'
##
produceToken1Cat <- function(category, p, entrenchment,ntrench, E, lambda){
  
  ## sample 1 exemplar from the category, weighted by activation
  ## weight ('strength', here)
  ##
  ## (NB: there's a bug in sample() -- can't
  ## sample from a vector of length 1 if it's a positive integer)
  if(nrow(category)==1){
    x <- category$cue[1]
  }
  else{
    x <- sample(category$cue, 1, prob = category$strength)
  }

  ## if we're doing entrenchment
  if(entrenchment){
    ## memory-weighted distance of each exemplar
    memWeightedDist <- abs(x - category$cue)*(1/category$strength)

    ## number of exemplars to average over (all exemplars if there are
    ## fewer than ntrench)
    num <- min(length(memWeightedDist), ntrench)

    ## calculate memory-weighted average of these num exemplars
    inds <- order(memWeightedDist)[1:num]
    weights <-  category$strength[inds] 
    cues <- category$cue[inds]
    x <- weighted.mean(cues, weights)
  }
  
  ## add noise (from uniform distribution) and lenition bias (constant), if any
  x <- x + runif(1, -E, E) + lambda
  
  return(x)
  
}

##
## score for token 'x' for category 'category'
##
classifyToken1Cat <- function(x, category, windowThresh){
  ## 1 if *any* exemplar is within window of x, else NA (0, for pierrehumbert)
  score <- ifelse( length(which(abs(x - category$cue) < windowThresh))>0, 1, NA)

  if(is.na(score)){cat("score=0: discarding exemplar\n")}

  return(score)
}



##
## variance of elements of x, weighted by w
##
weighted.var <- function(x,w){
  wmean <- weighted.mean(x, w)
  as.numeric(w%*%(x - wmean)^2)/sum(w)
}

## df: output of pierrehumbertOneCat
##
## plots evolution of category distribution at times 'times'
##
## times: a vector of integers at which summaryDf was stored (Ex: c(5000, 25000, 75000, 100000)). If times = c() (default), just plots at the times for which category structure was stored.
##
## facet = TRUE : distribution at each time in a separate panel
## facet = FALSE : distributions at different times overplotted in a single panel
oneCatEvoPlot <- function(df, facet=TRUE, times=c()){
  summaryDf <- df$summaryDf
  if(length(times)>0){
    summaryDf <- subset(summaryDf, t %in% times)
  }
  hackBuffer <- 0.15*(max(summaryDf$cue) -  min(summaryDf$cue))
  if(facet){
    g <- ggplot(aes(x=cue), data=summaryDf) + stat_density(aes(weight=overallStrength, color=as.factor(t)),size=1, geom="line", position="identity") + facet_grid(~t) + theme(legend.position='none') + xlab("Phonetic character") + xlim(min(summaryDf$cue)-hackBuffer, max(summaryDf$cue) + hackBuffer)
  }
  else{
    g <- ggplot(aes(x=cue), data=summaryDf) + stat_density(aes(weight=overallStrength,color=as.factor(t)),size=1, geom="line", position="identity") + scale_color_discrete(name="steps") + xlab("Phonetic character") + xlim(min(summaryDf$cue), max(summaryDf$cue))
  }
  return(g)
}

## df: output of pierrehumbertOneCat
## plots (weighted) category variance vs time
##
oneCatVariancePlot <- function(df){
  #temp <- ddply(df, .(t), function(x){data.frame(var=weighted.var(x$cue,x$overallStrength))})
  qplot(t, var, data=df$meanVarDf, ylab="Variance",xlab="Steps", geom="line") + geom_line()
}

## df: output of pierrehumbertOneCat
## plots (weighted) category mean vs time
##
oneCatMeanPlot <- function(df){
  #temp <- ddply(df, .(t), function(x){data.frame(mean=weighted.mean(x$cue,x$overallStrength))})
  qplot(t, mean, data=df$meanVarDf, ylab="Mean",xlab="Steps", geom="line") + geom_line()
}


