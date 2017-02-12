library(plyr)
library(ggplot2)


## two-category pierrehumbert simulation
##
## arguments
## -------
## tau: memory decay time (positive integer; P: 2000)
##
## windowThresh: window threshold for pierrehumbert decision rule (positive number; P: 0.05)
##
## lambda: lenition bias for each category ( list, e.g. list(c1=0,c2=0) )
##
## E: production noise from uniform dist between -E and E (positive number; P: 0.1)
##
## n: number of iterations (positive integer)
##
## p: production probability p (num between 0 and 1)
##
## entrenchment: do entrenchment? (TRUE or FALSE)
##
## ntrench: num of closest exemplars to consider (positive integer; P: 500)
##
## cutoff: num of most recent exemplars in each category to keep. set to -1 to keep all
##
## exemplars, otherwise a positive integer (P: -1)
##
## initialValues: initial value of the 'phonetic character' for each category (vector, i.e. c(0,1))
##
## storeIvl: store category structure every storeIvl steps, for plotting later (positive int)



pierrehumbertTwoCat <- function(tau=2000, windowThresh=0.05, lambda=list(c1=0,c2=0),
                                E=0.1, n=10000, p=0.5, entrenchment=FALSE, ntrench=500,
                                cutoff=10000, initialValues=c(0,1),
                                storeIvl=2000){
  ## check arguments
  try(stopifnot(is.numeric(tau) && tau>0))
  try(stopifnot(is.numeric(windowThresh) && windowThresh>0))
  try(stopifnot(is.list(lambda) && all(as.logical(lapply(lambda, is.numeric)))))
  try(stopifnot(is.numeric(E) && E>0))
  try(stopifnot(is.numeric(n) && n>0))
  try(stopifnot(is.numeric(p) && p>=0 && p<=1))
  try(stopifnot(is.logical(entrenchment)))
  try(stopifnot(is.numeric(ntrench) && ntrench>0))
  try(stopifnot(is.numeric(cutoff) && (cutoff==-1 || cutoff>0)))
  try(stopifnot(is.vector(initialValues) && all(as.logical(lapply(initialValues, is.numeric)))))
  try(stopifnot(is.numeric(storeIvl) && storeIvl>0))
  

  ## figure out which timesteps to store category structure at
  storeIts <- seq(from=storeIvl, to=n, by=storeIvl)

  if(!(n %in% storeIts)){ storeIts <- c(storeIts, n) }
  
  ## set up data frame where simulation results will be kept
  summaryDf <- data.frame(category = character(0), cue = numeric(0), strength = numeric(0), weight = numeric(0),  t = numeric(0))
  
  ## categories start with single values. strength = exp(-(t - t_i)/tau), = 1 because t = t_i = 1
  categories <- list(c1 = data.frame(cue = initialValues[1], T_i = 1, strength = 1),
                     c2 = data.frame(cue = initialValues[2], T_i = 1, strength = 1))
  
  initialCategories <- categories

  ## mean/var df
  meanVarDf <- data.frame(t = numeric(0), mean = numeric(0), var = numeric(0))

  ## agent talks to herself n times
  for(i in 2:n){
    if(i %% 100 == 0) cat(i, "\n")
    
    ## update strength of all exemplars currently in list
    categories <- lapply(categories, function(x){x$strength <- x$strength*exp(-1/tau); return(x)})

    
    ## store mean and variance
    df <- ldply(categories, .fun=function(x){data.frame(mean = weighted.mean(x$cue,x$strength),
                                             var=weighted.var(x$cue, x$strength), t=i)})
    names(df)[1] <- 'category'

    meanVarDf <- rbind(meanVarDf, df)

    
    ## if we're using a cutoff, check if there are >cutoff exemplars
    ## in each category, and discard oldest exemplar if so
    if(cutoff != -1){
      for(C in names(categories)){
        category <- categories[[C]]
        if(nrow(category) > cutoff){
          category <- category[2:nrow(category),]
          stopifnot(nrow(category)==cutoff)
          rownames(category) <- seq(1:cutoff)
          categories[[C]] <- category
        }
      }
    }
    
    
    ## pick category to produce from
    C <- ifelse(runif(1) < p, 'c1', 'c2')
    category <- categories[[C]]

    ## produce an exemplar
    x <- produceToken2Cat(category, entrenchment, ntrench, E, lambda[[C]])

    ## classify the exemplar
    C <- classifyToken2Cat(x, categories,windowThresh)

    ## if there's a tie (C = NA), don't store the utterance
    if(!is.na(C)){
      categories[[C]] <- rbind(categories[[C]], data.frame(cue = x, T_i = i, strength=1))
    }

    ## store current category (all exemplars) in summaryDf, every so often
    if(i %in% storeIts){
      ## catStrength: relative activation among exemplars in *same* category
      df <- ldply(categories, .fun=function(x){x$catStrength <- x$strength/sum(x$strength); return(x)})
      

      ## overallStrength: relative activation among *all* exemplars in
      ## all categories
      df$overallStrength <- df$strength / sum(df$strength) 

      names(df)[1] <- 'category'
      df$T_i <- NULL
      df$t <- i
      
      summaryDf <- rbind(summaryDf, df)
    }
    
  }

  return(list(summaryDf=summaryDf, meanVarDf=meanVarDf))
  
}





## produce a token from this category at time t
produceToken2Cat <- function(category, entrenchment, ntrench, E, lenition){

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
    x <- as.numeric(cues %*% weights/sum(weights)) 
  }

  ## add noise (from uniform distribution) and lenition bias (constant), if any
  x <- x + runif(1, -E, E) + lenition
  
  return(x)
  
}

##
## input: token 'x' to be categorized, among 'categories', using window function with threshold 'windowThresh'
## output: label of category with highest score, or NA if scores tie (including 0 for both)
classifyToken2Cat <- function(x, categories,windowThresh){

  ## sum of activation of all exemplars within x's window
  pierrehumbertScore <- function(category){sum(ifelse(abs(x - category$cue) < windowThresh, 1, 0) * category$strength)}
  
  scores <- as.numeric(lapply(categories, pierrehumbertScore))

  if(scores[1] == scores[2]){
    return(NA)
  }
  else{
    return(names(categories)[which(scores==max(scores))])
  }
}

## df: output of pierrehumbertOneCat
##
## plots evolution of distribution of each category at times 'times'
##
## times: a vector of integers at which summaryDf was stored (Ex: c(5000, 25000, 75000, 100000))
## facet = TRUE : distribution at each time in a separate panel
## facet = FALSE : distributions at different times overplotted in a single panel
twoCatEvoPlot <- function(df, facet=TRUE, times=c()){
  summaryDf <- df$summaryDf
  if(length(times)>0){
    summaryDf <- subset(summaryDf, t %in% times)
  }
  if(facet){
    g <- ggplot(aes(x=cue), data=summaryDf) + geom_density(aes(color=category, fill=category, weight=overallStrength),size=0,alpha=0.7) + facet_wrap(~t) + xlab("Phonetic character")
  }
  else{
    g <- ggplot(aes(x=cue), data=summaryDf) + geom_density(aes(color=as.factor(t), lty=category, weight=overallStrength),size=1) + xlab("Phonetic character")
  }
  return(g)
}

## df: output of pierrehumbertTwoCat
##
## plots evolution of (weighted) mean of c1 minus (weighted) mean of c2
##
twoCatMeanDiffPlot <- function(df){
  temp <- ddply(df$meanVarDf, .(t), function(x){data.frame(diff=subset(x, category=='c1')$mean - subset(x, category=='c2')$mean)})
  g <- qplot(t, diff, data=temp, geom="line",xlab="Steps", ylab="c1 mean - c2 mean") + geom_line()
  return(g)
}


## df: output of pierrehumbertTwoCat
## plots (weighted) mean of each category vs time
##
twoCatMeanPlot <- function(df){
  qplot(t, mean, color=category, data=df$meanVarDf, ylab="Mean",xlab="Steps", geom="line") + geom_line()
}

## df: output of pierrehumbertTwoCat
## plots (weighted) variance of each category vs time
##
twoCatVariancePlot <- function(df){
  qplot(t, var, color=category, data=df$meanVarDf, ylab="Variance",xlab="Steps", geom="line") + geom_line()
}





