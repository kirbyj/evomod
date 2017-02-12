## Functions to simulate one-category and two-category Pierrehumbert (2001) exemplar models
##
## LI 511
## 6/2013
##

##
## EXAMPLES
## --------
##
## 0. First, run these files in R (either using source(), or Source File in the File window
source("pierrehumbertOne.R")
source("pierrehumbertTwo.R")

##
##
## SINGLE-CATEGORY
## ---------------
##
## 1. This code runs a simulation with tau=2000, window width=0.05, 
## noise=0.1 (as in Pierrehumbert single-category simulations), no
## entrenchment, no lenition, for 10000 steps. The results are stored
## in a data frame named outputDf.

## the value for storeIvl controls at what timestep a set of results 
## is stored. in this example, where n=10000 and storeIvl=2000, the
## memory state is stored at 5 timesteps: 2000, 4000, 6000, 8000, and
## 10000.

outputDf <- pierrehumbertOneCat(tau=2000, windowThresh=0.05, lambda=0, E=0.1, n=10000, entrenchment=FALSE, storeIvl=2000)

## plot the distribution of the phonetic character over time.
## setting facet=TRUE (the default) will plot each curve in its own
## window; setting facet=FALSE (as in the example below) produces 
## an overplot.

oneCatEvoPlot(outputDf, facet=F)

## the functions oneCatMeanPlot and oneCatVariancePlot are 
## wrappers to visualize the evolution of the (memory-weighted)
## mean and variance over time.
oneCatMeanPlot(outputDf)
oneCatVariancePlot(outputDf)


## 2. This code runs a simulation with tau=2000, window width=0.05, 
## noise=0.1 (as in Pierrehumbert single-category simulations), 
## lenition of 0.01, and entrenchment over 100 exemplars, for 10000 steps. 
## The results are stored in a data frame named outputDf2. Again, the
## full exemplar list is stored every 2000 timesteps. Note that in this
## example, where entrenchment=TRUE, an additional parameter must be 
## specified: ntrench, the number of exemplars closest to the production 
## target which are averaged over to produce a new target. Here, we
## set ntrench=100.

outputDf2 <- pierrehumbertOneCat(tau=2000, windowThresh=0.05, lambda=0.01, E=0.1, n=10000, entrenchment=TRUE, ntrench=100,  storeIvl=2000)

oneCatEvoPlot(outputDf2, facet=F)
oneCatMeanPlot(outputDf2)
oneCatVariancePlot(outputDf2)


## 3. A simulation with tau=100, window width=0.05, noise=0.1,
## lenition of 0.01, and entrenchment over 100 exemplars, for 10000 steps. 
## The results are stored in a data frame named outputDf3.

outputDf3 <- pierrehumbertOneCat(tau=100, windowThresh=0.05, lambda=0.01, E=0.1, n=10000, entrenchment=TRUE, ntrench=100,  storeIvl=2000)

oneCatEvoPlot(outputDf3, facet=F)
oneCatMeanPlot(outputDf3)
oneCatVariancePlot(outputDf3)


## 4. For comparison with 3: a simulation with tau=50, 
## window width=0.05, noise=0.1, lenition of 0.01, 
## and entrenchment over 100 exemplars, for 10000 steps. 
## The results are stored in a data frame named outputDf4.

outputDf4 <- pierrehumbertOneCat(tau=50, windowThresh=0.05, lambda=0.01, E=0.1, n=10000, entrenchment=TRUE, ntrench=100,  storeIvl=2000)

oneCatEvoPlot(outputDf4, facet=F)
oneCatMeanPlot(outputDf4)
oneCatVariancePlot(outputDf4)


## note that the output data frames (outputDf, outputDf2, etc.) are actually 
## lists containing more than one data frame. If you want to examine the 
## contents of the individual data frames, you can use the $ syntax, e.g.
summary(outputDf$summaryDf)
summary(outputDf$meanVarDf)


##
## TWO CATEGORIES 
## ---------------

## 5. This code runs a simulation for two categories with means at 0
## and 1 (initialValues = c(0,1)) and equal frequency (p=0.5), with
## tau=2000, window width=0.05, noise=0.1, no entrenchment, no
## lenition (lambda = list(c1=0, c2=0)), for 10000 steps, with the
## memory state stored at 2000, 4000, 6000, 8000, and 10000 steps
outputDf5 <- pierrehumbertTwoCat(tau=2000, windowThresh=0.05, n=10000, lambda=list(c1=0, c2=0), p=0.5, E=0.1, initialValues=c(0,1), entrenchment=FALSE, storeIvl=2000)

##
##

## plot the distribution of the phonetic character for each category
## at 2000, 4000, 6000, and 10000 steps.  setting facet=TRUE (the default, below) will plot one panel
## for each time point; setting facet=FALSE will overplot all
## distributions in a single panel
twoCatEvoPlot(outputDf5, times=c(2000, 4000, 6000, 10000))

## twoCatMeanDiffPlot plots the evolution of the difference in
## (memory-weighted) category means over time
twoCatMeanDiffPlot(outputDf5)

## twoCatMeanPlot and twoCatVariancePlot plot the evolution of the
## (weighted) mean and variance, respectively, of the categories over
## time
twoCatMeanPlot(outputDf5)
twoCatVariancePlot(outputDf5)


## 6. Replication of the neutralization example in Pierrehumbert
## (2001): unmarked category 3x more likely than the marked category
## (p=0.75), lenition = -0.1 for the marked category, entrenchment
## with ntrench=500, for 15000 steps, with the memory state stored
## every 2000 steps. All other parameters as in #5.
##

outputDf6 <- pierrehumbertTwoCat(tau=2000, windowThresh=0.05, n=15000, lambda=list(c1=0, c2=-0.1), p=0.75, E=0.1, initialValues=c(0,1), cutoff=10000, entrenchment=TRUE, ntrench=500, storeIvl=2000)

## plot the distribution of the phonetic character for each category
## at 2000, 4000, 10000, and 15000 steps.  
twoCatEvoPlot(outputDf6, times=c(2000, 4000, 10000, 15000))

twoCatMeanDiffPlot(outputDf6)
twoCatMeanPlot(outputDf6)
twoCatVariancePlot(outputDf6)


