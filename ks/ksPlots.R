library(ggplot2)

## functions to visualize results of runSim()
##
##

## resultDf: list resulting from runSim()
## 
meanSdPlot <- function(resultDf){
  g <- ggplot(aes(x = t, y=mean, ymin=mean-1.96*sd, ymax=mean+1.96*sd), data=resultDf$summaryStatsDf) +
    geom_ribbon(alpha=0.2) + geom_line(color='blue',size=1) +
      xlab("Generation") + ylab("p (mean +- 2 SD)") 
  return(g)
}

## resultDf: list resulting from runSim()
##
spreadPlot <- function(resultDf){
  g <- ggplot(aes(x = t, y=q50, ymin=q05, ymax=q95), data=resultDf$summaryStatsDf) +
    geom_ribbon(alpha=0.2) + geom_line(color='blue',size=1) +
      xlab("Generation") + ylab("p (5-95%)")
  return(g)
}

## resultDf: list resulting from runSim()
##
densityPlot <- function(resultDf, type='PDF'){
  
  ## threshold PDF for cleaner plot
  lowerPThresh <- 0.0001
  upperPThresh <- 1
  resultDf$densityDf$probDist2 <- resultDf$densityDf$probDist
  resultDf$densityDf[which(resultDf$densityDf$probDist<lowerPThresh), 'probDist2'] <- lowerPThresh
  resultDf$densityDf[which(resultDf$densityDf$probDist>upperPThresh), 'probDist2'] <- upperPThresh

  g <- ggplot(aes(x=t,y=p,z=cumDist), data=resultDf$densityDf)
  
  if(type=='CDF'){
    g <- g + geom_raster(aes(fill=cumDist))
  }
  else if(type=='PDF'){
    g <- g + geom_raster(aes(fill=probDist2))
  }
  else{
    stop()
  }
  g <- g + xlab("Generation (t)") + ylab("Degree of coarticulation (c)") +
     scale_fill_gradient2(name=expression(log(f[C^t](c)))) +
     guides(fill=guide_legend(title.theme = element_text(size=12, angle=0)))
    #scale_fill_continuous(low='blue',high='red',guide=guide_legend(title=expression(Prob(p))))
  return(g)
}



