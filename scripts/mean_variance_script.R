#####
## plotMeanVariance.R
#####
## PLOT VARIANCE AS A FUNCTION OF MEAN TO ASSESS MEAN-VARIANCE RELATIONSHIP
#####
## GIVEN A MATRIX OF COUNTS counts WITH:
##   - ROWS PER FEATURE
##   - COLUMNS PER SAMPLE
## TREATMENT VECTOR treat FOR EACH SAMPLE
#####


plotMeanVariance <- function(counts, treat){
  
  require(edgeR)
  
  ## SCALE BY 75%ile
  p75 <- apply(counts, 2, function(x){ quantile(x, 0.75) })
  p75scaled <- p75 / mean(p75) * mean(colSums(counts))
  countsScaled <- t(apply(counts, 1, function(x){ x / p75scaled * mean(p75scaled) }))
  
  ## GET EDGER ESTIMATE OF THETA
  d <- DGEList(counts=counts, group=treat, lib.size=p75scaled)
  d <- estimateCommonDisp(d)
  theta <- 1/d$common.dispersion
  
  plotThis <- matrix(nrow=0, ncol=2)
  for(i in unique(treat)){
    tmpDat <- countsScaled[, treat == i]
    plotThis <- rbind(plotThis, cbind(rowMeans(tmpDat), apply(tmpDat, 1, var)))
  }
  
  ## FIND K FOR THE OD POISSON MODEL
  find.k <- lm(plotThis[, 2] ~ -1 + plotThis[, 1])
  k <- find.k$coefficients
  
  ## PLOT SQRT OF MEAN VERSUS SD
  plot(sqrt(as.numeric(plotThis[, 1])), sqrt(as.numeric(plotThis[, 2])), xlab = "sqrt(mean counts)", ylab = "standard deviation")
  sortx<-sort(as.numeric(plotThis[, 1]))
  y.nb<-sqrt(sortx+(sortx^2)/theta)
  lines(sqrt(sortx),y.nb,col="red",lwd=2)
  y.odp<-sqrt(sortx*k)
  lines(sqrt(sortx),y.odp,col="blue",lwd=2)
  y.p<-sqrt(sortx) 
  lines(sqrt(sortx),y.p,col="green",lwd=2)
  legend("topleft",legend = c("Poisson assumptions","over-dispersed Poisson (OD)","negativ binomial (NB) global Theta from edgeR"),col = c("green","blue","red"),lwd = 2,pch=NA)
  title("Mean-Variance plot to predict best variance model for Data")
  
}
