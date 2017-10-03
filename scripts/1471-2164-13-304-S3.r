#####
## fitNBedgeR.R
#####
## FUNCTION TO FIT A NEGATIVE BINOMIAL MODEL
## LOCAL ESTIMATES OF THETA
## SCALING TO 75%ILE OF COUNTS PER SAMPLE
## CREATE QQ PLOT OF PEARSON GOF STATISTICS
#####
fitNBedgeR <- function(counts, treat){
  
  require(edgeR)
  require(MASS)
  
  ## SCALE BY 75%ile
  p75 <- apply(counts, 2, function(x){ quantile(x, 0.75) })
  p75scaled <- p75 / mean(p75) * mean(colSums(counts))
  
  ## GET EDGE R ESTIMATE
  d <- DGEList(counts=counts, group=treat, lib.size=p75scaled)
  d <- estimateCommonDisp(d)
  ## LOOK FOR LOCAL ESTIMATES OF DISPERSION
  loc <- estimateTagwiseDisp(d, prior.n=3)
  
  theta <- 1/loc$tagwise.dispersion
  
  ## RUN THE MODELS PER GENE
  res <- lapply(as.list(1:nrow(counts)), function(x){
    try(glm(as.numeric(counts[x, ]) ~ factor(treat) + offset(log(p75scaled)), family=negative.binomial(theta[x])))
  })
  
  ## COMPUTE PEARSON RESIDUALS
  pearsonResid <- sapply(as.list(1:length(res)), function(this){
    x <- res[[this]]
    try(sum(((x$y-x$fitted.values)^2)/(x$fitted.values*(1+x$fitted.values/theta[this]))))
  })
  
  return(list(models = res,
              pearsonResid = pearsonResid))
  
}

## FUNCTION TO CREATE A QQ PLOT
myQQplot <- function(chis, df, ...){
  
  px <- qchisq(ppoints(chis), df)
  py <- sort(chis)
  
  plot(px, py, type="n", xlab="Theoretical Quantiles", ylab="Observed Quantiles", ...)
  points(px[1:round(length(px)*.95,0)], py[1:round(length(px)*.95,0)], col="black")
  points(px[(round(length(px)*.95,0)+1):round(length(px)*.99,0)], py[(round(length(px)*.95,0)+1):round(length(px)*.99,0)], col="blue")
  points(px[(round(length(px)*.99,0)+1):length(px)], py[(round(length(px)*.99,0)+1):length(px)], col="red")
  abline(0, 1)
  rug(quantile(px, 0:10/10), side=3)
  
}





#####
## EXAMPLE OF HOW TO EXTRACT INFORMATION FROM MODELS
## AND PRODUCE QQ PLOT
#####
## GIVEN A MATRIX OF COUNTS counts WITH:
##   - ROWS PER FEATURE
##   - COLUMNS PER SAMPLE
## TREATMENT VECTOR treat FOR EACH SAMPLE
#####

myOutput <- fitNBedgeR(counts, treat)
dfResid <- myOutput[["models"]][[1]]$df.residual

myQQplot(myOutput[["pearsonResid"]], dfResid, main = "My QQ Plot")

