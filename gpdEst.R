#gpdEst function definition
#Jon Nielsen - 2016

gpdEst <- function(data, percentileThreshold){
  threshold  <- quantile(data, percentileThreshold)
  tail  <- data[data>threshold]
  
  library(gPdtest)
  
  #Test to see if negative of positive shape parameters. Positive -> use "amle", else use "combined"
  gpd.test(tail) 
  fit <- gpd.fit(tail, "amle")
  
  fit[3] <- threshold
  
  return(fit)
  
}
