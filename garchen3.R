
#Jon Nielsen - 2016


garchen3 <- function(window, lastVol, parameters){
  
  c = parameters[1]
  g1 = parameters[2]
  g2 = parameters[3]
  
  mu = mean(window)
  
  v = c + g1*lastVol + g2*(window[1] - mu)^2
  
  return(v)
}