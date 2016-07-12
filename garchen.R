#Translation of .m-file "garchen2"
#Jon Nielsen - 2016

#Hur ska den här funktionen köras egentligen? Bör bytas ut mot en vectorized computation istället. 
#Vad är logiken bakom siffrorna? rita upp det. 
#I orginalskripten används v som vektor med ett steg längre än winL.


garchen <- function(window, c, g1, g2, bl=1){
  #If bl is set to 1 the function gives you the complete vector with variances for the whole window.
  
  m = mean(window)
  winL = length(window)
  v <- c()
  
  #Se över den här delen!
  
  if (bl == 1){
    v[winL+1] = var(window)
    
    for (i in 0:(winL-1)){
      v[length(v)-i-1] = c + g1 * v[length(v)-i] + g2 * (window[length(window)-i]-m)^2
    }
  }
  
  
  if (bl != 1) {
    v[1] <- var(window)
    
    for (i in 0:(winL-1)){
      v[1] = c + g1*v[1] + g2*(window[length(window)-i]-m)^2
    }
  }
  return(v)
}