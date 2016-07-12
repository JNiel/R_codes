#EVT.m script
#Jon Nielsen - 2016

#KONTROLLERA KODEN, KOMBINERAR FUNKTIONER MED MATRISINDEX

#Extreme value theory
#Parametric updating conditional
#Semi-parametric

#Parameter estimation settings

#gpdUpdateFreq <- 30

#alpha <- 99
#winMax <- 1460
#win <- 1460

evt_run <- function(potU=1, potC=1, paramCheck, gpdUpdateFreq=30, alpha=99, winMax=1460, win=1460, tailPercentileV=c(0.96, 0.98), emwa_par=c(0.85, 0.9, 0.95)){


  violPotUandC <-  matrix(data=NA, nrow = (N-winMax), ncol=10)
  violPotU <- matrix(data=NA, nrow = (N-winMax), ncol=2)
  diffPotU <- matrix(data=NA, nrow = (N-winMax), ncol=2)
  violPotC <- matrix(data=NA, nrow = (N-winMax), ncol=8)
  diffPotC <- matrix(data=NA, nrow = (N-winMax), ncol=8)
  potUVar <- matrix(data=NA, nrow = 2, ncol = (N-winMax))
  potCVar <- matrix(data=NA, nrow = 8, ncol = (N-winMax))
  potUEs <- matrix(data=NA, nrow = 2, ncol = (N-winMax))
  potCEs <- matrix(data=NA, nrow = 8, ncol = (N-winMax))
  
  #Defining what percentile the tail should be made of
  #tailPercentileV <- c(0.96, 0.98) 
  
  #---------------------------------------------------------------------------------
  #Initiate outer loop. 
  #--------------------------------------------------------------------------------- 
  
  for(k in 1:length(tailPercentileV)){
    tailPercentile <- tailPercentileV[k]
    
    #---------------------------------------------------------------------------------
    #POT - unconditional with updating parameters
    #---------------------------------------------------------------------------------
    if(potU==1){
      for(i in 0:(N-winMax)){
        
        window <- fx[(N-winMax+1-i):N, cur]
    
        if(i%%gpdUpdateFreq == 0){
          #Estimating parameters and threshold for every loop
           v <- gpdEst(window, tailPercentile) 
           shape <- v[1]
           scale <- v[2]
           u <- v[3]
        }
        
        beta <- scale
        xi <- shape
        M <- length(window)
        Nu <- sum(window>u)
        
        potUVar[k, (N - winMax - i)] <- u + beta/xi * ((M/Nu*(1-alpha/100))^xi-1)
        
        potUEs[k, (N - winMax - i)] <- potUVar[k, (N - winMax - i)]/(1-xi) + (beta-xi*u)/(1-xi)
      }
      
      violPotU[, k] <- fx[1:(N-winMax), cur]>t(potUVar[k, 1:dim(potUVar)[2]])
      diffPotU[, k] <- fx[1:(N-winMax), cur] - t(potUVar[k, 1:dim(potUVar)[2]])
      percentage <- sum(violPotU[, k]) / (N - winMax)
  }
  
  
  
  #---------------------------------------------------------------------------------
  #POT - conditional with parameter updating
  #---------------------------------------------------------------------------------
  if(potC==1){  
    quantiles <- rep(0, N-winMax + 1)
    dev <- rep(0, N - winMax + 1)
      
    reEstimate <- 0 #paramter updating toggle
    
    #start values when GARCH(1, 1)
    for (reEstimate in 1:4){
            
      for (i in 0:(N -winMax -1)){
        
        #Using increasing window
        window <- fx[(N - winMax + 1 - i):N, cur] 
        
        #Updating parameters for every step
        
        if(reEstimate==4){
          c <- paramCheck[(i+1), 1]
          g1 <- paramCheck[(i+1), 2]
          g2 <- paramCheck[(i+1), 3]
        }
        
        if(reEstimate<4){
          #g1 <- c(0.85, 0.9, 0.95)
          c <- 0
          g1 <- emwa_par[reEstimate]
          g2 <- 1-g1
        }
        
        #Creating vector with time series volatility per each loss in sample and for t+1
        stdWindow <- sqrt(garchen(window, c, g1, g2, 1)) 
        
        #rescaling with volatility
        scaledWindow <- (window - mean(window))/t(stdWindow[1:length(stdWindow)]) 
        
        if (i==0 || i%%gpdUpdateFreq==0){
          #estimating parameters and threshold according to update freq
          v <- gpdEst(scaledWindow, tailPercentile) 
          shape <- v[1]
          scale <- v[2]
          u <- v[3]
          
          beta <- scale
          xi <- shape
          
          #estimating parameters and threshold for every loop
          w <- gpdEst(window, tailPercentile) 
          wshape <- w[1]
          wscale <- w[2]
          wu <- w[3]
          
          wbeta <- wscale
          wxi <- wshape
        }
        
        M <- length(scaledWindow)
        Nu <- sum(scaledWindow>u)
        
        quantiles[(N-winMax-i+1)] <- u + beta/xi * ((M/Nu*(1-alpha/100))^-xi - 1)
        
        potCVar[((k-1)*4+reEstimate), (N-winMax-i)] <- mean(window) + quantiles[(N-winMax-i+1)]*stdWindow[1]
        
        potCEs[((k-1)*4+reEstimate), (N-winMax-i)] <- potCVar[((k-1)*4+reEstimate), (N-winMax-i)]/(1-wxi) + (wbeta - wu*wxi)/(1-wxi)
      }
    
      violPotC[, ((k-1)*4+reEstimate)] <- fx[(1:(N-winMax)), cur]>t(potCVar[((k-1)*4+reEstimate), 1:dim(potCVar)[2]])
      diffPotC[, ((k-1)*4+reEstimate)] <- fx[(1:(N-winMax)), cur]-t(potCVar[((k-1)*4+reEstimate), 1:dim(potCVar)[2]])
      
      percentage <- sum(violPotC) / (N-winMax)
    }
  }
}

#This step should be in the loop.
violPotUandC[, 1:2] <- violPotU
violPotUandC[, 3:10] <- violPotC

sum(violPotC)

diffPotUandC <- c(violPotU, violPotC)


}


#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#-------------------------------DU ÄR HÄR-----------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

