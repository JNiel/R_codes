#Historical simulation - translated from .m-file.
#Jon Nielsen - 2016
#Import FX data
#Assuming all short positions

#READ DATA AS 'fxMat'
#Name the columns of fxMat


#Historical simulattion global variables

#CLEAR vialBhs violAwhs violVwhs

#Global variable 'N' is set in global_variables.R

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

windows = c(365, 750, 1095, 1460) #Window sizes
winMax = max(windows)

bhsVar = matrix(0, ncol=length(windows), nrow = (N-winMax))
bhsEs =  matrix(0, ncol=length(windows), nrow = (N-winMax))
violBhs = matrix(NA, nrow=(N-winMax), ncol=length(windows))#NA's 
diffBhs = matrix(NA, nrow=(N-winMax), ncol=length(windows))#NA's 
  
for (k in 1:length(windows)){
  winSize = windows[k]
  wVector = rep(1, winSize) / winSize
  
  for (i in 0:(N - winMax - 1)){
    
    window = fx[((N-winMax+1-i):(N-winMax+winSize-i)), cur]
    window = cbind(window, wVector)
    #Sort window with weights in descending order. 
    sortedWindow = window[order(window[,1], decreasing=TRUE),] 
    
    for (j in 1:dim(sortedWindow)[1]){
      if (sum(sortedWindow[1:j,2]) > (1-alpha/100)){
        bhsVar[(N - winMax - i), k] = sortedWindow[j]
        bhsEs[(N - winMax - i), k] = mean(sortedWindow[1:j-1])
        break
      }
    }
  }  
  violBhs[,k] = fx[1:(N-winMax), cur] > bhsVar[, k]
  diffBhs[,k] = fx[1:(N-winMax), cur] - bhsVar[, k]
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#--------------------Age weighted historical simulation---------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

lamb <- c(0.99, 0.995, 0.999)
violAwhs <- matrix(NA, ncol=(length(windows)*length(lamb)), nrow=(N-winMax))
diffAwhs <- matrix(NA, ncol=(length(windows)*length(lamb)), nrow=(N-winMax))
varAwhs <- matrix(NA, ncol=(length(windows)*length(lamb)), nrow=(N-winMax))
esAwhs <- matrix(NA, ncol=(length(windows)*length(lamb)), nrow=(N-winMax))


for (l in 1:length(lamb)){
  #Reset variables
  winSize = c()
  wVector = c()
  awhsVar <- matrix(NA, ncol=(N-winMax), nrow=length(windows))
  awhsEs <- matrix(NA, ncol=(N-winMax), nrow=length(windows))
  
  for (k in 1:length(windows)){
    
    winSize = windows[k]
    #set the first weight
    wVector[1] <- (1-lamb[l])/(1-lamb[l]^winSize)
    
    #Set rest of weight vector
    for (i in 2:winSize) {
      wVector[i] <- wVector[i-1]*lamb[l]
    }
    
    #THIS BIT IS EQUAL IN ALL STEPS? WRITE AS FUNCTION AND USE APPLY!!!
    #THIS BIT IS EQUAL IN ALL STEPS? WRITE AS FUNCTION AND USE APPLY!!!
    #THIS BIT IS EQUAL IN ALL STEPS? WRITE AS FUNCTION AND USE APPLY!!!
    
    #Outer loop for creating simulation vector
    for (i in 0:(N-winMax-1)){
      #Creating rolling windows with losses
      window = fx[((N-winMax+1-i):(N-winMax+winSize-i)), cur]
      window = cbind(window, wVector)
      #Sort window with weights in descending order.
      sortedWindow = window[order(window[,1], decreasing=TRUE),] 
      
      for (j in 1:dim(sortedWindow)[1]){
        if (sum(sortedWindow[1:j, 2])>(1-alpha/100)){
          awhsVar[k, (N-winMax-i)] = sortedWindow[j]
          awhsEs[k, (N-winMax-i)] =mean(sortedWindow[1:j])
          break
        }
      }
    }
   
    violAwhs[,((l-1)*4+k)] = fx[(1:(N-winMax)), cur] > t(awhsVar[k,])
    diffAwhs[,((l-1)*4+k)] = fx[(1:(N-winMax)), cur] - t(awhsVar[k,])
    
    varAwhs[,((l-1)*4+k)] = awhsVar[k,]
    esAwhs[,((l-1)*4+k)] = awhsEs[k,]
  }
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#--------------------Volatility weighted historical simulation--------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#Rescale losses

#update frequency of GARHC parameters is set in param_update.R (currently 30 days)

wVector <- c()
violVwhs <- matrix(NA, ncol=(length(windows)*4), nrow=(N-winMax))
diffVwhs <- matrix(NA, ncol=(length(windows)*4), nrow=(N-winMax))
varVwhs <- matrix(NA, ncol=(length(windows)*4), nrow=(N-winMax))
esVwhs <- matrix(NA, ncol=(length(windows)*4), nrow=(N-winMax))

vwhsVar <- matrix(NA, ncol=(N-winMax), nrow=length(windows))
vwhsEs <- matrix(NA, ncol=(N-winMax), nrow=length(windows))

for(reEstimate in 1:4){
  
  for (k in 1:length(windows)){
  
    winSize = windows[k]
    wVector = rep(1, winSize)/winSize
  
    for(i in 0:(N-winMax-1)){
      #Rolling window with losses
      window = fx[(N-winMax+1-i):(N-winMax+winSize-i), cur] 
      
      if(reEstimate != 4){
        g1 <- c(0.85, 0.9, 0.94)
        c = 0
        g1 = g1[reEstimate]
        g2 = 1 - g1
      }
      
      if(reEstimate == 4){
        c  = paramCheck[(i+1), 1]
        g1 = paramCheck[(i+1), 2]
        g2 = paramCheck[(i+1), 3]
      }
      
      #stdWindow = sqrt(garchen(fx[((N-winMax+1-i):N), cur], c, g1, g2))
      stdWindow = sqrt(garchen(window, c, g1, g2))
      
      #Rescaling losses:
      scaledWin = window*(stdWindow[1]/stdWindow[2:(length(window)+1)])
      
      scaledWin = cbind(scaledWin, wVector)
      sortedWindow = scaledWin[order(scaledWin[,1], decreasing=TRUE), ]
    
      for(j in 1:dim(sortedWindow)[1]){
        if (sum(sortedWindow[1:j, 2]) > (1-alpha/100)){
          vwhsVar[k, (N-winMax-i)] = sortedWindow[j]
          vwhsEs[k, (N-winMax-i)] = mean(sortedWindow[1:(j-1)])
          break
        }
      }
    }
    violVwhs[, ((reEstimate-1)*4+k)] = fx[(1:(N-winMax)), cur] > vwhsVar[k,]
    diffVwhs[, ((reEstimate-1)*4+k)] = fx[(1:(N-winMax)), cur] - vwhsVar[k,]
    
    varVwhs[, ((reEstimate-1)*4+k)] = vwhsVar[k,]
    esVwhs[, ((reEstimate-1)*4+k)] = vwhsEs[k,]
  }
}







