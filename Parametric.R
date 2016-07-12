#Parametric
#Jon Nielsen - 2016

winMax = 1460

#----------------------------------------------------------------------------------------------------
#---------------------------Normal distribution - with GARCH-----------------------------------------
#----------------------------------------------------------------------------------------------------

violParamN <- matrix(NA, ncol=4, nrow=(N-winMax))
diffParamN <- matrix(NA, ncol=4, nrow=(N-winMax))
violParamT <- matrix()
paramN <- matrix(NA, ncol=4, nrow=(N-winMax+1))
paramN_Es <- matrix(NA, ncol=4, nrow=(N-winMax+1))
gStd <- rep(0, (N-winMax+1))
normCrit <- qnorm((alpha/100), 0, 1)
m <- rep(0, (N-winMax+1))

#Using GARCH otherwise it is EMWA
updateFreq = 30;
parameters = c(c, g1, g2)



for (i in 0:(N-winMax)){
  m[(length(m)-i)] = mean(fx[((N-winMax+1-i):(N-winMax+winSize-i)), cur])
}

for (k in 1:4){
    
  if (k==4){
    
    window = fx[((N-winMax+1):N) ,cur]
    c <- paramCheck[1,1]
    g1 <- paramCheck[1,2]
    g2 <- paramCheck[1,3]
    gStd[N-winMax+1] <- sqrt(garchen(window, c, g1, g2, 0))
    
    for(i in 1:(N-winMax-1)){
      
      window = fx[((N-winMax+1-i):(N-winMax+winSize-i)), cur]
      parameters = paramCheck[(i+1), ]
      
      gStd[length(gStd)-i] = sqrt(garchen3(window, (gStd[length(gStd)-i+1])^2, parameters))
      
    }
    
    dev = normCrit * gStd
    
    paramN[, k] = m + dev[1:(N-winMax+1)]
    
    paramN_Es[, k] = m + gStd[1:(N-winMax+1)] * pnorm(normCrit, 0, 1)/(1-alpha/100)
    
    violParamN[, k] = fx[(1:(N-winMax)), cur] > paramN[(2:dim(paramN)[1]), k]
    diffParamN[, k] = fx[(1:(N-winMax)), cur] - paramN[(2:dim(paramN)[1]), k]
    
  }
  
  if (k!=4) {
    
    g1 = c(0.85, 0.9, 0.95)
    c = 0
    g1 = g1[k]
    g2 = 1-g1
    window = fx[(N-winMax+1):N, cur]
    gStd[N-winMax+1] = sqrt(garchen(window, c, g1,g2, 0))
    
    for (i in 0:(N-winMax-1)){
      window = fx[(N-winMax+1-i):N, cur]
      #Calculating std for time series using EWMA param.
      gStd[length(gStd)-i+1] = sqrt(garchen(window, c, g1, g2, 0)) 
    }
    
    dev = normCrit * gStd
    paramN[, k] = m + dev[1:(N-winMax+1)]
    paramN_Es[, k] = m + gStd[1:(N-winMax+1)] * pnorm(normCrit, 0, 1)/(1-alpha/100)
    
    violParamN[, k] = fx[(1:(N-winMax)), cur] > paramN[(2:dim(paramN)[1]), k]
    diffParamN[, k] = fx[(1:(N-winMax)), cur] - paramN[(2:dim(paramN)[1]), k]
    
  }
}

#Remove estimate for tomorrow
paramN = paramN[2:dim(paramN)[1], ]
paramN_Es = paramN_Es[2:dim(paramN_Es)[1], ]


#----------------------------------------------------------------------------------------------------
#---------------------------Student-t distribution - with GARCH--------------------------------------
#----------------------------------------------------------------------------------------------------

m = rep(0, (N-winMax+1))
kurt = rep(0, (N-winMax+1))
rm(dgf)

violParamT = matrix(NA, ncol = 4, nrow=(N-winMax))
diffParamT = matrix(NA, ncol = 4, nrow=(N-winMax))
paramT = matrix(NA, ncol = 4, nrow=(N-winMax+1))
paramT_Es = matrix(NA, ncol = 4, nrow=(N-winMax+1))

gStd <- c()
require(moments)

for (i in 0:(N-winMax)){
  m[length(m)-i] = mean(fx[((N-winMax+1-i):(N-winMax+winSize-i)), cur])
  kurt[length(kurt)-i] = kurtosis(fx[((N-winMax+1-i):(N-winMax+winSize-i)), cur])
}

#Degrees of freedom vector
dgf = (4*kurt-6)/(kurt-3)

for (k in 1:4){
      
  if (k==4){
    
    window = fx[(N-winMax+1-i):N, cur]
    
    c = paramCheck[1,1]
    g1 = paramCheck[1,2]
    g2 = paramCheck[1,3]
    
    gStd[N-winMax+1] = sqrt(garchen(window, c, g1, g2,0))
    
    for (i in 1:(N-winMax-1)){
      
      window = fx[(N-winMax+1-i):N, cur]
      
      parameters = paramCheck[i+1, ]
      
      gStd[length(gStd)-i] = sqrt(garchen3(window, gStd[length(gStd)-i+1]^2, parameters))
    }
    
    tCrit = qt(alpha/100, dgf)
    devT = sqrt((dgf-2)/dgf)*tCrit*gStd[1:length(tCrit)]
    
    paramT[, k] = m + devT[1:(N-winMax+1)]
    paramT_Es[, k] = m + sqrt((dgf-2)/dgf) * gStd * pt(tCrit, dgf)/(1-alpha/100) * (dgf+tCrit^2)/(dgf-1)
    
  }
  if (k!=4){
    
    g1 = c(0.85, 0.9, 0.95)
    c = 0
    g1 = g1[k]
    g2 = 1-g1
    window = fx[(N-winMax+1):N, cur]
    gStd[(N-winMax+1)] = sqrt(garchen(window, c, g1, g2, 0))
    
    for (i in 0:(N-winMax-1)){
      
      window = fx[(N-winMax+1-i):N, cur]
      gStd[length(gStd)-i-1] = sqrt(garchen(window, c, g1, g2, 0))
      
    }
    
    tCrit = qt(alpha/100, dgf)
    devT = sqrt((dgf-2)/dgf) * tCrit * gStd[1:length(gStd)]
    
    paramT[, k] = m + devT[1:(N-winMax+1)]
    paramT_Es[, k] = m + sqrt((dgf-2)/dgf) * gStd[1:length(tCrit)] * pt(tCrit, dgf)/(1-alpha/100) * (dgf+tCrit^2)/(dgf-1)
    
  }


  violParamT[, k] = fx[(1:(N-winMax)), cur] > paramT[2:dim(paramT)[1], k]
  diffParamT[, k] = fx[(1:(N-winMax)), cur] - paramT[2:dim(paramT)[1], k]
}


