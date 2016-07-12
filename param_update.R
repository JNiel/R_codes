#Parameter updating
#Jon Nielsen 2016

garch_par_upd <- function(updateFreq=30, cur=1, winMax=1460){

  require(rugarch)

  paramCheck <- matrix(data = NA, nrow=(N-winMax+1), ncol=3)

  #GARCH(1,1) model parameter update. 
  Mdl = ugarchspec(variance.model = list(model='sGARCH', garchOrder=c(1,1)), distribution='std')
  for (i in 0:(N-winMax)) {
    if (i%%updateFreq == 0){
      #fx is now set as a matrix. Pick currency (row number)
      estMdl = ugarchfit(Mdl, fx[(N-winMax+1 - i) : N, cur], solver='hybrid')
      parameters = c(coef(estMdl)[4], coef(estMdl)[6], coef(estMdl)[5])
    }
    paramCheck[(i+1),] <- parameters
  }
  return(paramCheck)
}