#Assign global variables used in functions.

N <- dim(fx)[1]
gpdUpdateFreq <- 30

#Choose currency rate:
#(1=EUR; 2=USD; 3=GBP; 4=DKK; 5=NOK; 6=AUD; 7=CHF)
cur <- 1 

#Confidence level
alpha <- 99

winMax <- 1460
win <- 1460
