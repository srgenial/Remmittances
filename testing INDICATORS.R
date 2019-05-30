#clean
rm(list = ls())
sessionInfo()

# Set the working directory
setwd("C:/Users/srgen/OneDrive/Documentos/Applied Economics/Econometrics/Econometrics II class/HW EcoII-2 Corona/Final")
getwd()

#Library setup
libs <- c("WDI", "forecast")
class(libs)
length(libs)
for (i in libs) {
  if(!is.element(i,.packages(all.available = TRUE))) {
    install.packages(i,repos="https://cran.revolutionanalytics.com/")
  }
  library(i,character.only = TRUE)
}

#library("WDI")
#World Bank Development Indicators for R  https://github.com/vincentarelbundock/WDI
#WDIsearch('gdp.*capita.*constant')
#import DataFrame
gdp<-WDI(indicator = 'NY.GDP.PCAP.KD', country=c('US'), start=1960, end=2016)
names(gdp)<-c("iso2c", "county", "GDPPerCap", "year")
head(gdp)
#reArrange
names(gdp)<-c("code", "C","GDPPC","yr")
gdp<-gdp[order(gdp$yr),]
head(gdp)

#IDENTIFICATION -Visualize
plot(GDPPC~yr,data=gdp, type="l", col=4, lwd=3)

#Dataframe to Time seriess
US<-ts(gdp$GDPPC, start=min(gdp$yr), end=max(gdp$yr))

#ESTIMATE ARIMA model
(USopt<-auto.arima(US))
#THe function doesnt work with MX data

#FORECASTING
GDPfore<-forecast(object = USopt,h=5)
plot(GDPfore)
