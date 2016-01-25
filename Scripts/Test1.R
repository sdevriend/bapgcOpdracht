dfairquality <- airquality
setwd("C:/Users/jesse/Documents/Bio-informatica/Jaar 3/Periode 2/Bapgc")
getwd()

save(dfairquality, file="dfairquality.RData")
load(dfairquality.RData)


write.csv2(dfairquality, file="dfairquality2.csv")

convertFahrenheitToCelsius <- function(Ftemp) {
  tempC <- (Ftemp-32)*(5/9)
  return(tempC)
}

dfairquality$tempC <- convertFahrenheitToCelsius(dfairquality$Temp)

sum(is.na(dfairquality$Solar.R))
sum(is.na(dfairquality$Wind))
sum(is.na(dfairquality$Ozone))
sum(is.na(dfairquality$Month))
sum(is.na(dfairquality$Day))

DATEPARSE('M', dfairquality$Month)
