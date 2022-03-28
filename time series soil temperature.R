require(zoo)
require(xts)
require(forecast)
palette(c("black","red3","green3","blue2","darkcyan","purple","sienna","gray67"))
# soiltemp <- read.table("clipboard", header = TRUE, sep = "\t",
#                        stringsAsFactors = FALSE)
# write.csv(soiltemp, file = "soiltemp.csv", row.names = FALSE)
soiltemp <- read.csv("soiltemp2.csv",
                     stringsAsFactors = FALSE)
#
# do some checks of the data
summary(soiltemp)
plot(soiltemp$Temp_15cm, type = "l", col = "sienna")
#
# convert data frame to time series object
soilT_zoo <- read.csv.zoo("soiltemp2.csv",
                        format = "%Y-%m-%d %H:%M:%S", 
                        tz = "Australia/Perth", 
                        index.column=1,
                        header = TRUE)
par(mar = c(4,4,1,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2)
plot(soilT_zoo, ylab = "Soil temperature at 15 cm (\u00B0C)",
     xlab = "Date (2014)", col = 7, lwd = 2, cex.lab = 1.4)
# note the difference in the x-axis!
#
# convert zoo object to xts (needed for some functions)
soilT_xts <- as.xts(soilT_zoo)
plot(soilT_xts, col = "chocolate") # just to check
#
# create a new time series from a moving average for each 24h
# to remove daily periodicity
soilT_movAv <- rollmean(soilT_zoo,24)
plot(soilT_movAv) # just to check
#
# inspect autocorrelation function acf()
plot(acf(soilT_xts))
#
# inspect partial autocorrelation function (pacf)
plot(pacf(soilT_xts))
#
# use Box test for autocorrelation
Box.test(soilT_xts)
#
# create a linear model of the time series...
# ...the residuals give us just the periodicity component
lm0 <- lm(coredata(soilT_xts) ~ index(soilT_xts))
#
# make a time series of the lm residuals
soilT_periodic <- zoo(resid(lm0), index(soilT_xts))
plot(soilT_movAv); abline(lm0, col = "coral") # just to check
#
# and plot everything to show the components
plot(cbind(soilT_zoo,soilT_periodic,soilT_movAv, soilT_lmfit),
  main = "Time series decomposition:\nSoil temperature at 15 cm (\u00B0C)", 
     cex.main = 1.5, yax.flip = TRUE, col = c(7,8,5,6))
#
# use the forecast:: package to run an ARIMA model
auto.arima(soilT_xts, max.p = 10, max.q = 10, seasonal = TRUE)
#
# use output of auto.arima() to run arima
am0 <- arima(x = soilT_xts, order = c(2, 1, 2),
       seasonal = list(order = c(0, 1, 1), period = 24))
summary(am0)
confint(am0)
checkresiduals(am0)
fc0 <- forecast(am0, h = 168)
require(ggplot2)
autoplot(fc0)+
  ylab("Soil temperature at 15 cm (\u00B0C)") +
  xlab("Time since 2014-04-01") +
  ggtitle("Forecasted soil temperature")+
  theme_bw()
ls(fc0)
summary(fc0$x)
