# Introduction to time series analysis for ENVT5503
# UPDATED FOR 2021!
# 
# Load the packages we need for specialized time series
# functions in R
#
library(zoo)
library(xts)
library(forecast)
library(tseries)
#
# (optional) make a better colour palette than the R default!
palette(c("black","red3","green3","blue2",
          "darkcyan","purple","sienna","gray60"))
# co2_temp <- read.table("clipboard", header = TRUE, sep = "\t",
#                        stringsAsFactors = FALSE)
# write.csv(co2_temp, file = "co2_temp.csv", row.names = FALSE)
co2_temp <- read.csv("co2_temp.csv",
                     stringsAsFactors = FALSE)
#
# do some checks of the data
summary(co2_temp)
plot(co2_temp$CO2_ppm, type = "l", col = "sienna")
#
# convert data frame to time series object
co2_zoo <- read.csv.zoo("co2.csv",
                        format = "%Y-%m-%d", 
                        tz = "HST", 
                        index.column=1,
                        header = TRUE)
par(mar = c(4,4,1,1), mgp = c(1.7,0.3,0), tcl = 0.3, 
    font.lab=2, cex = 1.2)
plot(co2_zoo[630:749], 
     xlab = "Date", ylab = expression(paste(CO[2]," (ppm)")), 
     col = 6, cex.lab = 1.2, yax.flip = T, lwd = 2)
mtext(expression(paste("Last 10 years' atmospheric CO"[2],
                       " at Mauna Loa Observatory, Hawaii")),
      3,-1.5, family = "serif", font = 2, cex = 1.4, col = 6)
abline(lm0, lwd = 2, lty = 2)
# note the difference in the x-axis!
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plot differencing for stationarity
par(mfrow = c(3,1), mar = c(0,4,0,1), oma = c(4,0,1,0), 
    cex.lab = 1.4)
plot(co2_zoo, ylab = "Raw data, no differencing",
     xlab="", xaxt="n")
lines(loess.smooth(index(co2_zoo),coredata(co2_zoo)), 
      col = 4, lwd = 2)
abline(lm(coredata(co2_zoo)~index(co2_zoo)), col = 2)
legend("topleft", 
       legend = c("Data", "Loess smoothing","Linear model"),
       cex = 1.8, col = c(1,4,2), lwd = c(1,2,1), bty = "n")
plot(diff(co2_zoo,1),
     ylab = "First differencing",
     xlab="", xaxt="n")
rect(par("usr")[1], par("usr")[3], 
     par("usr")[2], par("usr")[4], 
     col = "#e8e8e8", border = 1)
lines(diff(co2_zoo,1), col = 1)
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(co2_zoo,1)),
                   coredata(diff(co2_zoo,1))), 
      col = 4, lwd = 2)
abline(lm(coredata(diff(co2_zoo,1))~index(diff(co2_zoo,1))), 
       col = 2)
plot(diff(diff(co2_zoo,1),1),
     ylab = "Second differencing",
     xlab="Date")
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(diff(co2_zoo,1),1)),
                   coredata(diff(diff(co2_zoo,1),1))), 
      col = 4, lwd = 2)
abline(lm(coredata(diff(diff(co2_zoo,1)))~
            index(diff(diff(co2_zoo,1)))), col = 2)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# testing for stationarity
# need tseries package for Augmented Dickey–Fuller (adf) Test
require(tseries) 
adf.test(co2_zoo)
adf.test(diff(co2_zoo,1))
adf.test(diff(diff(co2_zoo,1),1))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# convert zoo object to xts 
# (needed for some functions eg. findfrequency())
co2_xts <- as.xts(co2_zoo)
plot(co2_xts, col = 7) # just to check
#
# create a new time series from a moving average for 
# each 24h to remove daily fluctuations
co2_1yr_Av <- rollmean(co2_zoo,12)
plot(co2_1yr_Av, col = 5, lwd = 2) # just to check
#
# inspect autocorrelation function acf()
plot(acf(co2_xts))
#
# inspect partial autocorrelation function (pacf)
plot(pacf(co2_xts))
#
# use Box test for autocorrelation
Box.test(co2_xts)
#
# create a linear model of the time series...
# ...the residuals give us just the periodicity component
lm0 <- lm(coredata(co2_zoo) ~ index(co2_zoo))
#
# make a time series of the lm residuals
co2_periodic <- zoo(resid(lm0), index(co2_xts))
# just to check with a plot
plot(co2_zoo); abline(lm0, col = "coral", lwd = 2) 
#
# and plot everything to show the components
co2_lmfit <- zoo(lm0$fitted.values, index(co2_zoo))
plot(cbind(co2_zoo, co2_periodic, co2_1yr_Av, co2_lmfit),
  main = expression(paste("Time series of atmospheric CO"[2],
                      " at Mauna Loa Observatory, Hawaii")),
  xlab = "Date", lty = c(1,1,1,2), 
     cex.main = 1.5, yax.flip = TRUE, col = c(4,8,5,6))
#
# use the forecast:: package to run an ARIMA model
auto.arima(co2_xts, seasonal = T)
#
# use output of auto.arima() to run arima
am0 <- arima(x = co2_xts, order = c(2, 1, 2))
summary(am0)
confint(am0)
checkresiduals(am0)
#
# but we know the time series has seasonality, so...
# ...we add in a seasonality term too with 
# correct periodicity...
# ...and we can reduce the order (trial and error!)
#
# first check the periodicity...
periodicity(co2_xts)
findfrequency(co2_xts)
# ...then run a modified ARIMA model
am1 <- arima(x = co2_xts, order = c(0, 1, 1),
             seasonal = list(order = c(0, 1, 1), 
                             period = 12))
summary(am1)
confint(am1)
checkresiduals(am1)
fc0 <- forecast(am0, h = 60)
fc1 <- forecast(am1, h = 60)
#
par(mfrow = c(2, 1))
plot(fc0,ylab = expression(paste("Atmospheric CO"[2], 
                                 " at Mauna Loa")),
     xlab="Date/Time code",
     main = expression(bold("Forecast for atmospheric CO"[2])))
mtext("ARIMA with no seasonality", 3, -1, adj = 0.1, font = 2)
mtext(am0$call, 3, -2.2, adj = 0.1, col = "blue3", cex = 0.9, 
      family="mono")
plot(fc1,ylab = expression(paste("Atmospheric CO"[2], 
                                 " at Mauna Loa")),
     xlab="Date/Time code",
     main = expression(bold("Forecast for atmospheric CO"[2])))
mtext("ARIMA with annual seasonality", 3, -1, 
      adj = 0.1, font = 2)
mtext(am1$call, 3, -2.2, adj = 0.1, col = "blue3", 
      cex = 0.9, family="mono")
par(mfrow = c(1, 1))
#
# or we can use ggplot for a slightly 'prettier' plot
require(ggplot2)
autoplot(fc1)+
  ylab(expression(paste("Atmospheric CO"[2]," at Mauna Loa")))+
  xlab("Date/Time code") +
  ggtitle(expression(bold("Forecast for atmospheric CO"[2])))+
  theme_bw()

#### Exponential smoothing models ####

# Sometimes ARIMA models may not be the best option;
# another commonly used method is exponential smoothing.

# check help(forecast::ets) to correctly specify the model type!
co2_ets <- ets(co2_xts, model = "MMN")
summary(co2_ets)
#
# default plot function can plot forecasts too...
plot(forecast(co2_ets, h=28))
plot(co2_ets) # decomposition plot
# (exponential smoothing is not so good for the CO2 data!)
# _._._._._._._._._._._._._._._._._._._._._._._._._._._._
#
# plot the observed and predicted values, and the residuals 
# ideally the residuals represent the random variations
#   (but only if we have the 'perfect' time series model!)
plot(coredata(co2_xts) ~ seq(1,749), type = "l", lwd = 5, 
     col = 8, ylim = c(280,425),
     xlab = "Observation number", 
     ylab = expression(bold(paste("Atmospheric CO"[2], 
                                  " at Mauna Loa (ppm)"))))
lines(coredata(co2_xts)+am0$residuals ~ seq(1,749), col = 4)
abline(h=290, col = "coral2", lty = 3)
lines((am0$residuals*2)+290 ~ seq(1,749), col = 2)
legend("topleft", 
       legend = c("Observed","Fitted",
                  paste0("Residual (offset,\ntrue mean = ",
                         signif(mean(am1$resid),2),")")),
       col = c(8,4,2), lwd = c(5,1,1),
       bty = "n", cex = 1.2, inset = 0.01)
# _._._._._._._._._._._._._._._._._._._._._._._._._
#
# same for exponential smoothing model
plot(coredata(co2_xts) ~ seq(1,749), pch = 16, cex = 0.75, 
     col = 8, ylim = c(280,425),
     xlab = "Observation number", 
     ylab = expression(bold(paste("Atmospheric CO"[2], 
                                  " at Mauna Loa (ppm)"))))
sf = 100
lines(co2_ets$fitted ~ seq(1,749), col = 2)
lines((co2_ets$res*sf)+290 ~ seq(1,749), col = 4)
legend("topleft", 
       legend = c("Observed","Fitted",
           paste0("Residual \u00D7 ",sf,
                  " (offset,\ntrue mean = ",
                         signif(mean(co2_ets$resid),2),")")),
       col = c(8,2,4), lwd = c(NA,1,1), pch = c(16,NA,NA),
       pt.cex = 0.75, bty = "n", cex = 1.2, inset = 0.01)
rm(sf)

# end code