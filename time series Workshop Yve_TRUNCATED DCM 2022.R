#### Introduction to time series analysis for ENVT5503 ####
# 
# Load the packages we need for specialized time series
# functions in R
#
library(zoo)
library(xts) # we need the xts "eXtended Time Series" format for some functions
library(forecast)
library(tseries)
library(car)
#
# (optional) make a better colour palette than the R default!
palette(c("black","red3","green3","blue2",
          "darkcyan","purple","sienna","gray67"))
# 
#### data input ####
# read the data into a data frame - this is how R often stores
# data - it's not the format we need but we'll use it for comparison

# Yvelines <- read.csv("Yvelines_DCM.csv")
# colnames(Yvelines) <- c("Date","DCM")
#
# do some checks of the data
summary(Yvelines) # simple summary of each column
str(Yvelines) # more detailed information about the R object ('str'=structure)
plot(Yvelines$DCM, type = "l", col = 7)
#
# we really want the data in a different type of R object!
# we use the 'zoo' package to read the data in as a time series
# this is of class 'zoo' which is much more flexible than the
# default time series object in R

# Yvelines_DCM <- na.omit(Yvelines[,c("Collect.Date", "DCM")])
# write.csv(Yvelines_DCM, file="Yvelines_DCM.csv", row.names = F)
# rm(Yvelines_DCM)
Yvelines_DCM_zoo <- read.csv.zoo("Yvelines_DCM.csv",
                               format = "%d/%m/%Y", 
                               tz = "CET", 
                               index.column=1,
                               header = TRUE)
# run next comment line if log-transformed variable preferred!
# coredata(Yvelines_DCM_zoo) <- log10(coredata(Yvelines_DCM_zoo))

# note that we specified how the date was specified in our data file,
# and the time zone. Our time column was column 1, so we tell R this too.

# do some quick checks of the new R object:
summary(Yvelines_DCM_zoo) 
str(Yvelines_DCM_zoo) # POSIXct in the output is a date-time format

#### begin exploratory data analysis of time series ####
#
# try a plot of the data in our time series object...
# first we change the default plotting parameters using par(...)
par(mar = c(4,4,1,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2)
#
# then plot the object with custom axis labels
plot(Yvelines_DCM_zoo, type="l", ylab = "DCM (\u00B5g/L)",
     xlab = "Date", col = 7, lwd = 2, cex.lab = 1.4)
# note the difference in the x-axis from plot(Yvelines$DCM) !
#

pt0 <- powerTransform(coredata(Yvelines_DCM_zoo))
Yvelines_DCM_zoo <- Yvelines_DCM_zoo^pt0$roundlam

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plot differencing for stationarity
par(mfrow = c(3,1), mar = c(0,4,0,1), oma = c(4,0,1,0), cex.lab = 1.4)
plot(Yvelines_DCM_zoo, ylab = "transformed data, no differencing",
     xlab="", xaxt="n", col = 8)
lines(loess.smooth(index(Yvelines_DCM_zoo),coredata(Yvelines_DCM_zoo)), 
      col = 4, lwd = 2)
abline(lm(coredata(Yvelines_DCM_zoo)~index(Yvelines_DCM_zoo)), col = 2)
legend("topright", legend = c("Yvelines_DCM_Data", "Loess smoothing","Linear model"),
       cex = 1.8, col = c(1,4,2), lwd = c(1,2,1), bty = "n")
plot(diff(Yvelines_DCM_zoo,1),
     ylab = "First differencing",
     xlab="", xaxt="n", col = 8)
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(Yvelines_DCM_zoo,1)),coredata(diff(Yvelines_DCM_zoo,1))), 
      col = 4, lwd = 2)
abline(lm(coredata(diff(Yvelines_DCM_zoo,1))~index(diff(Yvelines_DCM_zoo,1))), col = 2)
plot(diff(diff(Yvelines_DCM_zoo,1),1),
     ylab = "Second differencing",
     xlab="Date", col = 8)
mtext("Date",side = 1, line = 2.2, font = 2)
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(diff(Yvelines_DCM_zoo,1),1)),
                   coredata(diff(diff(Yvelines_DCM_zoo,1),1))), 
      col = 4, lwd = 2)
abline(lm(coredata(diff(diff(Yvelines_DCM_zoo,1)))~index(diff(diff(Yvelines_DCM_zoo,1)))), col = 2)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# testing for stationarity
require(tseries) # needed for Augmented Dickey–Fuller (adf) Test
adf.test(Yvelines_DCM_zoo)
adf.test(diff(Yvelines_DCM_zoo,1))
adf.test(diff(diff(Yvelines_DCM_zoo,1),1))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# often we need to use another time series data format,
# xts (eXtended Time Series) which allows us to use more functions
#
# make an xts object from our zoo object:

Yvelines_DCM_xts <- as.xts(Yvelines_DCM_zoo)
plot(Yvelines_DCM_xts, col = 6, type = "b") # just to check
str(Yvelines_DCM_xts) # just to check

#### Finding the trend ####
## 1. using a moving average ___________________
# create a new time series from a moving average for each 12mo
# to remove daily periodicity using the zoo::rollmean function
Yvelines_DCM_movAv <- rollmean(Yvelines_DCM_zoo, 28)
plot(Yvelines_DCM_zoo, col=8) # original data
lines(Yvelines_DCM_movAv, lwd = 2) # just to check
#
# optionally add some other lines to this plot for comparison
lines(rollmean(Yvelines_DCM_zoo,6), col=2, lty = 2, lwd = 2) # 6 mo is incorrect periodicity
#
## 2. using a linear (regression) model ____________________
# create a linear model of the time series...
lm0 <- lm(coredata(Yvelines_DCM_zoo) ~ index(Yvelines_DCM_zoo))
summary(lm0)

# use a plot to look at the linear fit
plot(Yvelines_DCM_zoo, col = 8)
abline(lm0, col = 2)
#
#### isolating the periodicity ####
# To model the periodicity we need to understand the autocorrelation
# structure of the time series. First we can do this graphically:
#
# plot the autocorrelation function acf()
plot(acf(Yvelines_DCM_xts))
# what does this tell you about autocorrelation in this time series?
#
# plot the partial autocorrelation function (pacf)
plot(pacf(Yvelines_DCM_xts))
pacf(diff(Yvelines_DCM_xts,1)[-1])
pacf(diff(diff(Yvelines_DCM_xts,1),1)[-c(1,2)])
# interpreting partial autocorrelations is more complicated -
# refer to Long & Teetor (2019, Section 14.15). Simplistically,
# partial autocorrelation allows us to identify which and how many
# autocorrelations will be needed to model the time series data.
#
# use Box test for autocorrelation
# the null hypothesis is that no autocorrelation exists
# at any lag distance (so p <= 0.05 'rejects' null):
Box.test(Yvelines_DCM_xts)
Box.test(diff(Yvelines_DCM_xts,1))
Box.test(diff(diff(Yvelines_DCM_xts,1),1))
#
# remember we made a linear model of the time series...
# ...the residuals can give us just the periodicity component
lm0 <- lm(coredata(Yvelines_DCM_xts) ~ index(Yvelines_DCM_xts))
#
# make a time series of the lm residuals
Yvelines_DCM_periodic <- zoo(resid(lm0), index(Yvelines_DCM_xts))
plot(Yvelines_DCM_movAv); abline(lm0, col = "coral") # just to check
#
# ...but there is an argument that the trend is shown better by 
# a moving average than by a linear model, so we can subtract the
# moving average from the original data to get the periodicity:
Yvelines_DCM_periodic2 <- Yvelines_DCM_zoo - Yvelines_DCM_movAv
# plot, setting y axis limits to similar scale to original data:
plot(Yvelines_DCM_periodic2)
#
#### (optional) look at the unexplained variation ####
Yvelines_DCM_lmfit <- zoo(lm0$fitted.values, index(Yvelines_DCM_xts))

Yvelines_DCM_err = Yvelines_DCM_zoo - (Yvelines_DCM_movAv + Yvelines_DCM_periodic2)
plot(Yvelines_DCM_err, ylim = c(-200 ,200))
# 
# and plot everything to show the components
Yvelines_DCM_lmfit <- zoo(lm0$fitted.values, index(Yvelines_DCM_zoo))
plot(cbind(Yvelines_DCM_zoo,Yvelines_DCM_periodic2,
           Yvelines_DCM_movAv, Yvelines_DCM_err),
  main = "Time series decomposition:\nGroundwater DCM (mg/L)", 
     cex.main = 1.5, yax.flip = TRUE, col = c(7,8,5,6))
# notes: we could also include Yvelines_DCM_movAv
#        cbind is to combine columns;
#        \n inserts a line break into a text string;
#        \u followed by a 4-character code inserts a Unicode 
#           character (e.g. \u00B0 gives the degrees symbol);
#        yax.flip alternates vertical axis labels
#
#### end exploratory data analysis +=+=+=+=+=+=+=+=+=+=+=+
# ...which leads us into ARIMA forecast modelling ...
# -=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-
#
#### modelling time series with ARIMA ####

# use the forecast:: package to run an ARIMA model
auto.arima(Yvelines_DCM_xts, max.p = 3, max.q = 3, max.d = 1, 
           seasonal = F)
#
# use output of auto.arima() to run arima
#   1. with no seasonality
am0 <- arima(x = Yvelines_DCM_xts, order = c(1,1,0))
summary(am0)
confint(am0)
#   2. with seasonality
findfrequency(Yvelines_DCM_xts)
am1 <- arima(x = Yvelines_DCM_xts, order = c(1,1,0),
             seasonal = list(order = c(1, 1, 0), period = 30))
summary(am1)
confint(am1)
# ls(am1)
#
# checking residuals is our best diagnostic tool...
#   1. residual plot (top) should look like white noise
#   2. residuals should not be autocorrelated (bottom left plot)
#       AND p-value from Ljung-Box test (R console) should be 
#       GREATER THAN 0.05 
#   3. residuals should be normally distributed (bottom right plot)
checkresiduals(am0)
#
checkresiduals(am1)
# 
# use the ARIMA model to produce a forecast
fc0 <- forecast(am0, h = 25)
fc1 <- forecast(am1, h = 25)
#
par(mfrow = c(2, 1), cex.main = 0.9, mar = c(3,3,1,1),
    mgp=c(1.6,0.3,0), tcl=0.2, font.lab=2)
plot(fc0,
     ylab = expression(bold(sqrt("DCM, \u00B5g/L"))),
     xlab="Date/Time code",
     main = "Forecast for DCM (\u00B5g/L)")
# lines(index(Yvelines_DCM_zoo), coredata(Yvelines_DCM_zoo), col = "red3", cex=5)
# xvals <- as.numeric(index(Yvelines_DCM_xts))-as.numeric(index(Yvelines_DCM_xts)[1])
# lines(coredata(Yvelines_DCM_xts) ~ 
#        (as.numeric(index(Yvelines_DCM_xts))-as.numeric(index(Yvelines_DCM_xts)[1])), 
#      col = "red3")
mtext("ARIMA with no seasonality", 3, -1, adj = 0.9, font = 2)
mtext(am0$call, 3, -2.2, adj = 0.9, col = "blue3", cex = 0.9, 
      family="mono")

plot(fc1,ylab = expression(paste(log[10],"(DCM, \u00B5g/L)")),
     xlab="Date/Time code",
     main = "Forecast for DCM (\u00B5g/L)")
# lines(DCM_log, col = "grey70")
mtext(paste("ARIMA with periodicity =",am1$arma[5]), 
      side = 3, line = -1, adj = 0.8, font = 2)
mtext(am1$call, 3, -2.2, adj = 0.8, col = "blue3", cex = 0.9, family="mono")
par(mfrow = c(1, 1))
#
# we can also make a slightly 'prettier' plot using ggplot2
require(ggplot2) # gives best results using autoplot(...)
autoplot(fc1)+
  ylab("DCM (\u00B5g/L)") +
  xlab("Time since 2014-04-01") +
  ggtitle("Forecasted DCM concentration (\u00B5g/L)") +
  theme_bw()
#
# standard R plot
par(mar = c(4,4,2,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2, cex.lab=1.2)
plot(fc1, 
     ylab="DCM (\u00B5g/L)",
     xlab="Date", 
     main="Forecasted DCM concentration (\u00B5g/L)",
     xaxt="n")
date1 <- as.POSIXct("1993-01-01") # earliest date on axis
date2 <- as.POSIXct("1995-09-16") # latest date on axis
ndays <- 365.25 # between Date axis tick marks
axis(1, at = seq(as.numeric(date1), as.numeric(date2), 86400*ndays), 
     labels = as.Date.character(seq(date1, date2, 86400*ndays)) )
#
# _._._._._._._._._._._._._._._._._._._._._._._._._._._._
#    ___             _    _                          __  
#  .'   `.          / |_ (_)                        [  | 
# /  .-.  \ _ .--. `| |-'__   .--.   _ .--.   ,--.   | | 
# | |   | |[ '/'`\ \| | [  |/ .'`\ \[ `.-. | `'_\ :  | | 
# \  `-'  / | \__/ || |, | || \__. | | | | | // | |, | | 
#  `.___.'  | ;.__/ \__/[___]'.__.' [___||__]\'-;__/[___]
#           [__|                                                                              
#         _              ___    ___  
#        / |_          .' ..] .' ..]    |
#  .--. `| |-'__   _  _| |_  _| |_      |
# ( (`\] | | [  | | |'-| |-''-| |-'     |
#  `'.'. | |, | \_/ |, | |    | |     \ | /
# [\__) )\__/ '.__.'_/[___]  [___]     \|/
#
# sometimes ARIMA models may not be the best option
# another commonly used method is exponential smoothing
#
#### try an exponential smoothing model ####
# check help(forecast::ets) to correctly specify the model type!
# parameters for lower & upper are c("alpha","beta","gamma","phi")
Yvelines_DCM_ets <- ets(Yvelines_DCM_xts, model = "ZZZ")
summary(Yvelines_DCM_ets)
#
# default plot function can plot forecasts too...
date1 <- as.POSIXct("1993-12-10 CEST") # start date
date2 <- as.POSIXct("1995-09-16 CEST") # end date
plot(forecast(Yvelines_DCM_ets, h=12), col=7, xlab = "time since start")
# fac <- (as.numeric(date2) - as.numeric(date1))/(par("usr")[2] - par("usr")[1])
# lines(Yvelines_DCM_ets$fitted ~ 
#         seq(0,(as.numeric(date2) - as.numeric(date1))/140.1, 
#             length.out = length(Yvelines_DCM_ets$fitted)),
#       col = "#B0B0B080", lwd = 3)
mtext(names(Yvelines_DCM_ets$par),3,seq(-1,-5,-1), adj = 0.7, col = 4, font = 3)
mtext(signif(Yvelines_DCM_ets$par,3),3,seq(-1,-5,-1), adj = 0.8, col = 4, font = 3)
# ndays <- 365.25 # between Date axis tick marks
# axis(1, at = seq(as.numeric(date1), as.numeric(date2), 86400*ndays),
#      labels = as.Date.character(seq(date1, date2, 86400*ndays)),
#      cex.axis = 1)
# mtext("Date", side = 1, line = 1.8, font = 2)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))
rm(list = c("date1", "date2", "ndays"))

plot(Yvelines_DCM_ets) # decomposition plot
# (exponential smoothing is OK for the DCM data!)
#
# end code
# 
#### REFERENCES ####
#
# Hyndman R, Athanasopoulos G, Bergmeir C, Caceres G, Chhay L, O'Hara-Wild
# M, Petropoulos F, Razbash S, Wang E, Yasmeen F (2020). forecast:
# Forecasting functions for time series and linear models. R package
# version 8.12, http://pkg.robjhyndman.com/forecast.
#
# Long, J.D., Teetor, P., 2019. Time series analysis. Chapter 14, 
# The R Cookbook, Second Edition https://rc2e.com/timeseriesanalysis.
#
# R Core Team (2019). R: A language and environment for statistical
# computing. R Foundation for Statistical Computing, Vienna, Austria.
# URL https://www.R-project.org/.
#
# Ryan, J.A. and Ulrich, J.M. (2018). xts: eXtensible Time
# Series. R package version 0.11-2.
# https://CRAN.R-project.org/package=xts
#
# Zeileis, A. and Grothendieck, G. (2005). zoo: S3 Infrastructure
# for Regular and Irregular Time Series. Journal of Statistical
# Software, 14(6), 1-27. https://doi.org/10.18637/jss.v014.i06
#
