#### Introduction to time series analysis for ENVT5503 ####
# 
# Load the packages we need for specialized time series
# functions in R
#
library(zoo)
library(xts) # we'll include xts just so you know about it
library(forecast)
#
# (optional) make a better colour palette than the R default!
palette(c("black","red3","green3","blue2",
          "darkcyan","purple","sienna","gray67"))
# 
#### data input ####
# read the data into a data frame - this is how R often stores
# data - it's not the format we need but we'll use it for comparison

soiltemp <- read.csv("soiltemp2.csv",
                     stringsAsFactors = FALSE)
#
# do some checks of the data
summary(soiltemp) # simple summary of each column
str(soiltemp) # more detailed information about the R object ('str'=structure)
plot(soiltemp$Temp_15cm, type = "l", col = 7)
#
# we really want the data in a different type of R object!
# we use the 'zoo' package to read the data in as a time series
# this is of class 'zoo' which is much more flexible than the
# default time series object in R

soilT_zoo <- read.csv.zoo("soiltemp2.csv",
                        format = "%Y-%m-%d %H:%M:%S", 
                        tz = "Australia/Perth", 
                        index.column=1,
                        header = TRUE)
# note that we specified how the date was specified in our data file,
# and the time zone. Our time column was column 1, so we tell R this too.

# do some quick checks of the new R object:
summary(soilT_zoo) 
str(soilT_zoo) # POSIXct in the output is a date-time format

#### begin exploratory data analysis of time series ####
#
# try a plot of the data in our time series object...
# first we change the default plotting parameters using par(...)
par(mar = c(4,4,1,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2)
#
# then plot the object with custom axis labels
plot(soilT_zoo, ylab = "Soil temperature at 15 cm (\u00B0C)",
     xlab = "Date (2014)", col = 7, lwd = 2, cex.lab = 1.4)
# note the difference in the x-axis!
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plot differencing for stationarity
par(mfrow = c(3,1), mar = c(0,4,0,1), oma = c(4,0,1,0), cex.lab = 1.4)
plot(soilT_zoo, ylab = "Raw data, no differencing",
     xlab="", xaxt="n", col = 8)
lines(loess.smooth(index(soilT_zoo),coredata(soilT_zoo),
                   span = 0.2), 
      col = 4, lwd = 2, lty = 2)
abline(lm(coredata(soilT_zoo)~index(soilT_zoo)), col = 2)
legend("topright", cex = 1.5, bty = "n", 
       legend = c("Data", "Loess smoothing","Linear model"),
       lty = c(1,2,1), col = c(1,4,2), lwd = c(1,2,1))
plot(diff(soilT_zoo,1),
     ylab = "First differencing",
     xlab="", xaxt="n", col = 8)
rect(par("usr")[1], par("usr")[3], 
     par("usr")[2], par("usr")[4], 
     col = "#e8e8e8", border = 1)
lines(diff(soilT_zoo,1), col = "grey50")
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(soilT_zoo,1)),
                   coredata(diff(soilT_zoo,1)),
                   span = 0.2), 
      col = 4, lwd = 2, lty = 2)
abline(lm(coredata(diff(soilT_zoo,1))~index(diff(soilT_zoo,1))), col = 2)
plot(diff(diff(soilT_zoo,1),1),
     ylab = "Second differencing",
     xlab="Date", col = 8)
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(diff(soilT_zoo,1),1)),
                   coredata(diff(diff(soilT_zoo,1),1)),
                   span = 0.2), 
      col = 4, lwd = 2, lty = 2)
abline(lm(coredata(diff(diff(soilT_zoo,1)))~index(diff(diff(soilT_zoo,1)))), col = 2)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# testing for stationarity
require(tseries) # needed for Augmented Dickey–Fuller (adf) Test
adf.test(soilT_zoo)
adf.test(diff(soilT_zoo,1))
adf.test(diff(diff(soilT_zoo,1),1))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# sometimes we need to use another time series data format,
# xts (eXtended Time Series) which allows us to use more functions
#
# make an xts object from our zoo object:

soilT_xts <- as.xts(soilT_zoo)
summary(soilT_xts)

plot(soilT_xts, col = 6, lwd = 1) # just to check
#
#### Finding the trend ####
## 1. using a moving average ___________________
# create a new time series from a moving average for each 24h
# to remove daily periodicity using the zoo::rollmean function
periodicity(soilT_xts)
findfrequency(soilT_xts)

soilT_movAv <- rollmean(soilT_zoo, 24)
plot(soilT_movAv, lwd = 2) # just to check
#
# optionally add some other lines to this plot for comparison
lines(soilT_zoo, col=8) # original data
lines(rollmean(soilT_zoo, 12), col=2) # 12 h is incorrect periodicity
#
## 2. using a linear (regression) model ____________________
# create a linear model of the time series...
lm0 <- lm(coredata(soilT_zoo) ~ index(soilT_zoo))

# use a plot to look at the linear fit
plot(soilT_zoo, col = 8)
abline(lm0, col = 2)


#### isolating the periodicity ####

# To model the periodicity we need to understand the 
# autocorrelation structure of the time series. 
# First we can do this graphically:
#
# plot the autocorrelation function acf()
plot(acf(soilT_xts))
# what does this tell you about autocorrelation in 
# this time series?
#
# plot the partial autocorrelation function (pacf)
plot(pacf(soilT_xts))
# interpreting partial autocorrelations is more complicated -
# refer to Long & Teetor (2019, Section 14.15). Simplistically,
# partial autocorrelation allows us to identify which and
# how many autocorrelations will be needed to model the 
# time series data.
#
# use Box test for autocorrelation
# the null hypothesis is that no autocorrelation exists
# at any lag distance (so p <= 0.05 'rejects' null):
Box.test(soilT_xts)
#
# remember we made a linear model of the time series...
# ...the residuals can give us just the periodicity component
lm0 <- lm(coredata(soilT_xts) ~ index(soilT_xts))
#
# make a time series of the lm residuals
soilT_periodic <- zoo(resid(lm0), index(soilT_xts))
plot(soilT_movAv); abline(lm0, col = "coral") # just to check
#
# ...but there is an argument that the trend is shown better 
# by a moving average than by a linear model, so we can
# subtract the moving average from the original data to get 
# the periodicity:
soilT_periodic2 <- soilT_zoo - soilT_movAv
# plot, setting y axis limits to similar scale to original data:
plot(soilT_periodic2, ylim = c(-6 ,6))
#
#### (optional) look at the unexplained variation ####
soilT_err = soilT_zoo - (soilT_lmfit + soilT_periodic2)
plot(soilT_err, ylim = c(-6 ,6))
# 
# and plot everything to show the components
soilT_lmfit <- zoo(lm0$fitted.values, index(soilT_zoo))
plot(cbind(soilT_zoo,soilT_periodic2,soilT_movAv, soilT_err),
  main = "Time series decomposition:\nSoil temperature at 15 cm (\u00B0C)", 
     ylim = c(-3,22), cex.main = 1.5, yax.flip = TRUE, col = c(7,8,5,6))
# notes: we could also include soilT_lmfit
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
auto.arima(soilT_zoo, max.p = 5, max.q = 5, seasonal = TRUE)
#
# use output of auto.arima() to run arima
#   1. with no seasonality
am0 <- arima(x = soilT_zoo, order = c(2,1,2))
summary(am0)
confint(am0)
#   2. with seasonality
findfrequency(soilT_xts)
am1 <- arima(x = soilT_zoo, order = c(1,1,1),
             seasonal = list(order = c(0, 1, 1), period = 24))
summary(am1)
confint(am1)
#
# checking residuals is our best diagnostic tool...
#   1. residual plot (top) should look like white noise
#   2. residuals should not be autocorrelated (bottom left plot)
#       AND p-value from Ljung-Box test (R console) should be 
#       GREATER THAN 0.05 
#   3. residuals should be normally distributed 
#     (bottom right plot)
checkresiduals(am0)
#
checkresiduals(am1)
# 
# use the ARIMA model to produce a forecast
fc0 <- forecast(am0, h = 168)
fc1 <- forecast(am1, h = 168)
#
par(mfrow = c(2, 1), cex.main = 0.9)
plot(fc0,ylab = "Soil temperature at 15 cm (\u00B0C)",
     xlab="Date/Time code",
     main = "Forecast for Soil temperature at 15 cm (\u00B0C)")
mtext("ARIMA with no seasonality", 3, -1, adj = 0.9, font = 2)
mtext(am0$call, 3, -2.2, adj = 0.9, col = "blue3", 
      cex = 0.9, family="mono")
plot(fc1,ylab = "Soil temperature at 15 cm (\u00B0C)",
     xlab="Date/Time code",
     main = "Forecast for Soil temperature at 15 cm (\u00B0C)")
mtext("ARIMA with 24 hour periodicity", 3, -1, 
      adj = 0.9, font = 2)
mtext(am1$call, 3, -2.2, adj = 0.9, col = "blue3", 
      cex = 0.9, family="mono")
par(mfrow = c(1, 1))
#
# we can also make a slightly 'prettier' plot using ggplot2
require(ggplot2) # gives best results using autoplot(...)
autoplot(fc1)+
  ylab("Soil temperature at 15 cm (\u00B0C)") +
  xlab("Time since 2014-04-01") +
  ggtitle("Forecasted soil temperature") +
  theme_bw()
#
# standard R plot
par(mar = c(4,4,2,1), mgp = c(1.7,0.3,0), tcl = 0.3, 
    font.lab=2, cex.lab=1.2)
plot(fc1, 
     ylab="Soil temperature at 15 cm (\u00B0C)",
     xlab="Date", 
     main="Forecasted soil temperature",
     xaxt="n")
date1 <- as.POSIXct("2014-04-01") # earliest date on axis
date2 <- as.POSIXct("2014-06-01") # latest date on axis
ndays <- 10 # between Date axis tick marks
axis(1, 
     at = seq(as.numeric(date1),as.numeric(date2),86400*ndays), 
     labels = as.Date.character(seq(date1, date2, 86400*ndays)))

#### Exponential smoothing models ####

# Sometimes ARIMA models may not be the best option;
# another commonly used method is exponential smoothing.

# try an exponential smoothing model...
# check help(forecast::ets) to correctly specify the model type!
# parameters for lower and upper options are 
# c("alpha","beta","gamma","phi")
soilT_ets <- ets(soilT_xts, model = "ZZM", 
                 lower = rep(0.001, 4), 
                 upper = c(0.05, 0.2, 0.9999, 0.8))
summary(soilT_ets)
#
# default plot function can plot forecasts too...
plot(forecast(soilT_ets, h=168), col=7)
# lines(soilT_ets$fitted ~ seq(0,5.2668e6, length.out = length(soilT_ets$fitted)),
#       col = "#80008080", lwd = 3)
mtext(names(soilT_ets$par),3,seq(-1,-5,-1), adj = 0.7, col = 4, font = 3)
mtext(signif(soilT_ets$par,3),3,seq(-1,-5,-1), adj = 0.8, col = 4, font = 3)

plot(soilT_ets) # decomposition plot - can also use autoplot()
# (exponential smoothing is not so good for the soilT data!)
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
