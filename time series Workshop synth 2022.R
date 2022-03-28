#### Introduction to time series analysis for ENVT5503 ####
# 
# Load the packages we need for specialized time series
# functions in R
#
library(zoo)
library(Kendall)
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

# assessTS <- read.csv("assign.csv")
# colnames(assessTS) <- c("Date","conc")
#
# do some checks of the data
summary(assessTS) # simple summary of each column
str(assessTS) # more detailed information about the R object ('str'=structure)
plot(assessTS$conc, type = "l", col = 6)
#
# we really want the data in a different type of R object!
# we use the 'zoo' package to read the data in as a time series
# this is of class 'zoo' which is much more flexible than the
# default time series object in R

# assessTS_conc <- na.omit(assessTS[,c("Collect.Date", "conc")])
# write.csv(assessTS_conc, file="assessTS_conc.csv", row.names = F)
# rm(assessTS_conc)
assessTS_conc_zoo <- read.csv.zoo("assign.csv",
                               format = "%d/%m/%Y", 
                               tz = "Australia/Perth", 
                               index.column=1,
                               header = TRUE)
# run next comment line if log-transformed variable preferred!
# coredata(assessTS_conc_zoo) <- log10(coredata(assessTS_conc_zoo))

# note that we specified how the date was specified in our data file,
# and the time zone. Our time column was column 1, so we tell R this too.

# do some quick checks of the new R object:
summary(assessTS_conc_zoo) 
str(assessTS_conc_zoo) # POSIXct in the output is a date-time format

#### begin exploratory data analysis of time series ####
#
# try a plot of the data in our time series object...
# first we change the default plotting parameters using par(...)
par(mfrow = c(3,1), mar = c(4,4,1,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2)
#
# then plot the object with custom axis labels
plot(assessTS_conc_zoo, ylab = "conc (\u00B5g/L)",
     xlab = "Date", col = 7, lwd = 2, cex.lab = 1.4)
plot(log10(assessTS_conc_zoo), 
     ylab = expression(bold(paste(log[10],"(conc, \u00B5g/L)"))),
     xlab = "Date", col = 6, lwd = 2, cex.lab = 1.4)
# note the difference in the x-axis from plot(assessTS$conc) !

pt0 <- powerTransform(coredata(assessTS_conc_zoo))

plot((assessTS_conc_zoo^pt0$lambda), ylab = "power-transformed (conc, \u00B5g/L)",
     xlab = "Date", col = 4, lwd = 2, cex.lab = 1.4)

# the R Cookbook suggests Box-Cox (power) transforming the variable
# to "stabilize the variance" (i.e. reduce heteroscedasticity, and
# )
# this looks like it will work, but can be better to use log10[conc]
# as the values are easier to interpret
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plot differencing for stationarity
par(mfrow = c(3,1), mar = c(0,4,0,1), oma = c(4,0,1,0), cex.lab = 1.4)
plot(assessTS_conc_zoo, ylab = "Raw data, no differencing",
     xlab="", xaxt="n", col = 8)
lines(loess.smooth(index(assessTS_conc_zoo),coredata(assessTS_conc_zoo)), 
      col = 4, lwd = 2)
abline(lm(coredata(assessTS_conc_zoo)~index(assessTS_conc_zoo)), col = 2)
legend("topright", legend = c("assessTS_conc_Data", "Loess smoothing","Linear model"),
       cex = 1.8, col = c(1,4,2), lwd = c(1,2,1), bty = "n")
plot(diff(assessTS_conc_zoo,1),
     ylab = "First differencing",
     xlab="", xaxt="n", col = 8)
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(assessTS_conc_zoo,1)),coredata(diff(assessTS_conc_zoo,1))), 
      col = 4, lwd = 2)
abline(lm(coredata(diff(assessTS_conc_zoo,1))~index(diff(assessTS_conc_zoo,1))), col = 2)
plot(diff(diff(assessTS_conc_zoo,1),1),
     ylab = "Second differencing",
     xlab="Date", col = 8)
mtext("Date",side = 1, line = 2.2, font = 2)
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(diff(assessTS_conc_zoo,1),1)),
                   coredata(diff(diff(assessTS_conc_zoo,1),1))), 
      col = 4, lwd = 2)
abline(lm(coredata(diff(diff(assessTS_conc_zoo,1)))~index(diff(diff(assessTS_conc_zoo,1)))), col = 2)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# testing for stationarity
require(tseries) # needed for Augmented Dickey–Fuller (adf) Test
adf.test(assessTS_conc_zoo)
adf.test(diff(assessTS_conc_zoo,1))
adf.test(diff(diff(assessTS_conc_zoo,1),1))

Kendall(coredata(assessTS_conc_zoo), as.numeric(index(assessTS_conc_zoo)))
# SeasonalMannKendall(assessTS_conc_zoo) # needs base R time series object
cor.test(coredata(assessTS_conc_zoo), as.numeric(index(assessTS_conc_zoo)), method="kendall")
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# sometimes we need to use another time series data format,
# xts (eXtended Time Series) which allows us to use more functions
#
# make an xts object from our zoo object:

assessTS_conc_xts <- as.xts(assessTS_conc_zoo)
plot(assessTS_conc_xts, col = 6, type = "l") # just to check
str(assessTS_conc_xts) # just to check

#### Finding the trend ####
## 1. using a moving average ___________________
# create a new time series from a moving average for each 12mo
# to remove daily periodicity using the zoo::rollmean function

findfrequency(assessTS_conc_xts) # an approximation!! check the data!

assessTS_conc_movAv <- rollmean(assessTS_conc_zoo, 52)
plot(assessTS_conc_zoo, col=8) # original data
lines(assessTS_conc_movAv, lwd = 2) # just to check
lines(loess.smooth(index(assessTS_conc_zoo), 
                   coredata(assessTS_conc_zoo), span = 0.2), 
      col = "red3")

#
# optionally add some other lines to this plot for comparison
lines(rollmean(assessTS_conc_zoo,4), col=6, lty = 2, lwd = 2) # 4 steps is incorrect periodicity
#
## 2. using a linear (regression) model ____________________
# create a linear model of the time series...
lm0 <- lm(coredata(assessTS_conc_zoo) ~ index(assessTS_conc_zoo))
summary(lm0)

# use a plot to look at the linear fit
plot(assessTS_conc_zoo, col = 8, type = "b")
abline(lm0, col = 2)

conc_trend <- loess.smooth(index(assessTS_conc_zoo), 
                           coredata(assessTS_conc_zoo), 
                           span = 0.8, evaluation = length(assessTS_conc_zoo))

assessTS_conc_trend <- zoo(conc_trend$y, index(assessTS_conc_zoo))
lines(assessTS_conc_trend, col = "blue2")
lm1 <-  lm(log10(coredata(assessTS_conc_zoo)) ~ index(assessTS_conc_zoo))
lines(zoo(10^lm1$fitted.values, index(assessTS_conc_zoo)), col = "purple", lty="dashed")
legend("top", bty = "n", inset = 0.02, cex = 1.5, 
       legend = c("actual data","linear model","log-linear","loess smoothed"), pch = c(1,NA,NA,NA),
       col = c(8,2,"purple","blue2"), lty=c(1,1,2,1))
#
#### isolating the periodicity ####
# To model the periodicity we need to understand the autocorrelation
# structure of the time series. First we can do this graphically:
#
# plot the autocorrelation function acf()
plot(acf(assessTS_conc_xts))
# what does this tell you about autocorrelation in this time series?
#
# plot the partial autocorrelation function (pacf)
plot(pacf(assessTS_conc_xts))
# interpreting partial autocorrelations is more complicated -
# refer to Long & Teetor (2019, Section 14.15). Simplistically,
# partial autocorrelation allows us to identify which and how many
# autocorrelations will be needed to model the time series data.
#
# use Box test for autocorrelation
# the null hypothesis is that no autocorrelation exists
# at any lag distance (so p <= 0.05 'rejects' null):
Box.test(assessTS_conc_xts)
Box.test(diff(assessTS_conc_xts,1)) # 1 difference
Box.test(diff(diff(assessTS_conc_xts,1),1)) # 2 differences
#
# remember we made a loess model of the time series...
# ...the residuals can give us the periodicity component 
# plus random variation
#
# make a time series of the loess residuals
assessTS_conc_periodic <- assessTS_conc_zoo - assessTS_conc_trend
plot(assessTS_conc_trend, ylim = c(-2,10))
lines(assessTS_conc_periodic, col = "coral", lwd = 2) # just to check
#
# we can use less smoothing in the loess function to retain periodicity...
# ...we adjust the 'span' option (lower values give less smoothing)
conc_LOESS2 <- loess.smooth(index(assessTS_conc_zoo), 
                           coredata(assessTS_conc_zoo), 
                           span = 0.2, evaluation = length(assessTS_conc_zoo))
# then use the new loess model to make a time series...
# ...which contains both periodic and trend information
assessTS_conc_LOESS2 <- zoo(conc_LOESS2$y, index(assessTS_conc_zoo))
lines(assessTS_conc_LOESS2, col = "blue2")

# the difference between the data and the less smoothed loess 
# should be just 'noise' or 'error'
assessTS_conc_err <- assessTS_conc_zoo - assessTS_conc_LOESS2
# plot, setting y axis limits to similar scale to original data:
lines(assessTS_conc_err, col = 4)

# ...and the periodicity should be the difference between the
# very smoother and less smoother loess
assessTS_conc_periodic2 <- assessTS_conc_LOESS2 - assessTS_conc_trend

# and plot everything to show the components
plot(cbind(assessTS_conc_zoo,assessTS_conc_periodic2,
     assessTS_conc_trend, assessTS_conc_err),
     main = "Time series decomposition:\nGroundwater solute conc (mg/L)", 
     cex.main = 1.5, yax.flip = TRUE, col = c(8,6,4,2), ylim = c(-2,12))
mtext("Data", 3, -1, adj=0.95)
# notes: we could also include assessTS_conc_lmfit
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
auto.arima(assessTS_conc_xts, max.p = 3, max.q = 3, max.d = 0, 
           seasonal = TRUE)
#
# use output of auto.arima() to run arima
#   1. with no seasonality
am0 <- arima(x = assessTS_conc_xts, order = c(2,0,2))
summary(am0)
confint(am0)
#   2. with seasonality
# findfrequency(assessTS_conc_xts)
am1 <- arima(x = assessTS_conc_xts, order = c(1,0,1),
             seasonal = list(order = c(1, 0, 1), period = 52))
summary(am1)
confint(am1)
#
# checking residuals is our best diagnostic tool...
#   1. residual plot (top) should look like white noise
#   2. residuals should not be autocorrelated (bottom left plot)
#       AND p-value from Ljung-Box test (R console) should be 
#       GREATER THAN 0.05 
#   3. residuals should be normally distributed (bottom right plot)
checkresiduals(am0)
#
# checkresiduals(am1)
# 
# use the ARIMA model to produce a forecast
fc0 <- forecast(am0, h = 30)
fc1 <- forecast(am1, h = 30)
#
par(mfrow = c(2, 1), cex.main = 0.9, mar = c(3,3,1,1),
    mgp=c(1.6,0.3,0), tcl=0.2, font.lab=2)
plot(fc0,ylab = "Conc, \u00B5g/L]",
     fcol = 4, xlab="Date/Time code", main = "")
lines(assessTS_conc_zoo, col = "grey70")
mtext("ARIMA with no seasonality", 1, -2.2, adj = 0.1, font = 2)
mtext(am0$call, 1, -1, adj = 0.1, col = 4, cex = 0.9, 
      family="mono", font = 2)

plot(fc1,ylab = "Conc, \u00B5g/L]",
     xlab="Date/Time code", main = "", fcol = 6)
lines(assessTS$conc, col = "grey70")
mtext(paste("ARIMA with",am1$arma[5],"day periodicity"),
      side = 1, line = -2.2, adj = 0.1, font = 2)
mtext(am1$call, 1, -1, adj = 0.1, col = 6, cex = 0.9, family="mono")
par(mfrow = c(1, 1))
#
# we can also make a slightly 'prettier' plot using ggplot2
require(ggplot2) # gives best results using autoplot(...)
autoplot(fc1)+
  ylab("conc (\u00B5g/L)") +
  xlab(paste("Time (s) since",index(assessTS_conc_zoo)[1])) +
  ggtitle("Forecasted conc (\u00B5g/L)") +
  theme_bw()
#
# standard R plot
par(mar = c(4,4,2,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2, cex.lab=1.2)
plot(fc1, 
     ylab="conc (\u00B5g/L)",
     xlab=paste("Time (s) since",index(assessTS_conc_zoo)[1]), 
     main="Forecasted conc  (\u00B5g/L)",
     xaxt="n")
date1 <- index(assessTS_conc_zoo[1]) # earliest date on axis
date2 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)]) # latest date on axis
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
assessTS_conc_ets <- ets(assessTS_conc_xts, model = "ZZZ")
summary(assessTS_conc_ets)
#
# default plot function can plot forecasts too...
date1 <- index(assessTS_conc_zoo[1]) # earliest date on axis
date2 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)]) # latest date on axis
plot(forecast(assessTS_conc_ets, h=30), col=6, 
     xlab = "", xaxt="n", ylab = "Conc.")
mtext(names(assessTS_conc_ets$par),3,seq(-1,-5,-1), adj = 0.7, col = 4, font = 3)
mtext(signif(assessTS_conc_ets$par,3),3,seq(-1,-5,-1), adj = 0.8, col = 4, font = 3)
date1 <- index(assessTS_conc_zoo[1]) # earliest date on axis
date2 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)],
                    origin = as.POSIXct("1970/1/1 08:00")) # latest date on axis
date3 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)] - as.numeric(date1),
                    origin = as.POSIXct("1970/1/1 08:00")) # latest date on axis
ndays <- 365.25 # between Date axis tick marks
axis(1, at = seq(0, as.numeric(date3), 86400*ndays),
     labels = as.Date.character(seq(date1, date2, 86400*ndays)),
     cex.axis = 1)
mtext("Date", side = 1, line = 1.8, font = 2)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))

plot(assessTS_conc_ets, col = 2, xaxt="n") # ETS decomposition plots
axis(1, at = seq(0, as.numeric(date3), 86400*ndays),
     labels = as.Date.character(seq(date1, date2, 86400*ndays)),
     cex.axis = 1)

# (exponential smoothing is OK for the conc data!)
#
rm(list = c("date1", "date2", "ndays"))
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


