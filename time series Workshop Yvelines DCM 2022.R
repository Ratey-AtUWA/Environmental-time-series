#### Introduction to time series analysis for ENVT5503 ####
# 
# Load the packages we need for specialized time series
# functions in R
#
library(zoo)      # for basic irregular time series functions
library(xts)      # we need the xts "eXtended Time Series" format for 
                  # some functions
library(Kendall)  # for trend analysis with Mann-Kendall test
library(trend)    # for trend analysis using the Sen slope
library(forecast) # for time series forecasting with ARIMA and 
                  # exponential smoothing 
library(tseries)  # for assessing stationarity using the 
                  # Augmented Dickey-Fuller test
library(lmtest)   # for Breusch-Pagan heteroscedasticity test etc.
library(car)      # for various commonly-used functions
library(ggplot2)  # alternative to base R plots
#
# (optional) make a better colour palette than the R default!
palette(c("black","red3","green3","blue2",
          "darkcyan","purple","sienna","gray67"))
# 
#### data input ####
# read the data into a data frame - this is how R often stores
# data - it's not the format we need but we'll use it for comparison

# Yvelines <- read.csv("Yvelines.csv")
# colnames(Yvelines) <- c("Date","DCM")
#
# do some checks of the data
summary(Yvelines) # simple summary of each column
str(Yvelines) # more detailed information about the R object ('str'=structure)
plot(Yvelines$DCM, type = "l", col = 6)
#
# we really want the data in a different type of R object!
# we use the 'zoo' package to read the data in as a time series
# this is of class 'zoo' which is much more flexible than the
# default time series object in R

# Yvelines_DCM <- na.omit(Yvelines[,c("Collect.Date", "DCM")])
# write.csv(Yvelines_DCM, file="Yvelines_DCM.csv", row.names = F)
# rm(Yvelines_DCM)
Yvelines_DCM_zoo <- read.csv.zoo("Yvelines.csv",
                               format = "%Y-%m-%d", 
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
par(mfrow = c(3,1), mar = c(4,4,1,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2)
#
# then plot the object with custom axis labels
plot(Yvelines_DCM_zoo, ylab = "DCM (\u00B5g/L)",
     xlab = "Date", col = 7, lwd = 2, cex.lab = 1.4)
plot(log10(Yvelines_DCM_zoo), 
     ylab = expression(bold(paste(log[10],"(DCM, \u00B5g/L)"))),
     xlab = "Date", col = 6, lwd = 2, cex.lab = 1.4)
# note the difference in the x-axis from plot(Yvelines$DCM) !

pt0 <- powerTransform(coredata(Yvelines_DCM_zoo))

plot((Yvelines_DCM_zoo^pt0$roundlam), ylab = "power-transf.(DCM, \u00B5g/L)",
     xlab = "Date", col = 4, lwd = 2, cex.lab = 1.4)

# the R Cookbook suggests Box-Cox (power) transforming the variable
# to "stabilise the variance" (i.e. reduce heteroscedasticity, and
# )
# this looks like it will work, but we suggest using log10[conc]
# as the values are easier to interpret
Yvelines_DCM_zoo <- Yvelines_DCM_zoo^pt0$roundlam
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plot differencing for stationarity
par(mfrow = c(3,1), mar = c(0,4,0,1), oma = c(4,0,1,0), cex.lab = 1.4)
plot(Yvelines_DCM_zoo, ylab = "Raw data, no differencing",
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

# sometimes we need to use another time series data format,
# xts (eXtended Time Series) which allows us to use more functions
#
# make an xts object from our zoo object:

Yvelines_DCM_xts <- as.xts(Yvelines_DCM_zoo)
par(mfrow=c(1,1))
plot(Yvelines_DCM_xts, col = 6, type = "b") # just to check
str(Yvelines_DCM_xts) # just to check

#### Finding the trend ####
## 1. using a moving average ___________________
# create a new time series from a moving average for each 12mo
# to remove daily periodicity using the zoo::rollmean function
Yvelines_DCM_movAv <- rollmean(Yvelines_DCM_zoo, 20)
plot(Yvelines_DCM_zoo, col=8) # original data
lines(Yvelines_DCM_movAv, lwd = 2) # just to check
#
# optionally add some other lines to this plot for comparison
lines(rollmean(Yvelines_DCM_zoo,5), col=2, lty = 2, lwd = 2) # 5 steps is incorrect periodicity
#
## 2. using a linear (regression) model ____________________
# create a linear model of the time series...
lm0 <- lm(coredata(Yvelines_DCM_zoo) ~ index(Yvelines_DCM_zoo))
summary(lm0)

# use a plot to look at the linear fit
plot(Yvelines_DCM_zoo, col = 8, type = "b")
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
# interpreting partial autocorrelations is more complicated -
# refer to Long & Teetor (2019, Section 14.15). Simplistically,
# partial autocorrelation allows us to identify which and how many
# autocorrelations will be needed to model the time series data.
#
# use Box test for autocorrelation
# the null hypothesis is that no autocorrelation exists
# at any lag distance (so p <= 0.05 'rejects' null):
Box.test(Yvelines_DCM_xts)
Box.test(diff(Yvelines_DCM_xts,1)) # 1 difference
Box.test(diff(diff(Yvelines_DCM_xts,1),1)) # 2 differences
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
plot(Yvelines_DCM_periodic2, col = 2)
#
#### (optional) look at the unexplained variation ####
Yvelines_DCM_lmfit <- zoo(lm0$fitted.values, index(Yvelines_DCM_xts))

Yvelines_DCM_err = Yvelines_DCM_zoo - (Yvelines_DCM_movAv + Yvelines_DCM_periodic2)
plot(Yvelines_DCM_err, ylim = c(-200 ,200))
# 
# and plot everything to show the components
Yvelines_DCM_lmfit <- zoo(lm0$fitted.values, index(Yvelines_DCM_zoo))
plot(cbind(Yvelines_DCM_zoo,Yvelines_DCM_periodic2,
     Yvelines_DCM_movAv, Yvelines_DCM_err), ylim = c(-100,650),
     main = expression(paste("Time series decomposition: Groundwater ",
                             "[DCM (mg/L)]"^0.5)),
   ylab = c("Power-transf conc.","Apparent periodicity","Smooth trend","Error"), 
     cex.main = 1.5, yax.flip = TRUE, col = c(8,6,4,2))
# notes: we could also include Yvelines_DCM_lmfit
#        cbind is to combine columns;
#        \n inserts a line break into a text string;
#        \u followed by a 4-character code inserts a Unicode 
#           character (e.g. \u00B0 gives the degrees symbol);
#        yax.flip alternates vertical axis labels

# -=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-
#### using loess instead of moving average ####

# we can visualise trends and decompose time series using loess (locally 
# estimated polynomial smoothing). This does not have the relationship to ARIMA
# that moving averaging does, but does allow us to separate periodicity from 
# random error in time series decomposition.

y_trend1 <- as.data.frame(loess.smooth(index(Yvelines_DCM_zoo), 
                        coredata(Yvelines_DCM_zoo), 
                        span = 0.15, evaluation = length(Yvelines_DCM_zoo)))
Yvelines_DCM_loess1 <- as.zoo(y_trend1$y, order.by=y_trend$x, frequency=1)
plot(Yvelines_DCM_zoo, col = 8, type = "l")
lines(Yvelines_DCM_lmfit, col = 2, lty = 2)
lines(Yvelines_DCM_loess1, col = "skyblue", lwd = 3)

y_trend2 <- loess.smooth(index(Yvelines_DCM_zoo), 
                                    coredata(Yvelines_DCM_zoo), 
                                    span = 0.07, 
                                    evaluation = length(Yvelines_DCM_zoo))
Yvelines_DCM_loess2 <- as.zoo(y_trend2$y, order.by=y_trend$x, frequency=0)
lines(Yvelines_DCM_loess2, col = "purple", lwd = 2)

legend("topright", bty = "n", inset = 0.02, cex = 1.25, 
       legend = c("actual data","linear model","loess coarse", "loess fine"), 
       col = c(8,2,"skyblue","purple"), lty=c(1,2,1,1), lwd = c(1,1,3,2))

Yvelines_DCM_pseudPer <- (Yvelines_DCM_loess2 - Yvelines_DCM_loess1)
Yvelines_DCM_random <- as.zoo(coredata(Yvelines_DCM_zoo) - coredata(Yvelines_DCM_loess1),
                 order.by = index(Yvelines_DCM_zoo), frequency = 0)
par(xpd=T)
plot(cbind(Yvelines_DCM_zoo,Yvelines_DCM_loess1,
           Yvelines_DCM_pseudPer, Yvelines_DCM_random), 
     type="b", pch = 3, main = "", xlab = "Date", 
     cex.main = 1.5, yax.flip = TRUE, col = c(8,6,4,2), ylim = c(-5e4,5e5),
     ylab = c("Data","Trend","Pseudo-periodicity","Random error"))
x1 <- (0.5*(par("usr")[2]-par("usr")[1]))+par("usr")[1]
y <- (c(0.2, 0.45, 0.73, 0.96)*(par("usr")[4]-par("usr")[3]))+par("usr")[3]
text(rep(x1,4), y, 
     labels = c("Unaccounted variation","Apparent periodicity","Trend","Data"), 
     col = c(2,4,6,1), font=2)  ;  par(xpd=F)


#### end exploratory data analysis +=+=+=+=+=+=+=+=+=+=+=+
# ...which leads us into ARIMA forecast modelling ...
# -=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-
#
#### modelling time series with ARIMA ####

# use the forecast:: package to run an ARIMA model
auto.arima(Yvelines_DCM_xts, max.p = 3, max.q = 3, max.d = 3, 
           seasonal = TRUE)
#
# use output of auto.arima() to run arima
#   1. with no seasonality
amY <- arima(x = Yvelines_DCM_xts, order = c(2,1,1))
summary(amY)
confint(amY)
#   2. with seasonality
# findfrequency(Yvelines_DCM_PT)
# am1 <- arima(x = DCM_log, order = c(1,1,1),
#              seasonal = list(order = c(1, 0, 0), period = 20))
# summary(am1)
# confint(am1)
# ls(am1)
#
# checking residuals is our best diagnostic tool...
#   1. residual plot (top) should look like white noise
#   2. residuals should not be autocorrelated (bottom left plot)
#       AND p-value from Ljung-Box test should be > 0.05 
#   3. residuals should be normally distributed (bottom right plot)
checkresiduals(amY)
# checkresiduals(am1)

# use the ARIMA model to produce a forecast ####

# first calculate mean interval of our time series
mean(diff(index(Yvelines_DCM_xts)))
# then e.g. for 50 days prediction:
fcast_len <- 50/as.numeric(mean(diff(index(Yvelines_DCM_xts))))

fcY <- forecast(amY, h = fcast_len) # h value from days/(mean interval of 
                                    # irregular series) = 50/5.72331
# fc1 <- forecast(am1, h = 50)
#
# standard R plot
par(mfrow = c(1,1), mar = c(3,4,2,2), mgp = c(1.7,0.3,0), tcl = 0.3, 
    font.lab=2, cex.lab=1.2)
plot(fcY, 
     ylab=expression(bold("(DCM, \u00B5g/L)"^0.5)),
     xlab="Date", 
     main="Forecasted DCM concentration (\u00B5g/L, power-transformed)",
     col.main = 4, xaxt="n")
# axis tick code to get proper dates on x-axis:
# slope = diff(y1,y2)/diff(x1,x2)
b <- diff(as.numeric(range(index(Yvelines_DCM_zoo))))/diff(range(index(fcY$x)))
# y=a+bx, so a = y-bx
a <- as.numeric(min(index(Yvelines_DCM_zoo)))-(b*min(index(fcY$x)))
axis(1, at = axTicks(1), 
     labels = substr(as.POSIXct(a + b*axTicks(1), origin="1970-01-01"),1,10) )
mtext("ARIMA with no seasonality", 1, -2.2, adj = 0.1, font = 2, col = 4)
mtext(amY$call, 1, -1, adj = 0.1, col = 4, cex = 0.9, 
      family="mono", font = 2)

# we predicted y (date) from x (number) y = a + bx (from axis tick code above)
# so x = (y - a)/b
# the end time is
fc_end <- max(index(fcY$lower))
# then using linear model in axis tick code above:
fc_end_date <- as.POSIXct(a + b * fc_end, origin="1970-01-01")

dates3 <- c(range(index(fcY$x)),fc_end)
abline(v = dates3, col = c("grey70","grey70",4), lty = 2)
del <- 7e3 # offset for text on vertical lines
text(dates3-del, rep(500,3), srt = 90, cex = 0.9, col = c("grey50","grey50",4), 
     labels=substr(as.POSIXct(a + (b*dates3), origin="1970-01-01"),1,10))

# check how dates are stored in arima and forecast objects! ####
# first the actual data
plot(index(Yvelines_DCM_zoo), coredata(Yvelines_DCM_zoo)+120,
     type="p", ylim = c(-120,790))
lines(index(Yvelines_DCM_zoo), fitted.values(fcY)+120, col = "dodgerblue")
lines(index(Yvelines_DCM_zoo),residuals(fcY),col="plum")
points(index(Yvelines_DCM_zoo),rep(par("usr")[3],length(index(Yvelines_DCM_zoo))), pch=3)

# then the arima model and forecast model objects:
summary(index(amY$residuals))
summary(index(fcY$x)) # identical to summary(index(amY$residuals))
summary(diff(index(amY$residuals)))
summary(diff(index(fcY$x))) # identical to summary(diff(index(amY$residuals)))
summary(as.numeric(index(Yvelines_DCM_xts)))
mean(diff(as.numeric(index(Yvelines_DCM_xts))))
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
date2 <- as.POSIXct("1997-05-05 CEST") # end date
plot(forecast(Yvelines_DCM_ets, h=12), col=6, xlab = "time since start", ylab = "Box-Cox transformed conc.")
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

plot(Yvelines_DCM_ets, col = 2) # decomposition plots
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


