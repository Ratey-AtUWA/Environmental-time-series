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

# Yvelines_DCM <- read.csv("incr2.csv")
# colnames(Yvelines_DCM) <- c("Date","conc")
#
# do some checks of the data
summary(Yvelines_DCM) # simple summary of each column
str(Yvelines_DCM) # more detailed information about the R object ('str'=structure)
plot(Yvelines_DCM$conc, type = "l", col = 6)
#
# we really want the data in a different type of R object!
# we use the 'zoo' package to read the data in as a time series
# this is of class 'zoo' which is much more flexible than the
# default time series object in R

# Yvelines_DCM_conc <- na.omit(Yvelines_DCM[,c("Collect.Date", "conc")])
# write.csv(Yvelines_DCM_conc, file="Yvelines_DCM_conc.csv", row.names = F)
# rm(Yvelines_DCM_conc)
Yvelines_DCM_zoo <- read.csv.zoo("Yvelines2.csv",
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
par(mfrow = c(3,1), mar = c(4,4,1,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2, xpd=FALSE)
#
# then plot the object with custom axis labels
plot(Yvelines_DCM_zoo, ylab = "Concentration (\u00B5g/L)",
     xlab = "Date", col = 7, lwd = 2, cex.lab = 1.4)
plot(log10(Yvelines_DCM_zoo), 
     ylab = expression(bold(paste(log[10],"(conc, \u00B5g/L)"))),
     xlab = "Date", col = 6, lwd = 2, cex.lab = 1.4)
# note the difference in the x-axis from plot(Yvelines_DCM$conc) !

pt0 <- powerTransform(coredata(Yvelines_DCM_zoo))
if(pt0$lambda<0) {
  plot(
    -1 * (Yvelines_DCM_zoo ^ pt0$lambda),
    ylab = "power-transf. (Conc., \u00B5g/L)",
    xlab = "Date",
    col = 4, lwd = 2, cex.lab = 1.4
  )
} else {
  plot((Yvelines_DCM_zoo ^ pt0$lambda),
       ylab = "power-transf. (Conc., \u00B5g/L)",
       xlab = "Date",
       col = 4, lwd = 2, cex.lab = 1.4
  )
}
rm(pt0)
# the R Cookbook suggests Box-Cox (power) transforming the variable
# to "stabilize the variance" (i.e. reduce heteroscedasticity, and
# reduce extreme ranges in the data)
# this does work, but can be better to use log10[conc]
# as the values are easier to interpret

# bit let's use power transform anyway...
Yvelines_DCM_zoo <- Yvelines_DCM_zoo ^ 0.5


#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plot differencing for stationarity
par(mfrow = c(3,1), mar = c(0,4,0,1), oma = c(4,0,1,0), cex.lab = 1.4)
plot(Yvelines_DCM_zoo, ylab = "Raw data, no differencing",
     xlab="", xaxt="n", col = 8)
lines(loess.smooth(index(Yvelines_DCM_zoo),coredata(Yvelines_DCM_zoo)), 
      col = 4, lwd = 2)
abline(lm(coredata(Yvelines_DCM_zoo)~index(Yvelines_DCM_zoo)), col = 2)
legend("top", legend = c("Yvelines_DCM_Data", "Loess smoothing","Linear model"),
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

Kendall(coredata(Yvelines_DCM_zoo), as.numeric(index(Yvelines_DCM_zoo)))
# SeasonalMannKendall(Yvelines_DCM_zoo) # needs base R time series object
cor.test(coredata(Yvelines_DCM_zoo), as.numeric(index(Yvelines_DCM_zoo)), method="kendall")
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# often we need to use another time series data format,
# xts (eXtended Time Series) which allows us to use more functions
#
# make an xts object from our zoo object:

Yvelines_DCM_xts <- as.xts(Yvelines_DCM_zoo)
plot(Yvelines_DCM_xts, col = 6, type = "l") # just to check
str(Yvelines_DCM_xts) # just to check

#### Finding the trend ####
## 1. using a moving average ___________________
# create a new time series from a moving average for each 12mo
# to remove daily periodicity using the zoo::rollmean function

findfrequency(Yvelines_DCM_xts) # an approximation!! check the data!

Yvelines_DCM_movAv <- rollmean(Yvelines_DCM_zoo, 52)
plot(Yvelines_DCM_zoo, col=8) # original data
lines(Yvelines_DCM_movAv, lwd = 2) # just to check
lines(loess.smooth(index(Yvelines_DCM_zoo), 
                   coredata(Yvelines_DCM_zoo), span = 0.2), 
      col = 2)
lines(loess.smooth(index(Yvelines_DCM_zoo), 
                   coredata(Yvelines_DCM_zoo), span = 0.7), 
      col = 4, lty = 2, lwd = 2)
legend("topleft", bty = "n", inset = 0.02, cex = 1.5, 
       legend = c("actual data","moving average (52)","LOESS, span = 0.2","LOESS, span = 0.7"), 
       pch = NA, col = c(8,1,2,4), lty=c(1,1,1,2), lwd = c(1,2,1,2))


#
# optionally add some other lines to this plot for comparison
lines(rollmean(Yvelines_DCM_zoo,6), col=6, lty = 2, lwd = 2) # 6 steps is incorrect periodicity
#
## 2. using a linear (regression) model ____________________
# create a linear model of the time series...
lm0 <- lm(coredata(Yvelines_DCM_zoo) ~ index(Yvelines_DCM_zoo))
summary(lm0)
Yvelines_DCM_lmfit <- zoo(lm0$fitted.values, index(Yvelines_DCM_zoo))
# use a plot to look at the linear fit
plot(Yvelines_DCM_zoo, col = 8, type = "l")
lines(Yvelines_DCM_lmfit, col = 2)

conc_trend <- loess.smooth(index(Yvelines_DCM_zoo), 
                           coredata(Yvelines_DCM_zoo), 
                           span = 0.8, evaluation = length(Yvelines_DCM_zoo))

Yvelines_DCM_trend <- zoo(conc_trend$y, index(Yvelines_DCM_zoo))
lines(Yvelines_DCM_trend, col = "blue2")
lm1 <-  lm(log10(coredata(Yvelines_DCM_zoo)) ~ index(Yvelines_DCM_zoo))
lines(zoo(10^lm1$fitted.values, index(Yvelines_DCM_zoo)), col = "purple", lty="dashed")
legend("top", bty = "n", inset = 0.02, cex = 1.5, 
       legend = c("actual data","linear model","log-linear","loess smoothed"), pch = c(1,NA,NA,NA),
       col = c(8,2,"purple","blue2"), lty=c(1,1,2,1))
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
# remember we made a loess model of the time series...
# ...the residuals can give us the periodicity component 
# plus random variation
#
# make a time series of the loess residuals
Yvelines_DCM_periodic <- Yvelines_DCM_zoo - Yvelines_DCM_trend
plot(Yvelines_DCM_trend, ylim = c(-150,660))
lines(Yvelines_DCM_periodic, col = "gold", lwd = 2) # just to check
#
# we can use less smoothing in the loess function to 
# retain periodicity... IF IT EXISTS! . . .
# ...we adjust the 'span' option (lower values give less smoothing)
conc_LOESS2 <- loess.smooth(index(Yvelines_DCM_zoo), 
                           coredata(Yvelines_DCM_zoo), 
                           span = 0.1, evaluation = length(Yvelines_DCM_zoo))
# then use the new loess model to make a time series...
# ...which contains both periodic and trend information
Yvelines_DCM_LOESS2 <- zoo(conc_LOESS2$y, index(Yvelines_DCM_zoo))
lines(Yvelines_DCM_LOESS2, col = "blue2")

# the difference between the data and the less smoothed loess 
# should be just 'noise' or 'error'
Yvelines_DCM_err <- Yvelines_DCM_zoo - Yvelines_DCM_LOESS2
# plot, setting y axis limits to similar scale to original data:
lines(Yvelines_DCM_err, col = 4)

# ...and the periodicity should be the difference between the
# very smoothed and less smoothed loess
Yvelines_DCM_periodic2 <- Yvelines_DCM_LOESS2 - Yvelines_DCM_trend

# and plot everything to show the components
par(xpd=TRUE)
plot(cbind(Yvelines_DCM_zoo,Yvelines_DCM_periodic2,
     Yvelines_DCM_trend, Yvelines_DCM_err),
     main = "Time series decomposition:\nGroundwater solute conc (mg/L)",
     xlab = "Date", 
     cex.main = 1.5, yax.flip = TRUE, col = c(1,6,4,2), ylim = c(-150,660))
x1 <- (0.5*(par("usr")[2]-par("usr")[1]))+par("usr")[1]
y <- (c(0.15, 0.4, 0.6, 0.85)*(par("usr")[4]-par("usr")[3]))+par("usr")[3]
text(rep(x1,4), y, labels = c("Random error","Trend","Periodicity","Data"), col = c(2,4,6,1))
# notes: we could also include linear model if we make it into a time series object
#        cbind is to combine columns;
#        \n inserts a line break into a text string;
#        \u followed by a 4-character code inserts a Unicode 
#           character (e.g. \u00B0 gives the degrees symbol);
#        yax.flip alternates vertical axis labels

# IF WE DON'T THINK THERE IS PERIODICITY run the next lines
Yvelines_DCM_periodic2 <- Yvelines_DCM_periodic2*0
Yvelines_DCM_err <- Yvelines_DCM_zoo - Yvelines_DCM_trend
plot(cbind(Yvelines_DCM_zoo,Yvelines_DCM_periodic2,
           Yvelines_DCM_trend, Yvelines_DCM_err),
     main = "Time series decomposition:\nGroundwater solute conc (mg/L)",
     xlab = "Date", 
     cex.main = 1.5, yax.flip = TRUE, col = c(1,6,4,2), ylim = c(-150,660))
x1 <- (0.5*(par("usr")[2]-par("usr")[1]))+par("usr")[1]
y <- (c(0.15, 0.4, 0.6, 0.8)*(par("usr")[4]-par("usr")[3]))+par("usr")[3]
text(rep(x1,4), y, labels = c("Random error","Trend","Periodicity","Data"), col = c(2,4,6,1))
#
#### end exploratory data analysis +=+=+=+=+=+=+=+=+=+=+=+
# ...which leads us into ARIMA forecast modelling ...
# -=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-=+=-
#
#### modelling time series with ARIMA ####

# use the forecast:: package to run an ARIMA model
auto.arima(Yvelines_DCM_xts, max.p = 3, max.q = 3, max.d = 0, 
           seasonal = TRUE)
#
# use output of auto.arima() to run arima
#   1. with no seasonality
am0 <- arima(x = Yvelines_DCM_xts, order = c(0,1,1))
summary(am0)
confint(am0)
#   2. with seasonality
# findfrequency(Yvelines_DCM_xts)
am1 <- arima(x = Yvelines_DCM_xts, order = c(1,1,0),
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
checkresiduals(am1)
# 
# use the ARIMA model to produce a forecast
fc0 <- forecast(am0, h = 30)
fc1 <- forecast(am1, h = 30)
#
par(mfrow = c(2, 1), cex.main = 0.9, mar = c(3,3,1,1),
    mgp=c(1.6,0.3,0), tcl=0.2, font.lab=2)
plot(fc0,ylab = "Conc, \u00B5g/L]",
     fcol = 4, xlab="Date/Time code", main = "")
lines(Yvelines_DCM_zoo, col = "grey70")
mtext("ARIMA with no seasonality", 3, -1, adj = 0.1, font = 2)
mtext(am0$call, 3, -2.2, adj = 0.1, col = 4, cex = 0.9, 
      family="mono", font = 2)

plot(fc1,ylab = "Conc, \u00B5g/L]",
     xlab="Date/Time code", main = "", fcol = 6)
lines(Yvelines_DCM$conc, col = "grey70")
mtext(paste("ARIMA with",am1$arma[5],"day periodicity"),
      side = 3, line = -1, adj = 0.1, font = 2)
mtext(am1$call, 3, -2.2, adj = 0.1, col = 6, cex = 0.9, family="mono")
par(mfrow = c(1, 1))
#
# we can also make a slightly 'prettier' plot using ggplot2
require(ggplot2) # gives best results using autoplot(...)
autoplot(fc1)+
  ylab("conc (\u00B5g/L)") +
  xlab(paste("Time (s) since",index(Yvelines_DCM_zoo)[1])) +
  ggtitle("Forecasted conc (\u00B5g/L)") +
  theme_bw()
#
# standard R plot
par(mar = c(4,4,2,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2, cex.lab=1.2)
plot(fc1, 
     ylab="conc (\u00B5g/L)",
     xlab=paste("Time (s) since",index(Yvelines_DCM_zoo)[1]), 
     main="Forecasted conc  (\u00B5g/L)",
     xaxt="n")
date1 <- first(index(Yvelines_DCM_xts)) # earliest date on axis
date2 <- last(index(Yvelines_DCM_xts)+365.25*86400,
                    origin = as.POSIXct("1970/1/1 08:00"))
ndays <- 365.25 # between Date axis tick marks
axis(1, at = seq(0, as.numeric(date2)-as.numeric(date1), 86400*ndays), 
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
date1 <- index(Yvelines_DCM_zoo[1]) # earliest date on axis
date2 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)]) # latest date on axis
plot(forecast(Yvelines_DCM_ets, h=30), col=6, 
     xlab = "", xaxt="n", ylab = "Conc.")
mtext(names(Yvelines_DCM_ets$par),3,seq(-1,-5,-1), adj = 0.7, col = 4, font = 3)
mtext(signif(Yvelines_DCM_ets$par,3),3,seq(-1,-5,-1), adj = 0.8, col = 4, font = 3)
date1 <- first(index(Yvelines_DCM_xts)) # earliest date on axis
# date2 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)],
#                     origin = as.POSIXct("1970/1/1 08:00")) # latest date on axis
# date3 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)] - as.numeric(date1),
#                     origin = as.POSIXct("1970/1/1 08:00")) # latest date on axis
date2 <- last(index(Yvelines_DCM_xts))
ndays <- 365.25 # between Date axis tick marks
axis(1, at = seq(0, as.numeric(date2) - as.numeric(date1), 86400*ndays),
     labels = as.Date.character(seq(date1, date2, 86400*ndays)),
     cex.axis = 1)
mtext("Date", side = 1, line = 1.8, font = 2)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))

plot(Yvelines_DCM_ets, col = 2, xaxt="n") # ETS decomposition plots
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


