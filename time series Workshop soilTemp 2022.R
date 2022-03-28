#### Introduction to time series analysis for ENVT5503 ####
#   _____  _                     _____              _            
#  |_   _|(_)                   /  ___|            (_)           
#    | |   _  _ __ ___    ___   \ `--.   ___  _ __  _   ___  ___ 
#    | |  | || '_ ` _ \  / _ \   `--. \ / _ \| '__|| | / _ \/ __|
#    | |  | || | | | | ||  __/  /\__/ /|  __/| |   | ||  __/\__ \
#    \_/  |_||_| |_| |_| \___|  \____/  \___||_|   |_| \___||___/
#            ___                 _              _
#           / _ \               | |            (_)
#          / /_\ \ _ __    __ _ | | _   _  ___  _  ___
#          |  _  || '_ \  / _` || || | | |/ __|| |/ __| 
#          | | | || | | || (_| || || |_| |\__ \| |\__ \ 
#          \_| |_/|_| |_| \__,_||_| \__, ||___/|_||___/ 
#                                    __/ |
#                                   |___/

# Load the packages we need for specialized time series and
# other functions in R
#
library(zoo)      # for basic irregular time series functions
library(xts)      # we need the xts "eXtended Time Series" format for some functions
library(Kendall)  # for trend analysis with Mann-Kendall test
library(trend)    # for trend analysis using the Sen slope
library(forecast) # for time series forecasting with ARIMA and exponential smoothing 
library(tseries)  # for assessing stationarity using Augmented Dickey-Fuller test
library(lmtest)   # for Breusch-Pagan heteroscedasticity test and others
library(car)      # for various commonly-used functions
#
# (optional) make a better colour palette than the R default!
palette(c("black","red3","green3","blue2",
          "darkcyan","purple","sienna","gray67"))
# 
#### data input ####
# read the data into a data frame - this is how R often stores
# data - it's not the format we need but we'll use it for comparison

# soiltemp <- read.csv("soiltemp2.csv")
# colnames(soiltemp) <- c("Date","temp")
#
# do some checks of the data
summary(soiltemp) # simple summary of each column
str(soiltemp) # more detailed information about the R object ('str'=structure)
plot(soiltemp$temp, type = "l", col = 6)
#
# we really want the data in a different type of R object!
# we use the 'zoo' package to read the data in as a time series
# this is of class 'zoo' which is much more flexible than the
# default time series object in R

# soiltemp_conc <- na.omit(soiltemp[,c("Collect.Date", "conc")])
# write.csv(soiltemp_conc, file="soiltemp_conc.csv", row.names = F)
# rm(soiltemp_conc)
soiltemp_T15_zoo <- read.csv.zoo("soiltemp2.csv",
                               format = "%Y-%m-%d %H:%M:%S", 
                               tz = "Australia/Perth", 
                               index.column=1,
                               header = TRUE)
# run next comment line if log-transformed variable preferred!
# coredata(soiltemp_T15_zoo) <- log10(coredata(soiltemp_T15_zoo))

# note that we specified how the date was specified in our data file,
# and the time zone. Our time column was column 1, so we tell R this too.

# do some quick checks of the new R object:
summary(soiltemp_T15_zoo) 
str(soiltemp_T15_zoo) # POSIXct in the output is a date-time format
plot(soiltemp_T15_zoo, col = 2) # just to check

# for seveal functions we need to use another time series data format,
# xts (eXtended Time Series) which allows us to use more functions
#
# make an xts object from our zoo object:

soiltemp_T15_xts <- as.xts(soiltemp_T15_zoo)
plot(soiltemp_T15_xts, col = 6) # just to check
str(soiltemp_T15_xts) # just to check

#### begin exploratory data analysis of time series ####
#
# try a plot of the data in our time series object...
# first we change the default plotting parameters using par(...)
par(mfrow = c(3,1), mar = c(4,4,1,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2)
#
# then plot the object with custom axis labels
plot(soiltemp_T15_zoo, ylab = "Temperature (\u00B0C)",
     xlab = "Date", col = 7, lwd = 2, cex.lab = 1.4)

# let's compare with transformed time series data - we may
# need to do this to meet later modelling assumptions
plot(log10(soiltemp_T15_zoo), 
     ylab = expression(bold(paste(log[10],"(Temperature, \u00B0C)"))),
     xlab = "Date", col = 6, lwd = 2, cex.lab = 1.4)

pt0 <- powerTransform(coredata(soiltemp_T15_zoo))
if(pt0$lambda<0) {
  plot(
    -1 * (soiltemp_T15_zoo ^ pt0$lambda),
    ylab = "power-transf. (Temp., \u00B0C)",
    xlab = "Date",
    col = 4, lwd = 2, cex.lab = 1.4
  )
} else {
  plot((soiltemp_T15_zoo ^ pt0$lambda),
       ylab = "power-transf. (Temp., \u00B0C)",
       xlab = "Date",
       col = 4, lwd = 2, cex.lab = 1.4
  )
}
rm(pt0)

## Which of these plots looks like it might be homoscedastic, 
## i.e., have constant variance regardless of time?

# the R Cookbook suggests Box-Cox (power) transforming the variable
# to "stabilize the variance" (i.e. reduce heteroscedasticity)
# this does work, but can be better to use log10[conc]
# as the values are easier to interpret
#
cat("#### If we think a transformation is needed, then do it now !*!*!*!*!\n")
# e.g.
# soiltemp_T15_zoo <- soiltemp_T15_zoo^0.5 
# **don't forget to do the same for the xts version!**
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#### Assessing if a time series variable is stationary ####
# A stationary variable's mean and variance are not dependent on time.
# In other words, for a stationary series, the mean and variance at 
# any time are representative of the whole series.

# If there is a trend for the value of the variable to 
# increase or decrease, or if there are periodic fluctuations,
# we don't have a stationary time series.

# Many useful statistical analyses and models for
# time series models need a stationary time series as input, or
# a time series that can be made stationary with 
# transformations or differencing.

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# testing for stationarity
require(tseries) # needed for Augmented Dickey–Fuller (adf) Test
d0 <- adf.test(soiltemp_T15_zoo); print(d0)
d1 <- adf.test(diff(soiltemp_T15_zoo,1)); print(d1)
d2 <- adf.test(diff(diff(soiltemp_T15_zoo,1),1)); print(d2)

# plot differencing for stationarity
par(mfrow = c(3,1), mar = c(0,4,0,1), oma = c(4,0,1,0), cex.lab = 1.4)
plot(soiltemp_T15_zoo, ylab = "Raw data, no differencing",
     xlab="", xaxt="n", col = 8)
lines(loess.smooth(index(soiltemp_T15_zoo),coredata(soiltemp_T15_zoo)), 
      col = 4, lwd = 2)
abline(lm(coredata(soiltemp_T15_zoo)~index(soiltemp_T15_zoo)), col = 2)
legend("topright", legend = c("soiltemp_T15_Data", "Loess smoothing","Linear model"),
       cex = 1.8, col = c(1,4,2), lwd = c(1,2,1), bty = "n")
mtext(paste("adf.test p value =",signif(d0$p.value,3)), 
      side = 1, line = -1.2, adj = 0.05)
plot(diff(soiltemp_T15_zoo,1),
     ylab = "First differencing",
     xlab="", xaxt="n", col = 8)
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(soiltemp_T15_zoo,1)),coredata(diff(soiltemp_T15_zoo,1))), 
      col = 4, lwd = 2)
abline(lm(coredata(diff(soiltemp_T15_zoo,1))~index(diff(soiltemp_T15_zoo,1))), col = 2)
mtext(paste("adf.test p value =",signif(d1$p.value,3)), 
      side = 1, line = -1.2, adj = 0.05)
plot(diff(diff(soiltemp_T15_zoo,1),1),
     ylab = "Second differencing",
     xlab="Date", col = 8)
mtext("Date",side = 1, line = 2.2, font = 2)
abline(h = 0, col = "grey", lty = 2)
lines(loess.smooth(index(diff(diff(soiltemp_T15_zoo,1),1)),
                   coredata(diff(diff(soiltemp_T15_zoo,1),1))), 
      col = 4, lwd = 2)
abline(lm(coredata(diff(diff(soiltemp_T15_zoo,1)))~index(diff(diff(soiltemp_T15_zoo,1)))), col = 2)
mtext(paste("adf.test p value =",signif(d2$p.value,3)), 
      side = 1, line = -1.2, adj = 0.05)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))

#### Finding the trend ####

# 1. Determine if there is actually any trend.
#    (see Department of Water 2015)
# 1a. apply the Mann-Kendall test from the 'trend' package
mk.test(coredata(soiltemp_T15_zoo))
#    or
Kendall(coredata(soiltemp_T15_zoo), as.numeric(index(soiltemp_T15_zoo)))
# SeasonalMannKendall(soiltemp_T15_zoo) # needs base R regular time series object
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# 1b. estimate the Sen's slope
sens.slope(coredata(soiltemp_T15_zoo))

## 2. using a moving average ___________________
# create a new time series from a moving average for each period
# to remove daily periodicity using the zoo::rollmean function

findfrequency(soiltemp_T15_xts) # an approximation!! check the data!

# We know these are hourly data with diurnal fluctuation,
# so a rolling mean for chunks of length 24 should be OK
soiltemp_T15_movAv <- rollmean(soiltemp_T15_zoo, 24)
plot(soiltemp_T15_zoo, col=8, type="l") # original data
# add the moving average
lines(soiltemp_T15_movAv, lwd = 2)
legend("topright", legend = c("Data","Moving average"), col = c(8,1), lwd = c(1,2))

# optionally add some other lines to this plot for comparison
lines(rollmean(soiltemp_T15_zoo,15), col=3, lty = 5, lwd = 2) # 15 steps is incorrect periodicity
legend("topright", legend = c("Data","Moving average","Wrong periodicity!"), 
       col = c(8,1,3), lwd = c(1,2,2), lty = c(1,1,5))
#
## 3. using a linear (regression) model ____________________
# create a linear model of the time series...
lm0 <- lm(coredata(soiltemp_T15_zoo) ~ index(soiltemp_T15_zoo))
summary(lm0)

# use a plot to look at the linear fit
plot(soiltemp_T15_zoo, col = 8, type = "l", pch = 3, cex = 0.5)
abline(lm0, col = 2, lty = 2)

## 4. using Locally Estimate Scatterplot Smoothing (loess) ____________________
#     (a form of locally weighted non-parametric regression
#      to fit smooth curves to data)
y_trend <- loess.smooth(index(soiltemp_T15_zoo), 
                           coredata(soiltemp_T15_zoo), 
                           span = 0.15, evaluation = length(soiltemp_T15_zoo))

soiltemp_T15_trend <- zoo(y_trend$y, index(soiltemp_T15_zoo))
lines(soiltemp_T15_trend, col = "skyblue", lwd = 3)
legend("top", bty = "n", inset = 0.02, cex = 1.25, 
       legend = c("actual data","linear model","loess smoothed"), 
       col = c(8,2,"skyblue"), lty=c(1,2,1), lwd = c(1,1,3))
#
# looking at data homo/heteroscedasticity - this is fudging it a LOT!
lm1 <- lm(coredata(soiltemp_T15_zoo) ~ y_trend$y)
bptest(lm1) # check homoscedasticity of ts residuals
shapiro.test(lm1$residuals) ; qqPlot(lm1$residuals)

#### isolating the periodicity ####

# NOTE THAT TIME SERIES DON'T ALWAYS HAVE PERIODICITY !

# To model the periodicity we need to understand the autocorrelation
# structure of the time series. First we can do this graphically:
#
# plot the autocorrelation function acf()
acf(soiltemp_T15_xts)
# what does this tell you about autocorrelation in this time series?
#
# plot the partial autocorrelation function (pacf)
pacf(soiltemp_T15_xts)
# interpreting partial autocorrelations is more complicated -
# refer to Long & Teetor (2019, Section 14.15). Simplistically,
# partial autocorrelation allows us to identify which and how many
# autocorrelations will be needed to model the time series data.
#
# use Box test for autocorrelation
# the null hypothesis is that no autocorrelation exists
# at any lag distance (so p <= 0.05 'rejects' null):
Box.test(soiltemp_T15_xts)
Box.test(diff(soiltemp_T15_xts,1)) # 1 difference
Box.test(diff(diff(soiltemp_T15_xts,1),1)) # 2 differences
#
# remember we made a loess model of the time series...
# ...the residuals can give us the periodicity component 
# plus random variation
#
# make a time series of the loess residuals
soiltemp_T15_periodic <- soiltemp_T15_zoo - soiltemp_T15_trend
plot(soiltemp_T15_trend, ylim = c(-2,20))
lines(soiltemp_T15_periodic, col = "coral", lwd = 2) # just to check
#
# we can use less smoothing in the loess function to retain periodicity...
# ...we adjust the 'span' option (lower values give less smoothing)
conc_LOESS2 <- loess.smooth(index(soiltemp_T15_zoo), 
                           coredata(soiltemp_T15_zoo), 
                           span = 0.012, evaluation = length(soiltemp_T15_zoo))
# then use the new loess model to make a time series...
# ...which contains both periodic and trend information
soiltemp_T15_LOESS2 <- zoo(conc_LOESS2$y, index(soiltemp_T15_zoo))
lines(soiltemp_T15_LOESS2, col = "blue2")

# the difference between the data and the less smoothed loess 
# should be just 'noise' or 'error'
soiltemp_T15_err <- soiltemp_T15_zoo - soiltemp_T15_LOESS2
# plot, setting y axis limits to similar scale to original data:
lines(soiltemp_T15_err, col = 4)

# ...and the periodicity should be the difference between the
# very smoothed and less smoothed loess
soiltemp_T15_periodic2 <- soiltemp_T15_LOESS2 - soiltemp_T15_trend

# and plot everything to show the components
plot(cbind(soiltemp_T15_zoo,soiltemp_T15_periodic2,
     soiltemp_T15_trend, soiltemp_T15_err),
     main = "Time series decomposition:\nSoil Temperature at 15cm (\u00B0C)", 
     cex.main = 1.5, yax.flip = TRUE, col = c(8,6,4,2), ylim = c(-2,22))
x1 <- (0.5*(par("usr")[2]-par("usr")[1]))+par("usr")[1]
y <- (c(0.21, 0.3, 0.6, 0.82)*(par("usr")[4]-par("usr")[3]))+par("usr")[3]
text(rep(x1,4), y, 
     labels = c("Unaccounted variation","Trend","Apparent periodicity","Data"), 
     col = c(2,4,6,1))
# notes: we could also include soiltemp_T15_lmfit
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
#    ___  ______  _____ ___  ___  ___                            _        _      
#   / _ \ | ___ \|_   _||  \/  | / _ \                          | |      | |     
#  / /_\ \| |_/ /  | |  | .  . |/ /_\ \    _ __ ___    ___    __| |  ___ | | ___ 
#  |  _  ||    /   | |  | |\/| ||  _  |   | '_ ` _ \  / _ \  / _` | / _ \| |/ __|
#  | | | || |\ \  _| |_ | |  | || | | |   | | | | | || (_) || (_| ||  __/| |\__ \
#  \_| |_/\_| \_| \___/ \_|  |_/\_| |_/   |_| |_| |_| \___/  \__,_| \___||_||___/
#                                                                                
# use the forecast:: package to run an ARIMA model
auto.arima(soiltemp_T15_zoo, max.p = 3, max.q = 3, max.d = 0, 
           seasonal = TRUE)
#
# the best ARIMA models will have the LOWEST AIC (or AICc) value
#
# use output of auto.arima() to run arima
#   1. with no seasonality
am0 <- arima(x = soiltemp_T15_xts, order = c(1,0,3))
summary(am0)
confint(am0)
#   2. with seasonality
# findfrequency(soiltemp_T15_xts)
am1 <- arima(x = soiltemp_T15_xts, order = c(1,0,2),
             seasonal = list(order = c(1, 0, 1), period = 24))
summary(am1)
confint(am1)
ls(am1)
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
fc0 <- forecast(am0, h = 168)
fc1 <- forecast(am1, h = 168)
#
par(mfrow = c(2, 1), cex.main = 0.9, mar = c(3,3,1,1),
    mgp=c(1.6,0.3,0), tcl=0.2, font.lab=2)
plot(fc0,ylab = "Temperature (\u00B0C)",
     fcol = 4, xlab="Date/Time code", main = "")
lines(soiltemp_T15_zoo, col = "grey70")
mtext("ARIMA with no seasonality", 1, -2.2, adj = 0.1, font = 2)
mtext(am0$call, 1, -1, adj = 0.1, col = 4, cex = 0.9, 
      family="mono", font = 2)

plot(fc1,ylab = "Temperature (\u00B0C)",
     xlab="Date/Time code", main = "", fcol = 6)
lines(soiltemp$conc, col = "grey70")
mtext(paste("ARIMA with",am1$arma[5],"h periodicity"),
      side = 1, line = -2.2, adj = 0.1, font = 2)
mtext(am1$call, 1, -1, adj = 0.1, col = 6, cex = 0.9, family="mono")
par(mfrow = c(1, 1))
#
# we can also make a slightly 'prettier' plot using ggplot2
require(ggplot2) # gives best results using autoplot(...)
autoplot(fc1)+
  ylab("Temperature (\u00B0C)") +
  xlab(paste("Time (s) since",index(soiltemp_T15_zoo)[1])) +
  ggtitle("Forecasted conc (\u00B5g/L)") +
  theme_bw()
#
# standard R plot
par(mar = c(4,4,2,1), mgp = c(1.7,0.3,0), tcl = 0.3, font.lab=2, cex.lab=1.2)
plot(fc1, 
     ylab="Temperature (\u00B0C)",
     xlab=paste("Time (s) since",index(soiltemp_T15_zoo)[1]), 
     main="Forecasted conc  (\u00B5g/L)")

#           |_    _|_       _  o _|_      
#   o  o  o |_)|_| |_   \^/(_| |  |_     
#                                   /      
#  _|_|_  _  __ _  /  _    __  _  __ _  | 
#   |_| |(/_ | (/_   _>    |||(_) | (/_ o   
#
#
#### Exponential Smoothing Models ####
# _._._._._._._._._._._._._._._._._._._._._._._._._._
#  _____                                   _   _       _ 
# | ____|_  ___ __   ___  _ __   ___ _ __ | |_(_) __ _| |
# |  _| \ \/ / '_ \ / _ \| '_ \ / _ \ '_ \| __| |/ _` | |
# | |___ >  <| |_) | (_) | | | |  __/ | | | |_| | (_| | |
# |_____/_/\_\ .__/ \___/|_| |_|\___|_| |_|\__|_|\__,_|_|
#   ____     |_|                _   _     _
#  / ___| _ _ _ _   ___   __ _ | |_| |__ (_)_ __   __ _   
#  \___ \| '_ ` _ \ / _ \ / _ \| __| '_ \| | '_ \ / _` |  
#   ___) | | | | | | (_) | (_) | |_| | | | | | | | (_| |  
#  |____/|_| |_| |_|\___/ \___/ \__|_| |_|_|_| |_|\__, |  
#                                                 |___/
# _._._._._._._._._._._._._._._._._._._._._._._._._._
#                                               #
# sometimes ARIMA models may not be the best option
# another commonly used method is exponential smoothing
#
#### try an exponential smoothing model ####
# check help(forecast::ets) to correctly specify the model type!
# parameters for lower & upper are c("alpha","beta","gamma","phi")
soiltemp_T15_ets <- ets(soiltemp_T15_xts, model = "ZZZ")
summary(soiltemp_T15_ets)
#
# default plot function can plot forecasts too...
date1 <- index(soiltemp_T15_zoo[1]) # earliest date on axis
date2 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)],
                    origin = as.POSIXct("1970/1/1 08:00")) # latest date on axis
date3 <- as.POSIXlt(index(fc0$lower)[NROW(fc0$lower)] - as.numeric(date1),
                    origin = as.POSIXct("1970/1/1 08:00")) # latest date on axis
plot(forecast(soiltemp_T15_ets, h=168), col=6, 
     xlab = "", xaxt="n", ylab = "Conc.")
mtext(names(soiltemp_T15_ets$par),3,seq(-1,-5,-1), adj = 0.7, col = 4, font = 3)
mtext(signif(soiltemp_T15_ets$par,3),3,seq(-1,-5,-1), adj = 0.8, col = 4, font = 3)
ndays <- 15 # between Date axis tick marks
axis(1, at = seq(0, as.numeric(date3), 86400*ndays),
     labels = as.Date.character(seq(date1, date2, 86400*ndays)),
     cex.axis = 1)
mtext("Date", side = 1, line = 1.8, font = 2)
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))

plot(soiltemp_T15_ets, col = 2, xaxt="n") # ETS decomposition plots
axis(1, at = seq(0, as.numeric(date3), 86400*ndays),
     labels = as.Date.character(seq(date1, date2, 86400*ndays)),
     cex.axis = 1)

# (exponential smoothing is maybe not so good for the soil temperature data!)
#
rm(list = c("date1", "date2", "date3", "ndays"))
# end code
# 
#### REFERENCES ####
#
# Department of Water (2015). Calculating trends in nutrient data, Government of
# Western Australia, Perth.
# https://water.wa.gov.au/__data/assets/pdf_file/0014/6800/Calculating-trends-in-nutrient-data-10-3-15.pdf
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


