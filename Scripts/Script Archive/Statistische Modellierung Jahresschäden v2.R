## ---------------------------------------------------------------------------------
## Statistische Modellierung der GVZ Jahresschaeden 
##
## 1. mit einer Lognormal Verteilung
## 2. als Zeitreihe
##
## Mirco Heidemann
## 04/2016
## ---------------------------------------------------------------------------------

## Arbeitsverzeichnis definieren
setwd('I:/R/Statistische Modellierung Jahresschäden/data/')

## Benoetigte Funktionen laden
require(tseries)
require(forecast)
require(fGarch)
require(MASS)
require(PerformanceAnalytics)
require(sfsmisc)
require(boot)
source('J:/Naturgefahren/Datenanalysen/Statistische Analysen/Rfunctions/logst.R')
source('J:/Naturgefahren/Datenanalysen/Statistische Analysen/Rfunctions/f.acf.R')
source('J:/Naturgefahren/Datenanalysen/Statistische Analysen/Rfunctions/lognormal.R')

## Daten einlesen:
data <- read.csv2('gvz jahreschadenstatistik 2015.csv')
## Nur ES
data <- data[, c(1,2,4)]
data$versSum <- as.numeric(as.character(data$versSum))
data$elementarschad <- as.numeric(as.character(data$elementarschad))

## schadensatz definieren:
## Indem nur die Anzahl der beschaedigten Gebaeude im Verhaeltnis zum
## Gesamtgebaeudebestand betrachtet wird, laesst sich die Schaden-
## entwicklung vollstaendig von Wert- bzw. Preisentwicklungen entkoppeln.
data$lr <- data$elementarschad/data$versSum
dat <- data$lr

## ----------------------------------------------------------------------
## 1. Ansatz: Geeignete Verteilungsfunktion der Schadensaetze/loss ratio
## finden und Parameter schaetzen.
## ----------------------------------------------------------------------

## Finde den besseren fit zwischen einer log-normal und einer gamma
## verteilung:
find.bestfit <- function(x){
  logN <- fitdistr(x, "log-normal")
  gam  <- fitdistr(x, "gamma")
  ans <- ifelse(AIC(logN) < AIC(gam), "logN", "gam")
  return(ans)
}
find.bestfit(dat)

## Paramterschaetzung der lognormal Verteilung
## Maximum-likelihood fitting, eingabe nicht logarithmiert
# mle <- fitdistr(dat, "lognormal")
# meanlog <- mle$estimate["meanlog"]
# sdlog <- mle$estimate["sdlog"]

## function 'lognormal' fits an observed distribution with respect
## to a lognormal model and computes p value
stat <- lognormal(dat, limit=100)
## auf log-normal skala
mu <- exp(as.numeric(stat$meanlog))
sigma <- exp(as.numeric(stat$sdlog))

## qq-plot, normal fit der log-es?
par(mfrow=c(1,2))
qqnorm(log(dat)); qqline(log(dat))
## histogram und dichte der lognormal es
hist(dat, density=20, ylim=c(0, 40), main="histogram and lognormal density")
curve(dlnorm(x, meanlog = stat$meanlog), add=T, col='darkblue', lwd=2)
par(mfrow=c(1,1))


## ----------------------------------------------------------------------
## 2. Ansatz: Schadensaetze/loss ratio als Zeitreihe betrachten und 
## Risikokenngroessen (VaR, ExSf) schaetzen
## -----------------------------------------------------------------------
## Ansatz eines autoregressiven models, schadensaetze als zeitreihe
## Define log-ES-loss ratio as time series of class ts
tsdat <- ts(data$lr , start=1960, frequency=1)
## frequency:
## monthly: frequency=12, deltat=1/12
## yearly:  frequency=1, deltat=1
## daily:   frequency=365, deltat=1/365

plot(tsdat,ylab="Schadensatz in Promille",main="TS-ES-Schadensatz")
# ## Plot the part of the time series from 2006 to 2015.
# plot(window(tsdat, start=1975, end=2015))


## GARCH Model der loss ratio Elementar

## Compute and display Log Returns (lgR), 
## d.h. relative aenderungen in  den loss ratios
lg.tsdat <- log(tsdat)
lgR <- diff(lg.tsdat)
plot(data$Jahr[-1], lgR, ylab="Log Return", xlab="Jahr", type="l",
     main="Log Returns Loss Ratio Elementar")
abline(h=0, col="gray")

## Analyze Log Returns
par(mfrow=c(2,2))
acf(lgR, main="ACF of Log Returns")
acf(lgR^2, main="ACF of Squared Log Returns")
## Check assumption of normal-distribution:
qqnorm(lgR, pch=20, main="Normal QQ Plot of Log Returns")
qqline(lgR)
hist(lgR, col="lightblue", freq=F, main=NA)
title("Histogram of Log Returns"); box()
xx <- seq(-4, 4, length=300)
yy <- dnorm(xx, mean(lgR), sd(lgR))
zz <- dnorm(xx, median(lgR), mad(lgR))
lines(xx, yy, col="red")
lines(xx, zz, col="blue")
legend("topright", lty=1, col=c("red", "blue"), legend=c("plain", "robust"))
par(mfrow=c(1,1))
## --> Log Returns der Loss Ratio Elementar sind normalverteilt!!
## Zusaetzlich testen mit..
#jarque.bera.test(lgR)
## H0: lgR sind gaussverteilt (verwerfen wenn p-value < 0.05)

# ## Determine three reasonable guesses each for a pure ARCH
# ## model, as well as for a GARCH model --> Analyze squared Log Returns
# par(mfrow=c(1,2))
# acf(lgR^2, ylim=c(-1,1), main="ACF of Squared Log Returns")
# pacf(lgR^2, ylim=c(-1,1),main="PACF of Squared Log Returns")
# par(mfrow=c(1,1))
# 
# ## Fit a suitable GARCH(p,q) model
# ## --> Find the order (p,q) with minimal AIC
# ## Choose p and q from 0 to 2. Both cannot be zero at the same time
# 
# mAIC <- matrix(rep(NA, 9), nrow=3)
# colnames(mAIC) <- c("ARCH0","ARCH1","ARCH2")
# rownames(mAIC) <- c("GARCH0","GARCH1","GARCH2")
# for (i in 0:2){ for (j in 0:2){
#   if(i!=0 |j!=0){
#     fit <- garch(lgR,order=c(i,j), trace=F)
#     mAIC[i+1,j+1] <- AIC(fit)
#   }}}
# mAIC
# min(mAIC, na.rm=T)
# ## The AIC suggests a GARCH(0,2) --> ARCH(2)
# ## Nach den ACF Plots auch ARCH(1) moeglich


## Fitting a GARCH(0,1)/ARCH(1) Modell with function fGarch
## ---------------------------------------------------------
## gaussion distribution
gfit <- garchFit(~garch(1,0), lgR, cond.dis="norm", trace=F)
summary(gfit)
## Output:
## if p-Value of mu signif diff from 0, than
## fittting (modelling) a mean is required.

## Jarque Bera Test does reject the H0 (sign. p), than
## Residuals are NOT normal distributed

## Ljung-Box test does not reject the H0, than
## --> No Autocorrelation in squared resids to lag1
## --> Volatility (not good?)

## Residual Analysis GARCH model, Autocorrelation of residuals
res <- residuals(gfit)
par(mfrow=c(2,2))
plot(res, type="l", main="Time Series of Residuals")
acf(res, na.action=na.pass, ylim=c(-1,1), main="ACF of squared Residuals")
pacf(res, na.action=na.pass, ylim=c(-1,1), main="PACF of squared Residuals")
#'QQ-Plot
qqnorm(as.numeric(res, ylim=c(-10, 10), pch=20, main='Residual normal plot'))
qqline(as.numeric(res), col="grey", lty=3)
par(mfrow=c(1,1))

## Risk Management: Value at risk (VaR), expected shortfall (ExSf)

## VaR from data alone, model-free
p <- 1-(99.5/100)
## Mit einer W'keit von 99.5% ist der loss ratio nicht unter 'VaR.emp'
VaR.emp <- quantile(lgR, p)

## VaR mit (Gauss-) Random Walk fuer k jahre
## Annahme: log(loss ratios) sind unabhaengig und gauss-verteilt
## k-jahr 99.5%-VaR
k <- 2
VaR.rw <- qnorm(p, k*mean(lgR), sqrt(k)*sd(lgR))
## VaR mit GARCH/ARCH Modell
k <- 2
VaR.gf <- qnorm(p, mean=k*coef(gfit)[1],
                sd=sqrt(sum(predict(gfit, n.ahead=k)$standardDeviation^2)))

## from log return to loss ratio
## simpleReturn = exp(lgR) - 1
## lr.t = simpleReturn * lr.t-1 + lr.t-1
t <- (exp(lgR)-1)*tsdat[-length(tsdat)] + tsdat[-length(tsdat)]
VaR.gf.smpR <- exp(VaR.gf)-1

## ExSf from data alone, model-free
ExSf.emp <- mean(lgR[lgR<quantile(lgR, p)])




## Entwicklung des GVZ Versicherungskapitals
ts.vs<- ts(data$versSum , start=1960, frequency=1)
plot(ts.vs,ylab="Versicherungskapital in Millionen",
     main="Entwicklung des GVZ Versicherungskapitals")
abline(reg=lm(versSum~Jahr, data), col='red')

# ## lineare regression und residual Analyse:
# lm <- lm(versSum~Jahr, data)
# summary(lm)
# par(mfrow=c(1,2))
# plot(lm, which=1:2)
# par(mfrow=c(1,1))

## Compute and display Log Returns (lgR), 
## d.h. relative aenderungen in  den loss ratios
lg.tsvs <- log(ts.vs)
lgR.vs <- diff(lg.tsvs)
plot(data$Jahr[-1], lgR.vs, ylab="Log Return", xlab="Jahr", type="l",
     main="Log Returns Versicherungskapital")
abline(h=0, col="gray")

## Analyze Log Returns
par(mfrow=c(2,2))
acf(lgR.vs, main="ACF of Log Returns")
acf(lgR.vs^2, main="ACF of Squared Log Returns")
## Check assumption of normal-distribution:
qqnorm(lgR.vs, pch=20, main="Normal QQ Plot of Log Returns")
qqline(lgR.vs)
hist(lgR.vs, col="lightblue", freq=F, main=NA)
title("Histogram of Log Returns"); box()
xx <- seq(-4, 4, length=300)
yy <- dnorm(xx, mean(lgR.vs), sd(lgR.vs))
zz <- dnorm(xx, median(lgR.vs), mad(lgR.vs))
lines(xx, yy, col="red")
lines(xx, zz, col="blue")
legend("topright", lty=1, col=c("red", "blue"), legend=c("plain", "robust"))
par(mfrow=c(1,1))
## --> Log Returns der VS sind NICHT normalverteilt!!
## Zusaetzlich testen mit..
#jarque.bera.test(lgR.vs)
## H0: lgR.vs sind gaussverteilt (verwerfen wenn p-value < 0.05)


