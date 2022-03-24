## ---------------------------------------------------------------------------------
## Statistische Modellierung der GVZ Jahresschaeden mit einer Lognormal Verteilung
##
## 1. Logarithmieren der Schadensaetze von 1960 bis 2015
## 2. Trendbereinigung
## 3. Parameterschaetzung der gesamten Zeitreihe mit der MQ-Methode
##    (Momentum-Quantil Methode)
##
## Mirco Heidemann
## 03/2016
## ---------------------------------------------------------------------------------

## Arbeitsverzeichnis definieren
setwd('I:/R/Statistische Modellierung Jahresschäden/data/')

## Benoetigte Funktionen laden
require(forecast)
require(MASS)
require(PerformanceAnalytics)
require(sfsmisc)
source('J:/Naturgefahren/Datenanalysen/Statistische Analysen/Rfunctions/logst.R')
source('J:/Naturgefahren/Datenanalysen/Statistische Analysen/Rfunctions/f.acf.R')

## Daten einlesen:
data <- read.csv2('gvz jahreschadenstatistik 2015.csv')
## Nur ES
dat <- data[, c(1,2,4)]
dat$versSum <- as.numeric(as.character(dat$versSum))
dat$elementarschad <- as.numeric(as.character(dat$elementarschad))

## schadensatz definieren:
## Indem nur die Anzahl der beschaedigten Gebaeude im Verhaeltnis zum
## Gesamtgebaeudebestand betrachtet wird, laesst sich die Schaden-
## entwicklung vollstaendig von Wert- bzw. Preisentwicklungen abkoppeln.
dat$lr <- dat$elementarschad/dat$versSum

## Logarithmieren der Schadensaetze von 1960 bis 2015
dat$log.lr <- log(dat$lr)

## Define log-ES-loss ratio as time series of class ts
tsdat <- ts(dat$log.lr , start=1960, frequency=1)
## frequency:
## monthly: frequency=12, deltat=1/12
## yearly:  frequency=1, deltat=1
## daily:   frequency=365, deltat=1/365

## Plot TS
plot(tsdat,ylab="Schadensatz in Promille",main="TS-ES-Schadensatz")
# ## Plot the part of the time series from 2006 to 2015.
# plot(window(tsdat, start=1975, end=2015))

## Trendschaetzung mit OLS
## yt=beta0+beta1t+residt

# fit.vs <- lm(essatz ~ Jahr, data=dat)
# summary(fit.vs)

tsfit.vs <- tslm(tsdat ~ trend)
summary(tsfit.vs)
trend <- as.numeric(tsfit.vs$coefficients[2])

f <- forecast(tsfit.vs, h=5,level=c(80,95))
plot(f, ylab="logarithmierte Schadensätze Elementar",
     xlab="t")
lines(fitted(tsfit.vs),col="blue")

## Residual autocorrelation
par(mfrow=c(1,2))
res.tsvs <- resid(tsfit.vs)
plot(res.tsvs,ylab="res (ts.lm)")
abline(0,0)
Acf(res.tsvs)
par(mfrow=c(1,1))

## Fit a suitable ARIMA(p,d,q) model --> correlogram
## p: Auto-Regressive model --> look at a partial autocorrelation graph of the data
## q: Moving Average model --> look at an autocorrelation graph of the data
## d: Differencing --> original ts, d=1, diff-data, d=0

## Is this a stationary time series?
## Draw a time series plot. If need be, transform the original data. Than do diff()
f.acf(tsdat) # Peaks? stationary?
f.acf(diff(tsdat)) # stationary now?

## Alternative:
## try Akaike's Information Criterion (AIC) on a set of models and investigate the models 
## with the lowest AIC values
## try the Schwartz Bayesian Information Criterion (BIC) and investigate the models with 
## the lowest BIC values
## Exemple: require(forecast); auto.arima(x, ic = "aic")

## Write down the model with its estimated coefficients:
## ARIMA(p,d,q): Yt=Xt-Xt-1 mit 
## Yt = alpha1.hat*Yt-1 +Et - beta1.hat*Et-1; sigma^2(Et)
## Et follows the standard normal distribution N(0;1)
ar.fit.es <- arima(tsdat, order=c(2,1,0)) ## or order=c(2,1,0) mit aic
# Coefficients alpha1.hat, beta1.hat and sigma:
ar.fit.es$coef[1]; ar.fit.es$coef[2]; ar.fit.es$sigma2 # or ar.fit

## Check an AR(3) process without differencing
ar.fit.es <- arima(tsdat, order=c(2,1,0))

## Check Residulas, must look like white noise:
ts.resid.es <- ts(resid(ar.fit.es))
f.acf(ts.resid.es)
## Check assumption of normal-distribution:
qqnorm(ts.resid.es, main="Normal QQ Plot of Resids-log Elementarschäden")
qqline(ts.resid.es)

## The loss ratio at year t be written as LRt:
## LRt=ar.fit.es$coef[1]*LRt-1 + ar.fit.es$coef[2]*LRt-2 + Et



# ## Paramter-Schaetzung der lognormal Verteilung. Da die
# ## schadensaetze logarithmiert sind, folgen sie einer
# ## Normalverteilung
# ## JAHR ALS FACTOR, SINNVOLL???
# m.lgn <- glm(log(lr) ~ Jahr, family=gaussian(link="identity"), data=dat)
# anova(m.lgn, test="Chisq")
# 
# m.lgn.sig <- summary(m.lgn)$dispersion
# m.lgn.pred <- exp(predict(m.lgn) + 0.5*m.lgn.sig)
# 
# ## E[y]=exp(alpha+beta*x+sigma/2)
# ## wheras exp(sigma/2) is called the "Volatility Adjustment Factor"
# 
# x <- 2020
# exp(coef(m.lgn)[1] + coef(m.lgn)[2]*x + 0.5 * m.lgn.sig)
# 
# ## GLM with a gamma error and the link log
# ## Schadensaetze Modellieren mit einer Gamma-Verteilung
# m.gamma <- glm(lr ~ Jahr, family=Gamma(link=log),data = dat)
# anova(m.gamma, test="Chisq")
# ## The Gamma model indicates that the mean of our response value y is about...
# # (gamma_coef <- exp(coef(m.gamma)[[1]]))
# (gamma_coef <- mean(exp(coef(m.gamma))))
# 
# ## Gamma Modell
# par(mfrow=c(2,2), mar = c(5, 3, 2, 2))
# ## Tukey-Anscombe Plot (Residual vs. Fitted): Zeigt Abweichung von 
# ## der Form der Regressionsfunktion
# TA.plot(m.gamma, res= residuals(m.gamma, type="deviance"), labels="*",
#         main='Tukey-Anscombe Plot - Devianz Residuen', show.call=F)
# TA.plot(m.gamma, res= residuals(m.gamma, type="working"), labels="*",
#         main='Tukey-Anscombe Plot - Arbeits Residuen', show.call=F)
# ## Partial-Residuals Plot: Plots regression terms against their predictors
# ## show the relationship between a given independent variable and the
# ## response variable given that other independent variables are also in the model
# termplot(m.gamma, partial.resid=TRUE, rug=F,
#          main='Partial-Residuals Plot')
# # ## Normal QQ Plot
# qqnorm(residuals(m.gamma, type="deviance"),
#        main="Normal Q-Q Plot der Devianz Residuen")
# qqline(residuals(m.gamma, type="deviance"))
# 
# par(mfrow=c(1,1))
# 
# # ## more diagnostic plots:
# # library("boot")
# # glm.diag.plots(m.gamma)
# # 
# # ## The estimation of Value at Risk and Expected Shortfall
# # VaR(tsdat, method="historical")
# # ES(spxret11, method="historical")
# 
# AIC(m.lgn, m.gamma)
