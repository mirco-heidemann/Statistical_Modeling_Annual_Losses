## ----
## Risk Modelling, Part I
## Schaetzung jaehrlicher Elementarschaeden mit einer Lognormalverteilung und 
## einem Paretomodell im tail
## Pro Gefahr einzel oder gleich als gesamte Elementarschaedenverteilung
##
## Modellierung anhand der Schadenerfahrung der GVZ
##
## Mirco Heidemann
## 12/2016, 10/2017
## ----

## relative pfade spezifizieren
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pth_data <- ("./data/")
pth_rdata <- ("./data/Rdata/")
pth_pres <- ("./presentation/")
pth_fun <- ("../R_functions/")
pth_port <- ("../GVZ_Portfolio/")

## Funktionen laden
library(MASS)
library(evir)
#library(evd)
library(QRM)
library(stringr)
library(dplyr)

## Funktion fuer die Aufbereitung der Elementarschaeden laden
source(paste0(pth_fun, 'Fun_SchadPrepAll.R'))

## Daten laden und aufbereiten
schad.file <- paste0(pth_port, 'schaeden_201801_georef.csv')
index.file <- paste0(pth_port, 'versicherungsindex_gvz.csv')

## indexierte elementar- und feuerschaeden
schad <- schadPrep(schad.file, index.file)

## elementarschaeden ab 1982, nur (indexierte) schaeden groesser 0
schad_elementar <- schad %>% filter(format(schad_datum, "%Y") > 1981,
                                    # schadendaten bis und mit letzten Jahres
                                    format(schad_datum, "%Y") < 2018,
                                    schad_art == "Elementar",
                                    schad_index > 0)

schad_hagel <- schad_elementar %>% filter(str_detect(schad_code, 'Hagel'))
schad_sturm <- schad_elementar %>% filter(str_detect(schad_code, 'Sturm'))
schad_wasser <- schad_elementar %>% filter(str_detect(schad_code, 'Hochwasser'))
## Restlichen gefahren zusammen: Schneedruck, Erdrutsch & Steinschlag, Lawinen
schad_rest <- schad_elementar %>% filter(!str_detect(schad_code,
                                                     paste(c('Sturm',
                                                             'Hochwasser',
                                                             'Hagel'), collapse="|")))

## Welche jahresschaeden sollen modelliert werden?
df <- schad_hagel

# ## schadendaten pro gefahr sturm, wasser, hagel und rest
# ## Sturmschaeden
# ind.sturm <- which(str_detect(schad$schad_code, 'Sturm'))
# ## Ueberschwemmungsschaeden
# ind.wasser <- which(str_detect(schad$schad_code, 'Hochwasser'))
# ## Hagelschaeden
# ind.hagel <- which(str_detect(schad$schad_code, 'Hagel'))
# ## Rest
# ind.rest <- which(!str_detect(schad$schad_code,
#                               paste(c('Sturm',
#                                       'Hochwasser',
#                                       'Hagel'), collapse="|")))
# 
# schadEvery <- schad
# ## waehle die gefahr, oder benutze gesamte schadenvertilung ueber
# ## alle Gefahren
# if(length(ind.hagel) > 0) schad <- schadEvery[ind.hagel,]

## yearly data
aggr.anz <- aggregate(df$schad_index,
                      by = list(as.numeric(format(df$schad_datum,'%Y'))),
                      'length')
aggr.sum <- aggregate(df$schad_index,
                      by = list(as.numeric(format(df$schad_datum,'%Y'))),
                      'sum')
yearly.dat <- as.data.frame(cbind(aggr.sum, aggr.anz$x))
names(yearly.dat) <- c('jahr', 'schadSum', 'schadAnz')

# ## The GPD fits a distribution to exceedences over a threshold.
# ## We therefore need to pick an appropriate threshold. A useful
# ## plot to assist in this choice is the mean excess plot.
# par(mfrow=c(1,2))
# meplot(df$schad_index) # From evir package
# ## mrlplot(df$schad_index) # From evd package
# shape(df$schad_index) # From evir package
# par(mfrow=c(1,1))

# threshold <- 1.75e5 ## alle schaeden zusammen
# threshold.sturm <- 8e4; threshold <- threshold.sturm
# threshold.wasser <- 1.25e5; threshold <- threshold.wasser
# threshold.hagel <- 8.75e4; threshold <- threshold.hagel
# threshold.rest <- 8.5e3; threshold <- threshold.rest

## 90% Quantile as an appropriate threshold for the GPD fit?
threshold <- quantile(df$schad_index, 0.9)

schadSum.body <- df %>% filter(schad_index <= threshold) %>% 
  select(schad_index)
schadSum.body <- as.matrix(schadSum.body)
# ind <- which(df$schad_index <= threshold)
# schadSum.body <- df$schad_index[ind]
schadSum.tail <- df %>% filter(schad_index > threshold) %>% 
  select(schad_index)
schadSum.tail <- as.matrix(schadSum.tail)
# ind <- which(df$schad_index > threshold)
# schadSum.tail <- df$schad_index[ind]

## A) Schadensumme der nicht Extremen modelieren mit einer lognormal Verteilung
# Paramterschaetzung mit Maximum-likelihood fitting, eingabe nicht logarithmiert

## lognorm Verteilung nur mit Werte > 0
# schadSum.body <- schadSum.body %>% filter(schad_index > 0)
# ind <- which(df$schad_index <= 0)
# if (length(ind) > 0) df <- df[-ind,]

fit.ln <- fitdistr(schadSum.body, "lognormal")
mulog <- fit.ln$estimate["meanlog"]
sdlog <- fit.ln$estimate["sdlog"]

## qq-plot, normal fit der log-es?
par(mfrow=c(1,2))
qqnorm(log(schadSum.body)); qqline(log(schadSum.body))
## histogram und dichte der lognormal es
truehist(schadSum.body, nbins = 50, col= 8, xlim = c(0, 6e4),
         main = "histogram and lognormal density")
curve(dlnorm(x, meanlog = mulog, sdlog = sdlog),
      add = T, col = 'darkblue', lwd = 2)
par(mfrow = c(1,1))

# Draw cumulative density functions
plot(ecdf(schadSum.body), cex = 0, xlim = c(0, 1e5))
curve(plnorm(x, meanlog = fit.ln$estimate[1],
             sdlog = fit.ln$estimate[2]),
      add = T, col = 'red')
legend('right', lwd = 1, col = c(1,2), legend = c("empirical CDF", "lognormal CDF"))

## B) Anzahl der Schaeden pro Jahr schaetzen
truehist(yearly.dat$schadAnz, col = 8, nbins = 10,
         main="Anzahl Schaeden pro Jahr")

## ... mit einer Poisson-Verteilung 
fit.pois <- fitdistr(yearly.dat$schadAnz, densfun = "Poisson")
## oder mean(yearly.dat$schadAnz)
## im Mittel lambda bezahlte schaeden pro jahr...
lambda <- fit.pois$estimate["lambda"]

## ... mit einer gamma-Verteilung, scale data first, than then scaling the rate
fit_gamma <- fitdistr(yearly.dat$schadAnz / 1000, densfun = "gamma") ## scaled the data by 0.1
rate_gamma <- coef(fit_gamma)["rate"] / 1000  ## scaled the rate by 0.1
shape_gamma <- coef(fit_gamma)["shape"]

## Anzahl der Schaeden mit ein log-Normalverteilung
fit.ln.anz <- fitdistr(yearly.dat$schadAnz, "lognormal")
mulog.anz <- fit.ln.anz$estimate["meanlog"]
sdlog.anz <- fit.ln.anz$estimate["sdlog"]
## varianz modification (not scientific...)
# sdlog.anz.mod <- sdlog.anz * 1.8 ## alle zusammen
# sdlog.anz.mod <- sdlog.anz * 2 ## wasser
# sdlog.anz.mod <- sdlog.anz * 1.8 ## Sturm
# sdlog.anz.mod <- sdlog.anz * 1.5 ## Hagel

## qq-plot, normal fit der log-Anzahl?
par(mfrow = c(1,2))
qqnorm(log(yearly.dat$schadAnz)); qqline(log(yearly.dat$schadAnz))
## histogram und dichte der lognormal es
truehist(yearly.dat$schadAnz, nbins = 50, col = 8, xlim = c(0, 6e4),
         main="histogram and lognormal density")
curve(dlnorm(x, meanlog = mulog.anz,sdlog = sdlog.anz),
      add = T, col = 'darkblue', lwd = 2)
curve(dlnorm(x, meanlog = mulog.anz, sdlog = sdlog.anz.mod),
      add = T, col = 'red', lwd = 2)
legend('topright',
       c('density with empirical variance', 'density with modificated variance'),
       lwd = 2, col = c('darkblue', 'red'), bty = 'n')
par(mfrow=c(1,1))

# Draw cumulative density functions
plot(ecdf(yearly.dat$schadAnz), cex = 0, xlim = c(0, 1e5), main = 'CDF Lognormal SchadAnz')
curve(plnorm(x, meanlog = fit.ln.anz$estimate[1],
             sdlog = fit.ln.anz$estimate[2]),
      add = T, col = 'darkblue')
curve(plnorm(x, meanlog = fit.ln.anz$estimate[1],
             sdlog = sdlog.anz.mod),
      add = T, col = 'red')
legend('right', lwd = 1, col = c('black', 'darkblue', 'red'),
       legend = c('empirical CDF', 'lognormal CDF', 'lognormal CDF modified varince'),
       bty = 'n')

## Wenn lambda gross, approximativ mit einer Normalverteilung...
fit.approx.norm <- fitdistr(yearly.dat$schadAnz, densfun = "normal")
mu.norm.approx <- as.numeric(fit.approx.norm$estimate["mean"])
sd.norm.approx <- as.numeric(fit.approx.norm$estimate["sd"])
mu.norm.approx <- round(rnorm(1,mean = lambda, sd = sqrt(lambda)))

## welche verteilung fuer die anzahl jahresschaeden?
plot(ecdf(yearly.dat$schadAnz), cex = 0, main = 'CDF', xlab = NA)
curve(ppois(x, lambda), add = T, col = '2')
lines(seq(-100, max(yearly.dat$schadAnz), 100),
      pnorm(seq(-100, max(yearly.dat$schadAnz), 100), mu.norm.approx, sd.norm.approx),
      type = "s", ylab = "F(x)", main = "Poisson CDF", col = '6')
curve(plnorm(x, meanlog=fit.ln.anz$estimate[1],
             sdlog=fit.ln.anz$estimate[2]), add = T, col = '4')
curve(plnorm(x, meanlog=fit.ln.anz$estimate[1],
             sdlog = sdlog.anz.mod), add = T, col = '3')
legend('right', lwd=1, col = c(1,2,6,4,3),
       legend = c('empirical CDF', 'poisson CDF', 'approx. normal CDF',
                'logNorm CDF', 'logNorm CDF modified'))

## C) The Generalised Pareto Distribution (GPD) describes exceedences over a
## threshold and is convenient for modelling insurance claims.
fit.gpd <- gpd(df$schad_index, threshold) # From evir package
xi <- fit.gpd$par.ests[1]
beta <- fit.gpd$par.ests[2]

## diagnostics mit dem QQ plot und der CDF
par(mfrow=c(1,2))
qplot(df$schad_index, xi = xi)
# CDF Plot
plot(ecdf(df$schad_index[which(df$schad_index>threshold)]),
     cex = 0, main = "CDF Plot")
lines(sort(df$schad_index),
      pgpd(sort(df$schad_index),xi = xi, beta = beta, mu = threshold), col = 2)
par(mfrow=c(1,1))
## plot of the tail of the underlying distribution of the data
gpd.q(tailplot(fit.gpd), 0.995)


## D) Convolution Method Obtain Overall Annual Loss Distribution

# Quantile function of lognormal-GPD severity distribution
qlnorm.gpd = function(p, theta, theta.gpd, u)
{
  ## Fu: W'keit, dass eine logNorm verteilte Zufallszahl kleiner ist als treshold
  Fu = plnorm(u, meanlog = theta[1], sdlog = theta[2])
  x = ifelse(p < Fu, ## fuer die kleinen w'keiten (aus uniformvert.), simuliere aus GPD
             qlnorm( p = p, meanlog = theta[1], sdlog = theta[2]),
             qgpd(p = (p - Fu) / (1 - Fu) , xi = theta.gpd[1],
                   mu = theta.gpd[2], beta = theta.gpd[3]))
  return(x)
}
# Random sampling function of lognormal-GPD severity distribution
rlnorm.gpd = function(n, theta, theta.gpd, u)
{ r = qlnorm.gpd(runif(n), theta, theta.gpd, u)} ## runif: n Zahlen zwischen 0 und 1
#set.seed(1000)

## nSim = 1e6 # Number of simulated annual losses
## fuer modellentwicklung
nSim = 1e4 # Number of simulated annual losses

H = threshold # Threshold body-tail
# lambda = lambda # Parameter of Poisson body
theta1 = mulog # Parameter mu of lognormal (body)
theta2 = sdlog # Parameter sigma of lognormal (body)
theta1.tail = xi # Shape parameter of GPD (tail)
theta2.tail = H # Location parameter of GPD (tail)
theta3.tail = beta # Scale parameter of GPD (tail)
sj = rep(0, nSim) # Annual loss distribution inizialization
# freq = rpois(nSim, lambda) # Random sampling from Poisson
freq = round(rgamma(nSim, shape_gamma, rate_gamma), 0) # Random sampling from gamma
# freq = abs(rnorm(nSim, mu.norm.approx, sd.norm.approx)) # Random sampling from normal distribution
# freq <- rlnorm(nSim, mulog.anz, sdlog.anz) # Random sampling from lognormal
# freq <- rlnorm(nSim, mulog.anz, sdlog.anz.mod) # Random sampling from lognormal, modified

## DAS DAUERT!
for(i in 1:nSim) # Convolution with Monte Carlo method
  sj[i] = sum(rlnorm.gpd(n = freq[i], theta = c(theta1, theta2),
                         theta.gpd = c(theta1.tail, theta2.tail,
                                     theta3.tail), u = H))

## Jahresschadenverteilung
truehist(sj/1e6, col = '#bdd7e7', prob = TRUE, xlim = c(0, 200),
         # main = "Simulierte jaehrliche Elementarschadenverteilung\n(ohne copula)",
         main = "Simulierte jährliche Hagelschadenverteilung",
         ylab = 'relative frequency density', xlab = 'Jährlicher Schaden [Mio. CHF]')

## Risikokennzahlen:
## Value at Risk 99.9% (1x in 1'000 Jahren)
VaR_1000 <- quantile(sj, 0.999)
## Expected shortfall 99.9% (durschnitt aller schaeden seltener als 1'000 Jahren)
ExS_1000 <- mean(sj[which(sj > vaR_1000)])

## Value at Risk 99.5% (1x in 200 Jahren)
VaR_200 <- quantile(sj, 0.995)
## Expected shortfall 99.5% (durschnitt aller schaeden seltener als 200 Jahren)
ExS_200 <- mean(sj[which(sj > vaR_200)])

## Expected loss
ExLoss <- quantile(sj, 0.5) 

legend('topright',
       legend = c(paste('Value at Risk zum Niveau 99.5% (1x in 200 Jahren):',
                        round(VaR_200 / 1e6, 0), 'Mio CHF'),
                  paste("Expected Shortfall zum Niveau 99.5%:",
                        round(ExS_200 / 1e6, 0), 'Mio CHF'),
                  paste('Expected loss:',
                        round(ExLoss / 1e6, 0), 'Mio CHF'),
                  paste('Rechnungsdatum:', format(Sys.time(), '%Y-%m'))),
       bty = "n")

# ## ----
# ## Save the simulated and the yealy data as R Objects for Risk Modelling
# ## Part II: Overal yearly loss with Copula
# saveRDS(sj, "sj_hagel.rds"); saveRDS(yearly.dat, "yearly_dat_hagel.rds")
# ## ----
