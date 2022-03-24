## ----
## Risk Modelling, Part I
## Schaetzung jaehrlicher Elementarschaeden mit einer Lognormalverteilung und 
## einem Paretomodell im tail
## Pro Gefahr einzel oder gleich als gesamte Elementarschaedenverteilung
##
## Modellierung anhand der Schadenerfahrung der GVZ
##
## Mirco Heidemann
## 05/2018
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
title_char <- "Hagelschadenverteilung"

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

## 95% Quantile as an appropriate threshold for the GPD fit?
threshold <- quantile(df$schad_index, 0.95)

schadSum.body <- df %>% filter(schad_index <= threshold) %>% 
  select(schad_index)
schadSum.body <- as.matrix(schadSum.body)
schadSum.tail <- df %>% filter(schad_index > threshold) %>% 
  select(schad_index)
schadSum.tail <- as.matrix(schadSum.tail)

## A) Schadensumme der nicht Extremen modelieren mit einer lognormal Verteilung
# Paramterschaetzung mit Maximum-likelihood fitting, eingabe nicht logarithmiert

## !! Estimate parameters of Frechet distribution? !!
## !! Frechet distribution provided in VGAM, RTDE, ReIns, extraDistr and evd packages 

## !! Estimate parameters of Gamma distribution? !!

## beruecksichtige nur ueber 33% Quantil
q_min <- quantile(schadSum.body, 0.33)
ind <- which(schadSum.body > q_min)

fit.ln <- fitdistr(schadSum.body[ind], "lognormal")
mulog <- fit.ln$estimate["meanlog"]
sdlog <- fit.ln$estimate["sdlog"]

# Draw cumulative density functions
plot(ecdf(schadSum.body), cex = 0, xlim = c(0, 1e5))
curve(plnorm(x, meanlog = fit.ln$estimate[1],
             sdlog = fit.ln$estimate[2]),
      add = T, col = 'red')
legend('right', lwd = 1, col = c(1,2), legend = c("empirical CDF", "lognormal CDF"))

## B) Anzahl der Schaeden pro Jahr schaetzen
## ... mit einer gamma-Verteilung, scale data first, than then scaling the rate

## beruecksichtige nur ueber 33% Quantil
q_min <- quantile(yearly.dat$schadAnz, 0.33)
ind <- which(yearly.dat$schadAnz > q_min)
## scaled the data by 0.1
fit_gamma <- fitdistr(yearly.dat$schadAnz[ind] / 1000, densfun = "gamma")
## scaled the rate by 0.1
rate_gamma <- coef(fit_gamma)["rate"] / 1000
shape_gamma <- coef(fit_gamma)["shape"]

## C) The Generalised Pareto Distribution (GPD) describes exceedences over a
## threshold and is convenient for modelling insurance claims.
fit.gpd <- gpd(df$schad_index, threshold) # From evir package
xi <- fit.gpd$par.ests[1] ## shape parameter
beta <- fit.gpd$par.ests[2] ## scale parameter

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
# set.seed(1000)

## nSim = 1e6 # Number of simulated annual losses
## fuer modellentwicklung
nSim = 1e4 # Number of simulated annual losses

H = threshold # Threshold body-tail
# lambda = lambda # Parameter of Poisson body
theta1 = mulog # Parameter mu of lognormal (body)
# theta2 = sdlog # Parameter sigma of lognormal (body)

## !! Parameter sigma of lognormal (body), manuell erhoet !!
## Wasser: theta2 = 3 * sdlog
## Hagel: theta2 = 1 * sdlog
## Sturm: theta2 = 3 * sdlog
theta2 = 1 * sdlog # Parameter sigma of lognormal (body)

## !! Shape parameter of GPD, manuell erhoet !!
# theta1.tail = xi # Shape parameter of GPD (tail)
## Wasser: theta1.tail = 1.3 * xi
## Hagel: theta1.tail = 1.15 * xi
## Sturm: theta1.tail = 1.3 * xi
theta1.tail = 1.15 * xi # Shape parameter of GPD (tail)

theta2.tail = H # Location parameter of GPD (tail)
theta3.tail = beta # Scale parameter of GPD (tail)
sj = rep(0, nSim) # Annual loss distribution inizialization

## !! Random sampling from gamma, manuell erhoet !!
## Wasser: 3
## Hagel: 4
## Sturm: 3
freq = 3 * round(rgamma(nSim, shape_gamma, rate_gamma), 0)
## Random sampling from gamma
# freq = round(rgamma(nSim, shape_gamma, rate_gamma), 0)

## DAS DAUERT!
for(i in 1:nSim) # Convolution with Monte Carlo method
  sj[i] = sum(rlnorm.gpd(n = freq[i], theta = c(theta1, theta2),
                         theta.gpd = c(theta1.tail, theta2.tail,
                                     theta3.tail), u = H))

## Risikokennzahlen:
## Value at Risk 99.5% (1x in 200 Jahren)
VaR_200 <- quantile(sj, 0.995)
## Expected shortfall 99.5% (durschnitt aller schaeden seltener als 200 Jahren)
ExS_200 <- mean(sj[which(sj > VaR_200)])

## Expected loss
ExLoss <- quantile(sj, 0.5) 

## Value at Risk 99.9% (1x in 1'000 Jahren)
VaR_1000 <- quantile(sj, 0.999)
## Expected shortfall 99.9% (durschnitt aller schaeden seltener als 1'000 Jahren)
ExS_1000 <- mean(sj[which(sj > VaR_1000)])

## Jahresschadenverteilung
truehist(sj/1e6, col = '#bdd7e7', prob = TRUE, xlim = c(0, 600),
         # main = "Simulierte jaehrliche Elementarschadenverteilung\n(ohne copula)",
         main = paste0("Simulierte jährliche ", title_char),
         ylab = 'relative frequency density', xlab = 'Jährlicher Schaden [Mio. CHF]')

legend('topright',
       legend = c(paste('VaR zum Niveau 99.5% (1x in 200 Jahren):',
                        round(VaR_200 / 1e6, 0), 'Mio CHF'),
                  paste("ExS zum Niveau 99.5%:",
                        round(ExS_200 / 1e6, 0), 'Mio CHF'),
                  paste('Expected loss:',
                        round(ExLoss / 1e6, 0), 'Mio CHF'),
                  paste('Rechnungsdatum:', format(Sys.time(), '%Y-%m'))),
       bty = "n", cex = 0.8)

## Verteilung der Simulierten nicht-extremen Schaeden:
summary(rlnorm(nSim, theta1, theta2))

## Verteilung der Simulierten extremen Schaeden:
summary(rgpd(nSim, beta = theta3.tail, xi = theta1.tail, mu = H))

## Verteilung der Simulierten Anzahl Schaeden in Tausend
quantile(freq, c(0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 1))/1e3

## simulierter schaden in Millionen
quantile(sj, c(0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 1))/1e6

# ## ----
# ## Save the simulated and the yealy data as R Objects for Risk Modelling
# ## Part II: Overal yearly loss with Copula
# saveRDS(sj, "sj.hagel.rds"); saveRDS(yearly.dat, "yearly.dat.hagel.rds")
# ## ----
