## ----
## Jahresschaden Modellierung, Part I
## Schaetzung jaehrlicher Elementarschaeden mit einer Lognormalverteilung 
## (besserer Fit als mit einer Gammaverteilung) und einem Paretomodell im Tail.
## Simulation pro Gefahr Überschwemmung, Sturm und Hagel separat.
##
## Modellierung anhand der Schadenerfahrung der GVZ
##
## Mirco Heidemann
## Mai 2018
## ----
## Funktionen laden
library(MASS)
library(evir)
# library(evd)
library(QRM)
library(stringr)
library(dplyr)

## --- MANUELLE ANPASSUNGEN
## Welche jahresschaeden sollen modelliert werden?
## hazard_char = "Hochwasser", "Sturm", "Hagel"
hazard_char <- "Hochwasser"
## Number of simulated years
n_sim_year = 1e5

# title_char <- "Überschwemmungsschadenverteilung"
title_char <- paste0(hazard_char, "schadenverteilung")

## Tail der simulierten schadenverteilung begrenzen mit groest moeglichem
## Schadenszenario.
max_loss <- case_when(
  hazard_char == "Hochwasser" ~ 5e9,
  hazard_char == "Sturm" ~ 2e9,
  hazard_char == "Hagel" ~ 550e6)
## anzahl schaeden nicht hoeher als gesamt anzahl gebaeude
fact_max_geb <- 3e5

## varianz der schadenfrequenz erhoehen
fact_sdfreq <- case_when(
  hazard_char == "Hochwasser" ~ 2,
  hazard_char == "Sturm" ~ 1.5,
  hazard_char == "Hagel" ~ 1)
## gesamt schadenfrequenz erhoehen
fact_freq <- case_when(
  hazard_char == "Hochwasser" ~ 5,
  hazard_char == "Sturm" ~ 2,
  hazard_char == "Hagel" ~ 1.5)

## varianz extremer einzelschaeden erhoehen (tail)
fact_tail_shape <- case_when(
  hazard_char == "Hochwasser" ~ 1.2,
  hazard_char == "Sturm" ~ 1.4,
  hazard_char == "Hagel" ~ 1)

## varianz nicht-extremer einzelschaeden erhoehen (body)
fact_body_shape <- case_when(
  hazard_char == "Hochwasser" ~ 2,
  hazard_char == "Sturm" ~ 1.5,
  hazard_char == "Hagel" ~ 1)
## ---

## relative pfade spezifizieren
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pth_data <- ("./data/")
pth_rdata <- ("./data/Rdata/")
pth_pres <- ("./presentation/")
pth_fun <- ("../R_functions/")
pth_port <- ("../GVZ_Portfolio/")

## Funktion fuer die Aufbereitung der einzel schaeden (Feuer, Elementar) laden
source(paste0(pth_fun, 'Fun_SchadPrepAll.R'))

## Daten laden, aufbereiten und indexieren
schad.file <- paste0(pth_port, 'schaeden_201801_georef.csv')
index.file <- paste0(pth_port, 'versicherungsindex_gvz.csv')
schad <- schadPrep(schad.file, index.file)

## nur elementarschaeden ab 1982 mit (indexierten) schaeden groesser 0
schad_elementar <- schad %>% filter(format(schad_datum, "%Y") > 1981,
                                    # schadendaten bis und mit letzten Jahres
                                    format(schad_datum, "%Y") < 2018,
                                    schad_art == "Elementar",
                                    schad_index > 0)
## welche gefahr soll modelliert werden?
df <- schad_elementar %>% filter(str_detect(schad_code, hazard_char))

## aggregiere die einzenlschaeden zu jahresschaden
aggr_anz <- aggregate(df$schad_index,
                      by = list(as.numeric(format(df$schad_datum,'%Y'))),
                      'length')
aggr_sum <- aggregate(df$schad_index,
                      by = list(as.numeric(format(df$schad_datum,'%Y'))),
                      'sum')
yearly_dat <- as.data.frame(cbind(aggr_sum, aggr_anz$x))
names(yearly_dat) <- c('jahr', 'schadSum', 'schadAnz')

# ## The GPD fits a distribution to exceedences over a threshold.
# ## We therefore need to pick an appropriate threshold. A useful
# ## plot to assist in this choice is the mean excess plot.

## 95% Quantile as an appropriate threshold for the GPD fit...
## "5% der Schaeden machen 50% der Schadensumme aus"
threshold <- quantile(df$schad_index, 0.95)

schadSum_body <- df %>% filter(schad_index <= threshold) %>% 
  select(schad_index)
# schadSum_tail <- df %>% filter(schad_index > threshold) %>% 
#   select(schad_index)

# par(mfrow=c(1,2))
# evir::meplot(df$schad_index) # From evir package
# abline(v = threshold, col = "red")
# evd::mrlplot(df$schad_index) # From evd package
# abline(v = threshold, col = "red")
# evir::shape(df$schad_index) # From evir package
# par(mfrow=c(1,1))

## A) Schadensumme der nicht-Extremen Einzelschaeden modelieren mit einer
##    lognormal Verteilung. Paramterschaetzung mit Maximum-likelihood fitting
##    (eingabe nicht logarithmiert).

## beruecksichtige nur schaeden ueber 25% Quantil fuer die parameter schaetzung
schadSum_body <- schadSum_body %>% filter(schad_index >
                                            quantile(schad_index, 0.25))

fit_loss_lnorm <- fitdistr(schadSum_body$schad_index, "lognormal")
meanlog_loss <- fit_loss_lnorm$estimate["meanlog"]
sdlog_loss <- fit_loss_lnorm$estimate["sdlog"]

# ## Diagnosis: Wie gut passt die gewaehlte Verteilung?
# ## qq-plot, normal fit der log-es?
# qqnorm(log(schadSum_body$schad_index)); qqline(log(schadSum_body$schad_index))
# ## histogram und dichte der lognormal es
# truehist(schadSum_body$schad_index, nbins = 50, col= 8, xlim = c(0, 6e4),
#          main = "histogram and lognormal density")
# curve(dlnorm(x, meanlog = meanlog_loss, sdlog = sdlog_loss),
#       add = T, col = 'darkblue', lwd = 2)
# 
# ## Draw cumulative density functions
# plot(ecdf(schadSum_body$schad_index), cex = 0, xlim = c(0, 5e4))
# curve(plnorm(x, meanlog = meanlog_loss,
#              sdlog = sdlog_loss),
#       add = T, col = 2)
# curve(plnorm(x, meanlog = meanlog_loss,
#              sdlog = (fact_body_shape * sdlog_loss)),
#       add = T, col = 3)
# legend('right', lwd = 1, col = c(1,2,3), legend = c("empirical CDF",
#                                                   "lognormal CDF",
#                                                   "modified CDF"))

## B) Anzahl der Schaeden pro Jahr schaetzen mit einer lognormal-Verteilung
##    (besser als mit einer Poisson Verteilung)

## beruecksichtige nur jaerhliche anzahl schaeden ueber 25% Quantil
yearly_dat <- yearly_dat %>% filter(schadAnz >
                                            quantile(schadAnz, 0.25))

fit_freq_lnorm <- fitdistr(yearly_dat$schadAnz, "lognormal")
meanlog_freq <- fit_freq_lnorm$estimate["meanlog"]
sdlog_freq <- fit_freq_lnorm$estimate["sdlog"]

# ## Diagnosis: Wie gut passt die gewaehlte Verteilung?
# ## qq-plot, normal fit der log-Anzahl?
# qqnorm(log(yearly_dat$schadAnz)); qqline(log(yearly_dat$schadAnz))
# ## histogram und dichte der lognormal es
# truehist(yearly_dat$schadAnz, nbins = 50, col = 8, xlim = c(0, 1e3),
#          main="histogram and lognormal density")
# curve(dlnorm(x, meanlog = meanlog_freq, sdlog = sdlog_freq),
#       add = T, col = 'darkblue', lwd = 2)
# curve(dlnorm(x, meanlog = meanlog_freq, sdlog = (sdlog_freq * fact_sdfreq)),
#       add = T, col = 'red', lwd = 2)
# legend('topright',
#        c('density with empirical variance',
#          'density with modificated variance'),
#        lwd = 2, col = c('darkblue', 'red'), bty = 'n')
# 
# # Draw cumulative density functions
# plot(ecdf(yearly_dat$schadAnz), cex = 0, xlim = c(0, 2e3),
#      main = 'CDF Lognormal SchadAnz')
# curve(plnorm(x, meanlog = meanlog_freq,
#              sdlog = sdlog_freq),
#       add = T, col = 'red')
# curve(plnorm(x, meanlog = meanlog_freq,
#              sdlog = (sdlog_freq * fact_sdfreq)),
#       add = T, col = 'green')
# legend('right', lwd = 1, col = c('black', 'red', 'green'),
#        legend = c('empirical CDF', 'lognormal CDF',
#                   'lognormal CDF modified varince'),
#        bty = 'n')

## C) The Generalised Pareto Distribution (GPD) describes exceedences over a
## threshold and is convenient for modelling insurance claims:
## W'keit eines Schadens zw. c (z.B. 50'000 oder ein anderer hoher Wert) und x,
## respektive W'keit dass ein Schaden kleiner oder hoechstens gleich x ist:
## P(X<=x) = 1 - (c/x)^alpha
## c (oder beta), scale parameter: kleinster, von der Paretoverteilung, in
##   Betracht gezogener Schaden
## alpha (oder xi), shape parameter: Paretoparameter (erfahrungsgemaess zw. 0  
##       und 1 fuer Elementarschaeden, zw. 1 und 2 fuer Feuerschaeden)

## Schaetzung des Paretoparameters nach Maximum Likelihood:
## alpha = n/sum(ln(xi/c))
# pareto.MLE <- function(X)
# {
#   n <- length(X)
#   m <- min(X)
#   a <- n/sum(log(X/m))
#   return( c(m,a) ) 
# }
# pareto.MLE(dat)
## oder mit der gpdf funktion aus dem evir package

fit_gpd <- gpd(df$schad_index, threshold) # From evir package
xi <- fit_gpd$par.ests[1] ## shape parameter
beta <- fit_gpd$par.ests[2] ## scale parameter

# ## diagnostics mit dem QQ plot und der CDF
# par(mfrow=c(1,2))
# qplot(df$schad_index, xi = xi)
# # CDF Plot
# plot(ecdf(df$schad_index[which(df$schad_index > threshold)]),
#      cex = 0, main = "CDF Plot", xlim = c(0, 1e6))
# lines(sort(df$schad_index),
#       evir::pgpd(sort(df$schad_index), xi = xi, beta = beta, mu = threshold),
#       col = "red")
# par(mfrow=c(1,1))
# ## plot of the tail of the underlying distribution of the data
# gpd.q(tailplot(fit_gpd), 0.995)

## D) Statistische Faltung fuer die jaehrliche schadensverteilung
##    (Convolution Method to obtain the annual loss distribution)

# Quantile function of lognormal-GPD severity distribution
qlnorm.gpd = function(p, theta, theta.gpd, u)
{
  # Fu: W'keit, dass eine logNorm verteilte Zufallszahl kleiner ist als treshold
  Fu = plnorm(u, meanlog = theta[1], sdlog = theta[2])
  ## fuer die kleinen w'keiten (aus uniformvert.), simuliere aus GPD
  x = ifelse(p < Fu,
             qlnorm(p = p, meanlog = theta[1], sdlog = theta[2]),
             evir::qgpd(p = (p - Fu) / (1 - Fu) , xi = theta.gpd[1],
                   mu = theta.gpd[2], beta = theta.gpd[3]))
  return(x)
}
# Random sampling function of lognormal-GPD severity distribution
rlnorm.gpd = function(n, theta, theta.gpd, u)
  ## runif: n Zahlen zwischen 0 und 1
{ r = qlnorm.gpd(runif(n), theta, theta.gpd, u)}
# set.seed(1000)

## Parameter definition:

## Number of simulated annual losses
n_sim = n_sim_year
## Threshold body-tail
H = threshold
## Parameter mu of lognormal (body)
theta1 = meanlog_loss
## Parameter sigma of lognormal (body), modifiziert mit "fact_body_shape" 
theta2 = fact_body_shape * sdlog_loss
## Shape parameter of GPD (tail), modifiziert mit "fact_tail_shape"
theta1.tail = fact_tail_shape * xi
## Location parameter of GPD (tail)
theta2.tail = H
## Scale parameter of GPD (tail)
theta3.tail = beta
## Vektor Inizialisierung fuer die Annual loss distribution
sj = rep(0, n_sim)

## Number of claims per year: Random sampling from lognormal
## modifiziert mit "fact_sdfreq"
freq = fact_freq * round(rlnorm(n_sim, meanlog_freq,
                                (fact_sdfreq * sdlog_freq)), 0)
## nicht hoeher als gesamt anzahl gebaeude
ind <- which(freq > fact_max_geb)
freq[ind] <- fact_max_geb

## DAS DAUERT!
for(i in 1:n_sim) # Convolution with Monte Carlo method
  sj[i] = sum(rlnorm.gpd(n = freq[i], theta = c(theta1, theta2),
                         theta.gpd = c(theta1.tail, theta2.tail,
                                     theta3.tail), u = H))

## Tail begrenzen mit den groesst moeglichen Schadensszenarien
ind <- which(sj > max_loss)
sj[ind] <- max_loss

## Risikokennzahlen:
## Value at Risk 99.5% (1x in 200 Jahren)
VaR_200 <- quantile(sj, 0.995)
## Expected shortfall 99.5% (Durschnitt aller Schaeden seltener als 200 Jahren)
ExS_200 <- mean(sj[which(sj > VaR_200)])

## Expected loss
ExLoss <- quantile(sj, 0.5) 

## Histogramm der simulierten Jahresschadenverteilung
truehist(sj/1e6, col = '#bdd7e7', prob = TRUE, xlim = c(0, 600),
         main = paste0("Simulierte jährliche ", title_char),
         ylab = 'relative frequency density',
         xlab = 'Jährlicher Schaden [Mio. CHF]')

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
summary(rlnorm(n_sim, theta1, theta2))

## Verteilung der Simulierten extremen Schaeden:
summary(evir::rgpd(n_sim, beta = theta3.tail, xi = theta1.tail, mu = H))

## Verteilung der Simulierten Anzahl Schaeden in Tausend
quantile(freq, c(0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 1))/1e3

## simulierter schaden in Millionen
quantile(sj, c(0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 1))/1e6

# ## ----
# ## Save the simulated and the yealy data as R Objects for loss modelling
# ## Part II: Overal yearly loss with Copula
# saveRDS(sj, paste0(pth_rdata, "sj_", hazard_char, ".rds"))
# saveRDS(yearly_dat, paste0(pth_rdata, "yearly_dat_", hazard_char, ".rds"))
# ## ----


## ----
## Jahresschaden Modellierung, Part II
## Simuliere gesamte Schadenverteilung - ueber alle Gefahren.
## Abhaengigkeit der simulierten Schadenverteilungen der
## einzelnen Gefahren (aus Part I) mit copula modelliert.
##
## Mirco Heidemann
## 05/2018
## ----

## Funktionen laden
# library(PerformanceAnalytics)

# ## Load simulated and the yealy data
# sj_sturm <- readRDS(paste0(pth_rdata, "sj_Sturm.rds"))
# yearly_dat_sturm <- readRDS(paste0(pth_rdata, "yearly_dat_Sturm.rds"))
# sj_hagel <- readRDS(paste0(pth_rdata, "sj_Hagel.rds"))
# yearly_dat_hagel <- readRDS(paste0(pth_rdata, "yearly_dat_Hagel.rds"))
# sj_wasser <- readRDS(paste0(pth_rdata, "sj_Hochwasser.rds"))
# yearly_dat_wasser <- readRDS(paste0(pth_rdata, "yearly_dat_Hochwasser.rds"))
# 
## Loss distribution risk class 1
# s1 = rlnorm(n = n_sim, meanlog = 4.5, sdlog = 2.3)
## Loss distribution risk class 2
# s2 = rlnorm(n = n_sim, meanlog = 5, sdlog = 2.5)
# s1 = sj_sturm # Loss distribution risk class sturm
# s2 = sj_hagel # Loss distribution risk class hagel
# s3 = sj_wasser # Loss distribution risk class wasser
# 
# # VaR.s1 = quantile(s1, 0.999) # VaR risk class 1
# # VaR.s2 = quantile(s2, 0.999) # VaR risk class 2
# # corr = 0.6 # Correlation among risk classes
# 
# ## beobachtete jahresschaeden summe pro gefahr in ein data frame
# yearlyDat <- data.frame(jahr = yearly_dat_sturm[,1],
#                         sturm = yearly_dat_sturm[,2],
#                         hagel = yearly_dat_hagel[,2],
#                         wasser = yearly_dat_wasser[,2])
# 
# ## korrelations matrix der einzelnen gefahren schaetzen
# ## keine negative Korrelation zulassen
# corr_matrix <- abs(cor(yearlyDat[,-1], method = "pearson"))
# # chart.Correlation((yearlyDat[,-1])) # from PerformanceAnalytics
# 
# # ## korrelations matrix mit zwei variablen
# # corr = cor(yearly_dat_Sturm$schadAnz, yearly_dat_hagel$schadAnz)
# # corr_matrix = matrix(data=c(1,corr,corr,1), nrow=2) # correlation matrix
# 
# dof <- 5 # degrees of freedom
# n_sim <- 1e6 # Number of simulated overall annual losses
# # Simulation from Student-t copula
# sim_copula_t <- rcopula.t(n = n_sim, df = dof, Sigma = corr_matrix)
# 
# ## !! Vielleicht besser Bernstein copula (nicht-parametrisch)? !!
# 
# ## overall annual loss distribution:
# ## sample quantiles for s1, s2, ..., with probabilities from sim_copula_t
# s <- quantile(s1, probs = sim_copula_t[,1]) +
#   quantile(s2, probs = sim_copula_t[,2]) +
#   quantile(s3, probs = sim_copula_t[,3])
# 
# ## Risikokennzahlen:
# ## Value at Risk 99.5% (1x in 200 Jahren)
# VaR_200 <- quantile(s, 0.995)
# ## Expected shortfall 99.5% (durschnitt aller schaeden seltener als 200 Jahre)
# ExS_200 <- mean(s[which(s > VaR_200)])
# 
# ## Expected loss
# ExLoss <- quantile(s, 0.5)
# 
# ## Value at Risk 99.9% (1x in 1'000 Jahren)
# VaR_1000 <- quantile(s, 0.999)
# ## Expected shortfall 99.9% (durschnitt aller schaeden seltener als 1'000 Jahre)
# ExS_1000 <- mean(s[which(s > VaR_1000)])
# 
# ## histogram of overall annual loss distribution
# truehist(s/1e6, col= '#bdd7e7', prob = TRUE, xlim=c(0, 600),
#          main="Simulierte jährliche Elementarschadenverteilung\n(Student-t copula)",
#          ylab='relative frequency density', xlab='jährlicher Schaden [Mio. CHF]')
# 
# legend('topright',
#        legend = c(paste('VaR zum Niveau 99.5% (1x in 200 Jahren):',
#                         round(VaR_200 / 1e6, 0), 'Mio CHF'),
#                   paste("ExS zum Niveau 99.5%:",
#                         round(ExS_200 / 1e6, 0), 'Mio CHF'),
#                   paste('Expected loss:',
#                         round(ExLoss / 1e6, 0), 'Mio CHF'),
#                   paste('Rechnungsdatum:', format(Sys.time(), '%Y-%m'))),
#        bty = "n", cex = 0.8)
# 
# ## simulierter schaden in Millionen
# quantile(s, c(0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 1))/1e6
