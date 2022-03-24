## ----
## Risk Modelling, Part II
## Simuliere gesamte Schadenverteilung - ueber alle Gefahren.
## Abhaengigkeit der simulierten Schadenverteilungen der
## einzelnen Gefahren (aus Part I) mit copula modelliert.
##
## Mirco Heidemann
## 10/2017
## ----

## relative pfade spezifizieren
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pth_data <- ("./data/")
pth_rdata <- ("./data/Rdata/")
pth_pres <- ("./presentation/")

## Funktionen laden
require(MASS)
# require(PerformanceAnalytics)
require(QRM)

## Load simulated and the yealy data
sj.sturm <- readRDS(paste0(pth_rdata, "sj.sturm.rds"))
yearly.dat.sturm <- readRDS(paste0(pth_rdata, "yearly.dat.sturm.rds"))
sj.hagel <- readRDS(paste0(pth_rdata, "sj.hagel.rds"))
yearly.dat.hagel <- readRDS(paste0(pth_rdata, "yearly.dat.hagel.rds"))
sj.wasser <- readRDS(paste0(pth_rdata, "sj.wasser.rds"))
yearly.dat.wasser <- readRDS(paste0(pth_rdata, "yearly.dat.wasser.rds"))

# s1 = rlnorm(n=nSim, meanlog=4.5, sdlog=2.3) # Loss distribution risk class 1
# s2 = rlnorm(n=nSim, meanlog=5, sdlog=2.5) # Loss distribution risk class 2
s1 = sj.sturm # Loss distribution risk class sturm
s2 = sj.hagel # Loss distribution risk class hagel
s3 = sj.wasser # Loss distribution risk class wasser

# VaR.s1 = quantile(s1, 0.999) # VaR risk class 1
# VaR.s2 = quantile(s2, 0.999) # VaR risk class 2
# corr = 0.6 # Correlation among risk classes

## beobachtete jahresschaeden summe pro gefahr in ein data frame
yearlyDat <- data.frame(jahr = yearly.dat.sturm[,1],
                        sturm = yearly.dat.sturm[,2],
                        hagel = yearly.dat.hagel[,2],
                        wasser = yearly.dat.wasser[,2])

## korrelations matrix der einzelnen gefahren schaetzen
corrMatrix <- cor(yearlyDat[,-1], method = "pearson")
# chart.Correlation(yearlyDat[,-1])

# ## korrelations matrix mit zwei variablen
# corr = cor(yearly.dat.sturm$schadAnz, yearly.dat.hagel$schadAnz)
# corrMatrix = matrix(data=c(1,corr,corr,1), nrow=2) # correlation matrix

dof <- 5 # degrees of freedom
nSim <- 1e6 # Number of simulated overall annual losses
# Simulation from Student-t copula
simCopulaT <- rcopula.t(n = nSim, df = dof, Sigma = corrMatrix)

## !! Vielleicht besser Bernstein copula (nicht-parametrisch)? !!

## overall annual loss distribution:
## sample quantiles for s1, s2, ..., with probabilities from simCopulaT
s <- quantile(s1, probs = simCopulaT[,1]) + 
  quantile(s2, probs = simCopulaT[,2]) +
  quantile(s3, probs = simCopulaT[,3])

## Risikokennzahlen:
## Value at Risk 99.5% (1x in 200 Jahren)
VaR_200 <- quantile(s, 0.995)
## Expected shortfall 99.5% (durschnitt aller schaeden seltener als 200 Jahren)
ExS_200 <- mean(s[which(s > VaR_200)])

## Expected loss
ExLoss <- quantile(s, 0.5) 

## Value at Risk 99.9% (1x in 1'000 Jahren)
VaR_1000 <- quantile(s, 0.999)
## Expected shortfall 99.9% (durschnitt aller schaeden seltener als 1'000 Jahren)
ExS_1000 <- mean(s[which(s > VaR_1000)])

## histogram of overall annual loss distribution
truehist(s/1e6, col= '#bdd7e7', prob = TRUE, xlim=c(0, 600),
         main="Simulierte jährliche Elementarschadenverteilung\n(Student-t copula)",
         ylab='relative frequency density', xlab='jährlicher Schaden [Mio. CHF]')

legend('topright',
       legend = c(paste('VaR zum Niveau 99.5% (1x in 200 Jahren):',
                        round(VaR_200 / 1e6, 0), 'Mio CHF'),
                  paste("ExS zum Niveau 99.5%:",
                        round(ExS_200 / 1e6, 0), 'Mio CHF'),
                  paste('Expected loss:',
                        round(ExLoss / 1e6, 0), 'Mio CHF'),
                  paste('Rechnungsdatum:', format(Sys.time(), '%Y-%m'))),
       bty = "n", cex = 0.8)


