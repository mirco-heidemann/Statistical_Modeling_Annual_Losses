## Risk Modelling
## Schaetzung jaehrlicher Elementarschaeden mit einem Paretomodell
## Frequenzextrapolation, erwartete Schadenlast, Schwankungszuschlag
##
## Mirco Heidemann
## 05/2016


## Arbeitsverzeichnis definieren
setwd('I:/R/Statistische Modellierung Jahresschäden/data/')
## Funktionen laden
require(MASS)
require(evd)
require(evir)

## --------------------BENUTZEREINGABEN---------------------------------------------
## treshold...
## - 6e4 fuer sturmschaeden
## - 7e4 fuer hagelschaeden
## - 9e4 fuer ueberschwemmungsschaeden
## ... und 8e4 fuer alle Elementarschaeden zusammen
tresh <- 6e4

## welche gefahr?
## sturm, wasser, hagel
ng <- 'sturm'

## Jaehrlichkeit des groessten, beobachtete ereignis (zB lothar fuer sturm)
RTP <- 80

## aktuelle GVZ Versicherungssumme in CHF
vs.gvz <- 473e9
## ---------------------------------------------------------------------------------

## Funktion fuer die Aufbereitung der Elementarschaeden laden
source('J:/Naturgefahren/FG2_Datenanalysen/Statistische Analysen/Rfunctions/f.schadPrep.R')

## Daten laden und aufbereiten
schad.file <- paste('J:/Naturgefahren/FG2_Datenanalysen/Portfolio gvz/gebäudebestand gvz/',
                    '2017/GemDat Export/Schaeden 2016.csv',
                    sep='')
index.file <- 'J:/Naturgefahren/FG2_Datenanalysen/Portfolio gvz/versicherungsindex.gvz.csv'

## Funktion fuer die Aufbereitung der Elementarschaeden
schad <- schadPrep(schad.file, index.file)

## schadensumme ueber NUll
ind <- which(schad$indexSchad > 0)
if(length(ind)>0) schad <- schad[ind,]

## schadendaten ab 1982
ind <- which(as.numeric(format(schad$schadat,'%Y')) > 1981)
if(length(ind)>0) schad <- schad[ind,]

# ## Finde die treshold pro Gefahr:
# ## The GPD fits a distribution to exceedences over a threshold.
# ## We therefore need to pick an appropriate threshold. A useful
# ## plot to assist in this choice is the mean excess plot.
# ind <- which(schad$artcode == '1 Sturm')
# sturm <- schad$indexSchad[ind]
# ind <- which(schad$artcode == '2 Hagel')
# hagel <- schad$indexSchad[ind]
# ind <- which(schad$artcode == '3 Hochwasser, ueberschwemmung')
# wasser <- schad$indexSchad[ind]
# 
# dat <- sturm
# quantile(dat, 0.995)
# 
# u <- 9e4
# ind <- which(dat>u)
# ## body
# truehist(dat[-ind])
# fit.ln <- fitdistr(dat[-ind], "lognormal")
# mulog <- fit.ln$estimate["meanlog"]
# sdlog <- fit.ln$estimate["sdlog"]
# curve(dlnorm(x, meanlog=mulog,sdlog=sdlog),
#       add=T, col='darkblue', lwd=2)
# ## tail
# truehist(dat[ind])
# 
# fit.gpd <- gpd(dat[ind], u) # From evir package
# xi <- fit.gpd$par.ests[1]
# beta <- fit.gpd$par.ests[2]
# curve(dgpd(x, xi = xi, mu = 0, beta = beta), add=T, col='darkblue', lwd=2)
# 
# shape(dat) # From evir package
# 
# meplot(dat, xlim=c(0,3e5)) # From evir package
# abline(v=u, col = 'red')
# mrlplot(dat, col=c("green","black","green"),nt=200,xlim=c(0,3e5)) # From evd package
# abline(v=u, col = 'red')


## Frequenzextraploation pro Gefahr, bis jetzt nur fuer Sturm
pareto.ES <- function(dat, c, Gefahr, rtp)
{
  #dat <- schad; c <- tresh; Gefahr <- ng; rtp <- RTP ## LOESCHEN nach Modellentwicklung!!!!
  
  ## Sturm
  if(Gefahr == 'sturm'){
    ind <- which(dat$artcode == '1 Sturm')
    if(length(ind)>0) dat <- dat[ind,]
    ## Lothar, groesster Sturmschaden
    ind <- which(format(dat$schadat, '%Y') == '1999')
    schad.lothar <- sum(dat$indexSchad[ind])
    schad.ereignis <- schad.lothar
  }
  # ## Ueberschwemmung
  # if(Gefahr == 2){
  #   ind <- which(dat$artcode == '3 Hochwasser, ueberschwemmung')
  #   if(length(ind)>0) dat <- dat[ind,]
  # }  
  # ## Hagel
  # if(Gefahr == 3){
  #   ind <- which(data$artcode == '2 Hagel')
  #   data <- data[ind,] 
  # }

  u <- c ## Treshold
  ind <- which(dat$indexSchad >= u)
  if(length(ind)>0) dat <- dat[ind,]
  dat <- dat$indexSchad
  ## histogram und dichte der lognormal es
  truehist(dat)
  # Draw cumulative density functions
  plot(ecdf(dat), cex=0, xlim=c(u, u*10))
  
  ## Paretoverteilung:
  ## W'keit eines Schadens zw. c (z.B. 50'000 oder ein anderer hoher Wert) und x,
  ## respektive W'keit dass ein Schaden kleiner oder hoechstens gleich x ist:
  ## P(X<=x) = 1 - (c/x)^alpha
  ## c: kleinster, von der Paretoverteilung, in Betracht gezogener Schaden
  ##    (scale parameter)
  ## alpha: Paretoparameter (erfahrungsgemaess zw. 0 und 1 fuer Elementarschaeden,
  ##        zw. 1 und 2 fuer Feuerschaeden)
  ##        (shape parameter)
  
  ## Schaetzung des Paretoparameters nach Maximum Likelihood: alpha = n/sum(ln(xi/c))
  # pareto.MLE <- function(X)
  # {
  #   n <- length(X)
  #   m <- min(X)
  #   a <- n/sum(log(X/m))
  #   return( c(m,a) ) 
  # }
  # pareto.MLE(dat)
  
  ## oder
  fit.gpd <- gpd(dat, u) # From evir package
  a <- as.numeric(fit.gpd$par.ests[1])
  b <- as.numeric(fit.gpd$par.ests[2])
  # plot(fit.gpd)
  # ## oder ein MLE Fitting mit dem package extRemes
  # require(extRemes)
  # pot.ext <- fevd(schad$schadSum, method = "MLE", type="GP", threshold=c)
  
  ## Lothar Schaden in promille der VS..
  schadRatio.ereignis <- schad.ereignis/vs.gvz
  ## wie haeufig trifft das groesste, beobachtete ereignis (zB lothar fuer sturm)
  ## den kanton zh?
  freq.ereignis <- 1/rtp
  
  ## wie oft ein sturmereignis welches 1 Promille der vs zerstoert?
  ## 1.
  ratio.hohePrio <- 1e-3
  hohe.schadprio <- vs.gvz*ratio.hohePrio ## 1 Promille der vs
  tiefe.schadprio <- schad.ereignis
  ratio <- hohe.schadprio/tiefe.schadprio
  freq.tiefePrio <- freq.ereignis
  
  ## 2.
  alpha <- a
  ## verhaeltnisse tiefer zu hoher schadenprioritaet
  x <- seq(1, 10, 0.01)
  y <- 1-(1/(x^alpha))
  xlab.string <- paste('hohe Schadenpriorität\n',
                       '--------------------------------\n',
                       'tiefe Schadenpriorität', sep="")
  title.string <- paste('Paretoverteilungsfunktion, alpha = ',
                        round(alpha,2), sep="")
  plot(x, y, typ='l', col=2, main=title.string,
       ylab=NA, xlab=NA)
  title(xlab = xlab.string, line = 4)
  
  ## aus tabelle mit alpha = a und verhaeltnis von hoher zu tiefer schadenprioritaet = ratio:
  ## w'keit, dass ein schaden, der die tiefe schadenprioritaet uebersteigt, die hohe
  ## schadenprioritaet hoechstens gerade noch erreicht...
  ind <- which(x==as.character(round(ratio, 2)))
  p <- y[ind]
  ## die komplementaere w'keit, dass ein schaden auch noch die hohe schadenprio uebersteigt:
  pk <- 1-p
  
  ## 3.
  ## freq der schaeden, welche die hohe prioritaet ueberschreitet ist...
  freq.hohePrio <- freq.tiefePrio * pk
  return(c(freq.hohePrio))
}


(schadFreq.hohePrio <- pareto.ES(schad, tresh, ng, RTP))
print(paste('Die Wiederkehrperiode eines Sturm-Ereignisses mit einer Schadensumme',
            ' von 1 Promille der Versicherungssumme (entspricht ',vs.gvz*1e-9,
            ' Mio CHF) ist ', round(1/schadFreq.hohePrio, 0), ' Jahre',sep=''))
  