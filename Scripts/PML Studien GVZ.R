## ----
## Risikoeinschaetzung anhand der PML Studien von 2014
## (Swiss Re, Munich Re, Aon Benfield)

## Mirco Heidemann, Oktober 2017
## ----

## elementarschaden daten, stuetzpunkte, aus der PML Studie
## (Presentation Ewa kozlowsky irv)

## Daten aus dem AON Benfield Modell (abgelesen)
dat <- data.frame(rtp=c(10, 20, 50, 100, 200, 300, 400,
                        500, 600, 700, 800, 900, 1000),
                  schaden=c(90, 190, 400, 700, 1200, 1600, 1700,
                            1750, 1800, 1900, 1950, 2000, 2100))

plot(schaden~rtp, data=dat, type='l', col='red', 
     xlab = "RTP", ylab = "Jahresschaden [Mio. CHF]",
     main='PML Studien 2014 (AON Benfield Modell) - Resultat GVZ')
abline(v=dat$rtp, col="lightgray", lty="dotted")
abline(h=dat$schaden, col="lightgray", lty="dotted")

rtp <- dat$rtp
zeit <- 1/rtp
schaden <- dat$schaden

zeitvalues <- seq(0, 0.1, 0.001)
## A) exponential model
exp.model <- lm(log(schaden)~ zeit)
predicted.exp <- exp(predict(exp.model,list(zeit=zeitvalues)))

## B) quadratic model
zeit2 <- zeit^2

quad.model <-lm(schaden ~ zeit + zeit2)
predicted.quad <- predict(quad.model,list(zeit=zeitvalues, zeit2=zeitvalues^2))

## C) Gamma GLM mit einem log link
gam.model <- glm(schaden ~ zeit, family=Gamma(link="log"))
predicted.gam <- exp(predict(gam.model,list(zeit=zeitvalues)))

## D) Log-normal

# ACHTUNG: fit <- glm(y ~ x, gaussian(link = "log"))
# does NOT fit a model with a lognormal response distribution.  It fits a
# non-linear regression model with an ordinary gaussian response
# distribution.  This model has constant variance, whereas the lognormal
# model (which you would fit by transforming the response) has constant
# coefficient of variation.  You would transform the response for two
# reasons, namely it should linearize the relationship between
# (transformed) response and predictors AND it should change a constant CV
# into homoscedasticity, or constant variance.  This latter property as
# important as the first, usually.  You should not think of a glm with log
# link as a kind of handy alternative to a log-transformed regression as
# they are in reality very different models.

lognorm.model <- lm(log(schaden)~log(zeit))
predicted.lognorm <- exp(predict(lognorm.model,list(zeit=zeitvalues)))

## Plotten der modelle
plot(zeit, schaden, pch=16, xlab = "1/RTP", ylab = "Jahresschaden [Mio. CHF]",
     cex.lab = 1, col = "black", main = 'PML Studien Elementar 2014 - GVZ')
lines(zeitvalues, predicted.exp, col = "#a6611a", lwd = 2,  lty = 1)
lines(zeitvalues, predicted.quad, col = "#dfc27d", lwd = 2, lty = 2)
lines(zeitvalues, predicted.gam, col = "#80cdc1", lwd = 2, lty = 3)
lines(zeitvalues, predicted.lognorm, col = "#018571", lwd = 2, lty = 4)
legend('topright', legend = c('Exponential', 'Quadratic', 'Gamma', 'Log-Normal'),
       lty = c(1:4), lwd = 2, col = c("#a6611a", "#dfc27d", "#80cdc1", "#018571"),
       bty = 'n')
