#############################################################################################################
### P o w e r   A n a l y s i s   f o r   a   S i m p l e   A P I M   u s i n g   C o r r e l a t i o n s ###
### Author: Thomas Ledermann                                                                              ###
### Date created: November 2020                                                                           ###
### Last update: April 5, 2021                                                                            ###
#############################################################################################################

sampleSize <- 100
alphaLevel <- .05		# significance level
nsim <- 10000		# number of iterations
setSeed <- 123

rx1x2 <- .10	# correlation between X1 and X2
rx1y1 <- .30	# correlation between X1 and Y1
rx1y2 <- .20	# correlation between X1 and Y2
rx2y1 <- .20	# correlation between X2 and Y1
rx2y2 <- .30	# correlation between X2 and Y2
ry1y2 <- .30	# correlation between Y1 and Y2

vX1 <- 1		# VAR(X1), can be set to 1
vX2 <- 1		# VAR(X2), can be set to 1
vY1 <- 1		# VAR(Y1), can be set to 1
vY2 <- 1		# VAR(Y2), can be set to 1


# install and load packages
if(!require("lavaan")) install.packages("lavaan")
if(!require("paramtest")) install.packages("paramtest")
if(!require("simsem")) install.packages("BinNor")
if(!require("dplyr")) install.packages("dplyr")

library(lavaan)
library(paramtest)
library(simsem)
library(dplyr)

# Covariance matrix
rMat <- matrix(NA, 4, 4)
rMat[lower.tri(rMat, diag = TRUE)] <- c(1, rx1x2, rx1y1, rx1y2, 1, rx2y1, rx2y2, 1, ry1y2, 1)
rMat[upper.tri(rMat)] <- t(rMat)[upper.tri(rMat)]
vNames <- c('x1', 'x2', 'y1', 'y2')                    # Variable names
dimnames(rMat) <- list(vNames, vNames)
sds <- sqrt(c(vX1, vX2, vY1, vY2))
covMat <- sds %*% t(sds) * rMat
covMat

#### Distinguishable members
## Estimate the parameters using the covariance matrix
APIM <- '	
	y1 ~ a1*x1 + p21*x2
	y2 ~ p12*x1 + a2*x2
	x1 ~~ vx1*x1 + cx*x2
	x2 ~~ vx2*x2
	y1 ~~ vy1*y1 + cy*y2
	y2 ~~ vy2*y2
	k1 := p21/a1
	k2 := p12/a2
'

fit <- sem(APIM, sample.cov = covMat, sample.nobs = 100000)
summary(fit, standardized = TRUE, rsquare = TRUE)

ests <- parameterEstimates(fit)
ests

# extract results from the lavaan output
AE1 <- ests[ests$label == 'a1', 'est']
AE2 <- ests[ests$label == 'a2', 'est']
PE12 <- ests[ests$label == 'p12', 'est']
PE21 <- ests[ests$label == 'p21', 'est']
covX <- ests[ests$label == 'cx', 'est']
covY <- ests[ests$label == 'cy', 'est']
varX1 <- ests[ests$label == 'vx1', 'est']
varX2 <- ests[ests$label == 'vx2', 'est']
varY1 <- ests[ests$label == 'vy1', 'est']
varY2 <- ests[ests$label == 'vy2', 'est']

popEst <- cbind.data.frame(AE1, AE2, PE21, PE12, covX, covY, varX1, varX2, varY1, varY2)
popEst

## Power Simulation
models <- popEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y1 ~ ', .$AE1, '*X1
		Y1 ~ ', .$PE21, '*X2
		Y2 ~ ', .$PE12, '*X1
		Y2 ~ ', .$AE2, '*X2
		X1 ~~ ', .$covX, '*X2
		Y1 ~~ ', .$covY, '*Y2
		X1 ~~ ', .$varX1, '*X1
		X2 ~~ ', .$varX2, '*X2
		Y1 ~~ ', .$varY1, '*Y1
		Y2 ~~ ', .$varY2, '*Y2
	')
	fitModel <-'
		Y1 ~ AE1*X1 + PE21*X2
		Y2 ~ PE12*X1 + AE2*X2
		X1 ~~ covX*X2
		Y1 ~~ covY*Y2
	'
	data.frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY,
		gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
models

powerSim <- models %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY,
	pSimEst = list(pSim))
})

# Point estimates of the population model
powerSim

round(summaryParam(powerSim$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)
CoverN <- getCoverage(powerSim$pSimEst[[1]], coverParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'))


pEst <- getPower(powerSim$pSimEst[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'))
pEst
pEst <- getPower(powerSim$pSimEst[[1]])
pEst


# Distinguishable members - No parameter names in the fitModel for additional statistics 
modelsNoNames <- popEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y1 ~ ', .$AE1, '*X1
		Y1 ~ ', .$PE21, '*X2
		Y2 ~ ', .$PE12, '*X1
		Y2 ~ ', .$AE2, '*X2
		X1 ~~ ', .$covX, '*X2
		Y1 ~~ ', .$covY, '*Y2
		X1 ~~ ', .$varX1, '*X1
		X2 ~~ ', .$varX2, '*X2
		Y1 ~~ ', .$varY1, '*Y1
		Y2 ~~ ', .$varY2, '*Y2
	')
	fitModel <-'
		Y1 ~ X1 + X2
		Y2 ~ X1 + X2
		X1 ~~ X2
		Y1 ~~ Y2
	'
	data.frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY,
		gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
modelsNoNames 

powerSimNoNames <- modelsNoNames %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY,
	pSimEst = list(pSim))
})

round(summaryParam(powerSimNoNames$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)
CoverN <- getCoverage(powerSimNoNames$pSimEst[[1]], coverParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'))
CoverN

# Power Estimates
powerEst <- getPower(powerSimNoNames$pSimEst[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'))
powerEst


## Find Sample Size
sampleSizes <- seq(25, 800, 25)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)
nRep <- 1000	# number of replications

SimNest <- models %>% do({
	SimN <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], 
	lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	SimEstN = list(SimN))
})
Nest <- getPower(SimNest$SimEstN[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'), nVal = lowerN:upperN)
sampleSizeEst <- getfindPower(Nest, iv = "N", power = 0.80)
sampleSizeEst


### Indistinguishable members

# Estimate the parameters using the covariance matrix
iAPIM <- '	
	y1 ~ a1*x1 + p21*x2
	y2 ~ p12*x1 + a2*x2
	x1 ~~ vx1*x1 + cx*x2
	x2 ~~ vx2*x2
	y1 ~~ vy1*y1 + cy*y2
	y2 ~~ vy2*y2

	a1 == a2
	p12 == p21
	vx1 == vx2
	vy1 == vy2
	
	k1 := p21/a1
	k2 := p12/a2
'

ifit <- sem(iAPIM, sample.cov = covMat, sample.nobs = 100000)
summary(ifit, standardized = TRUE, rsquare = TRUE)

iests <- parameterEstimates(ifit)
iests

# extract results from the lavaan output
AE1 <- iests[iests$label == 'a1', 'est']
AE2 <- iests[iests$label == 'a2', 'est']
PE12 <- iests[iests$label == 'p12', 'est']
PE21 <- iests[iests$label == 'p21', 'est']
covX <- iests[iests$label == 'cx', 'est']
covY <- iests[iests$label == 'cy', 'est']
varY1 <- iests[iests$label == 'vy1', 'est']
varY2 <- iests[iests$label == 'vy2', 'est']

ipopEst <- cbind.data.frame(AE1, AE2, PE21, PE12, covX, covY, varY1, varY2)
ipopEst

imodels <- ipopEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y1 ~ ', .$AE1, '*X1
		Y1 ~ ', .$PE21, '*X2
		Y2 ~ ', .$PE12, '*X1
		Y2 ~ ', .$AE2, '*X2
		X1 ~~ ', .$covX, '*X2
		Y1 ~~ ', .$covY, '*Y2
		Y1 ~~ ', .$varY1, '*Y1
		Y2 ~~ ', .$varY2, '*Y2
	')
	fitModel <-'
		Y1 ~ AE1*X1 + PE21*X2
		Y2 ~ PE12*X1 + AE2*X2
		X1 ~~ varX1*X1 + covX*X2
		X2 ~~ varX2*X2
		Y1 ~~ varY1*Y1 + covY*Y2
		Y2 ~~ varY2*Y2

		# equality constraints
		AE1 == AE2
		PE12 == PE21
		varX1 == varX2
		varY1 == varY2
	'
	data.frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2, gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
imodels

ipowerSim <- imodels %>% do({
	ipSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	ipSimEst = list(ipSim))
})

# Point estimates of the model
ipowerSim

# Power Estimates
ipowerEst <- getPower(ipowerSim$ipSimEst[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'))
ipowerEst

round(summaryParam(ipowerSim$ipSimEst[[1]], detail = TRUE, alpha = 0.05), 3)


## Find Sample Size
sampleSizes <- seq(25, 800, 25)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)
nRep <- 1000	# number of replications

iSimNest <- imodels %>% do({
	iSimN <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	iSimEstN = list(iSimN))
})

iNest <- getPower(iSimNest$iSimEstN[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'), nVal = lowerN:upperN)
isampleSizeEst <- findPower(iNest, iv = "N", power = 0.80)
isampleSizeEst


## Skewness: skewness = (0, 0, 3, 3)
distSk <- bindDist(skewness = c(0, 0, 3, 3), kurtosis = c(0, 0, 0, 0))
ipowerSimSk <- imodels %>% do({
	ipSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], indDist = distSk, lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	ipSimEst = list(ipSim))
})

# Point estimates of the model
ipowerSimSk

# Power Estimates
ipowerEstSk <- getPower(ipowerSimSk$ipSimEst[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'))
ipowerEstSk

round(summaryParam(ipowerSimSk$ipSimEst[[1]], detail = TRUE, alpha = 0.05), 3)

# Find sample size
sampleSizes <- seq(75, 3200, 25)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)

iSimNestSk <- imodels %>% do({
	iSimNSk <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], indDist = distsk, lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	iSimEstNSk = list(iSimNSk))
})
iNestSk <- getPower(iSimNestSk$iSimEstNSk[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'), nVal = lowerN:upperN)
isampleSizeEstSk <-findPower(iNestSk, iv = "N", power = 0.80)
isampleSizeEstSk


## Kurtosis: skewness = (0, 0, 0, 0), kurtosis = (0, 0, 3, 3)
distKu <- bindDist(skewness = c(0, 0, 0, 0), kurtosis = c(0, 0, 3, 3))
ipowerSimKu <- imodels %>% do({
	ipSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], indDist = distKu, lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	ipSimEst = list(ipSim))
})

# Point estimates of the model
ipowerSimKu

# Power Estimates
ipowerEstKu <- getPower(ipowerSimKu$ipSimEst[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'))
ipowerEstKu

round(summaryParam(ipowerSimNku$ipSimEst[[1]], detail = TRUE, alpha = 0.05), 3)

# Find sample size
sampleSizes <- seq(25, 800, 25)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)

iSimNestKu <- imodels %>% do({
	iSimNKu <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], indDist = distku, lavaanfun = "sem")
	data_frame(AE1 = .$AE1, AE2 = .$AE2, PE21 = .$PE21, PE12 = .$PE12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	iSimEstNKu = list(iSimNKu))
})

iNestKu <- getPower(iSimNestKu$iSimEstNKu[[1]], powerParam = c('AE1', 'AE2', 'PE21', 'PE12', 'covX', 'covY'), nVal = lowerN:upperN)
isampleSizeEstKu <- findPower(iNestKu, iv = "N", power = 0.80)
isampleSizeEstKu

