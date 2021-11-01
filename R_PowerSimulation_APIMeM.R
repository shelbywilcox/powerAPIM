#######################################################################################################################
### P o w e r   A n a l y s i s   f o r   t h e   M e d i a t i o n   A P I M   u s i n g   C o r r e l a t i o n s ###
### Author: Thomas Ledermann                                                                                        ###
### Created: November 22, 2020                                                                                      ###
### Last update: November 1, 2021                                                                                   ###
#######################################################################################################################

sampleSize <- 100
alphaLevel <- .05		# significance level
nsim <- 10000		# number of iterations
setSeed <- 123

rx1x2 <- .10	# correlation between X1 and X2
rx1m1 <- .30	# correlation between X1 and M1
rx1m2 <- .05	# correlation between X1 and M2
rx1y1 <- .10	# correlation between X1 and Y1
rx1y2 <- .05	# correlation between X1 and Y2
rx2m1 <- .05	# correlation between X2 and M1
rx2m2 <- .30	# correlation between X2 and M2
rx2y1 <- .05	# correlation between X2 and Y1
rx2y2 <- .10	# correlation between X2 and Y2
rm1m2 <- .20	# correlation between M1 and M2
rm1y1 <- .50	# correlation between M1 and Y1
rm1y2 <- .40	# correlation between M1 and Y2
rm2y1 <- .50	# correlation between M2 and Y1
rm2y2 <- .40	# correlation between M2 and Y2
ry1y2 <- .50	# correlation between Y1 and Y2

vX1 <- 1		# VAR(X1), can be set to 1
vX2 <- 1		# VAR(X2), can be set to 1
vM1 <- 1		# VAR(M1), can be set to 1
vM2 <- 1		# VAR(M2), can be set to 1
vY1 <- 1		# VAR(Y1), can be set to 1
vY2 <- 1		# VAR(Y2), can be set to 1


# install and load package
if(!require("lavaan")) install.packages("lavaan")
if(!require("paramtest")) install.packages("paramtest")
if(!require("simsem")) install.packages("simsem")
if(!require("dplyr")) install.packages("dplyr")

library(lavaan)
library(paramtest)
library(simsem)
library(dplyr)

# Covariance matrix
rMat <- matrix(NA, 6, 6)
rMat[lower.tri(rMat, diag = TRUE)] <- c(1, rx1x2, rx1m1, rx1m2, rx1y1, rx1y2, 1, rx2m1, rx2m2, rx2y1, rx2y2, 1, rm1m2, rm1y1, rm1y2, 1, rm2y1, rm2y2, 1, ry1y2, 1)
rMat[upper.tri(rMat)] <- t(rMat)[upper.tri(rMat)]
vNames <- c('x1', 'x2', 'm1', 'm2', 'y1', 'y2')                    # Variable names
dimnames(rMat) <- list(vNames, vNames)
sds <- sqrt(c(vX1, vX2, vM1, vM2, vY1, vY2))
covMat <- sds %*% t(sds) * rMat
covMat

### Distinguishable members
## Estimate the parameters using the covariance matrix
APIMeM <- '	
	y1 ~ ca1*x1 + cp21*x2 + ba1*m1 + bp21*m2
	y2 ~ ca2*x2 + cp12*x1 + ba2*m2 + bp12*m1

	m1 ~ aa1*x1 + ap21*x2
	m2 ~ aa2*x2 + ap12*x1

	y1 ~~ vy1*y1 + cy*y2
	y2 ~~ vy2*y2
	m1 ~~ vm1*m1 + cm*m2
	m2 ~~ vm2*m2
	x1 ~~ vx1*x1 + cx*x2
	x2 ~~ vx2*x2

	# k
	ka1 := ap21/aa1
	ka2 := ap12/aa2
	kb1 := bp21/ba1
	kb2 := bp12/ba2
	kc1 := cp21/ca1
	kc2 := cp12/ca2

	# IE
	ie.a1 := aa1*ba1
	ie.a2 := aa2*ba2
	ie.a1p12 := aa1*bp12
	ie.a2p21 := aa2*bp21
	ie.p12a2 := ap12*ba2
	ie.p21a1 := ap21*ba1
	ie.p12p21 := ap12*bp21
	ie.p21p12 := ap21*bp12

	# Total IE
	tie11 := aa1*ba1 + ap12*bp21
	tie22 := aa2*ba2 + ap21*bp12
	tie12 := aa1*bp12 + ap12*ba2
	tie21 := aa2*bp21 + ap21*ba1

	# Total
	t11 := aa1*ba1 + ap12*bp21 + ca1
	t22 := aa2*ba2 + ap21*bp12 + ca2
	t12 := aa1*bp12 + ap12*ba2 + cp12
	t21 := aa2*bp21 + ap21*ba1 + cp21
'
fit <- sem(APIMeM, sample.cov = covMat, sample.nobs = 100000)
summary(fit, standardized = TRUE, rsquare = TRUE)

ests <- parameterEstimates(fit)
ests

# extract results from the lavaan output
aA1 <- ests[ests$label == 'aa1', 'est']
aA2 <- ests[ests$label == 'aa2', 'est']
aP12 <- ests[ests$label == 'ap12', 'est']
aP21 <- ests[ests$label == 'ap21', 'est']
bA1 <- ests[ests$label == 'ba1', 'est']
bA2 <- ests[ests$label == 'ba2', 'est']
bP12 <- ests[ests$label == 'bp12', 'est']
bP21 <- ests[ests$label == 'bp21', 'est']
cA1 <- ests[ests$label == 'ca1', 'est']
cA2 <- ests[ests$label == 'ca2', 'est']
cP12 <- ests[ests$label == 'cp12', 'est']
cP21 <- ests[ests$label == 'cp21', 'est']
covX <- ests[ests$label == 'cx', 'est']
covM <- ests[ests$label == 'cm', 'est']
covY <- ests[ests$label == 'cy', 'est']
varX1 <- ests[ests$label == 'vx1', 'est']
varX2 <- ests[ests$label == 'vx2', 'est']
varM1 <- ests[ests$label == 'vm1', 'est']
varM2 <- ests[ests$label == 'vm2', 'est']
varY1 <- ests[ests$label == 'vy1', 'est']
varY2 <- ests[ests$label == 'vy2', 'est']

popEst <- cbind.data.frame(aA1, aA2, aP21, aP12, bA1, bA2, bP21, bP12, cA1, cA2, cP21, cP12, covX, covM, covY, varX1, varX2, varM1, varM2, varY1, varY2)
popEst

## MC Simulation
models <- popEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y1 ~ ', .$bA1, '*M1
		Y1 ~ ', .$bP21, '*M2
		Y1 ~ ', .$cA1, '*X1
		Y1 ~ ', .$cP21, '*X2
		Y2 ~ ', .$bP12, '*M1
		Y2 ~ ', .$bA2, '*M2
		Y2 ~ ', .$cP12, '*X1
		Y2 ~ ', .$cA2, '*X2
		M1 ~ ', .$aA1, '*X1
		M1 ~ ', .$aP21, '*X2
		M2 ~ ', .$aP12, '*X1
		M2 ~ ', .$aA2, '*X2
		X1 ~~ ', .$covX, '*X2
		M1 ~~ ', .$covM, '*M2
		Y1 ~~ ', .$covY, '*Y2
		X1 ~~ ', .$varX1, '*X1
		X2 ~~ ', .$varX2, '*X2
		M1 ~~ ', .$varM1, '*M1
		M2 ~~ ', .$varM2, '*M2
		Y1 ~~ ', .$varY1, '*Y1
		Y2 ~~ ', .$varY2, '*Y2
	')
	fitModel <-'	
		Y1 ~ cA1*X1 + cP21*X2 + bA1*M1 + bP21*M2
		Y2 ~ cA2*X2 + cP12*X1 + bA2*M2 + bP12*M1

		M1 ~ aA1*X1 + aP21*X2
		M2 ~ aA2*X2 + aP12*X1

		X1 ~~ varX1*X1 + covX*X2
		X2 ~~ varX2*X2
		M1 ~~ varM1*M1 + covM*M2
		M2 ~~ varM2*M2
		Y1 ~~ varY1*Y1 + covY*Y2
		Y2 ~~ varY2*Y2

		A1A1 := aA1*bA1
		A2A2 := aA2*bA2
		A1P12 := aA1*bP12
		A2P21 := aA2*bP21
		P12A2 := aP12*bA2
		P21A1 := aP21*bA1
		P12P21 := aP12*bP21
		P21P12 := aP21*bP12

		tIE11 := aA1*bA1 + aP12*bP21
		tIE22 := aA2*bA2 + aP21*bP12
		tIE12 := aA1*bP12 + aP12*bA2
		tIE21 := aA2*bP21 + aP21*bA1

		TE11 := aA1*bA1 + aP12*bP21 + cA1
		TE22 := aA2*bA2 + aP21*bP12 + cA2
		TE12 := aA1*bP12 + aP12*bA2 + cP12
		TE21 := aA2*bP21 + aP21*bA1 + cP21
	'
	data.frame(aA1 = .$aA1, aA2 = .$aA2, bA1 = .$bA1, bA2 = .$bA2, cA1 = .$cA1, cA2 = .$cA2, aP21 = .$aP21, aP12 = .$aP12, bP21 = .$bP21, bP12 = .$bP12, cP21 = .$cP21, cP12 = .$cP12, covX = .$covX, covM = .$covM, covY = .$covY,
		     gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
models

powerSim <- models %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(aA1 = .$aA1, aA2 = .$aA2, bA1 = .$bA1, bA2 = .$bA2, cA1 = .$cA1, cA2 = .$cA2, aP21 = .$aP21, aP12 = .$aP12, bP21 = .$bP21, bP12 = .$bP12, cP21 = .$cP21, cP12 = .$cP12, covX = .$covX, covM = .$covM, covY = .$covY,
	pSimEst = list(pSim))
})

# Point estimates of the population model
powerSim

# Power estimates and other statistics
round(summaryParam(powerSim$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)

# Distinguishable members - No names for additional statistics 
modelsNoNames <- popEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y1 ~ ', .$bA1, '*M1
		Y1 ~ ', .$bP21, '*M2
		Y1 ~ ', .$cA1, '*X1
		Y1 ~ ', .$cP21, '*X2
		Y2 ~ ', .$bP12, '*M1
		Y2 ~ ', .$bA2, '*M2
		Y2 ~ ', .$cP12, '*X1
		Y2 ~ ', .$cA2, '*X2
		M1 ~ ', .$aA1, '*X1
		M1 ~ ', .$aP21, '*X2
		M2 ~ ', .$aP12, '*X1
		M2 ~ ', .$aA2, '*X2
		X1 ~~ ', .$covX, '*X2
		M1 ~~ ', .$covM, '*M2
		Y1 ~~ ', .$covY, '*Y2
		X1 ~~ ', .$varX1, '*X1
		X2 ~~ ', .$varX2, '*X2
		M1 ~~ ', .$varM1, '*M1
		M2 ~~ ', .$varM2, '*M2
		Y1 ~~ ', .$varY1, '*Y1
		Y2 ~~ ', .$varY2, '*Y2
	')
	fitModel <-'
		Y1 ~ X1 + X2 + M1 + M2
		Y2 ~ X2 + X1 + M2 + M1

		M1 ~ X1 + X2
		M2 ~ X2 + X1

		X1 ~~ X2
		M1 ~~ M2
		Y1 ~~ Y2
	'
	data.frame(aA1 = .$aA1, aA2 = .$aA2, bA1 = .$bA1, bA2 = .$bA2, cA1 = .$cA1, cA2 = .$cA2, aP21 = .$aP21, aP12 = .$aP12, bP21 = .$bP21, bP12 = .$bP12, cP21 = .$cP21, cP12 = .$cP12, covX = .$covX, covM = .$covM, covY = .$covY,
		gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
modelsNoNames 

powerSimNoNames <- modelsNoNames %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(aA1 = .$aA1, aA2 = .$aA2, bA1 = .$bA1, bA2 = .$bA2, cA1 = .$cA1, cA2 = .$cA2, aP21 = .$aP21, aP12 = .$aP12, bP21 = .$bP21, bP12 = .$bP12, cP21 = .$cP21, cP12 = .$cP12, covX = .$covX, covM = .$covM, covY = .$covY,
	pSimEst = list(pSim))
})

# Power estimates and other statistics
round(summaryParam(powerSimNoNames$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)


## Find Sample Size
sampleSizes <- seq(25, 800, 25)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)
nRep <- 1000	# number of replications

SimNest <- models %>% do({
	SimN <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], lavaanfun = "sem")
	data_frame(aA1 = .$aA1, aA2 = .$aA2, bA1 = .$bA1, bA2 = .$bA2, cA1 = .$cA1, cA2 = .$cA2, aP21 = .$aP21, aP12 = .$aP12, bP21 = .$bP21, bP12 = .$bP12, cP21 = .$cP21, cP12 = .$cP12, covX = .$covX, covM = .$covM, covY = .$covY,
	SimEstN = list(SimN))
})

Nest <- getPower(SimNest$SimEstN[[1]])
sampleSizeEst <- findPower(Nest, iv = "N", power = 0.80)
sampleSizeEst


### Indistinguishable members

## Estimate the parameters using the covariance matrix
iAPIMeM <- '	
	y1 ~ ca1*x1 + cp21*x2 + ba1*m1 + bp21*m2
	y2 ~ ca2*x2 + cp12*x1 + ba2*m2 + bp12*m1

	m1 ~ aa1*x1 + ap21*x2
	m2 ~ aa2*x2 + ap12*x1

	y1 ~~ vy1*y1 + cy*y2
	y2 ~~ vy2*y2
	m1 ~~ vm1*m1 + cm*m2
	m2 ~~ vm2*m2
	x1 ~~ vx1*x1 + cx*x2
	x2 ~~ vx2*x2

	aa1 == aa2
	ba1 == ba2
	ca1 == ca2

	ap12 == ap21
	bp12 == bp21
	cp12 == cp21

	vx1 == vx2
	vm1 == vm2
	vy1 == vy2

	# k
	ka1 := ap21/aa1
	ka2 := ap12/aa2
	kb1 := bp21/ba1
	kb2 := bp12/ba2
	kc1 := cp21/ca1
	kc2 := cp12/ca2

	# IE
	ie.a1 := aa1*ba1
	ie.a2 := aa2*ba2
	ie.a1.p12 := aa1*bp12
	ie.a2.p21 := aa2*bp21
	ie.p12.a2 := ap12*ba2
	ie.p21.a1 := ap21*ba1
	ie.p12.p21 := ap12*bp21
	ie.p21.p12 := ap21*bp12

	# Total IE
	tie11 := aa1*ba1 + ap12*bp21
	tie22 := aa2*ba2 + ap21*bp12
	tie12 := aa1*bp12 + ap12*ba2
	tie21 := aa2*bp21 + ap21*ba1

	# Total
	t11 := aa1*ba1 + ap12*bp21 + ca1
	t22 := aa2*ba2 + ap21*bp12 + ca2
	t12 := aa1*bp12 + ap12*ba2 + cp12
	t21 := aa2*bp21 + ap21*ba1 + cp21
'

ifit <- sem(iAPIMeM, sample.cov = covMat, sample.nobs = 100000)
summary(ifit, standardized = TRUE, rsquare = TRUE)

iests <- parameterEstimates(ifit)
iests

# extract results from the lavaan output
aA1 <- iests[iests$label == 'aa1', 'est']
aA2 <- iests[iests$label == 'aa2', 'est']
aP12 <- iests[iests$label == 'ap12', 'est']
aP21 <- iests[iests$label == 'ap21', 'est']
bA1 <- iests[iests$label == 'ba1', 'est']
bA2 <- iests[iests$label == 'ba2', 'est']
bP12 <- iests[iests$label == 'bp12', 'est']
bP21 <- iests[iests$label == 'bp21', 'est']
cA1 <- iests[iests$label == 'ca1', 'est']
cA2 <- iests[iests$label == 'ca2', 'est']
cP12 <- iests[iests$label == 'cp12', 'est']
cP21 <- iests[iests$label == 'cp21', 'est']
covX <- iests[iests$label == 'cx', 'est']
covM <- iests[iests$label == 'cm', 'est']
covY <- iests[iests$label == 'cy', 'est']
varX1 <- iests[iests$label == 'vx1', 'est']
varX2 <- iests[iests$label == 'vx2', 'est']
varM1 <- iests[iests$label == 'vm1', 'est']
varM2 <- iests[iests$label == 'vm2', 'est']
varY1 <- iests[iests$label == 'vy1', 'est']
varY2 <- iests[iests$label == 'vy2', 'est']

ipopEst <- cbind.data.frame(aA1, aA2, aP21, aP12, bA1, bA2, bP21, bP12, cA1, cA2, cP21, cP12, covX, covM, covY, varX1, varX2, varM1, varM2, varY1, varY2)
ipopEst

imodels <- ipopEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y1 ~ ', .$bA1, '*M1
		Y1 ~ ', .$bP21, '*M2
		Y1 ~ ', .$cA1, '*X1
		Y1 ~ ', .$cP21, '*X2
		Y2 ~ ', .$bP12, '*M1
		Y2 ~ ', .$bA2, '*M2
		Y2 ~ ', .$cP12, '*X1
		Y2 ~ ', .$cA2, '*X2
		M1 ~ ', .$aA1, '*X1
		M1 ~ ', .$aP21, '*X2
		M2 ~ ', .$aP12, '*X1
		M2 ~ ', .$aA2, '*X2
		X1 ~~ ', .$covX, '*X2
		M1 ~~ ', .$covM, '*M2
		Y1 ~~ ', .$covY, '*Y2
		X1 ~~ ', .$varX1, '*X1
		X2 ~~ ', .$varX2, '*X2
		M1 ~~ ', .$varM1, '*M1
		M2 ~~ ', .$varM2, '*M2
		Y1 ~~ ', .$varY1, '*Y1
		Y2 ~~ ', .$varY2, '*Y2
	')
	fitModel <-'	
		Y1 ~ cA1*X1 + cP21*X2 + bA1*M1 + bP21*M2
		Y2 ~ cA2*X2 + cP12*X1 + bA2*M2 + bP12*M1

		M1 ~ aA1*X1 + aP21*X2
		M2 ~ aA2*X2 + aP12*X1

		X1 ~~ varX1*X1 + covX*X2
		X2 ~~ varX2*X2
		M1 ~~ varM1*M1 + covM*M2
		M2 ~~ varM2*M2
		Y1 ~~ varY1*Y1 + covY*Y2
		Y2 ~~ varY2*Y2

		aA1 == aA2
		bA1 == bA2
		cA1 == cA2
		aP12 == aP21
		bP12 == bP21
		cP12 == cP21
		varX1 == varX2
		varM1 == varM2
		varY1 == varY2

		A1A1 := aA1*bA1
		A2A2 := aA2*bA2
		A1P12 := aA1*bP12
		A2P21 := aA2*bP21
		P12A2 := aP12*bA2
		P21A1 := aP21*bA1
		P12P21 := aP12*bP21
		P21P12 := aP21*bP12

		tIE11 := aA1*bA1 + aP12*bP21
		tIE22 := aA2*bA2 + aP21*bP12
		tIE12 := aA1*bP12 + aP12*bA2
		tIE21 := aA2*bP21 + aP21*bA1

		TE11 := aA1*bA1 + aP12*bP21 + cA1
		TE22 := aA2*bA2 + aP21*bP12 + cA2
		TE12 := aA1*bP12 + aP12*bA2 + cP12
		TE21 := aA2*bP21 + aP21*bA1 + cP21
	'
	data.frame(aA1 = .$aA1, aA2 = .$aA2, bA1 = .$bA1, bA2 = .$bA2, cA1 = .$cA1, cA2 = .$cA2, aP21 = .$aP21, aP12 = .$aP12, bP21 = .$bP21, bP12 = .$bP12, cP21 = .$cP21, cP12 = .$cP12, covX = .$covX, covM = .$covM, covY = .$covY,
		     gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
imodels

ipowerSim <- imodels %>% do({
	ipSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(aA1 = .$aA1, aA2 = .$aA2, bA1 = .$bA1, bA2 = .$bA2, cA1 = .$cA1, cA2 = .$cA2, aP21 = .$aP21, aP12 = .$aP12, bP21 = .$bP21, bP12 = .$bP12, cP21 = .$cP21, cP12 = .$cP12, covX = .$covX, covM = .$covM, covY = .$covY,
	ipSimEst = list(ipSim))
})


# Point estimates of the population model
ipowerSim

# Power estimates and other statistics
round(summaryParam(ipowerSim$ipSimEst[[1]], detail = TRUE, alpha = 0.05), 3)


## Find Sample Size
sampleSizes <- seq(25, 800, 25)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)
nRep <- 1000	# number of replications

iSimNest <- imodels %>% do({
	iSimN <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], lavaanfun = "sem")
	data_frame(aA1 = .$aA1, aA2 = .$aA2, bA1 = .$bA1, bA2 = .$bA2, cA1 = .$cA1, cA2 = .$cA2, aP21 = .$aP21, aP12 = .$aP12, bP21 = .$bP21, bP12 = .$bP12, cP21 = .$cP21, cP12 = .$cP12, covX = .$covX, covM = .$covM, covY = .$covY,
#		     A1A1 = .$A1A1, A2A2 = .$A2A2, A1P12 = .$A1P12, A2P21 = .$A2P21, P12A2 = .$P12A2, P21A1 = .$P21A1, P12P21 = .$P12P21, P21P12 = .$P21P12,
#		     tIE11 = .$tIE11, tIE22 = .$tIE22, tIE12 = .$tIE12, tIE21 = .$tIE21, TE11 = .$TE11, TE22 = .$TE22, TE12 = .$TE12, TE21 = .$TE21, covX = .$covX, covM = .$covM, covY = .$covY,
	iSimEstN = list(iSimN))
})

iNest <- getPower(iSimNest$iSimEstN[[1]])
isampleSizeEst <- findPower(iNest, iv = "N", power = 0.80)
isampleSizeEst
