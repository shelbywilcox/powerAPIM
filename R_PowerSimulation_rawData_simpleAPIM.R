#########################################################################################################
### P o w e r   S i m u l a t i o n   f o r   a   S i m p l e   A P I M   u s i n g   R a w   D a t a ###
### Author: Thomas Ledermann                                                                          ###
### Date: January 2022                                                                                ###
#########################################################################################################

sampleSize <- 100
alphaLevel <- .05		# significance level
nsim <- 10000		# number of iterations
setSeed <- 123

# Install and load packages
if(!require("lavaan")) install.packages("lavaan")
library(lavaan)

power <- function(raw.data, alpha = 0.05, reps = 1000, sample.size) {
	id <- 1:nrow(raw.data)
	results  <- sapply(1:reps, function(x) {
		index <- sample(id, size = sample.size, replace = TRUE) 
		b.data <- raw.data[index, ]
		mod <-'
			Y1 ~ AE1*X1 + PE21*X2
			Y2 ~ PE12*X1 + AE2*X2
			X1 ~~ X2
			Y1 ~~ Y2
	'
	fit <- sem(mod, b.data, missing = "ML")
	ests <- parameterEstimates(fit)
	pAE1 <- ests[ests$label == 'AE1', 'pvalue']
	pAE2 <- ests[ests$label == 'AE2', 'pvalue']
	pPE21 <- ests[ests$label == 'PE21', 'pvalue']
	pPE12 <- ests[ests$label == 'PE12', 'pvalue']
	resMatrix <- cbind(pAE1, pAE2, pPE21, pPE12)
	resMatrix 
    })
	powerEst <- rowSums(results < alpha)/reps
	names(powerEst) <- c("AE1", "AE2", "PE21", "PE12")
	powerEst
}
power(raw.data = dat, alpha = alphaLevel, reps = nsim, sample.size = sampleSize)

# Generate example data
n <- 100
id <- 1:n
X1 <- round(rnorm(n, 5, 2), 0)
X2 <- abs(round(X1 + c(NA, NA, rnorm(n-2, 0, 3)), 0))
Y1 <- abs(round(X1 + 0.5*X2 + c(rnorm(n-3, 0, 3), NA, NA, NA), 0))
Y2 <- abs(round(X2 + 0.5*X1 + rnorm(n, 0, 3), 0))

dat <- cbind.data.frame(X1, X2, Y1, Y2)
dat
round(cor(dat, use = "pairwise.complete.obs"), 2)

# Estimate the model
mod <-'
			Y1 ~ AE1*X1 + PE21*X2
			Y2 ~ PE12*X1 + AE2*X2
			X1 ~~ X2
			Y1 ~~ Y2
	'
	fit <- sem(mod, dat, missing = "ML")
summary(fit)

ests <- parameterEstimates(fit)
ests
pAE1 <- ests[ests$label == 'AE1', 'pvalue']
pAE1
