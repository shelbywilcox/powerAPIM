###############################################################################################################
### P o w e r   S i m u l a t i o n   f o r   a   M e d i a t i o n   A P I M   u s i n g   R a w   D a t a ###
### Author: Thomas Ledermann                                                                                ###
### Date: October 2022                                                                                      ###
###############################################################################################################

sampleSize <- 500
alphaLevel <- .05		# significance level
nsim <- 10		# number of iterations
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
			Y1 ~ ca1*X1 + cp21*X2 + ba1*M1 + bp21*M2
			Y2 ~ ca2*X2 + cp12*X1 + ba2*M2 + bp12*M1
			M1 ~ aa1*X1 + ap21*X2
			M2 ~ aa2*X2 + ap12*X1

			X1 ~~ X2
			M1 ~~ M2
			Y1 ~~ Y2

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
	fit <- sem(mod, b.data, missing = "ML")
	ests <- parameterEstimates(fit)
	paa1 <- ests[ests$label == 'aa1', 'pvalue']
	paa2 <- ests[ests$label == 'aa2', 'pvalue']
	pba1 <- ests[ests$label == 'ba1', 'pvalue']
	pba2 <- ests[ests$label == 'ba2', 'pvalue']
	pca1 <- ests[ests$label == 'ca1', 'pvalue']
	pca2 <- ests[ests$label == 'ca2', 'pvalue']
	pap21 <- ests[ests$label == 'ap21', 'pvalue']
	pap12 <- ests[ests$label == 'ap12', 'pvalue']
	pbp21 <- ests[ests$label == 'bp21', 'pvalue']
	pbp12 <- ests[ests$label == 'bp12', 'pvalue']
	pcp21 <- ests[ests$label == 'cp21', 'pvalue']
	pcp12 <- ests[ests$label == 'cp12', 'pvalue']

	piea1 <- ests[ests$label == 'ie.a1', 'pvalue']
	piea2 <- ests[ests$label == 'ie.a2', 'pvalue']
	piea1p12 <- ests[ests$label == 'ie.a1p12', 'pvalue']
	piea2p21 <- ests[ests$label == 'ie.a2p21', 'pvalue']

	piep12a2 <- ests[ests$label == 'ie.p12a2', 'pvalue']
	piep21a1 <- ests[ests$label == 'ie.p21a1', 'pvalue']
	piep12p21 <- ests[ests$label == 'ie.p12p21', 'pvalue']
	piep21p12 <- ests[ests$label == 'ie.p21p12', 'pvalue']

	ptie11 <- ests[ests$label == 'tie11', 'pvalue']
	ptie22 <- ests[ests$label == 'tie22', 'pvalue']
	ptie12 <- ests[ests$label == 'tie12', 'pvalue']
	ptie21 <- ests[ests$label == 'tie21', 'pvalue']

	pt11 <- ests[ests$label == 't11', 'pvalue']
	pt22 <- ests[ests$label == 't22', 'pvalue']
	pt12 <- ests[ests$label == 't12', 'pvalue']
	pt21 <- ests[ests$label == 't21', 'pvalue']
	resMatrix <- cbind(paa1, paa2, pba1, pba2, pca1, pca2, pap21, pap12, pbp21, pbp12, pcp21, pcp12,
			   piea1, piea2, piea1p12, piea2p21, piep12a2, piep21a1, piep12p21, piep21p12,
			   ptie11, ptie22, ptie12, ptie21, pt11, pt22, pt12, pt21)
	resMatrix 
    })
	powerEst <- rowSums(results < alpha)/reps
	names(powerEst) <- c("aa1", "aa2", "ba1", "ba2", "ca1", "ca2", "ap21", "ap12", "bp21", "bp12", "cp21", "cp12",
			     "piea1", "piea2", "piea1p12", "piea2p21", "piep12a2", "piep21a1", "piep12p21", "piep21p12",
			     "ptie11", "ptie22", "ptie12", "ptie21", "pt11", "pt22", "pt12", "pt21")
	powerEst
}
power(raw.data = dat, alpha = alphaLevel, reps = nsim, sample.size = sampleSize)

# Generate example data
n <- 100
id <- 1:n
X1 <- round(rnorm(n, 5, 2), 0)
X2 <- abs(round(X1 + c(NA, NA, rnorm(n-2, 0, 3)), 0))
M1 <- abs(round(X1 + 0.5*X2 + c(rnorm(n-3, 0, 3), NA, NA, NA), 0))
M2 <- abs(round(X2 + 0.5*X1 + rnorm(n, 0, 3), 0))
Y1 <- abs(round(X1 + 0.5*M2 + c(rnorm(n-3, 0, 3), NA, NA, NA), 0))
Y2 <- abs(round(X2 + 0.5*M1 + rnorm(n, 0, 3), 0))

dat <- cbind.data.frame(X1, X2, M1, M2, Y1, Y2)
dat
round(cor(dat, use = "pairwise.complete.obs"), 2)

# Estimate the model
mod <-'	Y1 ~ ca1*X1 + cp21*X2 + ba1*M1 + bp21*M2
	Y2 ~ ca2*X2 + cp12*X1 + ba2*M2 + bp12*M1
	M1 ~ aa1*X1 + ap21*X2
	M2 ~ aa2*X2 + ap12*X1

	X1 ~~ X2
	M1 ~~ M2
	Y1 ~~ Y2

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
fit <- sem(mod, dat, missing = "ML")
summary(fit)
