# Description
#   Calculate the probability density function (pdf) for a (truncated)
#   two-component Gaussian mixture.
#
# Example calls(s)
#
#   bayDem_calcGaussMixPdf <- function(th,y) {
# 
# Input(s)
#   Name    Type           Description
#   th      vector-like    A vectore-like object with the following entries:
#                          (1) z1    -- Weight of the first mixture
#                          (2) z2    -- Weight of the second mixture
#                          (3) mu1   -- Mean of the first mixture
#                          (4) sig1  -- Standard deviation of the first mixture
#                          (5) mu2   -- Mean of the second mixture
#                          (6) sig2  -- Standard deviation of the second mixture
#                          (7) ymin  -- Minimum value
#                          (8) ymax  -- Maximum value
#                                       [Samples are truncated to the interval
#                                        ymin to ymax]
#   y       vector         Vector at which to calculate pdf (calendar dates)
#
# Output(s)
#   Name    Type           Description
#   f       vector         Output pdf


bayDem_calcGaussMixPdf <- function(th,y) {
	numSamp <- dim(th)[2]
	f <- rep(0,length(y))

	# Pre-calculate the normalization vector that accounts for truncation
	normVect1 <- pnorm(th['ymax',],th['mu1',],th['sig1',]) - pnorm(th['ymin',],th['mu1',],th['sig1',])
	normVect2 <- pnorm(th['ymax',],th['mu2',],th['sig2',]) - pnorm(th['ymin',],th['mu2',],th['sig2',])

	for(ii in 1:length(f)) {
		f1 <- rep(0,numSamp)
		f2 <- rep(0,numSamp)
		# B is true for non-truncated samples
		B <- th['ymin',] <= y[ii] & y[ii] <= th['ymax',]
		f1[B] <- dnorm(y[ii],th['mu1',B],th['sig1',B]) / normVect1[B]
		f2[B] <- dnorm(y[ii],th['mu2',B],th['sig2',B]) / normVect2[B]
		f[ii] <- mean(th['z1',] * f1 + th['z2',] * f2)
	}
	return(f)
}
