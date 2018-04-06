# Description
#   Sample from a (truncated) two-component Gaussian mixture. This provides
#   calendar dates of radiocarbon samples from the demographic model specified
#   by the two-component Gaussian mixture.
#
# Example calls(s)
#
#   samp <- bayDem_sampleGaussMix(N,th)
# 
# Input(s)
#   Name    Type      Description
#   ygrid   vector    The locations at which to calculate the likelihood
#                                        ymin to ymax]
#
# Output(s)
#   Name    Type           Description
#   samp    vector         The samples (length = N)

bayDem_calcObsMatrix <- function(ygrid,sampRcMeas,calibDf) {
	# ygrid and TH_mc are in AD
	ygrid_BP <- 1950 - ygrid
	# sampRcMeas is radiocarbon measurements referenced to BP
	y_curve     <- rev(calibDf$yearBP)
	mu_k_curve  <- exp(-rev(calibDf$uncalYearBP)/8033)
	sig_k_curve <- rev(calibDf$uncalYearBPError) * mu_k_curve / 8033

	# Interpolate curves at ygrid_BP to yield mu_k and sig_k
	mu_k  <- approx(y_curve,mu_k_curve,ygrid_BP)
	mu_k  <- mu_k$y
	sig_k <- approx(y_curve,sig_k_curve,ygrid_BP)
	sig_k <- sig_k$y

	PHI_m <- replicate(length(ygrid_BP),sampRcMeas$phi_m)
	SIG_m <- replicate(length(ygrid_BP),sampRcMeas$sig_m)

	MU_k  <- t(replicate(length(sampRcMeas$phi_m),mu_k))
	SIG_k <- t(replicate(length(sampRcMeas$sig_m),sig_k))

	SIG_sq <- SIG_m^2 + SIG_k^2

	obsMat <- exp(-(PHI_m - MU_k)^2 / (SIG_sq) / 2) / SIG_sq / sqrt(2*pi)
}

