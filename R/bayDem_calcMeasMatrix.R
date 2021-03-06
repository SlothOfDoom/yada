# Description
#   Calculate the measurement matrix, which is the likelihood of the radiocarbon
#   measurements in sampRcMeas calculated at the calendar dates in the vector
#   ygrid:
#
#  M = [p(phi_m,1|y_1)  p(phi_m,1:y_2) ... ]
#    = [p(phi_m,2|y_1)       ---       ... ]
#    = [     ...             ---       ... ]
#
#   The total uncertainty of the measurement comes from measurement error
#   (SIG_M, calculated using the measurement error for each measurement) and
#   the calibration curve error (SIK_k, calculated using the uncertainty for
#   the calibration curve at each grid point). These uncertainties (and the
#   associated measurements) should already be "projected" to 1950 equivalents
#   (e.g., in bayDem_dateSampToC14Samp).
#
# Example calls(s)
#
#   M <- bayDem_calcMeasMatrix(ygrid,sampRcMeas,calibDf)
# 
# Input(s)
#   Name          Type      Description
#   ygrid         vector    The locations at which to calculate the likelihood
#   sampRcMeas    list      Radiocarbon samples (see bayDem_dateSampToC14samp)
#   calibDf       dframe    Calibration curve (see bayDem_loadCalibCurve)
#
# Output(s)
#   Name          Type      Description
#   M             matrix    [Nmeas x Ngrid] The measurement matrix -- that is,
#                           the likelihood of this measurement calculated for
#                           the calendar dates in ygrid.

bayDem_calcMeasMatrix <- function(ygrid,sampRcMeas,calibDf) {
	# ygrid is in AD
	ygrid_BP <- 1950 - ygrid

	# calibration curve
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

	M <- exp(-(PHI_m - MU_k)^2 / (SIG_sq) / 2) / sqrt(SIG_sq) / sqrt(2*pi)
	return(M)
}
