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
#   Name    Type           Description
#   N       integer        The number of samples
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
#
# Output(s)
#   Name    Type           Description
#   samp    vector         The samples (length = N)

bayDem_sampleGaussMix <- function(N,th) {
    # Use distr to sample from a two-component, truncated Gaussian mixture
    normMix <- distr::UnivarMixingDistribution(dist::Norm(mean=th[3],sd=th[4]),distr::Norm(mean=th[5],sd=th[6]),mixCoeff=th[1:2])
    normMixTrunc <- distr::Truncate(normMix,th[7],th[8])
    samp <- distr::r(normMixTrunc)(N)
    return(samp)
}

