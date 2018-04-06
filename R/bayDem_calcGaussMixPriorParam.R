# Description
#   Calculate the parameterization for the prior probability, p(th|alpha), for
#   the (truncated) Gaussian mixture model.
#
#   (1) The mean (mu) is drawn uniformly from ymin to ymax
#
#   (2) The scale (sigma) is drawn from an offset gamma distribution:
#
#       sigma ~ startOffset + gamma(shape=alpha,rate=beta)
#
#       For a gamma distribution, the following holds for the mode (mode0 is
#       b/c momentarily ignoring offset):
#
#       (a) mode0 = (alpha-1)/beta --> beta = (alpha-1)/mode0
#
#       The input location of the mode, modeOffset, is relative to the origin, so
#
#       mode0 = modeOffset - startOffset
#
#       Hence:
#
#       beta = (alpha-1)/(modeOffset - startOffset)
#
#   (3) The two-component mixture is drawn from the dirichlet distribution
#       with parameter dirichParam
#
# Example calls(s)
#
#   gmPriorParam <- bayDem_calcGaussMixPriorParam(ymin,ymax,startOffset,modeOffset,gammaAlpha,dirichParam)
# 
# Input(s)
#   Name           Type      Description
#   ymin           scalar    Minimum date
#   ymax           scalar    Maximum date
#   startOffset    scalar    Constant added to gamma draw for the scale (sigma)
#   modeOffset     scalar    Location of the mode for gamma draw for the scale
#                            [see above]
#   gammaAlpha     scalar    alpha for the gamma draw for the scale
#   dirichParam    vector    dirichlet parameter vector for the mixture
#                            proportions
#
# Output(s)
#   Name            Type     Description
#   gmPriorParam    list     A list that parametrizes the model parameter
#                            theta, with the following components:
#                            ymin -- Mininum calendar date
#                            ymax -- Maximum calendar date
#                            startOffset -- constant added to gamma draw
#                            modeOffset  -- location of mode for gamma draw
#                            gammaAlpha  -- alpha parameter for gamma draw
#                            gammaBeta   -- beta parameter for gamma draw
#                            dirichParam -- dirichlet parameters for z

bayDem_calcGaussMixPriorParam <- function(ymin,ymax,startOffset,modeOffset,gammaAlpha,dirichParam) {
    dy <- ymax - ymin
    gmPriorParam <- list(ymin=ymin,ymax=ymax,startOffset=startOffset,modeOffset=modeOffset,gammaAlpha=gammaAlpha,dirichParam=dirichParam)
    gmPriorParam$gammaBeta  <- (gammaAlpha-1)/(modeOffset-startOffset)
    return(gmPriorParam)
}

