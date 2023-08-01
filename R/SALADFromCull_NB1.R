
predSALADFromCULL <- function( cull_density, pred_type="PI", SALADBottomTime=106.9333, intLevel=0.95){
  #cull_density contains the cull data that is to be predicted from.
  #It is the #number of cots per minute
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)
  #SALADBottomTime is the total dive time (over all divers) for the cull.  If not specified the average
  #of the D2 fieldtrips is used.  At site level.  Either one value per site, or a single numeric (to
  #be replicated)
  
  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * log( cull_density+0.01) + log( SALADBottomTime) 
    #not including reef effect -- conditioned on 'average reef'
    # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
    mu <- exp( lps)
    p <- xx["r"] / (xx["r"] + mu)
    ys <- NULL
    if( type=="PI")
      ys <- stats::rnbinom( length(mu), size = xx["r"], prob=p)#rpois( length( mu), lambda = mu)
    if( type=="CI")
      ys <- mu
    if( is.null( ys))
      stop("type must be PI or CI")
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_SALADFromCull.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_SALADFromCull.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)#, newdatty=predDat)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  tmpDens <- tmp / SALADBottomTime
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( cull_density, rowMeans( tmp), rep( SALADBottomTime, nrow( tmp)), t( apply( tmp, 1, stats::quantile, prob=limmies)), 
                   rowMeans( tmpDens), t( apply( tmpDens, 1, stats::quantile, prob=limmies)))
  
  colnames( tmptmp) <- c( "cull_density","SALAD_nCOTS", "SALAD_BottomTime","nCOTS_lower","nCOTS_median","nCOTS_upper","SALAD_density","density_lower","density_median","density_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
  
}
