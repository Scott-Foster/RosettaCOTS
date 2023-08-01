
predCullFromSalad <- function( SALAD_density, pred_type="PI", cullBottomTime=473.6154, intLevel=0.95){
  #SALAD_density contains the SALAD data that is to be predicted from.
  #It is the #number of cots or #scars per minute
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)
  #cullBottomTime is the total dive time (over all divers) for the cull.  If not specified the average
  #of the D2 fieldtrips is used.  At site level.  Either one value per site, or a single numeric (to
  #be replicated)
  
  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * log( SALAD_density) + log( cullBottomTime)
    #not including reef effect -- conditioned on 'average reef'
    # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
    mu <- exp( lps)
    p <- xx["r"] / ( mu + xx["r"])
    ys <- NULL
    if( type=="PI")
      ys <- stats::rnbinom( length( mu), prob=p, size=xx["r"])#prpois( length( mu), lambda = mu)
    if( type=="CI")
      ys <- mu
    if( is.null( ys))
      stop("type must be PI or CI")
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_cullFromSALAD.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_cullFromSALAD.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)#, newdatty=predDat)
  if( !inherits( tmp, "array"))
    tmp <- matrix( tmp, nrow=1)
  tmpDens <- tmp / cullBottomTime  #done twice (see offset in prediction), but who cares...
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( SALAD_density, rowMeans( tmp), rep( cullBottomTime, nrow( tmp)), t( apply( tmp, 1, stats::quantile, prob=limmies)), 
                   rowMeans( tmpDens), t( apply( tmpDens, 1, stats::quantile, prob=limmies)))
  
  colnames( tmptmp) <- c( "SALAD_density","CULL_nCOTS", "CULL_BottomTime","nCOTS_lower","nCOTS_median","nCOTS_upper","CULL_density","density_lower","density_median","density_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
  
}
