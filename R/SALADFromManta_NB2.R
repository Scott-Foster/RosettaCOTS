
predSALADFromManta <- function( manta_scars, manta.nTows=4, pred_type="PI", SALADarea=10193.33, intLevel=0.95){
  #manta_scars contains the manta data that is to be predicted from (number of tows with scars seen).
  #manta.nTows is the number of tows in zone.
  SALADdist <- SALADarea
  #SALADarea is the total area searched for the SALAD dive. If not specified the average
  #of the D2 fieldtrips is used.  At site level.  Either one value per site, or a single numeric (to
  #be replicated)
  #SALADdist is the total dive distance (over all divers) for the dive. Naming for historical reasons.
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)
  
  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * log( ( manta_scars / manta.nTows) + 0.01) + log( SALADdist) 
    #not including reef effect -- conditioned on 'average reef'
    # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
    mu <- exp( lps)
    p <- xx['r'] / (xx['r'] + mu)
    ys <- NULL
    if( type=="PI")
      ys <- stats::rnbinom( length(mu), prob=p, size=xx['r'])#pois( length( mu), lambda = mu)
    if( type=="CI")
      ys <- mu
    if( is.null( ys))
      stop("type must be PI or CI")
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_SALADFromManta.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_SALADFromManta.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)#, newdatty=predDat)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  tmpDens <- tmp / SALADdist
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( (manta_scars/manta.nTows), rowMeans( tmp), rep( SALADdist, nrow( tmp)), t( apply( tmp, 1, stats::quantile, prob=limmies)), 
                   rowMeans( tmpDens), t( apply( tmpDens, 1, stats::quantile, prob=c(0.025,0.5,0.975))))
  
  colnames( tmptmp) <- c( "scaled_manta_scars","SALAD_Tot", "SALAD_dist","nCOTS_lower","nCOTS_median","nCOTS_upper","SALAD_density","density_lower","density_median","density_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
  
}
