
predTransomFromTransect <- function( transect_density=newDensity, pred_type="PI", intLevel=0.95){
  #transect_density contains the diver transect data that is to be predicted from.
  #It is the #number of locations with corals (any type soft or hard) / #number of hard substrate
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)

  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * transect_density
    #not including reef effect -- conditioned on 'average reef'
    # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
    p <- exp( lps) / (1+exp( lps))
    ys <- NULL
    if( type=="PI")
      ys <- stats::rbeta( length( transect_density), p*xx["phi"], (1-p)*xx["phi"])
    if( type=="CI")
      ys <- p
    if( is.null( ys))
      stop("type must be PI or CI")
    
    return( ys)
  }

  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_TransomFromTransects.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_TransomFromTransects.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( transect_density, rowMeans( tmp), t( apply( tmp, 1, stats::quantile, prob=limmies)))
  
  colnames( tmptmp) <- c( "transect_density","TransomFromTransect_coral","coral_lower","coral_median","coral_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
}
