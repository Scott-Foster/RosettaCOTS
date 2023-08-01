
predSALADFromeDNAbin <- function( av_PCR_binary, pred_type="PI", SALADarea=10193.33 / 10000, intLevel=0.90){
  #av_PCR_binary contains the PCR average binary data that are to be predicted from.  
  #averaging is performed within samples (between reps) and then within sites (between samples)
  #Averaged over technical reps and then samples within a site
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)
  #SALADarea is the total area searched (over all divers) for the transect. Units are ha  If not specified the average
  #of the D2 fieldtrips is used.  At site level.  Either one value per site, or a single numeric (to
  #be replicated)
  
  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * log( av_PCR_binary + 0.01) + log( SALADarea) 
    #not including reef effect -- conditioned on 'average reef'
    # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
    mu <- exp( lps)
    p <- xx['r'] / ( xx['r'] + mu)
    ys <- NULL
    if( type=="PI")
      ys <- stats::rnbinom( length( mu), prob=p, size=xx['r'])#rpois( length( mu), lambda = mu)
    if( type=="CI")
      ys <- mu
    if( is.null( ys))
      stop("type must be PI or CI")
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_saladFromeDNAbinary.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_saladFromeDNAbinary.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)#, newdatty=predDat)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  tmpDens <- tmp / SALADarea
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( av_PCR_binary, rowMeans( tmp), rep( SALADarea, nrow( tmp)), t( apply( tmp, 1, stats::quantile, prob=limmies)), 
                   rowMeans( tmpDens), t( apply( tmpDens, 1, stats::quantile, prob=limmies)))
  
  colnames( tmptmp) <- c( "av_PCR_binary","SALAD_total", "SALAD_SweptArea","nCOTS_lower","nCOTS_median","nCOTS_upper","SALAD_density","density_lower","density_median","density_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
  
}

