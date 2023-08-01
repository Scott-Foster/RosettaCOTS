
predCullFromeDNAoonc <- function( eDNA_PCR_conc, pred_type="PI", cullBottomTime=473.6154, intLevel=0.95){
  #eDNA_PCR_conc contains the PCR concentration data that is to be predicted from.  
  #Averaged over technical reps and then samples within a site
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)
  #cullBottomTime is the total dive time (over all divers) for the cull.  If not specified the average
  #of the D2 fieldtrips is used.  At site level.  Either one value per site, or a single numeric (to
  #be replicated)
  
  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * log( eDNA_PCR_conc + 0.01) + log( cullBottomTime) 
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
  fm.post <- readRDS( system.file("extdata", "samples_cullFromeDNA.RDS", package="RosettaCOTS"))
  #fm.post <- readRDS( samples_cullFromeDNA.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)#, newdatty=predDat)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  tmpDens <- tmp / cullBottomTime
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( eDNA_PCR_conc, rowMeans( tmp), rep( cullBottomTime, nrow( tmp)), t( apply( tmp, 1, stats::quantile, prob=limmies)), 
                   rowMeans( tmpDens), t( apply( tmpDens, 1, stats::quantile, prob=limmies)))
  
  colnames( tmptmp) <- c( "eDNA_PCR_conc","CULL_nCOTS", "CULL_BottomTime","nCOTS_lower","nCOTS_median","nCOTS_upper","CULL_density","density_lower","density_median","density_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
  
}

