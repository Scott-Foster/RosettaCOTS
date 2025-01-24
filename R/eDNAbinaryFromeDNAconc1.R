
predeDNAbinFromeDNAconc <- function( eDNA_PCR_conc, pred_type="PI", nreps=12, intLevel=0.95){

  #cull_density contains the cull data that is to be predicted from.
  #It is the #number of cots or #scars per minute
  #4 is the mode of the number of tows
  #814. is the average total tow length
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)
  #mantaBottomDist is the total trasnect distance (over all divers/segments) for the zone.  If not specified the average
  #of the D2 fieldtrips is used.  At site level.  Either one value per site, or a single numeric (to
  #be replicated)
  
  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * log( eDNA_PCR_conc + 0.01) 
    #not including reef effect -- conditioned on 'average reef'
    # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
    #not including the site effect -- conditioned on average site.
    p <- exp( lps) / (1+exp( lps))  #inverse logit
    
    ys <- NULL
    if( type=="PI")
      ys <- stats::rbinom( n = length( p), prob = p, size = nreps)
    if( type=="CI")
      ys <- p
    if( is.null( ys))
      stop("type must be PI or CI")
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_eDNAbinaryFromeDNAconc.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_eDNAbinaryFromeDNAconc.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)#, newdatty=predDat)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( eDNA_PCR_conc, rowMeans( tmp), rep( nreps, nrow( tmp)), t( apply( tmp, 1, stats::quantile, prob=limmies)))
  colnames( tmptmp) <- c( "av_eDNA_conc","Prob_eDNA", "nReps","Prob_lower", "Prob_median","Prob_upper")
  if( pred_type=="PI")
    tmptmp[,-1] <- tmptmp[,-1] / nreps
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
  
}
