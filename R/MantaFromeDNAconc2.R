
predMantaFromeDNA <- function( eDNA_PCR_conc, pred_type="PI", mantaBottomDist=200.3393, manta.nTows=4, intLevel=0.95){
  #eDNA_PCR_conc contains the PCR concentration data that is to be predicted from.  
  #Averaged over technical reps and then samples within a site
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
      xx["beta"] * log( eDNA_PCR_conc + 0.01) + log( mantaBottomDist) +
      #not including reef effect -- conditioned on 'average reef'
      # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
      rep( rnorm( length( eDNA_PCR_conc), mean=0, sd=sqrt( xx["sigma.S2"])), each=manta.nTows) #SALAD_density, not the constructed df
    
    p <- 1-exp( -exp( lps)) #cloglog
    
    ys <- NULL
    if( type=="PI")
      ys <- stats::rbinom( n = length( p), prob = p, size = 1)
    if( type=="CI")
      ys <- p
    if( is.null( ys))
      stop("type must be PI or CI")
    
    ys <- tapply( ys, eDNA_PCR_conc1$site, sum)
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_MantaFromeDNA.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_mantaFromeDNA.RDS")
  
  #rearranging for multiple mantaTows per site.
  eDNA_PCR_conc1 <- rep( eDNA_PCR_conc, each=manta.nTows)
  eDNA_PCR_conc1 <- data.frame( eDNA_PCR_conc=eDNA_PCR_conc1, site=rep( 1:length( eDNA_PCR_conc), each=manta.nTows))
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)#, newdatty=predDat)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  
  tmpDens <- tmp / manta.nTows
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( eDNA_PCR_conc, rowMeans( tmp), rep( mantaBottomDist, nrow( tmp)), t( apply( tmp, 1, stats::quantile, prob=limmies)), 
                   rowMeans( tmpDens), t( apply( tmpDens, 1, stats::quantile, prob=limmies)))
  colnames( tmptmp) <- c( "eDNA_PCR_conc","manta_nScars", "manta_Distance","nScars_lower","nScars_median","nScars_upper","manta_density","density_lower","density_median","density_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
}
