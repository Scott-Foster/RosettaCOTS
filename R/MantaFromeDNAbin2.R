
predMantaFromeDNAbin <- function( av_PCR_bin, pred_type="PI", mantaBottomDist=200.3393, manta.nTows=4, intLevel=0.90){
  #av_PCR_bin contains the PCR binary records that are averaged of reps and then samples  
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
      xx["beta"] * log( av_PCR_bin + 0.01) + log( mantaBottomDist)  +
      #not including reef effect -- conditioned on 'average reef'
      # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
      rep( rnorm( length( av_PCR_bin), mean=0, sd=sqrt( xx["sigma.S2"])), each=manta.nTows) #SALAD_density, not the constructed df
    
    p <- 1-exp( -exp( lps)) #cloglog
    
    ys <- NULL
    if( type=="PI")
      ys <- stats::rbinom( n = length( p), prob = p, size = 1)
    if( type=="CI")
      ys <- p
    if( is.null( ys))
      stop("type must be PI or CI")
    
    ys <- tapply( ys, av_PCR_bin1$site, sum)
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_MantaFromeDNABin.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_mantaFromeDNABin.RDS")
  
  #rearranging for multiple mantaTows per site.
  av_PCR_bin1 <- rep( av_PCR_bin, each=manta.nTows)
  av_PCR_bin1 <- data.frame( av_PCR_bin=av_PCR_bin1, site=rep( 1:length( av_PCR_bin), each=manta.nTows))
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)#, newdatty=predDat)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  
  tmpDens <- tmp / manta.nTows
  
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( av_PCR_bin, rowMeans( tmp), rep( mantaBottomDist, nrow( tmp)), t( apply( tmp, 1, stats::quantile, prob=limmies)), 
                   rowMeans( tmpDens), t( apply( tmpDens, 1, stats::quantile, prob=limmies)))
  colnames( tmptmp) <- c( "av_PCR_bin","manta_nScars", "manta_Distance","nScars_lower","nScars_median","nScars_upper","manta_density","density_lower","density_median","density_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
}
