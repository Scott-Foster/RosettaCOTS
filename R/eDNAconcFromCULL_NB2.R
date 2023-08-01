
predEDNAconcFromCULL <- function( cull_density, pred_type="PI", includeOutliers=FALSE, intLevel=0.95){
  #cull_density contains the cull data that is to be predicted from.
  #It is the #number of cots per minute dive
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)
  
  #includeOutliers is a boolean flag to indicate if the two outlier points should be included in the analysis.
  
  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * log( cull_density+0.01)
    mu <- exp( lps)
    #not including reef effect -- conditioned on 'average reef'
    #    mu <- exp( lps + xx["sigma2"]/2)
    p <- xx["r"] / (xx["r"] + mu)
    ys <- NULL
    if( type=="PI")
      ys <- stats::rnbinom( n=length( mu), size=xx["r"], prob=p)#rpois( n=length( mu), lambda=mu)#rlnorm( n=length( lps), meanlog=lps, sdlog=sqrt(xx["sigma2"])))#stats::rbinom( n = length( p), prob = p, size = newData$nTows)
    if( type=="CI")
      ys <- exp( lps)#+0.5*xx["sigma2"])
    if( is.null( ys))
      stop("type must be PI or CI")
    
    #average them up to the site level
    #all balanced so averages of averages is just the average of all.
    preds <- ys #preds <- tapply( ys, cull_density, mean)
    # preds <- cbind( as.numeric( names( preds)), preds)
    return( preds)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_eDNAconc_cleanFromCULL.RDS", package="RosettaCOTS"))
#  if( includeOutliers)
#    fm.post <- readRDS( samples_eDNAconcFromCULL.RDS")
#  else
#    fm.post <- readRDS( samples_eDNAconc_cleanFromCULL.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( cull_density, rowMeans( tmp), t( apply( tmp, 1, stats::quantile, prob=limmies)))
  
  colnames( tmptmp) <- c( "cull_density","eDNA_conc", "eDNA_conc_lower","eDNA_conc_median","eDNA_conc_upper")
  tmptmp <- as.data.frame( tmptmp)
  rownames( tmptmp) <- NULL
  
  return( tmptmp)
  
}
