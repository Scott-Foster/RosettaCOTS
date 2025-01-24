
predMantaCoralFromTransom <- function( transom_perc=newDensity, pred_type="PI", intLevel=0.95){
  #transom_perc contains the transomscan data that is to be predicted from.
  #It is the percentage of HARD corals
  
  #I'm hoping that the argument names are self-explanatory
  #pred_type is either "PI" (prediction interval) or "CI" (confidence interval)
  
  #function for predicting from each posterior sample.  Nested
  #as it only really applies to this model.
  pred.from.samp <- function(xx, type){
    lps <- xx["alpha"] + 
      xx["beta"] * transom_perc
    #not including reef effect -- conditioned on 'average reef'
    # + rnorm( nrow( SALAD_dat), mean=0, sd=sqrt( xx["sigma2"]))
    p <- exp( lps) / (1+exp( lps))
    ys <- NULL
    if( type=="PI")
      ys <- stats::rbeta( length( transom_perc), p*xx["phi"], (1-p)*xx["phi"])
    if( type=="CI")
      ys <- p
    if( is.null( ys))
      stop("type must be PI or CI")
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_MantaCoralFromTransom.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_MantaCoralFromTransom.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( transom_perc, rowMeans( tmp), t( apply( tmp, 1, stats::quantile, prob=limmies)))
  
  colnames( tmptmp) <- c( "transom_perc","MantaFromTransom_coral","coral_lower","coral_median","coral_upper")
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
}
