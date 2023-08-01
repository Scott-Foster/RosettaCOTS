
predTransectFromTransom <- function( transom_perc, numHard=400, pred_type="PI", intLevel=0.95){
  #transom_perc contains the diver transect data that is to be predicted from.
  #It is the ML predicted proportion of hard and soft corals
  #numHard is the number of points sampled on hard substrate.  Default is 400, the target for transect sampling
  
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
      ys <- stats::rbinom( n=length( transom_perc), size=numHard, prob=p)# / numHard
    if( type=="CI")
      ys <- p
    if( is.null( ys))
      stop("type must be PI or CI")
    
    return( ys)
  }
  
  #Loading the saved samples
  fm.post <- readRDS( system.file("extdata", "samples_TransectsFromTransom.RDS", package="RosettaCOTS"))
#  fm.post <- readRDS( samples_TransectsFromTransom.RDS")
  
  samps <- as.matrix( fm.post, iters=TRUE, chains=TRUE)
  tmp <- apply( samps, 1, pred.from.samp, type=pred_type)
  if( !inherits( tmp, "matrix"))
    tmp <- matrix( tmp, nrow=1)
  limmies <- c((1-intLevel)/2, 0.5, 1-(1-intLevel)/2)
  tmptmp <- cbind( transom_perc, rowMeans( tmp), t( apply( tmp, 1, stats::quantile, prob=limmies)))
  colnames( tmptmp) <- c( "transom_perc","TransectsFromTransom_coral","coral_lower","coral_median","coral_upper")
  if( pred_type=="PI")
    tmptmp[,-1] <- tmptmp[,-1] / numHard
  tmptmp <- as.data.frame( tmptmp)
  
  return( tmptmp)
}
