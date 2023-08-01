
COTS_calibrate <- function( from, to, from.metric, from.effort=NULL, to.effort=NULL, to.nTows=NULL, sites=NULL, intLevel=0.95){
  ####  elementary checking 
  #from.metric and from.effort
  if( is.null( from.metric))
    stop( "from.metric must be specified")
  if( is.null( from.effort) & !tolower(from) %in% tolower( c("Manta","eDNA_conc","eDNA_prop","coral_transom","coral_deep","coral_manta")))
    stop( "from.effort must be specified for your type of input data (specified in 'from' argument)")
  #sites
  if( is.null( sites)){
    warning("No sites argument supplied.  Assuming that all data comes from a single site (average from.metric and from.effort).")
    sites <- rep( 1, length( from.metric))
  }
  if( length( from.metric) != length( sites))
    stop( "Number of values in from.metric is different to the number of labels in sites. Please check.")
  #to.effort and to.nTows
  if( !is.null( to.effort) & length( to.effort) > 1)
    stop( "Multiple values of to.effort supplied. Please provide a single scalar.")
  if( tolower( to) %in% tolower( c("Manta","MANTA","manta")) & !is.null( to.nTows) & length( to.nTows)>1)
    stop( "Multiple values of to.nTows supplied. Please provide a single scalar.")
  #messages
  if( tolower( to) == "manta"){
    if( is.null( to.nTows))
      message( "Predictions are made using 4 tows (mode of the CCIP-D2 data).")
    if( is.null( to.effort))
      message( "Predictions are made using 200.5431 metre transects (mean of the CCIP-D2 data).")
  }
  if( tolower( to) == "cull" & is.null( to.effort))
    message( "Predictions are made using 473.6154 minutes of dive effort (mean of the CCIP-D2 data).")
  if( tolower(to) == "edna_prop" & is.null( to.effort))
    message( "Predictions are made using 12 samples per site (number taken at each site in the CCIP-D2 data).")
  if( tolower( to)=="salad" & is.null( to.effort))
    message( "Predictions are made using 106.9333 minutes search time OR 2048.571 metres OR 1.019333 hectares (means of CCIP-D2 data). Choice to match effort units of 'from.metric'.")
  if( tolower(to)=="coral_transect" & is.null( to.effort))
    message( "Predictions are made using 400 observations (over hard-substrate). This was the goal in the CCIP-D2 data.")
  ############
  if( tolower(from) == "cull"){
    if( is.null( from.effort))
      stop( "Please supply from.effort argument. It should represent the minutes diving for each dive (cull data).")
    if( length( from.metric) != length( from.effort))
      stop( "Number of values in from.metric is different to the number of values in from.effort. Please check.")
    sites <- as.factor( sites)
    #This loop and accumulating results is inefficient but I can't imagine that there is going to be millions of sites (where the ineffiency matters)
    res <- list()  
    for( ii in levels( sites))
      suppressWarnings( res[[ii]] <- predictSingleMeasurement( from=tolower(from), to=tolower(to), from.metric=from.metric[sites==ii], from.effort=from.effort[sites==ii], 
                                             to.effort=to.effort, to.nTows=to.nTows, intLevel=intLevel))
    message( attributes( res[[1]])$message)
    res <- do.call( "rbind",res)
  }
  ############
  if( tolower(from) %in% c("manta","edna_conc","edna_prop")){
    sites <- as.factor( sites)
    #This loop and accumulating results is inefficient but I can't imagine that there is going to be millions of sites (where the ineffiency matters)
    res <- list()  
    for( ii in levels( sites))
      suppressWarnings( res[[ii]] <- predictSingleMeasurement( from=tolower(from), to=tolower(to), from.metric=from.metric[sites==ii], from.effort=NULL,
                                             to.effort=to.effort, to.nTows=to.nTows, intLevel=intLevel))
    message( attributes( res[[1]])$message)     
    res <- do.call( "rbind",res)
  }
  ############
  if( tolower(from) == "salad"){
    if( is.null( from.effort))
      stop( "Please supply from.effort argument. For predicting Cull data, it should be SALAD dive time (minutes). For predicting Manta Scar and eDNA data, it should be SALAD dive distance (m). For eDNA, it should be SALAD search area (Ha).")
    if( length( from.metric) != length( from.effort))
      stop( "Number of values in from.metric is different to the number of values in from.effort.  Please check.")
    sites <- as.factor( sites)
    #This loop and accumulating results is inefficient but I can't imagine that there is going to be millions of sites (where the ineffiency matters)
    res <- list()  
    for( ii in levels( sites))
      suppressWarnings( res[[ii]] <- predictSingleMeasurement( from="SALAD", to=to, from.metric=from.metric[sites==ii], from.effort=from.effort[sites==ii], 
                                             to.effort=to.effort, to.nTows=to.nTows, intLevel=intLevel))
    message( attributes( res[[1]])$message)     
    res <- do.call( "rbind",res)
  }
  ############
  if( tolower(from) == "coral_transect"){
    if( is.null( from.effort))
      stop( "Please supply from.effort argument. For predicting from coral_transect data, it should be the number of hard substrate points in the transect.")
    if( length( from.metric) != length( from.effort))
      stop( "Number of values in from.metric is different to the number of values in from.effort.  Please check.")
    sites <- as.factor( sites)
    #This loop and accumulating results is inefficient but I can't imagine that there is going to be millions of sites (where the ineffiency matters)
    res <- list()  
    for( ii in levels( sites))
      suppressWarnings( res[[ii]] <- predictSingleMeasurement( from="coral_transect", to=to, from.metric=from.metric[sites==ii], from.effort=from.effort[sites==ii], 
                                             to.effort=to.effort, to.nTows=to.nTows, intLevel=intLevel))
    message( attributes( res[[1]])$message)     
    res <- do.call( "rbind",res)
  }
  ############
  if( tolower(from) %in% c("coral_transom","coral_deep","coral_manta")){
    sites <- as.factor( sites)
    #This loop and accumulating results is inefficient but I can't imagine that there is going to be millions of sites (where the ineffiency matters)
    res <- list()  
    for( ii in levels( sites))
      suppressWarnings( res[[ii]] <- predictSingleMeasurement( from=from, to=to, from.metric=from.metric[sites==ii], from.effort=from.effort[sites==ii], 
                                             to.effort=to.effort, to.nTows=to.nTows, intLevel=intLevel))
    message( attributes( res[[1]])$message)     
    res <- do.call( "rbind",res)
  }
  
  return( res)
}

# COTS_calibrate(from="SALAD", to="manta", from.metric=D2data$SALAD$No._Scars,
#               from.effort=D2data$SALAD$Distance, to.effort=NULL, to.nTows=NULL, sites=D2data$SALAD$Site, intLevel=0.95)

#predictSingleMeasurement <- function( from, to, from.metric=NULL, from.effort=NULL, from.nReps=NULL, to.effort=NULL, to.nTows=NULL, intLevel=0.95){
  
