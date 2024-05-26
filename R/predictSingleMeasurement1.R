
predictSingleMeasurement <- function( from, to, from.metric=NULL, from.effort=NULL, to.effort=NULL, to.nTows=NULL, intLevel=0.95){
  #from is text (e.g. "CULL","SALAD", "Manta", "Manta_coral")
  #to is text, see from.
  #from.metric is a scalar or vector.  It's precise value will depend on what from is.
  #from.effort is a scalar or vector. 
  #to.effort is a scalar -- effort for a site.
  #to.nTows is a scalar -- number of manta tows at a site.
  
  #################
  if( tolower(from)=="salad"){
    #assume mutliple dive records per site.  Sum #COTS over those dives (and dive time)
    if( is.null( from.metric))
      stop("SALAD nCOTS not supplied. Please do so (in argument from.metric)")
    from.metric <- sum( from.metric)
    if( is.null( from.effort))
      stop("SALAD minutes searched not supplied. Please do so (in argument from.effort)")
    from.effort <- sum( from.effort)
    ####
    if( tolower(to)=="cull"){
      if( is.null( to.effort)){
        warning( "Cull search time not supplied.  Assuming 473.6154 minutes (the average of the CCIP-D2 project data)")
        to.effort <- 473.6154 #mean of D2 data per site
      }
      predVal.PI <- predCullFromSalad( SALAD_density=from.metric / from.effort, pred_type="PI", cullBottomTime=to.effort, intLevel=intLevel)
      predVal.CI <- predCullFromSalad( SALAD_density=from.metric / from.effort, pred_type="CI", cullBottomTime=to.effort, intLevel=intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, predVal.CI[,c("CULL_density","density_lower","density_upper")], predVal.PI[,c("density_lower","density_upper")]))
      names( predVal) <- c("SALAD_COTS","SALAD_EFFORT", "SALAD_density", "CULL_BottomTime", "mean_CULL", "lowerCI_CULL", "upperCI_CULL","lowerPI_CULL","upperPI_CULL")
      attr( predVal, "message") <- "Predicting Cull from SALAD. SALAD measured in nCOTS per minute. Cull also measured in nCOTS per minute"
    }
    ####
    if( tolower(to)=="manta"){
      if( is.null( to.nTows)){
        warning( "Manta nTows not supplied.  Assuming to.nTows=4 (the mode of the CCIP-D2 project data)")
        to.nTows <- 4
      }
      if( is.null( to.effort)){
        warning( "Manta tow distance not supplied (as to.effort argument). Assuming 200.5431 metres (the average of the CCIP-D2 project data)")
        to.effort <- 200.5431
      }
      predVal.PI <- predMantaFromSalad( SALAD_density=from.metric / from.effort, pred_type="PI", mantaBottomDist=to.effort, manta.nTows=to.nTows, intLevel=intLevel)
      predVal.CI <- predMantaFromSalad( SALAD_density=from.metric / from.effort, pred_type="CI", mantaBottomDist=to.effort, manta.nTows=to.nTows, intLevel=intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, to.nTows, predVal.CI[, c( "manta_density","density_lower","density_upper")], predVal.PI[, c( "density_lower","density_upper")]))
      names( predVal) <- c("SALAD_COTS","SALAD_EFFORT","SALAD_density","manta.avTowLength","manta.nTows","Prob_Manta","lowerCI_Manta","upperCI_Manta","lowerPI_Manta","upperPI_Manta")
      attr( predVal,"message") <- "Predicting Manta from SALAD. SALAD measured in nCOTS per metre. Prediction is probability of scar presence (given tow distance)"
    }
    ####
    if( tolower(to)=="salad"){
      tmp <- from.metric / from.effort
      attr( tmp, "message") <- "You are trying to predict SALAD from SALAD.  Please reconsider:-)  Returning input data."
      return( tmp)
    }
    ####
    if( tolower(to)=="edna_conc"){
      predVal.PI <- predEDNAconcFromSalad( SALAD_density=from.metric / from.effort, pred_type="PI", intLevel=intLevel)
      predVal.CI <- predEDNAconcFromSalad( SALAD_density=from.metric / from.effort, pred_type="CI", intLevel=intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, predVal.CI[, c( "eDNA_conc","eDNA_conc_lower","eDNA_conc_upper")], predVal.PI[, c( "eDNA_conc_lower","eDNA_conc_upper")]))
      names( predVal) <- c("SALAD_COTS","SALAD_EFFORT","SALAD_density","eDNA_conc","lowerCI_eDNA_conc","upperCI_eDNA_conc","lowerPI_eDNA_conc","upperPI_eDNA_conc")
      attr( predVal, "message") <- "Predicting eDNA concentration from SALAD. SALAD measured in nCOTS per hectare. Prediction is concentration per sample (sum of technical reps)."
    }
    ####
    if( tolower(to)=="edna_prop"){
      if( is.null( to.effort)){
        warning( "nSamples not supplied (as to.effort argument).  Assuming 12 (the number taken during the CCIP-D2 project data)")
        to.effort <- 12
      }
      predVal.PI <- predeDNAbinFromSALAD( SALAD_density = from.metric / from.effort, nreps=to.effort, pred_type="PI")
      predVal.CI <- predeDNAbinFromSALAD( SALAD_density = from.metric / from.effort, nreps=to.effort, pred_type="CI")
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, predVal.CI[,c("Prob_eDNA","Prob_lower","Prob_upper")], predVal.PI[,c("Prob_lower","Prob_upper")]))
      names( predVal) <- c("SALAD_COTS","SALAD_EFFORT","SALAD_density", "nSamples", "eDNA_prop","lowerCI_eDNA_prop","upperCI_eDNA_prop","lowerPI_eDNA_prop","upperPI_eDNA_prop")
      attr( predVal, "message") <- "Predicting eDNA proportion positive from SALAD. SALAD measured in nCOTS per hectare. Prediction is proportion of samples at a site that are positive."
    }
  }
  #################
  if( from=="cull"){
    if( is.null( from.metric))
      stop( "CULL nCOTS not supplied.  Please do so (in argument from.metric)")
    from.metric <- sum( from.metric)
    if( is.null( from.effort))
      stop( "CULL mintues searched not supplied.  Please do so (in argument from.effort)")
    from.effort <- sum( from.effort)
    ####
    if( tolower(to)=="cull"){
      tmp <- from.metric / from.effort
      attr( tmp, "message") <- "You are trying to predict CULL from CULL.  Please reconsider:-)  Returning input data."
      return( tmp)
    }
    ####
    if( tolower(to)=="manta"){
      if( is.null( to.nTows)){
        warning( "Manta nTows not supplied.  Assuming to.nTows=4 (the mode of the CCIP-D2 project data)")
        to.nTows <- 4
      }
      if( is.null( to.effort)){
        warning( "Manta tow distance not supplied (as to.effort argument). Assuming 200.5431 metres (the average of CCIP-D2 project data)")
        to.effort <- 200.5431
      }
      predVal.PI <- predMantaFromCull( cull_density=from.metric / from.effort, pred_type="PI", mantaBottomDist=to.effort, manta.nTows=to.nTows, intLevel=intLevel)
      predVal.CI <- predMantaFromCull( cull_density=from.metric / from.effort, pred_type="CI", mantaBottomDist=to.effort, manta.nTows=to.nTows, intLevel=intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, to.nTows, predVal.CI[, c( "manta_density","density_lower","density_upper")], predVal.PI[, c( "density_lower","density_upper")]))
      names( predVal) <- c("CULL_COTS","CULL_EFFORT","CULL_density","manta.avTowLength","manta.nTows","Prob_Manta","lowerCI_Manta","upperCI_Manta","lowerPI_Manta","upperPI_Manta")
      attr( predVal, "message") <- "Predicting Manta from CULL. CULL measured in nCOTS per minute. Prediction is probability of scar presence (given tow distance)"
    }
    ####
    if( tolower(to)=="salad"){
      if( is.null( to.effort)){
        warning( "SALAD swept area not supplied.  Assuming 10193.33 (the average of the CCIP-D2 project data")
        to.effort <- 10193.33
      }
      predVal.PI <- predSALADFromCULL( cull_density = from.metric / from.effort, SALADBottomTime = to.effort, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predSALADFromCULL( cull_density = from.metric / from.effort, SALADBottomTime = to.effort, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, predVal.CI[, c( "SALAD_density","density_lower","density_upper")], predVal.PI[, c( "density_lower","density_upper")]))
      names( predVal) <- c("CULL_COTS","CULL_EFFORT","CULL_density","SALAD_EFFORT","SALAD_density","lowerCI_SALAD","upperCI_SALAD","lowerPI_SALAD","upperPI_SALAD")
      attr( predVal, "message") <- "Predicting SALAD from CULL. CULL measured in nCots per minute. SALAD also measured in nCOTS per minute."
    }
    ####
    if( tolower(to)=="edna_conc"){
      predVal.PI <- predEDNAconcFromCULL( cull_density=from.metric / from.effort, pred_type="PI", intLevel=intLevel)
      predVal.CI <- predEDNAconcFromCULL( cull_density=from.metric / from.effort, pred_type="CI", intLevel=intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, predVal.CI[, c( "eDNA_conc","eDNA_conc_lower","eDNA_conc_upper")], predVal.PI[, c( "eDNA_conc_lower","eDNA_conc_upper")]))
      names( predVal) <- c("CULL_COTS","CULL_EFFORT","CULL_density","eDNA_conc","lowerCI_eDNA_conc","upperCI_eDNA_conc","lowerPI_eDNA_conc","upperPI_eDNA_conc")
      attr( predVal, "message") <- "Predicting eDNA concentration from CULL. CULL measured in nCOTS per minute. Prediction is concentration per sample (sum of technical reps)."
    }
    ####
    if( tolower(to)=="edna_prop"){
      if( is.null( to.effort)){
        warning( "nSamples not supplied (as to.effort argument).  Assuming 12 (the number taken during the CCIP-D2 project data)")
        to.effort <- 12
      }
      predVal.PI <- predeDNAbinFromCULL( CULL_density=from.metric / from.effort, nreps=to.effort, pred_type="PI", intLevel = intLevel)
      predVal.CI <- predeDNAbinFromCULL( CULL_density=from.metric / from.effort, nreps=to.effort, pred_type="CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, predVal.CI[,c("Prob_eDNA","Prob_lower","Prob_upper")], predVal.PI[,c("Prob_lower","Prob_upper")]))
      names( predVal) <- c("CULL_COTS","CULL_EFFORT","CULL_density", "nSamples", "eDNA_prop","lowerCI_eDNA_prop","upperCI_eDNA_prop","lowerPI_eDNA_prop","upperPI_eDNA_prop")
      attr( predVal, "message") <- "Predicting eDNA proportion positive from CULL. CULL measured in nCOTS per minute. Prediction is proportion of samples at a site that are positive."
    }
  }
  #################
  if( from=="manta"){
    from.effort <- length( from.metric)
    if( is.null( from.metric))
      stop( "Manta nScars not supplied.  Please do so (in argument from.metric)")
    from.metric <- sum( 1-(from.metric=="a"))
    ####
    if( tolower(to)=="cull"){
      if( is.null( to.effort)){
        warning( "Cull search time not supplied.  Assuming 473.6154 minutes (the average of the CCIP-D2 project data)")
        to.effort <- 473.6154 #mean of D2 data per site
      }
      predVal.PI <- predCullFromManta( manta_scars=from.metric, manta.nTows=from.effort, pred_type="PI", cullBottomTime=to.effort, intLevel=intLevel)
      predVal.CI <- predCullFromManta( manta_scars=from.metric, manta.nTows=from.effort, pred_type="CI", cullBottomTime=to.effort, intLevel=intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, predVal.CI[,c("CULL_density","density_lower","density_upper")], predVal.PI[,c("density_lower","density_upper")]))
      names( predVal) <- c("Manta_NSCARS","MANTA_NTOWS", "MANTA_density", "CULL_BottomTime", "mean_CULL", "lowerCI_CULL", "upperCI_CULL","lowerPI_CULL","upperPI_CULL")
      attr( predVal, "message") <-"Predicting Cull from Manta scars. Manta scars measured in presence or absence of scars in a tow. Cull measured in nCOTS per minute. Using the number of Manta nScars data points supplied as the number of tows"
    }
    ####
    if( tolower(to)=="manta"){
      tmp <- from.metric / from.effort
      attr( tmp, "message") <- "You are trying to predict Manta from Manta.  Please reconsider:-)  Returning input data."
      return( tmp) 
    }
    ####
    if( tolower(to)=="salad"){
      if( is.null( to.effort)){
        warning( "SALAD swept area not supplied.  Assuming 10193.33 (the average of the CCIP-D2 project data")
        to.effort <- 10193.33
      }
      predVal.PI <- predSALADFromManta( manta_scars = from.metric, manta.nTows=from.effort, SALADdist = to.effort, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predSALADFromManta( manta_scars = from.metric, manta.nTows=from.effort, SALADdist = to.effort, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, predVal.CI[, c( "SALAD_density","density_lower","density_upper")], predVal.PI[, c( "density_lower","density_upper")]))
      names( predVal) <- c("Manta_scars","Manta_nTows","Manta_prop","SALAD_EFFORT","SALAD_density","lowerCI_SALAD","upperCI_SALAD","lowerPI_SALAD","upperPI_SALAD")
      attr( predVal, "message") <- "Predicting SALAD from Manta measured in proportion of tows seen scars (input is just scar presence-absence). SALAD measured in nCOTS per metre."
    }
    ####
    if( tolower(to)=="edna_conc"){
      predVal.PI <- predEDNAconcFromManta( manta_scars=from.metric, manta.nTows=from.effort, pred_type="PI", intLevel=intLevel)
      predVal.CI <- predEDNAconcFromManta( manta_scars=from.metric, manta.nTows=from.effort, pred_type="CI", intLevel=intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, predVal.CI[, c( "eDNA_conc","eDNA_conc_lower","eDNA_conc_upper")], predVal.PI[, c( "eDNA_conc_lower","eDNA_conc_upper")]))
      names( predVal) <- c("Manta_scars","Manta_nTows","Manta_prop","eDNA_conc","lowerCI_eDNA_conc","upperCI_eDNA_conc","lowerPI_eDNA_conc","upperPI_eDNA_conc")
      attr( predVal, "message") <- "Predicting eDNA concentration from Manta. Manta measured in proportion of tows that saw scars (input is just scar presence-absence). Prediction is concentration per sample (sum of technical reps)."
    }
    ####
    if( tolower(to)=="edna_prop"){
      if( is.null( to.effort)){
        warning( "nSamples not supplied (as to.effort argument).  Assuming 12 (the number taken during the CCIP-D2 project data)")
        to.effort <- 12
      }
      predVal.PI <- predeDNAbinFromManta( manta_scars=from.metric, manta.nTows=from.effort, nreps=to.effort, pred_type="PI", intLevel = intLevel)
      predVal.CI <- predeDNAbinFromManta( manta_scars=from.metric, manta.nTows=from.effort, nreps=to.effort, pred_type="CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, from.effort, from.metric/from.effort, to.effort, predVal.CI[,c("Prob_eDNA","Prob_lower","Prob_upper")], predVal.PI[,c("Prob_lower","Prob_upper")]))
      names( predVal) <- c("Manta_scars","Manta_nTows","Manta_prop", "nSamples", "eDNA_prop","lowerCI_eDNA_prop","upperCI_eDNA_prop","lowerPI_eDNA_prop","upperPI_eDNA_prop")
      attr( predVal, "message") <- "Predicting eDNA proportion positive from Manta scars. Manta measured in proportion of tows that saw scars (input here is scar presence-absence). Prediction is proportion of samples at a site that are positive."
    }
  }
  #################
  if( from=="edna_conc"){
    if( is.null( from.metric))
      stop( "eDNA_conc not supplied.  Please do so (in argument from.metric).")
    if( length( from.metric) > 1)
      warning( "More than one concentration given. Using their average (assume that they are technical reps)")
    from.metric <- mean( from.metric, na.rm=TRUE)
    ####
    if( tolower(to)=="cull"){
      if( is.null( to.effort)){
        warning( "Cull search time not supplied.  Assuming 473.6154 minutes (the average of the CCIP-D2 project data)")
        to.effort <- 473.6154 #mean of D2 data per site
      }
      predVal.PI <- predCullFromeDNAoonc( eDNA_PCR_conc=from.metric, pred_type="PI", cullBottomTime=to.effort, intLevel=intLevel)
      predVal.CI <- predCullFromeDNAoonc( eDNA_PCR_conc=from.metric, pred_type="CI", cullBottomTime=to.effort, intLevel=intLevel)
      predVal <- unlist( c( from.metric, to.effort, predVal.CI[,c("CULL_density","density_lower","density_upper")], predVal.PI[,c("density_lower","density_upper")]))
      names( predVal) <- c("eDNA_conc", "CULL_effort","CULL_density", "lowerCI_CULL", "upperCI_CULL","lowerPI_CULL","upperPI_CULL")
      attr( predVal, "message") <- "Predicting Cull from eDNA_conc. Cull measured in nCOTS per minute"
    }
    ####
    if( tolower(to)=="manta"){
      if( is.null( to.nTows)){
        warning( "Manta nTows not supplied.  Assuming to.nTows=4 (the mode of the CCIP-D2 project data)")
        to.nTows <- 4
      }
      if( is.null( to.effort)){
        warning( "Manta tow distance not supplied (as to.effort argument). Assuming 200.5431 metres (the average of CCIP-D2 project data)")
        to.effort <- 200.5431
      }
      predVal.PI <- predMantaFromeDNA( eDNA_PCR_conc=from.metric, pred_type="PI", mantaBottomDist=to.effort, manta.nTows=to.nTows, intLevel=intLevel)
      predVal.CI <- predMantaFromeDNA( eDNA_PCR_conc=from.metric, pred_type="CI", mantaBottomDist=to.effort, manta.nTows=to.nTows, intLevel=intLevel)
      predVal <- unlist( c( from.metric, to.effort, to.nTows, predVal.CI[, c( "manta_density","density_lower","density_upper")], predVal.PI[, c( "density_lower","density_upper")]))
      names( predVal) <- c("eDNA_conc","manta.avTowLength","manta.nTows","Prob_Manta","lowerCI_Manta","upperCI_Manta","lowerPI_Manta","upperPI_Manta")
      attr( predVal, "message") <- "Predicting Manta from eDNA_conc. Prediction is probability of scar presence (given tow distance)"
    }
    ####
    if( tolower(to)=="salad"){
      if( is.null( to.effort)){
        warning( "SALAD swept area not supplied.  Assuming 10193.33 (the average of the CCIP-D2 project data")
        to.effort <- 10193.33
      }
      predVal.PI <- predSALADFromeDNAconc( eDNA_PCR_conc=from.metric, SALADarea=to.effort, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predSALADFromeDNAconc( eDNA_PCR_conc=from.metric, SALADarea=to.effort, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, to.effort, predVal.CI[, c( "SALAD_density","density_lower","density_upper")], predVal.PI[, c( "density_lower","density_upper")]))
      names( predVal) <- c("eDNA_conc","SALAD_EFFORT","SALAD_density","lowerCI_SALAD","upperCI_SALAD","lowerPI_SALAD","upperPI_SALAD")
      attr( predVal, "message") <- "Predicting SALAD from eDNA_conc. SALAD measured in nCOTS per hectare"
    }
    ####
    if( tolower(to)=="edna_conc"){
      attr( from.metric, "message") <- "You are trying to predict eDNA_conc from eDNA_conc.  Please reconsider:-)  Returning input data."
      return( from.metric)
    }
    ####
    if( tolower(to)=="edna_prop"){
      if( is.null( to.effort)){
        warning( "nSamples not supplied (as to.effort argument).  Assuming 12 (the number taken during the CCIP-D2 project data)")
        to.effort <- 12
      }
      predVal.PI <- predeDNAbinFromeDNAconc( eDNA_PCR_conc=from.metric, nreps=to.effort, pred_type="PI", intLevel = intLevel)
      predVal.CI <- predeDNAbinFromeDNAconc( eDNA_PCR_conc=from.metric, nreps=to.effort, pred_type="CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, to.effort, predVal.CI[,c("Prob_eDNA","Prob_lower","Prob_upper")], predVal.PI[,c("Prob_lower","Prob_upper")]))
      names( predVal) <- c("eDNA_conc", "nSamples", "eDNA_prop","lowerCI_eDNA_prop","upperCI_eDNA_prop","lowerPI_eDNA_prop","upperPI_eDNA_prop")
      attr( predVal, "message") <- "Predicting eDNA proportion positive from eDNA_conc. Prediction is proportion of samples at a site that are positive."      
    }
  }
  #################
  if( from=="edna_prop"){
    if( is.null( from.metric))
      stop( "eDNA_prop not supplied.  Please do so (in argument from.metric).")
    if( length( from.metric) > 1)
      warning( "More than one proportion given. Using their average -- this may not be appropriate...?")
    from.metric <- mean( from.metric, na.rm=TRUE)
    ####
    if( tolower(to)=="cull"){
      if( is.null( to.effort)){
        warning( "Cull search time not supplied.  Assuming 473.6154 minutes (the average of the CCIP-D2 project data)")
        to.effort <- 473.6154 #mean of D2 data per site
      }
      predVal.PI <- predCullFromeDNAbinary( av_PCR_binary=from.metric, pred_type="PI", cullBottomTime=to.effort, intLevel=intLevel)
      predVal.CI <- predCullFromeDNAbinary( av_PCR_binary=from.metric, pred_type="CI", cullBottomTime=to.effort, intLevel=intLevel)
      predVal <- unlist( c( from.metric, to.effort, predVal.CI[,c("CULL_density","density_lower","density_upper")], predVal.PI[,c("density_lower","density_upper")]))
      names( predVal) <- c("eDNA_prop", "CULL_effort","CULL_density", "lowerCI_CULL", "upperCI_CULL","lowerPI_CULL","upperPI_CULL")
      attr( predVal, "message") <- "Predicting Cull from eDNA_prop. Cull measured in nCOTS per minute."
    }
    ####
    if( tolower(to)=="manta"){
      if( is.null( to.nTows)){
        warning( "Manta nTows not supplied.  Assuming to.nTows=4 (the mode of the CCIP-D2 project data)")
        to.nTows <- 4
      }
      if( is.null( to.effort)){
        warning( "Manta tow distance not supplied (as to.effort argument). Assuming 200.5431 metres (the average of CCIP-D2 project data)")
        to.effort <- 200.5431
      }
      predVal.PI <- predMantaFromeDNAbin( av_PCR_bin=from.metric, pred_type="PI", mantaBottomDist=to.effort, manta.nTows=to.nTows, intLevel=intLevel)
      predVal.CI <- predMantaFromeDNAbin( av_PCR_bin=from.metric, pred_type="CI", mantaBottomDist=to.effort, manta.nTows=to.nTows, intLevel=intLevel)
      predVal <- unlist( c( from.metric, to.effort, to.nTows, predVal.CI[, c( "manta_density","density_lower","density_upper")], predVal.PI[, c( "density_lower","density_upper")]))
      names( predVal) <- c("eDNA_prop","manta.avTowLength","manta.nTows","Prob_Manta","lowerCI_Manta","upperCI_Manta","lowerPI_Manta","upperPI_Manta")
      attr( predVal, "message") <- "Predicting Manta from eDNA_prop. Prediction is probability of scar presence (given tow distance)"
    }
    ####
    if( tolower(to)=="salad"){
      if( is.null( to.effort)){
        warning( "SALAD swept area not supplied.  Assuming 10193.33 (the average of the CCIP-D2 project data")
        to.effort <- 10193.33
      }
      predVal.PI <- predSALADFromeDNAbin( av_PCR_binary=from.metric, SALADarea=to.effort, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predSALADFromeDNAbin( av_PCR_binary=from.metric, SALADarea=to.effort, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, to.effort, predVal.CI[, c( "SALAD_density","density_lower","density_upper")], predVal.PI[, c( "density_lower","density_upper")]))
      names( predVal) <- c("eDNA_prop","SALAD_EFFORT","SALAD_density","lowerCI_SALAD","upperCI_SALAD","lowerPI_SALAD","upperPI_SALAD")
      attr( predVal, "message") <- "Predicting SALAD from eDNA_prop. SALAD measured in nCOTS per hectare."
    }
    ####
    if( tolower(to)=="edna_conc"){
      predVal.PI <- predEDNAconcFromeDNAbin( av_PCR_bin=from.metric, pred_type="PI", intLevel=intLevel)
      predVal.CI <- predEDNAconcFromeDNAbin( av_PCR_bin=from.metric, pred_type="CI", intLevel=intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[, c( "eDNA_conc","eDNA_conc_lower","eDNA_conc_upper")], predVal.PI[, c( "eDNA_conc_lower","eDNA_conc_upper")]))
      names( predVal) <- c("eDNA_prop","eDNA_conc","lowerCI_eDNA_conc","upperCI_eDNA_conc","lowerPI_eDNA_conc","upperPI_eDNA_conc")
      attr( predVal, "message") <- "Predicting eDNA concentration from eDNA_prop. Prediction is concentration per sample (sum of technical reps)."
    }
    ####
    if( tolower(to)=="edna_prop"){
      attr( from.metric, "message") <- "You are trying to predict eDNA_prop from eDNA_prop.  Please reconsider:-)  Returning input data."
      return( from.metric)
    }
  }
  
  #################
  ####  CORAL
  
  #################
  if( from=="coral_transect"){
    if( is.null( from.metric))
      stop( "Coral Transect data not supplied (count of corals). Please do so (in argument from.metric).")
    if( is.null( from.effort)){
      warning( "Coral Transect effort data not supplied (count of hard corals). Assuming 100 in each transect (target used in CCIP-D2 project).")
      from.effort <- rep( 100, length=length(from.metric))
    }
    from.metric <- sum( from.metric, na.rm=TRUE) / sum( from.effort, na.rm=TRUE)
    ####
    if( tolower(to)=="coral_transect"){
      attr( from.metric, "message") <- "You are trying to predict transect data from transect data. Please reconsider:-) Return input data."
      return( from.metric)
    }
    ####
    if( tolower(to)=="coral_transom"){
      predVal.PI <- predTransomFromTransect(transect_density = from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predTransomFromTransect(transect_density = from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("TransomFromTransect_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Transect","Coral_Transom","lowerCI_Coral_Transom", "upperCI_Coral_Transom", "lowerPI_Coral_Transom", "upperPI_Coral_Transom")
      attr( predVal, "message") <- "Predicting coral_transom from coral_transect. Prediction is coral coverage proportion."
    }
    ####
    if( tolower(to)=="coral_deep"){
      predVal.PI <- predDeepFromTransect(transect_density = from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predDeepFromTransect(transect_density = from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("DeepFromTransect_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Transect","Coral_Deep","lowerCI_Coral_Deep", "upperCI_Coral_Deep", "lowerPI_Coral_Deep", "upperPI_Coral_Deep")
      attr( predVal, "message") <- "Predicting coral_deep from coral_transect. Prediction is coral coverage proportion."
    }
    ####
    if( tolower(to)=="coral_manta"){
      predVal.PI <- predMantaCoralFromTransect(transect_density = from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predMantaCoralFromTransect(transect_density = from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("MantaFromTransect_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Transect","Coral_Manta","lowerCI_Coral_Manta", "upperCI_Coral_Manta", "lowerPI_Coral_Manta", "upperPI_Coral_Manta")
      attr( predVal, "message") <- "Predicting coral_manta from coral_transect. Prediction is coral coverage proportion (enumerated from diver scores)."
    }
  }
  if( from=="coral_transom"){
    if( is.null( from.metric))
      stop( "Coral Transom data not supplied (proportion of corals). Please do so (in argument from.metric).")
    ####
    if( tolower( to)=="coral_transect"){
      if( is.null( to.effort)){
        warning( "Coral transect effort not supplied (number of measured locations with hard substrate). Assuming that it is 400 (as was the target in the CCIP-D2 project).")
        to.effort <- 400
      }
      predVal.PI <- predTransectFromTransom(transom_perc=from.metric, numHard=to.effort, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predTransectFromTransom(transom_perc=from.metric, numHard=to.effort, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("TransectsFromTransom_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Transom","Coral_Transect","lowerCI_Coral_Transect", "upperCI_Coral_Transect", "lowerPI_Coral_Transect", "upperPI_Coral_Transect")
      attr( predVal, "message") <- "Predicting coral_transect from coral_transom. Prediction is the probability of a measured location being coral (and not bare hard substrate)."
    }
    ####
    if( tolower(to)=="coral_transom"){
      attr( from.metric, "message") <- "You are trying to predict transom data from transom data. Please reconsider:-) Return input data."
      return( from.metric)
    }
    ####
    if( tolower( to)=="coral_deep"){
      predVal.PI <- predDeepFromTransom(transom_perc=from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predDeepFromTransom(transom_perc=from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("DeepFromTransom_perc","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Transom","Coral_Deep","lowerCI_Coral_Deep", "upperCI_Coral_Deep", "lowerPI_Coral_Deep", "upperPI_Coral_Deep")
      attr( predVal, "message") <- "Predicting coral_deep from coral_transom. Prediction is the proportion of image that is coral."
    }
    ####
    if( tolower(to)=="coral_manta"){
      predVal.PI <- predMantaCoralFromTransom(transom_perc = from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predMantaCoralFromTransom(transom_perc = from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("MantaFromTransom_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Transom","Coral_Manta","lowerCI_Coral_Manta", "upperCI_Coral_Manta", "lowerPI_Coral_Manta", "upperPI_Coral_Manta")
      attr( predVal, "message") <- "Predicting coral_manta from coral_transom. Prediction is coral coverage proportion (enumerated from diver scores)."
    }
  }
  if( from=="coral_deep"){
    if( is.null( from.metric))
      stop( "Coral deep data not supplied (proportion of corals). Please do so (in argument from.metric).")
    ####
    if( tolower( to)=="coral_transect"){
      if( is.null( to.effort)){
        warning( "Coral transect effort not supplied (number of measured locations with hard substrate). Assuming that it is 400 (as was the target in the CCIP-D2 project).")
        to.effort <- 400
      }
      predVal.PI <- predTransectFromDeep(deep_perc = from.metric, numHard=to.effort, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predTransectFromDeep(deep_perc = from.metric, numHard=to.effort, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("TransectsFromDeep_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Deep","Coral_Transect","lowerCI_Coral_Transect", "upperCI_Coral_Transect", "lowerPI_Coral_Transect", "upperPI_Coral_Transect")
      attr( predVal, "message") <- "Predicting coral_transect from coral_deep. Prediction is the probability of a measured location being coral (and not bare hard substrate)."
    }
    ####
    if( tolower( to)=="coral_transom"){
      predVal.PI <- predTransomFromDeep(deep_perc=from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predTransomFromDeep(deep_perc=from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("TransomFromDeep_perc","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Deep","Coral_Transom","lowerCI_Coral_Transom", "upperCI_Coral_Transom", "lowerPI_Coral_Transom", "upperPI_Coral_Transom")
      attr( predVal, "message") <- "Predicting coral_transom from coral_deep. Prediction is the proportion of image that is coral."
    }
    ####
    if( tolower(to)=="coral_deep"){
      attr( from.metric, "message") <- "You are trying to predict deep data from deep data. Please reconsider:-) Return input data."
      return( from.metric)
    }
    ####
    if( tolower(to)=="coral_manta"){
      predVal.PI <- predMantaCoralFromDeep(deep_perc=from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predMantaCoralFromDeep(deep_perc = from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("MantaFromDeep_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Deep","Coral_Manta","lowerCI_Coral_Manta", "upperCI_Coral_Manta", "lowerPI_Coral_Manta", "upperPI_Coral_Manta")
      attr( predVal, "message") <- "Predicting coral_manta from coral_deep. Prediction is coral coverage proportion (enumerated from diver scores)."
    }
  }
  ####
  if( from=="coral_manta"){
    if( is.null( from.metric))
      stop( "Coral manta data not supplied (number of measurements coral). Please do so (in argument from.metric).")
    from.metric <- mean( from.metric, na.rm=TRUE)
    ####
    if( tolower( to)=="coral_transect"){
      if( is.null( to.effort)){
        warning( "Coral transect effort not supplied (number of measured locations with hard substrate). Assuming that it is 400 (as was the target in the CCIP-D2 project).")
        to.effort <- 400
      }
      predVal.PI <- predTransectFromManta(manta_perc=from.metric, numHard=to.effort, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predTransectFromManta(manta_perc=from.metric, numHard=to.effort, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("TransectsFromManta_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Manta","Coral_Transect","lowerCI_Coral_Transect", "upperCI_Coral_Transect", "lowerPI_Coral_Transect", "upperPI_Coral_Transect")
      attr( predVal, "message") <- "Predicting coral_transect from coral_manta. Prediction is the probability of a measured location being coral (and not bare hard substrate)."
    }
    ####
    if( tolower( to)=="coral_transom"){
      predVal.PI <- predTransomFromManta(manta_perc=from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predTransomFromManta(manta_perc=from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("TransomFromManta_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Manta","Coral_Transom","lowerCI_Coral_Transom", "upperCI_Coral_Transom", "lowerPI_Coral_Transom", "upperPI_Coral_Transom")
      attr( predVal, "message") <- "Predicting coral_transom from coral_manta. Prediction is the proportion of image that is coral."
    }
    ####
    if( tolower(to)=="coral_deep"){
      predVal.PI <- predDeepFromManta(manta_perc=from.metric, pred_type = "PI", intLevel = intLevel)
      predVal.CI <- predDeepFromManta(manta_perc=from.metric, pred_type = "CI", intLevel = intLevel)
      predVal <- unlist( c( from.metric, predVal.CI[,c("DeepFromManta_coral","coral_lower","coral_upper")], predVal.PI[,c("coral_lower","coral_upper")]))
      names( predVal) <- c( "Coral_Manta","Coral_Transom","lowerCI_Coral_Transom", "upperCI_Coral_Transom", "lowerPI_Coral_Transom", "upperPI_Coral_Transom")
      attr( predVal, "message") <- "Predicting coral_deep from coral_manta. Prediction is the proportion of image that is coral."
    }
    ####
    if( tolower(to)=="coral_manta"){
      attr( from.metric, "message") <- "You are trying to predict manta data from manta data. Please reconsider:-) Return input data."
      return( from.metric)
    }
  }
  
  return( predVal)
}


#predictSingleMeasurement(from = "coral_manta", to="coral_manta", from.metric = D2data$mantaTow$percHard[1])
