COTS_calibrate <- function( from, to, from.metric, from.effort=NULL, to.effort=NULL, to.nTows=NULL, sites=NULL, intLevel=0.95)

\name{COTS_calibrate}
\alias{COTS_calibrate}
\title{Predicts measurements from one measurement tool in another measurement too.}
\description{This function predicts one measurement from another, for COTS related measurement tools. The COTS tools considered, which directly or indirectly measure the number of COTS, are: SALAD, CULL, eDNA (concentration), eDNA (proportion positives) and Manta Tow (scars). The tools that measure coral cover are: Manta Tow, Coral Transects (part of SALAD), ReefScan-Transom and ReefScan-Deep. The predictions are made from a conditional modelling exercise (GLMM) for all two-way combinations of measurement tools, within COTS and Coral measurements.}
\usage{
 COTS_calibrate( from, to, 
	from.metric, from.effort=NULL, 
	to.effort=NULL, to.nTows=NULL, 
	sites=NULL, 
	intLevel=0.95)
}
\arguments{
\item{from}{A string giving the type of measurement to be converted. Options are: "Manta","eDNA_conc","eDNA_prop","coral_transom","coral_deep","coral_manta","SALAD","Cull". Capitalisation should not matter, but the presence of underscores will. See Details below.}
\item{to}{A string giving the type of the output measurement. Options are the same as \code{from}. If \code{from} and \code{to} are the same, then the \code{from} data are returned with a suitable warning.}
\item{from.metric}{The observation values for the \code{from} measurement. These are, as far as possible, recorded in the same format as is usually done. See Details below for some further information.}
\item{from.effort}{The amount of effort for the \code{from.metric} observation. These are, as far as possible, recorded in the same format as it routinely done. The units will change depending on the measuring too. For example, effort for cull dives are measured in minutes, whilst effort for Manta Tows is measured in metres. If NULL (default), AND effort is needed for prediction, then the average of the CCIP-D2 field trip is inserted.}
\item{to.effort}{The amount of effort that is assumed for the prediction. If NULL (default), AND it is needed for prediction, then the average of the CCIP-D2 field trip is inserted.}
\item{to.nTows}{The number of tows used in the prediction. Only used for manta tows whose measurement is at a tow-level. Default is NULL, which if manta predictions are needed, then 4 is inserted (the average of the CCIP-D2 field trip).}
\item{sites}{A character or factor indicating which observations belong to which sites. This is an important argument. It tells the function which observations should be treated as a group and which ones as different groups. When dealing with data from multiple reefs, make sure that the site labels are unique (e.g. a site "Z_02" may appear at both reefs and has the potential to be 'lumped' together inadvertently).}
\item{intLevel}{A scalar between 0.5 and 1, usually 0.95, or 0.9. This gives the limit of the confidence and prediction intervals from the prediction.}
}
\details{The prediction can only occur within COTS measurements and within Coral measurements. There is no scope for predicting COTS from Coral, and vice-versa.

The \code{from.metric} and \code{from.effrot} should be (for the different \code{from} sampling tools)
\describe{
\item{Cull}{\code{from.metric} is the number of COTS recorded at the site's cull dive. \code{from.effort} is the number of minutes. Note that all values that share a site label will be summed}
\item{SALAD}{Whilst there is a number of different metrics possible for the \code{from.metric}, we use swept area in this package (unit is metre-squared, not hectare). The \code{from.metric} is the total number of COTS or Scars recorded. These choices were made to match the desires of managers and researchers at the time of prodiction of RosettaCOTS.}
\item{Manta}{A vector of character representing each tow at a site. A value of "a" is given when no scars were seen. A value of "p" or "c" is used when there are scars.}
\item{eDNA_conc}{The concentration values for a eDNA measurement. When multiple values are given for each site, then the average is taken.}
\item{eDNA_prop}{The proportion of samples at a site with a positive reading. Alternatively, a vector of all samples at a site with a (0,1) reading. The proportion is then taken internally.}
\item{coral_manta}{A numeric giving the percentage of coral cover at a site. Note that this requires the ordered observations to be quantified by replacing the observations with their category's midpoint. If more than one observation is provided per site (multiple tows), then these are averaged.}
\item{coral_deep}{A numeric giving the proportion of coral cover (hard and soft).}
\item{coral_transom}{A numeric giving the proporion of coral cover (hard and soft).}
}

The \code{to.effort} argument specifies how much search effort is assumed in the prediction. If NULL (default), then it is assumed to be a sensible value chosen on the basis of the CCIP-D2 field trip.
}
\value{A data.frame that contains the relevant summary of the input data, at the site level, and the predicted output data. This includes a point prediction (for output) as well as confidence intervals (how far the mean prediction may vary), and prediction intervals (how much a new observation may vary). It is envisaged that the prediction intervals will be a lot more useful than the confidence intervals, in imagined analyses.
}

\examples{
#predicting Manta from SALAD (using default number of manta tows and their length).
SALADdat <- as.data.frame( readxl::read_excel( system.file("extdata", "CompiledData1_standardised.xlsx", package="RosettaCOTS"), sheet = "SALAD", na="NA"))
#make sure site names are unique
SALADdat$reefSite <- paste0(SALADdat$Reef_ID,":",SALADdat$Site)
predManta_fromSALAD <- COTS_calibrate( from="Salad", to="Manta", from.metric=SALADdat$Total, from.effort=SALADdat$Distance, to.effort=NULL, to.nTows=NULL, sites=SALADdat$reefSite)
#having a little look at the output -- each row is a site. In this example, all sites have the same manta tow length and number of tows per site.
utils::head( predManta_fromSALAD)

#predicting Cull rate from Manta scars
MantaDat <- as.data.frame( readxl::read_excel( system.file("extdata", "CompiledData1_standardised.xlsx", package="RosettaCOTS"), sheet = "MantaTow", na="NA"))
#make sure site names are unique
MantaDat$reefSite <- paste0(MantaDat$Reef_ID,":",MantaDat$Site)
#Get an idea of the Manta Tow scar data
table( MantaDat$reefSite, MantaDat$Feeding_Scars)
predCull_fromManta <- COTS_calibrate( from="Manta", to="CULL", from.metric=MantaDat$Feeding_Scars, from.effort=NULL, to.effort=NULL, to.nTows=NULL, sites=MantaDat$reefSite)
#having a little look at the output -- each row is a site. In this example, all sites have the same cull tow length.
utils::head( predCull_fromManta)
}



\author{Scott D. Foster}
