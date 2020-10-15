### DELPHINIUM EXALTATUM ECOLOGIAL NICHE MODELING
### Adam B. Smith | Missouri Botanical Garden | adamDOTsmithATmobotDOTorg | 2018-07
###
### source('C:/Ecology/Drive/Research/Delphinium exaltatum (Christy Edwards)/delphiniumExaltatum/Delphinium 00 Pre-submission.r')

	rm(list=ls())
	memory.limit(memory.limit() * 2^30)
	gc()

	makePlot <- TRUE # make all plots
	# makePlot <- FALSE # make some plots
	
### CONTENTS ###
### libraries and functions ###
### variables and options ###
### 00 collate species data ###
### 01 clean species data ###
### 02 incorporate presences into county shapefile, assign population, and define accessible regions ###
### calculate common climate PCA ###
### crop states/provinces to study region ###
### 03 match current environmental data with counties ###
### 04 assign geo-folds ###
### 05 calculate county weights ###
### 06 ENMs ###
### 07 visualize model performance ###
### create extractions of paleoclimate data ###
### 08a get maximum predicted value across all time ###
### 08b make maps of present ###
### 09 make maps of past ###
### 10 animate maps ###
### 11 make response curves ###

###############################
### libraries and functions ###
###############################

	library(taxize)
	library(raster)
	library(dismo)
	library(rgeos)
	
	library(scales)
	library(fpCompare)
	library(apng)

	library(omnibus)
	library(enmSdm)
	library(legendary)

	plotSimple <- function(w) {
	
		# w weights used for alpha in color (in range of [0, 1])
	
		if (class(w) == 'character') w <- delEx@data[ , w]
		if (anyNA(w)) w[is.na(w)] <- 0
	
		col <- alpha('darkblue', w)
		plot(delEx[delEx$accessible_popSpecies_extentBroad, ], border=NA)
		plot(delEx, col=alpha('darkblue', w), border=NA, add=TRUE)
		plot(delEx[delEx$delExObs == 1, ], add=TRUE)
		
	}
	
	plotNice <- function(v, title, leg, rast=NULL, stretch=FALSE, legMax=NULL, color='deeppink4', ext='broad', cexMult=1) {

		# v			values of response or name in delEx with values (from broad background)
		# title		title
		# leg		name of response (for legend) 
		# rast		raster to plot *under* the delEx polygon; ignored if NULL
		# stretch	Logical, if TRUE then stretch values to [lower, upper]
		# legMax	maximum value for scale of color legend, not of values shown on map
		# color		color
		# ext		'broad', 'all'
		# cexMult	multipler for cex

		say('   Plotting: ', title)
		
		load('./Other/Canada Admin Level 1 Albers Croppped to Study Region.Rdata')
		load('./Other/USA Admin Level 1 Albers Croppped to Study Region.Rdata')
		load('./Other/Mexico Admin Level 1 Albers Croppped to Study Region.Rdata')
	
		thisDelEx <- sp::spTransform(delEx, getCRS('albersNA', TRUE))
	
		if (class(v) == 'character') v <- thisDelEx@data[ , v]
	
		if (stretch) {
			w <- stretchMinMax(v, na.rm=TRUE)
			if (is.null(legMax)) legMax <- 1
		} else {
			w <- v
			if (is.null(legMax)) max(v, na.rm=TRUE)
		}
	
		if (anyNA(v)) w[is.na(v)] <- 0
	
		par(mar=mar)

		if (!is.null(rast)) {

			rast <- projectRaster(rast, crs=getCRS('albersNA'))
			breaks <- seq(0, 1, by=0.01)
			colFx <- grDevices::colorRampPalette(c('gray90', color))
			cols <- colFx(100)
			
			thisExt <- if (ext == 'broad') {
				extent(thisDelEx[which(thisDelEx$accessible_popSpecies_extentBroad), ])
			} else {
				extent(thisDelEx)
			}
			
			plot(rast, ext=thisExt, breaks=breaks, col=cols, legend=FALSE, box=FALSE, axes=FALSE)
			plot(thisDelEx, col='white', border='white', add=TRUE)
			
		} else {
		
			if (ext == 'broad') {
				plot(thisDelEx[which(thisDelEx$accessible_popSpecies_extentBroad), ], border=NA)
			} else {
				plot(thisDelEx, border=NA)
			}
		
		}
		
		colFx <- grDevices::colorRampPalette(c('gray90', color))
		cols <- colFx(100)
		w <- round(w * 100)
		if (any(w == 0)) w[w==0] <- 1
		cols <- cols[w]
		plot(thisDelEx, col=cols, border='gray', lwd=lwdCounties, add=TRUE)
		
		plot(canadaAdmin1AlbersNACropped, add=TRUE, lwd=lwdStates)
		plot(usaAdmin1AlbersNACropped, add=TRUE, lwd=lwdStates)
		plot(mexicoAdmin1AlbersNACropped, add=TRUE, lwd=lwdStates)

		plot(thisDelEx[thisDelEx$delExObs == 1, ], lwd=lwdCounties, add=TRUE)
		
		legendGrad('bottomright', inset=c(0.1, 0), height=0.6, width=0.1, col=c('white', color), labels=pretty(c(0, legMax)), title=leg, xpd=NA, boxBorder=NA, gradAdjY=c(0.175, 0.91), cex=0.7 * cexMult * cexLegend, swatches=list(list(swatchAdjY=c(0.93, 1), col='white', border='black', labels='Presence')))
		
		title(main=title, cex.main=cexMult * cexMain, sub=date(), cex.sub=cexSub * cexMult)
		
	}

	# equalizes train/test weights
	equalizeWeights <- function(wPresLowland, wBgLowland, wPresHighland, wBgHighland) {
	
		# wPresLowland, wPresHighland weights vectors for lowland and highland presences
		# wBgLowland, wBgHighland weights vectors for lowland and highland BG sites
	
		sumPresLowland <- sum(wPresLowland)
		sumBgLowland <- sum(wBgLowland)
	
		if (sumPresLowland > sumBgLowland) {
			wBgLowland <- wBgLowland * (sumPresLowland / sumBgLowland)
		} else {
			wPresLowland <- wPresLowland * (sumBgLowland / sumPresLowland)
		}
	
		sumPresHighland <- sum(wPresHighland)
		sumBgHighland <- sum(wBgHighland)
	
		if (sumPresHighland > sumBgHighland) {
			wBgHighland <- wBgHighland * (sumPresHighland / sumBgHighland)
		} else {
			wPresHighland <- wPresHighland * (sumBgHighland / sumPresHighland)
		}
	
		if (sum(wPresLowland) %!=% sum(wBgLowland) | sum(wPresHighland) %!=% sum(wBgHighland)) stop('Sum of presence weights != sum of background weights.')
		
		wLowland <- sum(c(wPresLowland, wBgLowland))
		wHighland <- sum(c(wPresHighland, wBgHighland))
		
		if (wLowland > wHighland) {
		
			wPresHighland <- wPresHighland * wLowland / wHighland
			wBgHighland <- wBgHighland * wLowland / wHighland
			
		} else {
			
			wPresLowland <- wPresLowland * wHighland / wLowland
			wBgLowland <- wBgLowland * wHighland / wLowland
			
		}
		
		out <- list()
		out$wPresLowland <- wPresLowland
		out$wBgLowland <- wBgLowland
		out$wPresHighland <- wPresHighland
		out$wBgHighland <- wBgHighland
		out
		
	}
	
	# calculate train/testCBI for species, lowland, highland
	testModel <- function(model, testData, ...) {
	
		predictions <- predict(model, testData, type='response')
		
		# species-level
		predPres <- predictions[testData$presBg == 1]
		predBg <- predictions[testData$presBg == 0]
		
		wPres <- testData$weight[testData$presBg == 1]
		wBg <- testData$weight[testData$presBg == 0]
		
		cbiSpecies <- contBoyce(pres=predPres, bg=predBg, presWeight=wPres, bgWeight=wBg)
		
		# lowland only
		predPres <- predictions[testData$presBg == 1 & testData$population == 'lowland']
		predBg <- predictions[testData$presBg == 0 & testData$population == 'lowland']
		
		wPres <- testData$weight[testData$presBg == 1 & testData$population == 'lowland']
		wBg <- testData$weight[testData$presBg == 0 & testData$population == 'lowland']
		
		cbiLowland <- contBoyce(pres=predPres, bg=predBg, presWeight=wPres, bgWeight=wBg)
		
		# highland only
		predPres <- predictions[testData$presBg == 1 & testData$population == 'highland']
		predBg <- predictions[testData$presBg == 0 & testData$population == 'highland']
		
		wPres <- testData$weight[testData$presBg == 1 & testData$population == 'highland']
		wBg <- testData$weight[testData$presBg == 0 & testData$population == 'highland']
		
		cbiHighland <- contBoyce(pres=predPres, bg=predBg, presWeight=wPres, bgWeight=wBg)
		
		out <- c(cbiSpecies, cbiLowland, cbiHighland)
		names(out) <- c('cbiSpecies', 'cbiLowland', 'cbiHighland')
		out
		
	}	

#############################
### variables and options ###
#############################

	setwd('C:/Ecology/Drive/Research/Delphinium exaltatum (Christy Edwards)')

	options(stringsAsFactors=FALSE)
	rasterOptions(format='GTiff', overwrite=TRUE)

	longLat <- c('longitude', 'latitude')
	
	# colors
	speciesCol <- 'magenta'
	lowlandCol <- 'chartreuse3'
	highlandCol <- 'firebrick1'
	contestedCol <- 'yellow'

	# plotting params, assuming res = 300
	mar <- c(2, 2, 3, 2) + 0.1
	cexMain <- 1
	cexSub <- 0.4
	cexLegend <- 1.3
	lwdCounties <- 1.2
	lwdStates <- 2
		
	# buffer widths to define regions
	narrow <- 300 # narrow accessible area (km)
	broad <- 600 # broad accessible area (km)
	entire <- 1200 # projection region
	
	# predictor variables (from Lorenz, D.J., Nieto-Lugilde, D., Blois, J.L., Fitzpatrick, M.C., and Williams, J.W.  2016.  Downscaled and debiased climate simulations for North America from 21,000 years ago to 2100AD.  Scientific Data 3:160048.)
	vars <- c('an_avg_TMAX', 'an_avg_TMIN', 'an_sd_TMAX', 'an_sd_TMIN', 'an_sum_GDD0', 'an_sd_GDD0', 'an_avg_ETR')
	varsNice <- c('Max Temp', 'Min Temp', 'Max Temp Var', 'Min Temp Var', 'GDD', 'GDD Var', 'Evapotrans')
	
	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries WGS84/Canada Admin Level 2.Rdata')
	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries WGS84/USA Admin Level 2.Rdata')
	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries WGS84/Mexico Admin Level 2.Rdata')

	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries Albers North America/Canada Admin Level 1.Rdata')
	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries Albers North America/USA Admin Level 1.Rdata')
	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries Albers North America/Mexico Admin Level 1.Rdata')

	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries Albers North America/Canada Admin Level 2.Rdata')
	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries Albers North America/USA Admin Level 2.Rdata')
	load('D:/Ecology/Political Geography/GADM/ver3pt6/Countries Albers North America/Mexico Admin Level 2.Rdata')

	dirCreate('./Analysis')
	dirCreate('./ENMs')
	
# say('###############################')
# say('### 00 collate species data ###')
# say('###############################')

	# bien <- read.csv('./Data/BIEN/00 Delphinium_exaltatum_observations.csv')
	# gbif <- read.csv('./Data/GBIF/02 occurrence Delphinium exaltatum.csv')
	# idigbio <- read.csv('./Data/iDigBio/01 Delphinium exaltatum/occurrence_raw.csv')
	# tropicos <- read.csv('./Data/TROPICOS/00 Delphinium_exaltatum.csv')
	# christy <- read.csv('./Data/Christy Edwards (Also in TROPICOS)/01 Vouchered Specimens.csv')
	
	# cw <- data.frame(
		# origSpecies=c('scrubbed_species_binomial', 'scientificName', 'dwc.scientificName', '_Delphinium exaltatum',  '_Delphinium exaltatum'),
		# origLong=c('longitude', 'decimalLongitude', 'dwc.decimalLongitude', 'Longitude', 'longitude'),
		# origLat=c('latitude', 'decimalLatitude', 'dwc.decimalLatitude', 'Latitude', 'latitude'),
		# origCoordUncer_m=c(NA, 'coordinateUncertaintyInMeters', 'dwc.coordinateUncertaintyInMeters', NA, '_5'),
		# datum_m=c(NA, NA, 'c("dwc.geodeticDatum", "dwc.verbatimCoordinateSystem")', NA, NA),
		# locality=c('locality_description', 'locality', 'dwc.locality', NA, 'siteCountyState'),
		# coordRemarks=c(NA, 'c("informationWithheld", "dataGeneralizations", "locationAccordingTo", "locationRemarks", "georeferencedBy", "georeferenceProtocol", "georeferenceSources", "georeferenceVerificationStatus", "georeferenceRemarks", "issue")', 'c("dwc.informationWithheld", "dwc.dataGeneralizations", "dwc.locationAccordingTo", "dwc.locationRemarks", "dwc.georeferencedBy", "dwc.georeferenceProtocol", "dwc.georeferenceSources", "dwc.georeferenceVerificationStatus", "dwc.georeferenceRemarks")', NA, 'code'),
		# georefIssuesBien=c('is_geovalid', NA, NA, NA, NA),
		# georefIssuesGbif=c(NA, 'hasGeospatialIssues', NA, NA, NA),
		# origElev=c('elevation_m', 'verbatimElevation', 'dwc.verbatimElevation', 'Elevation', NA),
		# admin0=c('country', 'countryCode', 'dwc.country', 'Country', '_United States'),
		# admin1=c('state_province', 'stateProvince', 'dwc.stateProvince', 'Upper', 'state'),
		# admin2=c('county', 'county', 'dwc.county', 'Lower', 'county'),
		# admin3=c(NA, 'municipality', 'dwc.municipality', NA, NA),
		# origDate=c('date_collected', 'eventDate', 'dwc.eventDate', 'Date', 'year'),
		# basis=c('observation_type', 'basisOfRecord', 'dwc.basisOfRecord', '_Vouchered specimen', '_Vouchered specimen'),
		# cultivated=c('is_cultivated_observation', NA, NA, NA, NA),
		# collectedBy=c('collector', NA, 'dwc.recordedBy', 'Collectors', 'voucher'),
		# idBy=c('identified_by', 'identifiedBy', 'dwc.identifiedBy', NA, NA),
		# notes=c(NA, 'c("occurrenceRemarks", "occurrenceStatus", "habitat", "fieldNotes")', 'c("dwc.occurrenceRemarks", "dwc.occurrenceStatus", "dwc.habitat", "dwc.fieldNotes")', NA, NA),
		# source=c('datasource', 'institutionCode', 'c("dcterms.rightsHolder", "dwc.datasetName", "dwc.institutionCode")', '_TROPICOS', '_Christy Edwards'),
		# downloadFrom=c('_BIEN', '_GBIF', '_iDigBio', '_TROPICOS', '_Christy Edwards')
	# )
		
	# records <- combineDf(bien, gbif, idigbio, tropicos, christy, crosswalk=cw, verbose=TRUE)
	# records$origCoordUncer_m <- as.numeric(records$origCoordUncer_m)
	# records$georefIssuesGbif <- as.logical(records$georefIssuesGbif)
	
	# # add presences to counties from NPS map
	# nps <- shapefile('./Data/NPS/countiesWithDelphiniumExaltatum')
	# for (i in 1:nrow(nps)) {
	
		# records <- rbind(
			# records,
			# data.frame(
				# origSpecies='Delphinium exaltatum',
				# origLong=NA,
				# origLat=NA,
				# origCoordUncer_m=NA,
				# datum_m=NA,
				# locality=NA,
				# coordRemarks=NA,
				# georefIssuesBien=NA,
				# georefIssuesGbif=NA,
				# origElev=NA,
				# admin0='United States',
				# admin1=nps$NAME_1[i],
				# admin2=nps$NAME_2[i],
				# admin3=NA,
				# origDate=NA,
				# basis='Observation (NPS)',
				# cultivated=FALSE,
				# collectedBy=NA,
				# idBy=NA,
				# notes=NA,
				# source='NPS',
				# downloadFrom='NPS'
			# )
		# )
		
	# }
	
	# write.csv(records, './Data/00 Delphinium exaltatum - All Occurrences.csv', row.names=FALSE)

# say('#############################')
# say('### 01 clean species data ###')
# say('#############################')

	# records <- read.csv('./Data/00 Delphinium exaltatum - All Occurrences.csv')
	
	# say('   ### species')
	# say('   ###########')
	
		# fill <- rep(NA, nrow(records))
		# records <- insertCol(data.frame(speciesClean=fill), into=records, at='origSpecies')
		# records$speciesClean[records$origSpecies == 'Delphinium exaltatum'] <- 'Delphinium exaltatum'
		# records$speciesClean[records$origSpecies == 'Delphinium exaltatum Ait.'] <- 'Delphinium exaltatum'
		# records$speciesClean[records$origSpecies == 'Delphinium urceolatum Jacq.'] <- 'Delphinium urceolatum'
		# records$speciesClean[records$origSpecies == 'Delphinium exaltatum Aiton'] <- 'Delphinium exaltatum'
		# records$speciesClean[records$origSpecies == 'Delphinium treleasi Bush ex K.C. Davis'] <- 'Delphinium treleasi'
		# records$speciesClean[records$origSpecies == 'Delphinium exaltatum W.T.Aiton'] <- 'Delphinium exaltatum'
		# records$speciesClean[records$origSpecies == 'Delphinium intermedium Aiton'] <- 'Delphinium intermedium'
		# records$speciesClean[records$origSpecies == 'Angelica venenosa (J. Greenway) Fernald'] <- 'Angelica venenosa'
		# records$speciesClean[records$origSpecies == 'Delphinium urceolatum'] <- 'Delphinium urceolatum'
		# records$speciesClean[records$origSpecies == 'Delphinium lilacinum'] <- 'Delphinium lilacinum'
		# records$speciesClean[records$origSpecies == 'Agrimonia rostellata Wallr.'] <- 'Agrimonia rostellata'
		# records$speciesClean[records$origSpecies == 'Delphinium californicum Torr. & A.Gray'] <- 'Delphinium californicum'
		# records$speciesClean[records$origSpecies == 'Delphinium hansenii subsp. kernense'] <- 'Delphinium hansenii'
		# records$speciesClean[records$origSpecies == 'Polygala senega L.'] <- 'Polygala senega'

		# # speciesTnrs <- tnrs(records$speciesClean)
		# # save(speciesTnrs, file='./Data/TNRS.RData')
		# load('./Data/TNRS.RData')

		# records <- insertCol(data.frame(species=fill), into=records, at='speciesClean')
		# records$species <- speciesTnrs$acceptedname
		
		# records <- records[records$species == 'Delphinium exaltatum', ]

	# say('   ### dates')
	# say('   #########')

		# years <- yearFromDate(x=records$origDate, yearLast=TRUE)
		# years[years == 0] <- NA
		
		# records <- insertCol(data.frame(year=years), into=records, at='origDate')

		# # # inspect
		# # write.csv(records, 'C:/ecology/!scratch/Records with Raw Dates.csv', row.names=FALSE)
		
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'Christopher T. Frye'] <- 2012
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9955 & records$collectedBy == 'Julian A. Steyermark'] <- 1955
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9955 & records$collectedBy == 'Clair L. Kucera'] <- 1955
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'Susan Broadrington & Randall J. Evans'] <- 2012
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9913 & records$collectedBy == 'Ginger Allington'] <- 2013
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'Grace Johnson'] <- 2012
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'George Yatskievych'] <- 2012
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9915 & records$collectedBy == 'George Yatskievych'] <- 2015
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9985 & records$collectedBy == 'Arthur Christ'] <- 1985
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'Yari Johnson'] <- 2012
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9915 & records$collectedBy == 'Rebekah Mohn'] <- 2015
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'S. P. Grund'] <- 2015
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'Rebecca Cook'] <- 2012
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'Abigail Hyduke & Thomas F. Wieboldt'] <- 2012
		# records$year[records$downloadFrom == 'TROPICOS' & records$year == 9912 & records$collectedBy == 'P.J. Harmon & Jeff Hajenga'] <- 2012
		
		# records$year[grepl(records$coordRemarks, pattern='RECORDED_DATE_INVALID')] <- NA
	
	# say('   ### other')
	# say('   #########')
	
		# records$cultivated <- as.logical(records$cultivated)
		# records$cultivated[grepl(records$notes, pattern='cultivated')] <- TRUE
		# records$cultivated[grepl(records$notes, pattern='garden')] <- TRUE
		# records$cultivated[records$locality == 'cult in Horto Upsaliensi'] <- TRUE
		# records$cultivated[records$locality == 'cult in Horto Hafniensi'] <- TRUE
		# records$cultivated[records$locality == 'Ex horto Hamburgensi [=presumably cultivated Botanical Gardens Hamburg].'] <- TRUE
		# records$cultivated[records$locality == 'Kiliae [=Presumably cultivated Botanical Gardens University of Kiel].'] <- TRUE
		# records$cultivated[records$locality == 'Ex horto meo. [Presumably cultivated in Europe].'] <- TRUE
		# records$cultivated[records$locality == 'Hort. bot. Halens [=presumably cultivated Botanical Gardens Halle].'] <- TRUE
		# records$cultivated[records$locality == 'Hort bot. Gryphic. [=presumably cultivated Greifswald Botanical Garden].'] <- TRUE
		# records$cultivated[records$locality == '[Sweden, Stockholm, Bergius Garden?]'] <- TRUE
		# records$cultivated[records$locality == 'MontrÃ©al, Mont Royal, jardin'] <- TRUE
		# records$cultivated[records$locality == 'Botanic Garden. | Cambridge, Mass.'] <- TRUE

	# say('   ### coordinates')
	# say('   ###############')
		
		# fill <- data.frame(coordApprox=rep(NA, nrow(records)), longitude=records$origLong, latitude=records$origLat)
		# records <- insertCol(fill, into=records, at='origLong')

		# tropicosRows <- which(records$downloadFrom == 'TROPICOS')

		# tropicosCoords <- convertTropicosCoords(
			# long=records$longitude[tropicosRows],
			# lat=records$latitude[tropicosRows],
			# returnApprox=TRUE
		# )
		
		# records$coordApprox[tropicosRows] <- tropicosCoords$approximate
		# records$longitude[tropicosRows] <- tropicosCoords$longitude
		# records$latitude[tropicosRows] <- tropicosCoords$latitude
		
		# records$longitude <- as.numeric(records$longitude)
		# records$latitude <- as.numeric(records$latitude)

	# say('   ### manual corrections')
	# say('   ######################')
		
		# # case 1
		# coords <- convertTropicosCoords("087°32′53'W", "37°52′56'N")
		# this <- which(records$locality == 'Audubon State Park' & records$admin2 == 'Henderson')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE
		
		# # case 2
		# coords <- convertTropicosCoords("087°32′53'W", "37°52′56'N")
		# this <- which(records$locality == 'Audubon State Park, Henderson.; John James Audubon State Park' & records$admin2 == 'Henderson')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE
		
		# # case 3
		# coords <- convertTropicosCoords("091°50′37'W", "36°38′45'N")
		# this <- which(records$locality == 'N-facing slope; 10 mi. S of West Plains heavy woods' & records$admin2 == 'Howell')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE
		
		# # case 4
		# coords <- convertTropicosCoords("078°52′24'W", "36°03′52'N")
		# this <- which(records$locality == 'North Carolina. Iredell soil, 3.5 miles north of Braggtown and 0.5 mile north of Eno river on co. rd. The Southeastern United States.')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE
		
		# # case 5
		# coords <- convertTropicosCoords("083°05′50'W", "35°28′19'N")
		# this <- which(records$locality == 'Crest of Oldfield Top, Waynesville')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE
		
		# # case 6
		# coords <- convertTropicosCoords("078°52′24'W", "36°02′55'N")
		# this <- which(records$locality == '1.7 miles North of Weaver (just North of Durham)')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE

		# # case 7
		# coords <- convertTropicosCoords("082°06′37'W", "36°06′13'N")
		# this <- which(records$locality == 'Tennessee-North Carolina border: Roan Mounatins')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE
		
		# # case 8
		# coords <- convertTropicosCoords("080°57′36'W", "39°54′20'N")
		# this <- which(records$locality == 'Captina Creek, Belmont County')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE

		# # case 9a
		# this <- which(records$locality == 'topeka' & records$collectedBy == 'c.p.t.')
		# records$longitude[this] <- -95.6752
		# records$latitude[this] <- 39.0473
		# records$coordApprox[this] <- TRUE
		
		# # case 9b
		# this <- which(records$locality == 'topeka' & records$collectedBy == 'scott, w.w.')
		# records$longitude[this] <- -95.6752
		# records$latitude[this] <- 39.0473
		# records$coordApprox[this] <- TRUE
		
		# # case 10
		# this <- which(records$locality == 'Statesville' & records$collectedBy == 'M. E. Hyams')
		# records$longitude[this] <- -80.8873
		# records$latitude[this] <- 35.7826
		# records$coordApprox[this] <- TRUE
		
		# # case 11
		# coords <- convertTropicosCoords("080°20′06'W", "37°10′18'N")
		# this <- which(records$locality == '4 mi. s.e. Ellett on shaley power line clearing')
		# records$longitude[this] <- coords$longitude
		# records$latitude[this] <- coords$latitude
		# records$coordApprox[this] <- TRUE
		
		# # case 12
		# this <- which(records$locality == 'Charlotte' & records$collectedBy == 'Collector unspecified')
		# records$longitude[this] <- -80.8431
		# records$latitude[this] <- 35.2271
		# records$coordApprox[this] <- TRUE
		
		# # case 13
		# this <- which(records$locality == 'Kirksville, Mo.' & records$collectedBy == 'Fisher GL')
		# records$longitude[this] <- -92.5832
		# records$latitude[this] <- 40.1948
		# records$coordApprox[this] <- TRUE
		
		# # case 14
		# this <- which(records$locality == 'Potrero' & records$collectedBy == 'A. Kellogg, W. G. W. Harford')
		# records$longitude[this] <- -116.6131
		# records$latitude[this] <- 32.6048
		# records$coordApprox[this] <- TRUE
		
		# # case 15
		# this <- which(records$locality == 'Junucta below Alexandria' & records$admin3 == 'Alexandria')
		# records$longitude[this] <- -78.0978
		# records$latitude[this] <- 40.5565
		# records$coordApprox[this] <- TRUE
		
		# # case 15a
		# this <- which(records$admin3 == 'Brodwell' & records$source == 'MPM')
		# records$admin3[this] <- 'Broadwell'
		# records$longitude[this] <- -89.4432
		# records$latitude[this] <- 40.0681
		# records$coordApprox[this] <- TRUE
		
		# # case 15b
		# this <- which(records$coordRemarks == 'precise location information available on request; location information withheld below municipality level; ; NA; ; ; ; ; ' & records$admin2 == 'Logan' & records$admin3 == 'Brodwell')
		# records$admin3[this] <- 'Broadwell'
		# records$longitude[this] <- -89.4432
		# records$latitude[this] <- 40.0681
		# records$coordApprox[this] <- TRUE

		# # case 16a
		# this <- which(records$coordRemarks == 'precise location information available on request; location information withheld below municipality level; ; NA; NA; ; ; NA; NA; RECORDED_DATE_INVALID' & records$admin2 == 'St. Clair' & records$admin3 == 'Belleville')
		# records$longitude[this] <- -89.9840
		# records$latitude[this] <- 38.5201
		# records$coordApprox[this] <- TRUE
		
		# # case 16b
		# this <- which(records$coordRemarks == 'precise location information available on request; location information withheld below municipality level; ; NA; ; ; ; ; ' & records$admin2 == 'St. Clair' & records$admin3 == 'Belleville')
		# records$longitude[this] <- -89.9840
		# records$latitude[this] <- 38.5201
		# records$coordApprox[this] <- TRUE
		
		# # case 17
		# this <- which(records$admin3 == 'Independence' & records$collectedBy == 'Frank Bush')
		# records$admin2[this] <- 'Jackson'
		# records$longitude[this] <- -94.4155
		# records$latitude[this] <- 39.0911
		# records$coordApprox[this] <- TRUE

		# records$coordApprox[records$georefIssuesBien == 1] <- TRUE
		# records$coordApprox[records$georefIssuesGbif] <- TRUE
		# records$coordApprox[records$origCoordUncer_m > 1000] <- TRUE
		# records$coordApprox[grepl(records$coordRemarks, pattern='COORDINATE_ROUNDED')] <- TRUE
		# records$coordApprox[grepl(records$coordRemarks, pattern='precise location information available on request')] <- TRUE
		# records$coordApprox[grepl(records$coordRemarks, pattern='; ; ; NA; michael smith (2017-09-26 07:20:19); ; georef batch tool 2017-09-26; GeoLocate; reviewed - high confidence; ')] <- FALSE

		# records$admin0[records$admin0 == 'U.S.A.'] <- 'United States'
		# records$admin0[records$admin0 == 'United States of America'] <- 'United States'
		# records$admin0[records$admin0 == 'US'] <- 'United States'
		# records$admin0[records$admin0 == 'usa'] <- 'United States'
		# records$admin0[records$admin0 == 'USA'] <- 'United States'
		
	# say('   ### assign county and coordinates based on county centroid if no coordinates')
	# say('   ############################################################################')
	
		# usa <- shapefile('D:/Ecology/Political Geography/GADM/ver2pt8/WGS84/USA_adm2')
		# usaAlbers <- sp::spTransform(usa, getCRS('albersNA', TRUE))
		
		# for (i in 1:nrow(records)) {
		
			# # get county for records with coordinates
			# if (!is.na(records$longitude[i])) {
			
				# gadmLoc <- extract(usa, SpatialPoints(records[i, longLat], proj4string=getCRS('wgs84', TRUE)))
				# records$admin0[i] <- gadmLoc$NAME_0
				# records$admin1[i] <- gadmLoc$NAME_1
				# records$admin2[i] <- gadmLoc$NAME_2
				# if (!is.na(records$origCoordUncer_m[i]) && records$origCoordUncer_m[i] <= 1000) records$coordApprox[i] <- FALSE
			
			# # assign county centroid if no coordinates but county is available
			# } else if (!is.na(records$admin1[i]) & !is.na(records$admin2[i]) & records$admin1[i] != '' & records$admin2[i] != '') {
			
				# if (!is.na(records$admin2[i])) {
				
					# if (grepl(records$admin2[i], pattern='St. ')) records$admin2[i] <- gsub(records$admin2[i], pattern='St. ', replacement='Saint ')
					# if (grepl(records$admin2[i], pattern='Ste. ')) records$admin2[i] <- gsub(records$admin2[i], pattern='Ste. ', replacement='Saint ')
					# if (grepl(records$admin2[i], pattern=' County')) records$admin2[i] <- gsub(records$admin2[i], pattern=' County', replacement='')
					# if (grepl(records$admin2[i], pattern=' Co.')) records$admin2[i] <- gsub(records$admin2[i], pattern=' Co.', replacement='')
					# if (grepl(records$admin2[i], pattern='hors MRC')) records$admin2[i] <- 'Montreal'
					
				# }
				
				# thisCo <- usaAlbers[usaAlbers$NAME_1 == records$admin1[i] & usaAlbers$NAME_2 == records$admin2[i], ]
				# cent <- gCentroid(thisCo)
				# cent <- sp::spTransform(cent, getCRS('wgs84', TRUE))
				# records$longitude[i] <- coordinates(cent)[ , 1]
				# records$latitude[i] <- coordinates(cent)[ , 2]
				# records$coordApprox[i] <- TRUE
				
			# }
		
		# }
	
		# records$admin0[records$admin0 == ''] <- NA
		# records$admin1[records$admin1 == ''] <- NA
		# records$admin2[records$admin2 == ''] <- NA
		# records$admin3[records$admin3 == ''] <- NA

		# # # inspect
		# # write.csv(records, 'C:/ecology/!Scratch/Records with Auto and Manually Cleaned Coordinates.csv', row.names=FALSE)
	
	# say('   ### remove erroneous records')
	# say('   ############################')
	
		# if (any(records$admin1 == 'California')) records <- records[records$admin1 != 'California', ]
	
	# say('   ### remove records that are cultivated, have no date, have no county, unverifiable record basis')
	# say('   ###############################################################################################')

		# records <- records[(!records$cultivated | is.na(records$cultivated)), ]

		# records <- records[!is.na(records$admin2), ]
		# records <- records[records$basis %in% c('living specimen', 'LIVING_SPECIMEN', 'Observation (NPS)', 'plot', 'Preserved Specimen', 'PRESERVED_SPECIMEN', 'preservedspecimen', 'PreservedSpecimen', 'specimen', 'Specimen', 'Vouchered specimen'), ]
		
	# say('   ### remove records that are redundant')
	# say('   #####################################')
	
		# say('      I am defining "redundant" records as those collected from the same county in the same year.')
		
		# nonRed <- data.frame()
		
		# # records
		# yearStateCounty <- paste(records$year, records$admin1, records$admin2)
		# yearStateCountyUni <- sort(unique(yearStateCounty))
	
		# # for each combination of collection year, state, and county
		# for (this in yearStateCountyUni) {
		
			# these <- records[this == yearStateCounty, ]
			
			# # prioritize
			# if (nrow(these) > 1) {
			
				# these <- rbind(
					# these[!these$coordApprox & !is.na(these$coordApprox), ],
					# these[these$coordApprox & !is.na(these$coordApprox), ],
					# these[is.na(these$coordApprox), ]
				# )

				# these <- rbind(
					# these[these$downloadFrom == 'Christy Edwards', ],
					# these[these$downloadFrom != 'Christy Edwards', ]
				# )
				
				# these <- rbind(
					# these[!is.na(these$year), ],
					# these[is.na(these$year), ]
				# )

				# these <- these[1, ]
				
			# }
			
			# nonRed <- rbind(nonRed, these)
		
		# }
		
		# write.csv(nonRed, './Data/01 Delphinium exaltatum - Cleaned & Non-Redundant (Date & County) Occurrences.csv', row.names=FALSE)

# say('########################################################################################################')
# say('### 02 incorporate presences into county shapefile, assign population, and define accessible regions ###')
# say('########################################################################################################')
	
	# records <- read.csv('./Data/01 Delphinium exaltatum - Cleaned & Non-Redundant (Date & County) Occurrences.csv')

	# delEx <- rbind(usaAdmin2, canadaAdmin2, mexicoAdmin2)
	
	# delEx <- delEx[ , which(names(delEx) %in% c('NAME_1', 'NAME_2'))]
	# names(delEx) <- c('admin1', 'admin2')
	
	# # remove "lake" delEx
	# delEx <- delEx[!(delEx$admin2 %in% c('Lake Superior', 'Lake Michigan', 'Lake Hurron', 'Lake Erie', 'Lake Ontario')), ]
	# delExEa <- sp::spTransform(delEx, getCRS('albersNA', TRUE))
	
	# # tally presences
	# delEx$delExObs <- 0
	# delEx$deNoYear <- 0
	
	# say('   ### collate records in counties')
	# say('   ###############################')

		# years <- sort(unique(records$year))
	
		# # add records with collection years to delEx
		# for (year in years) {

			# delEx$DUMMY <- 0
			# names(delEx)[ncol(delEx)] <- paste0('deYr', year)
			
			# theseRecords <- records[records$year == year & !is.na(records$year), ]
			
			# if (nrow(theseRecords) > 0) {
			
				# for (i in 1:nrow(theseRecords)) {
			
					# delEx@data[delEx$admin1 == theseRecords$admin1[i] & delEx$admin2 == theseRecords$admin2[i], which(names(delEx) == paste0('deYr', year))] <- 1
					
				# }
					
			# }
			
		# }
		
		# # collate number of presences in each county (temporary)
		# delEx$delExObs <- as.integer(rowSums(delEx@data[ , which(names(delEx) %in% paste0('deYr', years))]) > 0)

		# # add records without years to delEx
		# theseRecords <- records[is.na(records$year), ]
		
		# if (nrow(theseRecords) > 0) {
		
			# for (i in 1:nrow(theseRecords)) {
		
				# this <- which(delEx$admin1 == theseRecords$admin1[i] & delEx$admin2 == theseRecords$admin2[i] & delEx$delExObs == 0)
				# if (length(this) > 0) delEx@data$deNoYear[this] <- 1
				
			# }
				
		# }

		# # collate number of presences in each county (final value)
		# delEx$delExObs <- as.integer(rowSums(delEx@data[ , c(which(names(delEx) %in% paste0('deYr', years)), which(names(delEx) == 'deNoYear'))]) > 0)

	# say('   ### define accessible/projection area for species')
	# say('   #################################################')

		# delEx$areaKm2 <- gArea(sp::spTransform(delEx, getCRS('albersNA', TRUE)), byid=TRUE) / 1000000
	
		# delEx$areaKm2_popSpecies_extentNarrow <- NA
		# delEx$areaKm2_popSpecies_extentBroad <- NA
	
		# delEx$accessible_popSpecies_extentNarrow <- FALSE
		# delEx$accessible_popSpecies_extentBroad <- FALSE
	
		# pres <- delEx[which(delEx$delExObs == 1), ]
		# presEa <- sp::spTransform(pres, getCRS('albersNA', TRUE))

		# # buffers (near accessible, far accessible, projection region)
		# buffNarrow <- gBuffer(presEa, width=1000 * narrow)
		# buffBroad <- gBuffer(presEa, width=1000 * broad)
		# buffEntire <- gBuffer(presEa, width=1000 * entire)
		
		# buffNarrow <- sp::spTransform(buffNarrow, getCRS('wgs84', TRUE))
		# buffBroad <- sp::spTransform(buffBroad, getCRS('wgs84', TRUE))
		# buffEntire <- sp::spTransform(buffEntire, getCRS('wgs84', TRUE))

		# # flag counties in each buffer region
		# buffNarrow <- crop(delEx, buffNarrow)
		# buffBroad <- crop(delEx, buffBroad)
		# buffEntire <- crop(delEx, buffEntire)
		
		# delEx <- crop(delEx, extent(buffEntire))

		# buffNarrowAdmin <- paste(buffNarrow$admin1, buffNarrow$admin2)
		# buffBroadAdmin <- paste(buffBroad$admin1, buffBroad$admin2)
		# delExAdmin <- paste(delEx$admin1, delEx$admin2)
		
		# delEx$accessible_popSpecies_extentNarrow[delExAdmin %in% buffNarrowAdmin] <- TRUE
		# delEx$accessible_popSpecies_extentBroad[delExAdmin %in% buffBroadAdmin] <- TRUE

		# # calculate county area (may be cropped by buffer)
		# buffNarrow <- sp::spTransform(buffNarrow, getCRS('albersNA', TRUE))
		# buffBroad <- sp::spTransform(buffBroad, getCRS('albersNA', TRUE))
		
		# buffNarrow$area_km2 <- gArea(buffNarrow, byid=TRUE) / 1000000
		# buffBroad$area_km2 <- gArea(buffBroad, byid=TRUE) / 1000000

		# for (this in buffNarrowAdmin) {
			# delEx$areaKm2_popSpecies_extentNarrow[delExAdmin == this] <- buffNarrow$area_km2[buffNarrowAdmin == this]
		# }
			
		# for (this in buffBroadAdmin) {
			# delEx$areaKm2_popSpecies_extentBroad[delExAdmin == this] <- buffBroad$area_km2[buffBroadAdmin == this]
		# }

	# say('   ### plot accessible area for species')
	# say('   ####################################')
	
		# delExEa <- sp::spTransform(delEx, getCRS('albersNA', TRUE))
	
		# if (makePlot) {
		
			# png('./Analysis/02a Delphinium exaltatum - Accessible Region for Species.png', width=2400, height=1800, res=300)
			
				# par(mar=mar)
			
				# plot(delExEa[which(delExEa$accessible_popSpecies_extentBroad), ], border='gray', main='Delphinium exaltatum Accessible Region for Species', cex.main=cexMain, col='gray80')
				# plot(delExEa[delExEa$accessible_popSpecies_extentNarrow, ], border='gray', col='gray40', add=TRUE)

				# plot(usaAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)
				# plot(canadaAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)
				# plot(mexicoAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)

				# plot(usaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
				# plot(canadaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
				# plot(mexicoAdmin1AlbersNA, lwd=lwdStates, add=TRUE)

				# plot(delExEa[delExEa$delExObs == 1, ], col=speciesCol, add=TRUE)

				# legend('bottomright', inset=c(0, 0.05), legend=c('Presence', paste0('Narrow (', narrow, ' km)'), paste0('Broad (', broad, ' km)')), fill=c(speciesCol, 'gray40', 'gray80'), cex=cexLegend)
				
				# title(sub=date(), cex.sub=cexSub, line=0)
				
			# dev.off()
			
		# }

	# say('   ### assign populations and extents')
	# say('   ##################################')

		# delEx$population_biasLowland <- rep(NA, nrow(delEx))
		# delEx$population_biasHighland <- rep(NA, nrow(delEx))
		
		# # uncontested lowland
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Alabama'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Kansas'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Missouri'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Illinois'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Wisconsin'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Kentucky'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Tennessee'] <- 'lowland'
		
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Alabama'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Kansas'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Missouri'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Illinois'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Wisconsin'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Kentucky'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Tennessee'] <- 'lowland'
		
		# # uncontested highland
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Maryland'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Virginia'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'West Virginia'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Ohio'] <- 'lowland'
		
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Maryland'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Virginia'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'West Virginia'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Ohio'] <- 'lowland'

		# # uncontested mixed
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Alleghany'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Ashe'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Watauga'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Mitchell'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'McDowell'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Haywood'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Jackson'] <- 'highland'

		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Alleghany'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Ashe'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Watauga'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Mitchell'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'McDowell'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Haywood'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'North Carolina' & delEx$admin2 == 'Jackson'] <- 'highland'

		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Massachusetts' & delEx$admin2 == 'Middlesex'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Massachusetts' & delEx$admin2 == 'Franklin'] <- 'highland'

		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Massachusetts' & delEx$admin2 == 'Middlesex'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Massachusetts' & delEx$admin2 == 'Franklin'] <- 'highland'

		# # Pennsylvania (contested)
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Westmoreland'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Fayette'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Schuylkill'] <- 'lowland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Huntingdon'] <- 'highland'
		# delEx$population_biasLowland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Bedford'] <- 'highland'

		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania'] <- 'lowland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Westmoreland'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Fayette'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Schuylkill'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Huntingdon'] <- 'highland'
		# delEx$population_biasHighland[delEx$delExObs == 1 & delEx$admin1 == 'Pennsylvania' & delEx$admin2 == 'Bedford'] <- 'highland'

		# # summarize
		# say('      LOWLAND-BIAS population assignation:', pre=1)
		# say('      There are ', sum(delEx$population_biasLowland == 'lowland', na.rm=TRUE), ' presences assigned to the LOWLAND population.')
		# say('      There are ', sum(delEx$population_biasLowland == 'highland', na.rm=TRUE), ' presences assigned to the HIGHLAND population.', post=2)
		
		# say('      HIGHLAND-BIAS population assignation:')
		# say('      There are ', sum(delEx$population_biasHighland == 'lowland', na.rm=TRUE), ' presences assigned to the LOWLAND population.')
		# say('      There are ', sum(delEx$population_biasHighland == 'highland', na.rm=TRUE), ' presences assigned to the HIGHLAND population.', post=2)

		# # plot
		# if (makePlot) {
			
			# delExEa <- sp::spTransform(delEx, getCRS('albersNA', TRUE))

			# for (bias in c('Lowland', 'Highland')) {
			
				# png(paste0('./Analysis/02b Delphinium exaltatum - Population Assignation with a Bias toward ', bias, '.png'), width=2400, height=1800, res=300)
				
					# par(mar=mar)

					# plot(delExEa[which(delExEa$accessible_popSpecies_extentBroad), ], border=NA, main=paste0('Delphinium exaltatum Population Assignation with a ', bias, ' Bias'), cex.main=cexMain)

					# plot(delExEa[which(delExEa@data[ , paste0('population_bias', bias)] == 'lowland'), ], col=lowlandCol, add=TRUE)
					# plot(delExEa[which(delExEa@data[ , paste0('population_bias', bias)] == 'highland'), ], col=highlandCol, add=TRUE)

					# plot(usaAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)
					# plot(canadaAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)
					# plot(mexicoAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)

					# plot(usaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
					# plot(canadaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
					# plot(mexicoAdmin1AlbersNA, lwd=lwdStates, add=TRUE)

					# legend('bottomright', inset=c(0, 0.05), legend=c('Lowland', 'Highland'), fill=c(lowlandCol, highlandCol), cex=cexLegend)
					
					# title(sub=date(), cex.sub=cexSub, line=0)
					
				# dev.off()
				
			# }
			
		# }

		# ### define accessible regions for populations
		# #############################################
		
		# delEx$accessible_popLowland_extentNarrow_biasLowland <- FALSE
		# delEx$accessible_popHighland_extentNarrow_biasLowland <- FALSE
		
		# delEx$accessible_popLowland_extentNarrow_biasHighland <- FALSE
		# delEx$accessible_popHighland_extentNarrow_biasHighland <- FALSE
		
		# delEx$accessible_popLowland_extentBroad_biasLowland <- FALSE
		# delEx$accessible_popHighland_extentBroad_biasLowland <- FALSE
		
		# delEx$accessible_popLowland_extentBroad_biasHighland <- FALSE
		# delEx$accessible_popHighland_extentBroad_biasHighland <- FALSE

		# delEx$areaKm2_popLowland_extentNarrow_biasLowland <- NA
		# delEx$areaKm2_popHighland_extentNarrow_biasLowland <- NA
		
		# delEx$areaKm2_popLowland_extentNarrow_biasHighland <- NA
		# delEx$areaKm2_popHighland_extentNarrow_biasHighland <- NA
		
		# delEx$areaKm2_popLowland_extentBroad_biasLowland <- NA
		# delEx$areaKm2_popHighland_extentBroad_biasLowland <- NA
		
		# delEx$areaKm2_popLowland_extentBroad_biasHighland <- NA
		# delEx$areaKm2_popHighland_extentBroad_biasHighland <- NA

		# for (pop in c('lowland', 'highland')) {
		
			# for (bias in c('Lowland', 'Highland')) {

				# say('      Accessible region for population: ', pop, ' with ', tolower(bias), ' bias')
			
				# pres <- delEx[which(delEx$delExObs == 1 & delEx@data[ , paste0('population_bias', bias)] == pop), ]
				# presEa <- sp::spTransform(pres, getCRS('albersNA', TRUE))

				# # buffers
				# buffNarrow <- gBuffer(presEa, width=1000 * narrow)
				# buffBroad <- gBuffer(presEa, width=1000 * broad)
				
				# buffNarrow <- sp::spTransform(buffNarrow, getCRS('wgs84', TRUE))
				# buffBroad <- sp::spTransform(buffBroad, getCRS('wgs84', TRUE))

				# # flag counties in each buffer region
				# buffNarrow <- crop(delEx, buffNarrow)
				# buffBroad <- crop(delEx, buffBroad)
				
				# buffNarrowAdmin <- paste(buffNarrow$admin1, buffNarrow$admin2)
				# buffBroadAdmin <- paste(buffBroad$admin1, buffBroad$admin2)
				# delExAdmin <- paste(delEx$admin1, delEx$admin2)

				# delEx@data[delExAdmin %in% buffNarrowAdmin, paste0('accessible_pop', capIt(pop), '_extentNarrow_bias', bias)] <- TRUE
				# delEx@data[delExAdmin %in% buffBroadAdmin, paste0('accessible_pop', capIt(pop), '_extentBroad_bias', bias)] <- TRUE
				
				# # calculate county area (may be cropped by buffer)
				# buffNarrow <- sp::spTransform(buffNarrow, getCRS('albersNA', TRUE))
				# buffBroad <- sp::spTransform(buffBroad, getCRS('albersNA', TRUE))
				
				# buffNarrow$area_km2 <- gArea(buffNarrow, byid=TRUE) / 1000000
				# buffBroad$area_km2 <- gArea(buffBroad, byid=TRUE) / 1000000

				# for (this in buffNarrowAdmin) {
					# delEx@data[delExAdmin == this, paste0('areaKm2_pop', capIt(pop), '_extentNarrow_bias', bias)] <- buffNarrow$area_km2[buffNarrowAdmin == this]
				# }

				# for (this in buffBroadAdmin) {
					# delEx@data[delExAdmin == this, paste0('areaKm2_pop', capIt(pop), '_extentBroad_bias', bias)] <- buffBroad$area_km2[buffBroadAdmin == this]
				# }

				# png(paste0('./Analysis/02c Delphinium exaltatum - Accessible Region for ', capIt(pop), ' with a Bias toward ', bias, '.png'), width=2400, height=1800, res=300)

					# par(mar=mar)
				
					# delExEa <- sp::spTransform(delEx, getCRS('albersNA', TRUE))
					# delExEaBroad <- delExEa[which(delExEa@data[ , paste0('accessible_pop', capIt(pop), '_extentBroad_bias', bias), ]), ]
					# delExEaNarrow <- delExEa[which(delExEa@data[ , paste0('accessible_pop', capIt(pop), '_extentNarrow_bias', bias), ]), ]
				
					# plot(delExEa[delExEa@data[['accessible_popSpecies_extentBroad']], ], border=NA, main=paste0('Delphinium exaltatum ', capIt(pop), ' Population\'s Accessible Region\n(Assigned with a Bias toward ', bias, ')'), cex.main=cexMain)
					# plot(delExEaBroad, border='gray', col='gray80', lwd=lwdCounties, add=TRUE)
					# plot(delExEaNarrow, border='gray', col='gray40', lwd=lwdCounties, add=TRUE)

					# plot(usaAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)
					# plot(canadaAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)
					# plot(mexicoAdmin2AlbersNA, lwd=lwdCounties, border='gray', add=TRUE)

					# plot(usaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
					# plot(canadaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
					# plot(mexicoAdmin1AlbersNA, lwd=lwdStates, add=TRUE)

					# plot(delExEa[which(delExEa$delExObs == 1 & delExEa@data[ , paste0('population_bias', bias)] == pop), ], lwd=lwdCounties, col=ifelse(pop == 'lowland', lowlandCol, highlandCol), add=TRUE)

					# legend('bottomright', inset=c(0, 0.05), legend=c('Presence', paste0('Narrow (', narrow, ' km)'), paste0('Broad (', broad, ' km)')), fill=c(ifelse(pop == 'lowland', lowlandCol, highlandCol), 'gray40', 'gray80'), cex=cexLegend)
					
					# title(sub=date(), cex.sub=cexSub, line=0)
						
				# dev.off()
				
			# } # next bias

		# } # next population
	
	# save(delEx, file='./Data/02 Delphinium exaltatum - Assigned Populations and Accessible Regions - Spatial Polygons.RData')

# say('####################################')
# say('### calculate common climate PCA ###')
# say('####################################')

	# say('I am calculating a common PCA across present (1950-2005) and past (21 Kybp) climates to minimize the effect of changing PCA rotation when projecting to past climates.')

	# load('./Data/02 Delphinium exaltatum - Assigned Populations and Accessible Regions - Spatial Polygons.RData')
	# delExBroad <- delEx[delEx$accessible_popSpecies_extentBroad, ]
	
	# say('Extracting climate for present.')
	# sqClim <- stack(paste0('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/!ENSEMBLE 1950-2005 across 12 ESMs/', vars, '.tif'))
	# sqClim <- crop(sqClim, delEx)
	
	# delExMask <- rasterize(delExBroad, sqClim[[1]])
	# delExMask <- delExMask * 0 + 1
	# sqClim <- sqClim * delExMask
	
	# sites <- randomPoints(sqClim, 11000)
	
	# env <- extract(sqClim, sites)
	# colnames(env) <- vars

	# nas <- naRows(env)
	# if (length(nas) > 0) env <- env[-nas, ]

	# pca <- princomp(env, cor=TRUE)

	# # png('./Analysis/PCA on Current (1950-2005) Climate.png', width=1000, height=1000)
	
		# # par(mar=c(5, 4, 4, 2) + 1)
	
		# # pc1var <- pca$sdev[1]^2 / sum(pca$sdev^2)
		# # pc2var <- pca$sdev[2]^2 / sum(pca$sdev^2)
		# # plot(pca$scores[ , 1:2], xlab=paste0('PC 1 (', round(100 * pc1var, 1), '%)'), ylab=paste0('PC 2 (', 100 * round(pc2var, 1), '%)'), pch=16, cex=1.6, cex.lab=1.8, cex.axis=1.8, col='darkorange')
		
		# # ex <- 4
		
		# # for (thisVar in vars) {
		
			# # arrows(x0=0, y0=0, x1=ex * pca$loadings[thisVar, 1], y1=ex * pca$loadings[thisVar, 2], col='black', lwd=2)
			# # text(x=ex * pca$loadings[thisVar, 1], y=ex * pca$loadings[thisVar, 2], labels=thisVar, xpd=NA, col='black', cex=1.8, pos=2)
			
		# # }
		
		# # title(sub=date(), cex.sub=0.8)
		
	# # dev.off()
	
	# save(pca, file='./Data/PCA across Present (1950-2005) from Broad Accessible Region.RData')

# say('#############################################')
# say('### crop states/provinces to study region ###')
# say('#############################################')
	
	# say('This is for making visually-appealing maps later.')
	
	# dirCreate('./Other')
	
	# load('./Data/02 Delphinium exaltatum - Assigned Populations and Accessible Regions - Spatial Polygons.RData')
	# delExEa <- sp::spTransform(delEx, getCRS('albersNA', TRUE))
	
	# canadaAdmin1AlbersNACropped <- crop(canadaAdmin1AlbersNA, delExEa)
	# usaAdmin1AlbersNACropped <- crop(usaAdmin1AlbersNA, delExEa)
	# mexicoAdmin1AlbersNACropped <- crop(mexicoAdmin1AlbersNA, delExEa)
	
	# save(canadaAdmin1AlbersNACropped, file='./Other/Canada Admin Level 1 Albers Croppped to Study Region.Rdata')
	# save(usaAdmin1AlbersNACropped, file='./Other/USA Admin Level 1 Albers Croppped to Study Region.Rdata')
	# save(mexicoAdmin1AlbersNACropped, file='./Other/Mexico Admin Level 1 Albers Croppped to Study Region.Rdata')
	
# say('#########################################################')
# say('### 03 match current environmental data with counties ###')
# say('#########################################################')

	# say('Extrapolating climate and elevation data to counties using mean of cells across counties weighted for area of cell covered by the county. In a few cases after extraction counties have NAs so I am assigning these equal to the values of the closest county without NAs.')

	# load('./Data/02 Delphinium exaltatum - Assigned Populations and Accessible Regions - Spatial Polygons.RData')

	# cents <- gCentroid(delEx)
	# pairDist <- pointDist(delEx)
	# diag(pairDist) <- NA
	
	# ### extract climate
	# ###################
	
		# climate <- stack(paste0('D:/ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/!ENSEMBLE 1950-2005 across 12 ESMs/', vars, '.tif'))
		# climate <- crop(climate, delEx)

		# extracted <- raster::extract(climate, delEx, weights=TRUE, normalizeWeight=TRUE)
		
		# clim <- matrix(NA, nrow=nrow(delEx), ncol=length(vars))
		# clim <- as.data.frame(clim)
		# names(clim) <- vars
		
		# # collate
		# for (i in 1:nrow(clim)) {
		
			# countyClim <- extracted[[i]][ , vars]
			# if (class(countyClim) == 'numeric') {
			
				# clim[i, ] <- matrix(countyClim, nrow=1, ncol=length(vars))
				
			# } else {
		
				# wgts <- extracted[[i]][ , 'weight']
				# wgts <- matrix(wgts, nrow=nrow(countyClim), ncol=ncol(countyClim))
			
				# countyClim <- countyClim * wgts
				# clim[i, ] <- colSums(countyClim)
				
			# }
			
		# }

		# # assign NAs
		# for (thisVar in vars) {
		
			# nas <- which(is.na(clim[ , thisVar]))
				
			# if (length(nas) > 0) {

				# for (thisNa in nas) {
				
					# thisPairDist <- pairDist[thisNa, ]
					
					# while (is.na(clim[thisNa, thisVar])) {
					
						# newClim <- which.min(thisPairDist)
						# clim[thisNa, thisVar] <- clim[newClim, thisVar]
						# thisPairDist[newClim] <- NA
					
					# }
				
				# }


			# }
			
		# }
			
	# ### extract elevation
	# #####################
		
		# elevation <- raster('D:/Ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/10 arcmin/Elevation - 10 arcmin/elevation.tif')
		# elevation <- crop(elevation, delEx)

		# extracted <- raster::extract(elevation, delEx, weights=TRUE, normalizeWeight=TRUE)
		
		# elev <- rep(NA, nrow(delEx))
		
		# # collate
		# for (i in seq_along(elev)) {
		
			# countyElev <- extracted[[i]][ , c('weight', 'value')]
			# countyElev <- if (class(countyElev) == 'numeric') {
				# prod(countyElev)
			# } else {
				# countyElev[ , 'weight'] * countyElev[ , 'value']
			# }
			# elev[i] <- sum(countyElev)
				
		# }
		
		# # assign NAs
		# nas <- which(is.na(elev))
			
		# if (length(nas) > 0) {

			# for (thisNa in nas) {
			
				# thisPairDist <- pairDist[thisNa, ]
				
				# while (is.na(elev[thisNa])) {
				
					# newClim <- which.min(thisPairDist)
					# elev[thisNa] <- elev[newClim]
					# thisPairDist[newClim] <- NA
				
				# }
			
			# }


		# }
		
		# elev <- data.frame(elevation_m=elev)
		
	# ### add PCA values
	# ##################
		
		# load('./Data/PCA across Present (1950-2005) from Broad Accessible Region.RData')
		# pcaPred <- predict(pca, clim)
		# colnames(pcaPred) <- paste0('pc', 1:ncol(pcaPred))
		
	# ### remember
	# ############
	
	# delEx@data <- cbind(delEx@data, elev, clim, pcaPred)

	# ### plot
	# ########

	# # PCs
	# for (pc in 1:2) {
		
		# png(paste0('./Analysis/03 Delphinium exaltatum - Climate - PC ', pc, '.png'), width=2400, height=1800, res=300)
			# plotNice(v=delEx@data[ , paste0('pc', pc)], title=paste('Principal Component', pc), leg='PC', stretch=TRUE, color='steelblue4')
		# dev.off()
		
	# }

	# # elevation
	# png(paste0('./Analysis/03 Delphinium exaltatum - Elevation.png'), width=2400, height=1800, res=300)
		# thisElev <- delEx$elevation_m
		# thisElev[!delEx$accessible_popSpecies_extentBroad] <- 0
		# thisElev <- stretchMinMax(thisElev)
		# plotNice(v=thisElev, title='Elevation', leg='Elevation (m)', color='black')
	# dev.off()

	# # PCA
	# ex <- 4.4
	
	# png(paste0('./Analysis/03 Delphinium exaltatum - PCA.png'), width=1200, height=1200, res=300)
	
		# par(mar=mar + 2)
		# pc1var <- pca$sdev[1]^2 / sum(pca$sdev^2)
		# pc2var <- pca$sdev[2]^2 / sum(pca$sdev^2)
		# pcScores <- predict(pca, delEx@data)
		# smoothScatter(pcScores[delEx$accessible_popSpecies_extentBroad, 1:2], xlab=paste0('PC 1 (', round(100 * pc1var, 1), '%)'), ylab=paste0('PC 2 (', 100 * round(pc2var, 1), '%)'), main='PCA', nrpoints=0)

		# for (i in seq_along(vars)) {
			# arrows(x0=0, y0=0, x1=ex * pca$loadings[vars[i], 1], y1=ex * pca$loadings[vars[i], 2], col='black', lwd=1, length=0.05)
			# text(x=ex * pca$loadings[vars[i], 1], y=ex * pca$loadings[vars[i], 2], labels=varsNice[i], xpd=NA, col='black', cex=0.6, pos=2, xpd=NA)
		# }
			
		# points(pcScores[delEx$population_biasLowland == 'lowland' & delEx$population_biasHighland == 'lowland', 1:2], pch=22, bg=alpha(lowlandCol, 0.5), cex=1)
		# points(pcScores[delEx$population_biasLowland == 'highland' & delEx$population_biasHighland == 'highland', 1:2], pch=24, bg=alpha(highlandCol, 0.5), cex=1)
		# points(pcScores[delEx$population_biasLowland == 'lowland' & delEx$population_biasHighland == 'highland', 1:2], pch=21, bg=alpha(contestedCol, 0.5), cex=1)
		
		# legend('topright', inset=0.01, legend=c('Background', 'Lowland', 'Highland', 'Ambiguous'), fill=c('cornflowerblue', NA, NA, NA), border=c('black', NA, NA, NA), pt.bg=c(NA, lowlandCol, highlandCol, contestedCol), pch=c(NA, 22, 24, 21), cex=cexLegend / 2, bg=alpha('white', 0.3))
	
	# dev.off()


	# png(paste0('./Analysis/03 Delphinium exaltatum - Elevation.png'), width=1800, height=1000, res=300)

		# pres <- delEx@data[delEx$delExObs == 1, ]
		# pres <- pres[order(pres$elevation_m), ]
		
		# cols <- ifelse(
			# pres$population_biasLowland == 'lowland' & pres$population_biasHighland == 'lowland', lowlandCol,
			# ifelse(
				# pres$population_biasLowland == 'highland' & pres$population_biasHighland == 'highland', highlandCol,
				# 'orange'
			# )
		# )
		
		# plot(pres$elevation_m, pch=21, bg=cols, xlab='Rank', ylab='Elevation (m)', main='Elevation of Counties with Presences', cex=0.6, cex.main=cexMain)
		# text(1:nrow(pres), 50 + pres$elevation_m, labels=paste(pres$admin1, '-', pres$admin2), srt=90, adj=c(0, 0.5), xpd=NA, col=cols, cex=0.3)
		
		# legend('bottomright', inset=c(0.01, 0.01), pch=21, pt.bg=c(lowlandCol, highlandCol, 'orange'), legend=c('Lowland', 'Highland', 'Ambiguous'), cex=0.6)
		
		# title(sub=date(), cex.sub=cexSub, line=0)
		
	# dev.off()
	
	# save(delEx, file='./Data/03 Delphinium exaltatum - Extracted Current Climate 1950-2005 from Lorenz et al 2016 Scientific Data - Spatial Polygons.RData')

# say('###########################')
# say('### 04 assign geo-folds ###')
# say('###########################')

	# load('./Data/03 Delphinium exaltatum - Extracted Current Climate 1950-2005 from Lorenz et al 2016 Scientific Data - Spatial Polygons.RData')

	# cents <- gCentroid(delEx, byid=TRUE)

	# ### species
	# ###########
		
		# # assign presence geofolds
		# presCents <- cents[delEx$delExObs > 0, ]
		# gFoldPres <- geoFold(presCents, k=4, minIn=floor(sum(delEx$delExObs) / 4) - 1)
		
		# # assign background sites geofolds based on geofold of nearest presence
		# dists <- pointDist(cents, presCents)
		# gFoldPresIndex <- apply(dists, 1, which.min)
		# gFolds <- gFoldPres[gFoldPresIndex]
		
		# delEx$gFoldSpecies <- gFolds

		# cols <- rep(alpha('white', 0), nrow(delEx))
		# cols[delEx$gFoldSpecies == 1 & delEx$delExObs == 0] <- alpha('darkblue', 0.5)
		# cols[delEx$gFoldSpecies == 1 & delEx$delExObs > 0] <- alpha('darkblue', 1)
		# cols[delEx$gFoldSpecies == 2 & delEx$delExObs == 0] <- alpha('darkred', 0.5)
		# cols[delEx$gFoldSpecies == 2 & delEx$delExObs > 0] <- alpha('darkred', 1)
		# cols[delEx$gFoldSpecies == 3 & delEx$delExObs == 0] <- alpha('darkorange', 0.5)
		# cols[delEx$gFoldSpecies == 3 & delEx$delExObs > 0] <- alpha('darkorange', 1)
		# cols[delEx$gFoldSpecies == 4 & delEx$delExObs == 0] <- alpha('darkgreen', 0.5)
		# cols[delEx$gFoldSpecies == 4 & delEx$delExObs > 0] <- alpha('darkgreen', 1)
		
		# delExEa <- sp::spTransform(delEx, getCRS('albersNA', TRUE))
			
		# png('./Analysis/04a Delphinium exaltatum - Geo-folds for Species.png', width=2400, height=1800, res=300)

			# par(mar=mar)
		
			# plot(delExEa[delExEa$accessible_popSpecies_extentBroad, ], main='Geo-folds', border='gray', cex.main=cexMain, col=cols[delExEa$accessible_popSpecies_extentBroad])
			
			# plot(usaAdmin2AlbersNA, border='gray', lwd=lwdCounties, add=TRUE)
			# plot(canadaAdmin2AlbersNA, border='gray', lwd=lwdCounties, add=TRUE)
			# plot(mexicoAdmin2AlbersNA, border='gray', lwd=lwdCounties, add=TRUE)

			# plot(usaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
			# plot(canadaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
			# plot(mexicoAdmin1AlbersNA, lwd=lwdStates, add=TRUE)

			# legend('bottomright', inset=c(0, 0.05), legend=paste('g =', 1:4), fill=c('darkblue', 'darkred', 'darkorange', 'darkgreen'), cex=cexLegend)
			
			# title(sub=date(), cex.sub=cexSub, line=0)
			
		# dev.off()

	# ### population/bias
	# ###################
		
	# for (pop in c('lowland', 'highland')) {
				
		# for (bias in c('Lowland', 'Highland')) {
				
			# say('   Assigning g-folds for ', pop, ' population with ', tolower(bias), ' bias...')
				
			# # assign presence geofolds
			# presCents <- cents[which(delEx@data[ , paste0('population_bias', bias)] == pop), ]
			# gFoldPres <- geoFold(presCents, k=3, minIn=floor(length(presCents) / 3) - 1)
			
			# # assign background sites geofolds based on geofold of nearest presence
			# dists <- pointDist(cents, presCents)
			# gFoldPresIndex <- apply(dists, 1, which.min)
			# gFolds <- gFoldPres[gFoldPresIndex]
			
			# delEx$DUMMY <- gFolds
			# names(delEx)[ncol(delEx)] <- paste0('gFold', capIt(pop), '_bias', bias)

			# cols <- rep(alpha('white', 0), nrow(delEx))
			# cols[delEx@data[ , paste0('gFold', capIt(pop), '_bias', bias)] == 1] <- alpha('darkblue', 0.5)
			# cols[delEx@data[ , paste0('gFold', capIt(pop), '_bias', bias)] == 1 & delEx@data[ , paste0('population_bias', bias)] == pop] <- alpha('darkblue', 1)
			# cols[delEx@data[ , paste0('gFold', capIt(pop), '_bias', bias)] == 2] <- alpha('darkred', 0.5)
			# cols[delEx@data[ , paste0('gFold', capIt(pop), '_bias', bias)] == 2 & delEx@data[ , paste0('population_bias', bias)] == pop] <- alpha('darkred', 1)
			# cols[delEx@data[ , paste0('gFold', capIt(pop), '_bias', bias)] == 3] <- alpha('darkorange', 0.5)
			# cols[delEx@data[ , paste0('gFold', capIt(pop), '_bias', bias)] == 3 & delEx@data[ , paste0('population_bias', bias)] == pop] <- alpha('darkorange', 1)
			
			# delExEa <- sp::spTransform(delEx, getCRS('albersNA', TRUE))
				
			# png(paste0('./Analysis/04b Delphinium exaltatum - Geo-folds for ', capIt(pop), ' Population Assigned with Bias for ', bias, '.png'), width=2400, height=1800, res=300)

				# par(mar=mar)
			
				# thisDelExEa <- delExEa[which(delExEa@data[[paste0('accessible_pop', capIt(pop), '_extentBroad_bias', bias)]]), ]
				# cols <- cols[which(delExEa@data[[paste0('accessible_pop', capIt(pop), '_extentBroad_bias', bias)]])]
			
				# plot(delExEa[delExEa$accessible_popSpecies_extentBroad, ], main=paste0('Geo-folds for ', capIt(pop), ' Population Assigned with Bias toward ', bias), border='gray', cex.main=cexMain)
				# plot(thisDelExEa, col=cols, border='gray', lwd=lwdCounties, add=TRUE)
				
				# plot(usaAdmin2AlbersNA, border='gray', lwd=lwdCounties, add=TRUE)
				# plot(canadaAdmin2AlbersNA, border='gray', lwd=lwdCounties, add=TRUE)
				# plot(mexicoAdmin2AlbersNA, border='gray', lwd=lwdCounties, add=TRUE)

				# plot(usaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
				# plot(canadaAdmin1AlbersNA, lwd=lwdStates, add=TRUE)
				# plot(mexicoAdmin1AlbersNA, lwd=lwdStates, add=TRUE)

				# legend('bottomright', inset=c(0, 0.05), legend=c(paste(pop, 'g =', 1:3)), fill=c('darkblue', 'darkred', 'darkorange'), cex=cexLegend)
				
				# title(sub=date(), cex.sub=cexSub, line=0)
				
			# dev.off()
				
		# } # next bias
		
	# } # next population

	# save(delEx, file='./Data/04 Delphinium exaltatum - Assigned Geo-folds - Spatial Polygons.RData')

# say('###################################')
# say('### 05 calculate county weights ###')
# say('###################################')

	# say('I am calculating several types of weights.  All are scaled to (0, 1]). Each has a "Narrow" and a "Broad" variant.')
	# say('* areaWeight: log10(county area in the respective accessible region)')
	# say('* kdeAreaWeight: kdeWeight * areaWeight (rescaled to (0, 1])')
	# say('* fabaceaeWeight: log10(number of Fabaceae records in county + 1)')
	# say('* distAreaWeight: distWeight * areaWeight (rescaled to (0, 1])')

	# load('./Data/04 Delphinium exaltatum - Assigned Geo-folds - Spatial Polygons.RData')

	# say('   ### area weighting for species')
	# say('   ##############################')

		# for (ext in c('Narrow', 'Broad')) {
	
			# wPres <- wBg <- delEx@data[ , paste0('areaKm2_popSpecies_extent', ext)]
	
			# wPres[delEx$delExObs == 0 | !delEx@data[[paste0('accessible_popSpecies_extent', ext)]]] <- NA
			# wBg[!delEx@data[[paste0('accessible_popSpecies_extent', ext)]]] <- NA

			# areas <- delEx@data$areaKm2
			# wPres <- wPres * wPres / areas
			# wBg <- wBg * wBg / areas
			
			# sumPres <- sum(wPres, na.rm=TRUE)
			# sumBg <- sum(wBg, na.rm=TRUE)
			
			# if (sumPres > sumBg) {
				# wBg <- wBg * (sumPres / sumBg)
			# } else {
				# wPres <- wPres * (sumBg / sumPres)
			# }
			
			# maxVal <- max(wPres, wBg, na.rm=TRUE)

			# wPres <- wPres / maxVal
			# wBg <- wBg / maxVal

			# delEx@data$DUMMY1 <- wPres
			# delEx@data$DUMMY2 <- wBg
			# names(delEx@data)[(ncol(delEx) - 1):ncol(delEx)] <- c(
				# paste0('areaWeightPres_popSpecies_extent', ext),
				# paste0('areaWeightBg_popSpecies_extent', ext)
			# )

			# # plot
			# if (makePlot) {
				# png(paste0('./Analysis/05a Delphinium exaltatum - Area Weighting for Species Using ', ext, ' Accessible Region.png'), width=2400, height=1800, res=300)
					# plotNice(v=wBg, stretch=TRUE, title=paste0('Area Weighting for Species Using ', ext, ' Accessible Region'), leg='Weight', col='maroon4')
				# dev.off()
			# }
				
		# } # next extent

	# say('   ### area weighting for populations')
	# say('   ##################################')

		# areas <- delEx@data$areaKm2

		# for (pop in c('lowland', 'highland')) {
		
			# for (bias in c('Lowland', 'Highland')) {
	
				# for (ext in c('Narrow', 'Broad')) {
			
					# wPres <- wBg <- delEx@data[ , paste0('areaKm2_pop', capIt(pop), '_extent', ext, '_bias', bias)]

					# presIndex <- which(delEx@data[[paste0('population_bias', bias)]] == pop)
					# notPresIndex <- seq_along(wPres)[-presIndex]
					# wPres[notPresIndex] <- NA

					# notBgIndex <- which(!delEx@data[[paste0('accessible_pop', capIt(pop), '_extent', ext, '_bias', bias)]])
					# wBg[notBgIndex] <- NA

					# wBg <- wBg * (wBg / areas)
					
					# sumPres <- sum(wPres, na.rm=TRUE)
					# sumBg <- sum(wBg, na.rm=TRUE)
					
					# if (sumPres > sumBg) {
						# wBg <- wBg * (sumPres / sumBg)
					# } else {
						# wPres <- wPres * (sumBg / sumPres)
					# }
					
					# maxVal <- max(wPres, wBg, na.rm=TRUE)
		
					# wPres <- wPres / maxVal
					# wBg <- wBg / maxVal

					# delEx@data$DUMMY1 <- wPres
					# delEx@data$DUMMY2 <- wBg
					# names(delEx@data)[(ncol(delEx) - 1):ncol(delEx)] <- c(
						# paste0('areaWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias),
						# paste0('areaWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
					# )
		
					# if (makePlot) {
						# png(paste0('./Analysis/05b Delphinium exaltatum - Area Weighting for ', capIt(pop), ' Using ', ext, ' Accessible Region with Bias for ', bias, '.png'), width=2400, height=1800, res=300)
							# plotNice(v=wBg, stretch=FALSE, title=paste0('Area Weighting for ', capIt(pop), ' Using ', ext, ' Accessible Region with Bias for ', bias), leg='Weight')
						# dev.off()
					# }
						
				# } # next extent
				
			# } # next bias

		# } # next population

	# say('   ### KDE weighting for species')
	# say('   #############################')
		
		# data <- data.frame(delExObs=delEx$delExObs)
		
		# cents <- coordinates(gCentroid(delEx, byid=TRUE))
		# colnames(cents) <- c('centLong', 'centLat')
		# delEx@data <- cbind(delEx@data, cents)
		
		# centsScaled <- as.data.frame(scale(cents))
		# names(centsScaled) <- c('centLongScaled', 'centLatScaled')
		# w <- pmax(delEx@data[ , 'areaWeightBg_popSpecies_extentBroad'], delEx@data[ , 'areaWeightPres_popSpecies_extentBroad'], na.rm=TRUE)
		# data <- cbind(data, centsScaled, w)
		# dataPared <- data[-naRows(data), ]
		
		# kde <- trainGam(data=dataPared, resp='delExObs', preds=c('centLongScaled', 'centLatScaled'), w='w', construct=FALSE, select=FALSE, gamma=0.1, verbose=FALSE)

		# kdeWeight <- 1 - predict(kde, data, type='response')
		
		# for (ext in c('Narrow', 'Broad')) {
		
			# areas <- delEx@data$areaKm2
			# areasBg <- delEx@data[ , paste0('areaKm2_popSpecies_extent', ext)]
			
			# wPres <- kdeWeight * areasBg / areas
			# wPres[delEx@data$delExObs == 0] <- NA
			
			# wBg <- kdeWeight * areasBg / areas
			
			# sumPres <- sum(wPres, na.rm=TRUE)
			# sumBg <- sum(wBg, na.rm=TRUE)
			
			# if (sumPres > sumBg) {
				# wBg <- wBg * (sumPres / sumBg)
			# } else {
				# wPres <- wPres * (sumBg / sumPres)
			# }
			
			# maxVal <- max(wPres, wBg, na.rm=TRUE)

			# wPres <- wPres / maxVal
			# wBg <- wBg / maxVal

			# delEx@data$DUMMY1 <- wPres
			# delEx@data$DUMMY2 <- wBg
			# names(delEx@data)[(ncol(delEx) - 1):ncol(delEx)] <- c(
				# paste0('kdeAreaWeightPres_popSpecies_extent', ext),
				# paste0('kdeAreaWeightBg_popSpecies_extent', ext)
			# )

			# # plot
			# if (makePlot) {
				# png(paste0('./Analysis/05b Delphinium exaltatum - KDE x Area Weighting for Species Using ', ext, ' Extent.png'), width=2400, height=1800, res=300)
					# plotNice(v=wBg, stretch=TRUE, title=paste0('KDE * Area Weighting for Species Using ', ext, ' Extent'), leg='Weight', col='maroon4')
				# dev.off()
			# }

		# } # next extent

	# say('   ### KDE weighting for populations')
	# say('   #################################')
		
		# areas <- delEx@data$areaKm2

		# for (pop in c('lowland', 'highland')) {
		
			# for (bias in c('Lowland', 'Highland')) {
				
				# for (ext in c('Narrow', 'Broad')) {

					# areasBg <- delEx@data[ , paste0('areaKm2_pop', capIt(pop), '_extent', ext, '_bias', bias)]
					
					# wPres <- kdeWeight * areasBg / areas
					# presIndex <- which(delEx@data[[paste0('population_bias', bias)]] == pop)
					# notPresIndex <- seq_along(wPres)[-presIndex]
					# wPres[notPresIndex] <- NA

					# wBg <- kdeWeight * areasBg / areas
					# notBgIndex <- which(!delEx@data[[paste0('accessible_pop', capIt(pop), '_extent', ext, '_bias', bias)]])
					# wBg[notBgIndex] <- NA
					
					# sumPres <- sum(wPres, na.rm=TRUE)
					# sumBg <- sum(wBg, na.rm=TRUE)
					
					# if (sumPres > sumBg) {
						# wBg <- wBg * (sumPres / sumBg)
					# } else {
						# wPres <- wPres * (sumBg / sumPres)
					# }
					
					# maxVal <- max(wPres, wBg, na.rm=TRUE)
		
					# wPres <- wPres / maxVal
					# wBg <- wBg / maxVal

					# delEx@data$DUMMY1 <- wPres
					# delEx@data$DUMMY2 <- wBg
					# names(delEx@data)[(ncol(delEx) - 1):ncol(delEx)] <- c(
						# paste0('kdeAreaWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias),
						# paste0('kdeAreaWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
					# )
							
					# # plot
					# if (makePlot) {
						# png(paste0('./Analysis/05c Delphinium exaltatum - KDE x Area Background Weighting for ', capIt(pop), ' Using ', ext, ' Extent with ', bias, ' Bias.png'), width=2400, height=1800, res=300)
							# plotNice(v=wBg, stretch=TRUE, title=paste0('KDE * Area Weighting for ', capIt(pop), ' Using ', ext, ' Extent with ', bias, ' Bias'), leg='Weight', col='maroon4')
						# dev.off()
					# }

				# } # next extent
				
			# } # next bias
			
		# } # next population

	# say('   ### target background for species')
	# say('   #################################')

		# say('Using all records of Fabaceae in each county as county-level weights.')
		# say('NB You can clearly see state-level biases in records of Fabaceae. Also, the number of Delphinium is really low in middle/eastern states in the study region and likely just reflects the distribution of the genus over search effort.')
		
		# x <- data.frame(numDelphinium=rep(0, nrow(delEx)), numFabaceae=rep(0, nrow(delEx)))
		# delEx@data <- insertCol(x, into=delEx@data, at='delExObs', before=FALSE)
		
		# fabaceae <- read.csv('./Data/GBIF/02 occurrence Fabaceae (Unneeded Columns Deleted).csv')

		# # count number of records with coordinates in each county
		# coords <- fabaceae[!is.na(fabaceae$decimalLongitude) & !is.na(fabaceae$decimalLatitude) & !fabaceae$hasGeospatialIssues, ]
		# coords <- coords[!is.na(coords$decimalLongitude), ]
		# coords <- coords[!is.na(coords$decimalLatitude), ]
		# coords <- SpatialPoints(cbind(coords$decimalLongitude, coords$decimalLatitude), getCRS('wgs84', TRUE))

		# xInDelEx <- over(coords, delEx)
		# xInDelEx <- xInDelEx[!is.na(xInDelEx$admin2), ]
	
		# for (i in 1:nrow(xInDelEx)) {

			# thisRow <- which(delEx$admin1 == xInDelEx$admin1[i] & delEx$admin2 == xInDelEx$admin2[i])
			# delEx@data$numFabaceae[thisRow] <- delEx@data$numFabaceae[thisRow] + 1
			
		# }
	
		# # count number of records without coordinates in each county
		# noCoord <- fabaceae[is.na(fabaceae$decimalLongitude) & is.na(fabaceae$decimalLatitude) & !fabaceae$hasGeospatialIssues & !is.na(fabaceae$county) & !is.na(fabaceae$stateProvince), ]

		# if (nrow(noCoord) > 0) {
			
			# noCoord$county <- gsub(noCoord$county, pattern=' County', replacement='')
			# noCoord$county <- gsub(noCoord$county, pattern=' Cty.', replacement='')
			# noCoord$county <- gsub(noCoord$county, pattern=' Cty', replacement='')
			# noCoord$county <- gsub(noCoord$county, pattern=' Co.', replacement='')
			# noCoord$county <- gsub(noCoord$county, pattern='Ste ', replacement='Saint ')
			# noCoord$county <- gsub(noCoord$county, pattern='St. ', replacement='Saint ')
		
			# for (i in 1:nrow(noCoord)) {

				# thisRow <- which(delEx$admin1 == noCoord$stateProvince[i] & delEx$admin2 == noCoord$county[i])
				# if (length(thisRow) > 0) delEx@data$numFabaceae[thisRow] <- delEx@data$numFabaceae[thisRow] + 1
				
			# }
			
		# }

		# ### by extent for species
		# #########################

		# numFab <- delEx$numFabaceae
		# numFab[numFab == 0 & delEx$delExObs == 1] <- 1
		
		# for (ext in c('Narrow', 'Broad')) {
		
			# areas <- delEx@data$areaKm2
			# areasBg <- delEx@data[ , paste0('areaKm2_popSpecies_extent', ext)]
			
			# wPres <- numFab * areasBg / areas
			# wPres[delEx@data$delExObs == 0] <- NA
			# wPres <- log10(1 + wPres)
			
			# wBg <- numFab * areasBg / areas
			# wBg <- log10(1 + wBg)
			
			# sumPres <- sum(wPres, na.rm=TRUE)
			# sumBg <- sum(wBg, na.rm=TRUE)
			
			# if (sumPres > sumBg) {
				# wBg <- wBg * (sumPres / sumBg)
			# } else {
				# wPres <- wPres * (sumBg / sumPres)
			# }
			
			# maxVal <- max(wPres, wBg, na.rm=TRUE)

			# wPres <- wPres / maxVal
			# wBg <- wBg / maxVal

			# delEx@data$DUMMY1 <- wPres
			# delEx@data$DUMMY2 <- wBg
			# names(delEx@data)[(ncol(delEx) - 1):ncol(delEx)] <- c(
				# paste0('kdeAreaWeightPres_popSpecies_extent', ext),
				# paste0('kdeAreaWeightBg_popSpecies_extent', ext)
			# )

			# # plot
			# if (makePlot) {
				# png(paste0('./Analysis/05c Delphinium exaltatum - Target Background Weighting for Species Using ', ext, ' Extent.png'), width=2400, height=1800, res=300)
					# plotNice(v=wBg, stretch=TRUE, title=paste0('Target Background Weighting for Species Using ', ext, ' Extent'), leg='Weight', col='maroon4')
				# dev.off()
			# }

		# } # next extent

		# say('   ### target background for populations')
		# say('   #####################################')

		# areas <- delEx@data$areaKm2

		# numFab <- delEx$numFabaceae
		# numFab[numFab == 0 & delEx$delExObs == 1] <- 1
		
		# for (pop in c('lowland', 'highland')) {
		
			# for (bias in c('Lowland', 'Highland')) {
				
				# for (ext in c('Narrow', 'Broad')) {

					# areasBg <- delEx@data[ , paste0('areaKm2_pop', capIt(pop), '_extent', ext, '_bias', bias)]
					
					# wPres <- numFab * areasBg / areas
					# wPres <- log10(1 + wPres)
					# presIndex <- which(delEx@data[[paste0('population_bias', bias)]] == pop)
					# notPresIndex <- seq_along(wPres)[-presIndex]
					# wPres[notPresIndex] <- NA
					
					# wBg <- numFab * areasBg / areas
					# wBg <- log10(1 + wBg)
					# notBgIndex <- which(!delEx@data[[paste0('accessible_pop', capIt(pop), '_extent', ext, '_bias', bias)]])
					# wBg[notBgIndex] <- NA
					
					# wPres <- stretchMinMax(wPres, nudgeUp=TRUE, na.rm=TRUE)
					# wBg <- stretchMinMax(wBg, nudgeUp=TRUE, na.rm=TRUE)
					
					# sumPres <- sum(wPres, na.rm=TRUE)
					# sumBg <- sum(wBg, na.rm=TRUE)
					
					# if (sumPres > sumBg) {
						# wBg <- wBg * (sumPres / sumBg)
					# } else {
						# wPres <- wPres * (sumBg / sumPres)
					# }
					
					# maxVal <- max(wPres, wBg, na.rm=TRUE)
		
					# wPres <- wPres / maxVal
					# wBg <- wBg / maxVal

					# delEx@data$DUMMY1 <- wPres
					# delEx@data$DUMMY2 <- wBg
					# names(delEx@data)[(ncol(delEx) - 1):ncol(delEx)] <- c(
						# paste0('targetWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias),
						# paste0('targetWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
					# )
					
					# if (makePlot) {
						# png(paste0('./Analysis/05c Delphinium exaltatum - Target Background Weighting for ', capIt(pop), ' Using ', ext, ' Extent with ', bias, ' Bias.png'), width=2400, height=1800, res=300)
							# plotNice(v=wBg, stretch=TRUE, title=paste0('Target Background Weighting for ', capIt(pop), ' Using ', ext, ' Extent with ', bias, ' Bias'), leg='Weight', col='maroon4')
						# dev.off()
					# }

				# } # next extent
				
			# } # next bias
			
		# } # next population

	# save(delEx, file='./Data/05 Delphinium exaltatum - Defined County Weights - Spatial Polygons.RData')

# say('###############')
# say('### 06 ENMs ###')
# say('###############')

	# load('./Data/05 Delphinium exaltatum - Defined County Weights - Spatial Polygons.RData')
	
	# # evaluation results
	# results <- data.frame()
	
	# for (ext in c('Narrow', 'Broad')) {
	
		# for (bias in c('Lowland', 'Highland')) {
		
			# for (weightType in c('Area', 'Target', 'KDE')) {

				# say('Modeling with ', toupper(ext), ' extent with bias toward ', toupper(bias), ' and ', toupper(weightType), ' weighting', level=2)
			
				# ### collate lowland population data
				# ###################################
				
				# pop <- 'lowland'
				
				# weightColLowlandPres <- if (weightType == 'Area') {
					# paste0('areaWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# } else if (weightType == 'Target') {
					# paste0('targetWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# } else if (weightType == 'KDE') {
					# paste0('kdeAreaWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# }
					
				# weightColLowlandBg <- if (weightType == 'Area') {
					# paste0('areaWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# } else if (weightType == 'Target') {
					# paste0('targetWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# } else if (weightType == 'KDE') {
					# paste0('kdeAreaWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# }
					
				# presIndex <- which(delEx@data[[paste0('population_bias', bias)]] == pop)
				# presLowland <- delEx@data[presIndex, c('delExObs', 'pc1', 'pc2', weightColLowlandPres, paste0('gFold', capIt(pop), '_bias', bias))]
				# presLowland <- cbind(data.frame(index=presIndex, presBg=1), presLowland)

				# bgIndex <- which(delEx@data[[paste0('accessible_pop', capIt(pop), '_extent', ext, '_bias', bias)]])
				# bgLowland <- delEx@data[bgIndex, c('delExObs', 'pc1', 'pc2', weightColLowlandBg, paste0('gFold', capIt(pop), '_bias', bias))]

				# names(presLowland)[(ncol(presLowland) - 1):ncol(presLowland)] <- c('weight', 'geoFold')
				# names(bgLowland)[(ncol(bgLowland) - 1):ncol(bgLowland)] <- c('weight', 'geoFold')

				# bgLowland <- cbind(data.frame(index=bgIndex, presBg=0), bgLowland)

				# ### collate highland population data
				# ####################################

				# pop <- 'highland'
				
				# weightColHighlandPres <- if (weightType == 'Area') {
					# paste0('areaWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# } else if (weightType == 'Target') {
					# paste0('targetWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# } else if (weightType == 'KDE') {
					# paste0('kdeAreaWeightPres_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# }
					
				# weightColHighlandBg <- if (weightType == 'Area') {
					# paste0('areaWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# } else if (weightType == 'Target') {
					# paste0('targetWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# } else if (weightType == 'KDE') {
					# paste0('kdeAreaWeightBg_pop', capIt(pop), '_extent', ext, '_bias', bias)
				# }
					
				# presIndex <- which(delEx@data[[paste0('population_bias', bias)]] == pop)
				# presHighland <- delEx@data[presIndex, c('delExObs', 'pc1', 'pc2', weightColHighlandPres, paste0('gFold', capIt(pop), '_bias', bias))]
				# presHighland <- cbind(data.frame(index=presIndex, presBg=1), presHighland)

				# bgIndex <- which(delEx@data[[paste0('accessible_pop', capIt(pop), '_extent', ext, '_bias', bias)]])
				# bgHighland <- delEx@data[bgIndex, c('delExObs', 'pc1', 'pc2', weightColHighlandBg, paste0('gFold', capIt(pop), '_bias', bias))]

				# names(presHighland)[(ncol(presHighland) - 1):ncol(presHighland)] <- c('weight', 'geoFold')
				# names(bgHighland)[(ncol(bgHighland) - 1):ncol(bgHighland)] <- c('weight', 'geoFold')

				# bgHighland <- cbind(data.frame(index=bgIndex, presBg=0), bgHighland)

				# ### train by geo-fold set
				# #########################
				
				# for (set in 1:3) {
					
					# # # hand-picking geo-folds so they overlap as LITTLE as possible in extents
					# # if (bias == 'Lowland') {
					
						# # geoFoldLowland <- c(1, 2, 3)[set]
						# # geoFoldHighland <- c(1, 2, 3)[set]
						
					# # } else if (bias == 'Highland') {
					
						# # geoFoldLowland <- c(1, 2, 3)[set]
						# # geoFoldHighland <- c(1, 2, 3)[set]
						
					# # }
					
					# # hand-picking geo-folds so they overlap as MUCH as possible in extents
					# if (bias == 'Lowland') {
					
						# geoFoldLowland <- c(1, 2, 3)[set]
						# geoFoldHighland <- c(2, 1, 3)[set]
						
					# } else if (bias == 'Highland') {
					
						# geoFoldLowland <- c(1, 2, 3)[set]
						# geoFoldHighland <- c(2, 1, 3)[set]
						
					# }
					
					# trainPresLowland <- presLowland[presLowland$geoFold != geoFoldLowland, ]
					# trainBgLowland <- bgLowland[bgLowland$geoFold != geoFoldLowland, ]
					
					# trainPresHighland <- presHighland[presHighland$geoFold != geoFoldHighland, ]
					# trainBgHighland <- bgHighland[bgHighland$geoFold != geoFoldHighland, ]
					
					# testPresLowland <- presLowland[presLowland$geoFold == geoFoldLowland, ]
					# testBgLowland <- bgLowland[bgLowland$geoFold == geoFoldLowland, ]
					
					# testPresHighland <- presHighland[presHighland$geoFold == geoFoldHighland, ]
					# testBgHighland <- bgHighland[bgHighland$geoFold == geoFoldHighland, ]
					
					# ### equalize weights for training data
					# ######################################
						
						# ### wanting sum of weights of presences of populations to be equal and sum of weights of bg's to be equal

						# wPresLowland <- trainPresLowland$weight
						# wPresHighland <- trainPresHighland$weight

						# wBgLowland <- trainBgLowland$weight
						# wBgHighland <- trainBgHighland$weight

						# trainWeights <- equalizeWeights(wPresLowland=wPresLowland, wBgLowland=wBgLowland, wPresHighland=wPresHighland, wBgHighland=wBgHighland)
						
						# trainPresLowland$weight <- trainWeights$wPresLowland
						# trainBgLowland$weight <- trainWeights$wBgLowland
						
						# trainPresHighland$weight <- trainWeights$wPresHighland
						# trainBgHighland$weight <- trainWeights$wBgHighland
						
						# trainLowland <- rbind(trainPresLowland, trainBgLowland)
						# trainHighland <- rbind(trainPresHighland, trainBgHighland)

						# trainLowland$population <- 'lowland'
						# trainHighland$population <- 'highland'
						
						# trainData <- rbind(trainLowland, trainHighland)
						# trainData$population <- as.factor(trainData$population)
						
						# rownames(trainData) <- 1:nrow(trainData)

					# ### equalize weights for test data
					# ##################################
						
						# ### wanting sum of weights of presences of populations to be equal and sum of weights of bg's to be equal

						# wPresLowland <- testPresLowland$weight
						# wPresHighland <- testPresHighland$weight

						# wBgLowland <- testBgLowland$weight
						# wBgHighland <- testBgHighland$weight

						# testWeights <- equalizeWeights(wPresLowland=wPresLowland, wBgLowland=wBgLowland, wPresHighland=wPresHighland, wBgHighland=wBgHighland)
						
						# testPresLowland$weight <- testWeights$wPresLowland
						# testBgLowland$weight <- testWeights$wBgLowland
						
						# testPresHighland$weight <- testWeights$wPresHighland
						# testBgHighland$weight <- testWeights$wBgHighland
						
						# testLowland <- rbind(testPresLowland, testBgLowland)
						# testHighland <- rbind(testPresHighland, testBgHighland)

						# testLowland$population <- 'lowland'
						# testHighland$population <- 'highland'
						
						# testData <- rbind(testLowland, testHighland)
						# testData$population <- as.factor(testData$population)
						
						# rownames(testData) <- 1:nrow(testData)

					# ### model
					# #########
							
						# # GLM
						# # model <- trainGlm(data=trainData, resp='presBg', preds=c('pc1', 'pc2', 'population'), presPerTermInitial=5, presPerTermFinal=10, initialTerms=20, cubic=FALSE, w='weight', verbose=FALSE)
						# initialModel <- glm(presBg ~ pc1 + pc2 + pc1:pc2 + I(pc1^2) + I(pc2^2) + population + population:pc1 + population:pc2 + population:I(pc1^2) + population:I(pc2^2), data=trainData, family='binomial', method='brglmFit', na.action=stats::na.fail)
						
						# evals <- dredge(initialModel)
						
						# model <- get.models(evals, subset=rownames(evals)[1])[[1]]

							# cbiTest <- testModel(model, testData)
							# cbiTrain <- testModel(model, trainData)
							
							# thisResults <- data.frame(
								# algo = 'glm',
								# extent = tolower(ext),
								# bias = tolower(bias),
								# weightType = tolower(weightType),
								# modelType = 'single',
								# cbiTestSpecies = cbiTest[['cbiSpecies']],
								# cbiTestLowland = cbiTest[['cbiLowland']],
								# cbiTestHighland = cbiTest[['cbiHighland']],
								# cbiTrainSpecies = cbiTrain[['cbiSpecies']],
								# cbiTrainLowland = cbiTrain[['cbiLowland']],
								# cbiTrainHighland = cbiTrain[['cbiHighland']]
							# )
							
							# results <- rbind(results, thisResults)

							# save(model, file=paste0('./ENMs/BOTH populations ', toupper(ext), ' region ', toupper(bias), ' bias ', toupper(weightType), ' weights GLM set ', set, '.Rdata'))
						
						# # NS
						# model <- trainNs(data=trainData, resp='presBg', preds=c('pc1', 'pc2', 'population'), w='weight', presPerTermInitial=5, presPerTermFinal=10, verbose=FALSE)

							# cbiTest <- testModel(model, testData)
							# cbiTrain <- testModel(model, trainData)
							
							# thisResults <- data.frame(
								# algo = 'ns',
								# extent = tolower(ext),
								# bias = tolower(bias),
								# weightType = tolower(weightType),
								# modelType = 'single',
								# cbiTestSpecies = cbiTest[['cbiSpecies']],
								# cbiTestLowland = cbiTest[['cbiLowland']],
								# cbiTestHighland = cbiTest[['cbiHighland']],
								# cbiTrainSpecies = cbiTrain[['cbiSpecies']],
								# cbiTrainLowland = cbiTrain[['cbiLowland']],
								# cbiTrainHighland = cbiTrain[['cbiHighland']]
							# )
							
							# results <- rbind(results, thisResults)
						
							# save(model, file=paste0('./ENMs/BOTH populations ', toupper(ext), ' region ', toupper(bias), ' bias ', toupper(weightType), ' weights NS set ', set, '.Rdata'))
							
				# } # next train/test set
					
				# ### all-sites models
				# ####################
					
				# trainPresLowland <- presLowland
				# trainBgLowland <- bgLowland
				
				# trainPresHighland <- presHighland
				# trainBgHighland <- bgHighland
				
				# ### equalize weights for training data
				# ######################################
					
					# ### wanting sum of weights of presences of populations to be equal and sum of weights of bg's to be equal

					# wPresLowland <- trainPresLowland$weight
					# wPresHighland <- trainPresHighland$weight

					# wBgLowland <- trainBgLowland$weight
					# wBgHighland <- trainBgHighland$weight

					# trainWeights <- equalizeWeights(wPresLowland=wPresLowland, wBgLowland=wBgLowland, wPresHighland=wPresHighland, wBgHighland=wBgHighland)
					
					# trainPresLowland$weight <- trainWeights$wPresLowland
					# trainBgLowland$weight <- trainWeights$wBgLowland
					
					# trainPresHighland$weight <- trainWeights$wPresHighland
					# trainBgHighland$weight <- trainWeights$wBgHighland
					
					# trainLowland <- rbind(trainPresLowland, trainBgLowland)
					# trainHighland <- rbind(trainPresHighland, trainBgHighland)

					# trainLowland$population <- 'lowland'
					# trainHighland$population <- 'highland'
					
					# trainData <- rbind(trainLowland, trainHighland)
					# trainData$population <- as.factor(trainData$population)
					
					# rownames(trainData) <- 1:nrow(trainData)

				# ### model
				# #########
						
					# # GLM
					# # model <- trainGlm(data=trainData, resp='presBg', preds=c('pc1', 'pc2', 'population'), presPerTermInitial=5, presPerTermFinal=5, initialTerms=20, cubic=FALSE, w='weight', verbose=TRUE)
					# initialModel <- glm(presBg ~ pc1 + pc2 + pc1:pc2 + I(pc1^2) + I(pc2^2) + population + population:pc1 + population:pc2 + population:I(pc1^2) + population:I(pc2^2), data=trainData, family='binomial', method='brglmFit', na.action=stats::na.fail)
					
					# evals <- dredge(initialModel)
					
					# model <- get.models(evals, subset=rownames(evals)[1])[[1]]
					
						# save(model, file=paste0('./ENMs/BOTH populations ', toupper(ext), ' region ', toupper(bias), ' bias ', toupper(weightType), ' weights GLM ALL Sites.Rdata'))
					
					# # NS
					# model <- trainNs(data=trainData, resp='presBg', preds=c('pc1', 'pc2', 'population'), w='weight', presPerTermInitial=5, presPerTermFinal=10, verbose=FALSE)

						# save(model, file=paste0('./ENMs/BOTH populations ', toupper(ext), ' region ', toupper(bias), ' bias ', toupper(weightType), ' weights NS ALL sites.Rdata'))
							
			# } # next weight type
			
		# } # next extent
		
	# } # next bias

	# results$trainMinusTestSpecies <- results$cbiTestSpecies - results$cbiTrainSpecies
	# results$trainMinusTestLowland <- results$cbiTestLowland - results$cbiTrainLowland
	# results$trainMinusTestHighland <- results$cbiTestHighland - results$cbiTrainHighland
	
	# write.csv(results, './Analysis/06 Model Evaluations - Single Model for BOTH Populations.csv', row.names=FALSE)
	
# say('######################################')
# say('### 07 visualize model performance ###')
# say('######################################')

	# thisTextCex <- 0.9

	# results <- read.csv('./Analysis/06 Model Evaluations - Single Model for BOTH Populations.csv', as.is=TRUE)
	
	# results <- aggregate(results, by=list(results$algo, results$extent, results$bias, results$weightType), FUN=mean, na.rm=TRUE)
	# results$algo <- results$extent <- results$bias <- results$weightType <- NULL
	# names(results)[1:4] <- c('algo', 'extent', 'bias', 'weightType')
	
	# nas <- naRows(results[ , grepl(names(results), pattern='cbi')])
	# if (length(nas) > 0) results <- results[-nas, ]
	
	# codes <- paste(substr(results$algo, 1, 1), substr(results$extent, 1, 1), substr(results$bias, 1, 1), substr(results$weightType, 1, 1), sep='') 
	# codes <- toupper(codes)
	
	# # cols <- rep(0, length(codes))
	# # uniCodes <- unique(codes)
	# # for (thisCode in uniCodes) cols[codes == thisCode] <- max(cols, na.rm=TRUE) + 1
	
	# cols <- 1:nrow(results)
	
	# png('./Analysis/07 Model Performance.png', width=2400, height=800, res=300)
		
		# par(mar=c(5, 4, 3, 2) + 0.1, mfrow=c(1, 3), pty='s')

		# # species
		# plot(1, col='white', xlab='Train CBI', ylab='Test CBI', xlim=c(-1, 1), ylim=c(-1, 1), main='Both Populations Together')
		# abline(0, 1, col='gray')
		# abline(v=0, col='gray')
		# abline(h=0, col='gray')
		# text(results$cbiTrainSpecies, results$cbiTestSpecies, labels=codes, col=cols, xpd=NA, cex=thisTextCex)
		# mtext('algorithm | extent | bias | weight', side=1, line=4, outer=FALSE, xpd=NA, cex=0.6)

		# # lowland
		# plot(1, col='white', xlab='Train CBI', ylab='Test CBI', xlim=c(-1, 1), ylim=c(-1, 1), main='Lowland Population')
		# abline(0, 1, col='gray')
		# abline(v=0, col='gray')
		# abline(h=0, col='gray')
		# text(results$cbiTrainLowland, results$cbiTestLowland, labels=codes, col=cols, xpd=NA, cex=thisTextCex)
		# mtext('algorithm | extent | bias | weight', side=1, line=4, outer=FALSE, xpd=NA, cex=0.6)

		# # highland
		# plot(1, col='white', xlab='Train CBI', ylab='Test CBI', xlim=c(-1, 1), ylim=c(-1, 1), main='Highland Population')
		# abline(0, 1, col='gray')
		# abline(v=0, col='gray')
		# abline(h=0, col='gray')
		# text(results$cbiTrainHighland, results$cbiTestHighland, labels=codes, col=cols, xpd=NA, cex=thisTextCex)
		# mtext('algorithm | extent | bias | weight', side=1, line=4, outer=FALSE, xpd=NA, cex=0.6)

	# dev.off()

# say('###############################################')
# say('### create extractions of paleoclimate data ###')
# say('###############################################')

	# say('Using Lorenz et al 2016 Scientific Data climatologies.')
	
	# load('./Data/PCA across Present (1950-2005) from Broad Accessible Region.RData')

	# # for (gcm in c('CCSM', 'ECBilt')) {
	# # for (gcm in c('CCSM')) {
	# for (gcm in c('ECBilt')) {
	
		# for (t in seq(0, 21000, by=500)) {
	
			# say(gcm, ' ', t, ' ybp')
	
			# load('./Data/04 Delphinium exaltatum - Assigned Geo-folds - Spatial Polygons.RData')

			# cents <- gCentroid(delEx)
			# pairDist <- pointDist(delEx)
			# diag(pairDist) <- NA

			# climate <- stack(paste0('D:/ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/', gcm, '/', t, 'BP/', vars, '.tif'))
			# climate <- crop(climate, delEx)

			# extracted <- raster::extract(climate, delEx, weights=TRUE, normalizeWeight=TRUE)
			
			# clim <- matrix(NA, nrow=nrow(delEx), ncol=length(vars))
			# clim <- as.data.frame(clim)
			# names(clim) <- vars
			
			# # collate
			# for (i in 1:nrow(clim)) {
			
				# countyClim <- extracted[[i]][ , vars]
				# if (class(countyClim) == 'numeric') {
				
					# clim[i, ] <- matrix(countyClim, nrow=1, ncol=length(vars))
					
				# } else {
			
					# wgts <- extracted[[i]][ , 'weight']
					# wgts <- matrix(wgts, nrow=nrow(countyClim), ncol=ncol(countyClim))
				
					# countyClim <- countyClim * wgts
					# clim[i, ] <- colSums(countyClim)
					
				# }
				
			# }

			# # assign NAs
			# for (thisVar in vars) {
			
				# nas <- which(is.na(clim[ , thisVar]))
					
				# if (length(nas) > 0) {

					# for (thisNa in nas) {
					
						# thisPairDist <- pairDist[thisNa, ]
						
						# while (is.na(clim[thisNa, thisVar])) {
						
							# newClim <- which.min(thisPairDist)
							# clim[thisNa, thisVar] <- clim[newClim, thisVar]
							# thisPairDist[newClim] <- NA
						
						# }
					
					# }


				# }
				
			# }
			
			# # add PCA values
			# pcaPred <- predict(pca, clim)
			# colnames(pcaPred) <- paste0('pc', 1:ncol(pcaPred))
			
			# clim <- cbind(clim, pcaPred)
			# for (thisCol in colnames(clim)) delEx@data[ , thisCol] <- clim[ , thisCol]
			
			# save(delEx, file=paste0('./Data/06 Delphinium exaltatum - Extracted Past Climate ', prefix(t, 5), ' ybp with ', gcm, ' ESM from Lorenz et al 2016 Scientific Data - Spatial Polygons.RData'))
			
			# rm(clim, pastClim, countyClim, delEx); gc()
			
		# }
		
	# }

# say('#######################################################')
# say('### 08a get maximum predicted value across all time ###')
# say('#######################################################')

	# say('Doing this to know how to rescale maps.')

	# maxPredSpecies <- maxPredLowland <- maxPredHighland <- -Inf

	# for (gcm in c('CCSM', 'ECBilt')) {
	
		# for (t in seq(0, 21000, by=500)) {
	
			# say(gcm, ' ', t, ' Kybp...')
	
			# load(paste0('./Data/06 Delphinium exaltatum - Extracted Past Climate ', prefix(t, 5), ' ybp with ', gcm, ' ESM from Lorenz et al 2016 Scientific Data - Spatial Polygons.RData'))
				
			# delEx@data$population <- 'lowland'
			# delEx@data$population[1] <- 'highland'
			# delEx@data$population <- as.factor(delEx@data$population)
			# delEx@data$population <- 'lowland'
			
			# pred <- predict(model, delEx@data, type='response')
			# maxPredLowland <- max(maxPredLowland, pred)

			# delEx@data$population <- 'highland'
			
			# pred <- predict(model, delEx@data, type='response')
			# maxPredHighland <- max(maxPredHighland, pred)

			# maxPredSpecies <- max(maxPredSpecies, maxPredLowland, maxPredHighland)
			
		# }
		
	# }

	# maxPred <- data.frame(
		# maxPredSpecies = maxPredSpecies,
		# maxPredLowland = maxPredLowland,
		# maxPredHighland = maxPredHighland
	# )
	
	# write.csv(maxPred, './Analysis/08 Maximum predicted value across time.csv', row.names=FALSE)
	
# say('################################')
# say('### 08b make maps of present ###')
# say('################################')

	# load('./Data/05 Delphinium exaltatum - Defined County Weights - Spatial Polygons.RData')
	
	# maxPred <- read.csv('./Analysis/08 Maximum predicted value across time.csv', as.is=TRUE)
	
	# maxPredPrettyLowland <- max(pretty(c(0, maxPred$maxPredLowland)))
	# maxPredPrettyHighland <- max(pretty(c(0, maxPred$maxPredHighland)))
	
	# # GLM using BROAD extent with LOWLAND bias and KDE weighting
	# load('./ENMs/BOTH populations BROAD region LOWLAND bias KDE weights GLM ALL Sites.Rdata')
	
	# sink('./Analysis/08 Best Model.txt', split=TRUE)
	
		# say('Best/most consistent model was GLM on BOTH populations with BROAD accessible region and KDE weights with ambigous populations assigned to LOWLAND')
	
		# print(summary(model))
		
	# sink()
	
	# png(paste0('./Analysis/08 Delphinium exaltatum - 1950-2005.png'), width=4800, height=1800)
	
		# par(mfrow=c(1, 2))
		
		# ### lowland
		# ###########
		
		# delEx@data$population <- 'lowland'
		# delEx@data$population[1] <- 'highland'
		# delEx@data$population <- as.factor(delEx@data$population)
		# delEx@data$population <- 'lowland'
	
		# prediction <- predict(model, delEx@data, type='response')
		# prediction <- prediction / maxPredPrettyLowland
	
		# plotNice(v=prediction, title='Lowland (1950-2005)', legMax=maxPredPrettyLowland, leg='Suitability', color='darkgreen', cexMult=5)
	
		# ### highland
		# ############
		
		# delEx@data$population <- 'highland'
	
		# prediction <- predict(model, delEx@data, type='response')
		# prediction <- prediction / maxPredPrettyHighland
	
		# plotNice(v=prediction, stretch=FALSE, title='Highland (1950-2005)', legMax=maxPredPrettyHighland, leg='Suitability', color='darkgreen', cexMult=5)
	
	# dev.off()	

# say('############################')
# say('### 09 make maps of past ###')
# say('############################')

	# # GLM using BROAD extent with LOWLAND bias and KDE weighting
	# load('./ENMs/BOTH populations BROAD region LOWLAND bias KDE weights GLM ALL Sites.Rdata')

	# # PCA
	# load('./Data/PCA across Present (1950-2005) from Broad Accessible Region.RData')

	# # for rescaling predictions
	# maxPred <- read.csv('./Analysis/08 Maximum predicted value across time.csv', as.is=TRUE)
	# maxPredPrettyLowland <- max(pretty(c(0, maxPred$maxPredLowland)))
	# maxPredPrettyHighland <- max(pretty(c(0, maxPred$maxPredHighland)))

	# for (gcm in c('CCSM', 'ECBilt')) {
	# # for (gcm in c('CCSM')) {
	# # for (gcm in c('ECBilt')) {
	
		# for (t in seq(0, 21000, by=500)) {
		# # for (t in seq(0, 10000, by=500)) {
		# # for (t in seq(10500, 21000, by=500)) {
	
			# say(gcm, ' ', t, ' Kybp...')
	
			# load(paste0('./Data/06 Delphinium exaltatum - Extracted Past Climate ', prefix(t, 5), ' ybp with ', gcm, ' ESM from Lorenz et al 2016 Scientific Data - Spatial Polygons.RData'))
				
			# climRast <- stack(paste0('D:/ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/', gcm, '/', t, 'BP/', vars, '.tif'))
			
			# climRast <- crop(climRast, delEx)
			# climMat <- as.matrix(climRast)
			
			# climPc <- predict(pca, climMat)[ , 1:2]
			
			# climRastPc1 <- raster(nrows=nrow(climRast), ncols=ncol(climRast), crs=projection(climRast), ext=extent(climRast))
			# climRastPc2 <- raster(nrows=nrow(climRast), ncols=ncol(climRast), crs=projection(climRast), ext=extent(climRast))
			
			# climRastPc1[] <- climPc[ , 1]
			# climRastPc2[] <- climPc[ , 2]
			
			# climPc <- stack(climRastPc1, climRastPc2)
			# names(climPc) <- c('pc1', 'pc2')
			# climPc <- as.data.frame(climPc)
			
			# png(paste0('./Analysis/09 Delphinium exaltatum - ', gcm, ' for ', prefix(t, 5), ' ybp.png'), width=4800, height=1800)
			
				# par(mfrow=c(1, 2))
				
				# ### lowland
				# ###########
				
				# delEx@data$population <- 'lowland'
				# delEx@data$population[1] <- 'highland'
				# delEx@data$population <- as.factor(delEx@data$population)
				# delEx@data$population <- 'lowland'
			
				# predPoly <- predict(model, delEx@data, type='response')
				# predPoly <- predPoly / maxPredPrettyLowland
				
				# climPc$population <- 'lowland'
				# climPc$population[1] <- 'highland'
				# climPc$population <- as.factor(climPc$population)
				# climPc$population <- 'lowland'
				
				# predRastMat <- predict(model, climPc, type='response')
				# predRast <- raster(nrows=nrow(climRast), ncols=ncol(climRast), crs=projection(climRast), ext=extent(climRast))
				# predRast[] <- predRastMat
				# predRast <- predRast / maxPredPrettyLowland
			
				# plotNice(v=predPoly, rast=predRast, title=paste0('Lowland (', t, ' ybp using ', gcm, ' ESM)'), legMax=maxPredPrettyLowland, leg='Suitability', color='darkgreen', ext='all', cexMult=5)
				
				# text(-829618, 1587407, labels=paste0('Lowland\n', t, ' ybp'), adj=0, cex=cexMain * 8)
			
				# ### highland
				# ############
				
				# delEx@data$population <- 'highland'
			
				# predPoly <- predict(model, delEx@data, type='response')
				# predPoly <- predPoly / maxPredPrettyHighland
				
				# climPc$population <- 'highland'
				
				# predRastMat <- predict(model, climPc, type='response')
				# predRast <- raster(nrows=nrow(climRast), ncols=ncol(climRast), crs=projection(climRast), ext=extent(climRast))
				# predRast[] <- predRastMat
				# predRast <- predRast / maxPredPrettyHighland
			
				# plotNice(v=predPoly, rast=predRast, title=paste0('Highland (', t, ' ybp using ', gcm, ' ESM)'), legMax=maxPredPrettyHighland, leg='Suitability', color='darkgreen', ext='all', cexMult=5)
				
				# text(-829618, 1587407, labels=paste0('Highland\n', t, ' ybp'), adj=0, cex=cexMain * 8)
			
			# dev.off()
			
		# } # next time slice
		
	# } # next GCM

# say('#######################')
# say('### 10 animate maps ###')
# say('#######################')

	# # generate file names
	# for (gcm in c('CCSM', 'ECBilt')) {
	
		# say(gcm)
	
		# pngFiles <- paste0('./Analysis/09 Delphinium exaltatum - ', gcm, ' for ', prefix(seq(21000, 0, by=-500), 5), ' ybp.png')

		# apng(
			# pngFiles,
			# paste0('./Analysis/10 Delphinium exaltatum range shift from 0 to 21000 ybp using ', gcm, ' ESM (APNG).png'),
			# delay_num = 1,
			# delay_den = 1
		# )	
	# }
	
# say('###############################')
# say('### 11 make response curves ###')
# say('###############################')

	# load('./Data/05 Delphinium exaltatum - Defined County Weights - Spatial Polygons.RData')
	# delEx <- delEx[delEx$accessible_popSpecies_extentBroad, ]
	
	# # GLM using BROAD extent with LOWLAND bias and KDE weighting
	# load('./ENMs/BOTH populations BROAD region LOWLAND bias KDE weights GLM ALL Sites.Rdata')
	
	# countyClim <- delEx@data[ , c('pc1', 'pc2')]
	# data <- data.frame(pc1 = seq(min(countyClim$pc1), max(countyClim$pc1), length.out=1000),
		# pc2 = seq(min(countyClim$pc2), max(countyClim$pc2), length.out=1000),
		# population = 'lowland')
	
	# data$population[1] <- 'highland'
	# data$population <- as.factor(data$population)
	# data$population <- 'lowland'
		
	# png(paste0('./Analysis/11 Delphinium exaltatum - Response Curves.png'), width=2400, height=1400, res=300)
	
		# par(mfrow=c(1, 2), cex.main=cexMain)
		
		# ### PC1
		# #######

			# thisData <- data
			# thisData$pc2 <- median(delEx@data$pc2[delEx$delExObs == 1])
		
			# ### lowland
			# thisData$population <- 'lowland'
			# prediction <- predict(model, thisData, type='response')
			# plot(thisData$pc1, prediction, type='l', main='Response to PC1', xlab='PC1', ylab='Suitability', col=lowlandCol, ylim=c(0, 1))
		
			# ### highland
			# thisData$population <- 'highland'
			# prediction <- predict(model, thisData, type='response')
			# lines(thisData$pc1, prediction, type='l', col=highlandCol)
			
			# legend('topright', inset=0.01, legend=c('Lowland', 'Highland'), lwd=1, col=c(lowlandCol, highlandCol))
			
		# ### PC2
		# #######

			# thisData <- data
			# thisData$pc1 <- median(delEx@data$pc1[delEx$delExObs == 1])
		
			# ### lowland
			# thisData$population <- 'lowland'
			# prediction <- predict(model, thisData, type='response')
			# plot(thisData$pc2, prediction, type='l', main='Response to PC2', xlab='PC2', ylab='Suitability', col=lowlandCol, ylim=c(0, 1))
		
			# ### highland
			# thisData$population <- 'highland'
			# prediction <- predict(model, thisData, type='response')
			# lines(thisData$pc2, prediction, type='l', col=highlandCol)
			
		# title(sub=date(), cex.sub=cexSub, line=-1, outer=TRUE)
		
	# dev.off()	

say('DONE!!!', level=1, deco='%%%')
