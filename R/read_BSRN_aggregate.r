#'  clean BSRN data sync time and potential radiation

#' @filename  read_BSRN_aggregate.r
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
#' @depends read_BSRN_2datatable.r

#' @version 2018-11-16 meta data from Matthias Schwarz 
#' @version 2018-11-19 differences in 1min or 30min aggregates of rq regression results
#'                      1min value under cloudy sky can yield larger SWD for short periods
#'                      this effect can be reduced using a 30 min avg. 
#' @version 2018-11-19 reducing the quantile from 0.99 to 0.9 improves the fit and gets it closer to 
#'                      the Long and Ackermann 2000 estimate 
#'                      and gets similar trends !
#' @version 2018-11-20 BSRN data is reported in UTC, not local time!
#'                      eases calc of potential radiation with TimeZone_h.n = 0
#'                      this may require to calc local time and shift data                       
#' @version 2018-11-20 First eval across all sites very promissing
#'                      almost all sites show a small positive offset with tau = 0.9
#' @version 2018-11-21 implement aggregation procedure and flagging after Roesch2011 and the Long and Dutton BSRN guides
#'                      output a 15min aggregate file for all sites 
#'                      key is to compute a 15min aggregate which requires just 3 non-NA values 
#'                      DATA contains bad values, see when testing for XIA 
#' @version 2018-11-22 very similar estimates for 30min and 60min aggregates 
#' 
#' @version 2018-11-23 tried to identify the time steps when Rsd obs are in the QR regression within interval +-slope 
#'                      this not as easy since the extraterrestrial radiation is not perfectly fitting the shape of Rsd at surface
#'                      
#' @version 2018-11-23 CHECK for higher attenuation at high SZA ####
#### the power law allows a more linear relation between potential and actual solar radiation 
#### this alters the regression slope (increasing) by more than the sd of the coef
#### thus the modified potential solar radiation allows for a more reasonable identification 
#### of clear sky conditions also under higher SZA
#### Under low sun angles not correcting for the Cosine Response leads to an overestimation by the potential 
#### with a underestimation at noon (with higher loads). 
#' @version 2018-11-27 The coefficient seems to be broadly good with b = 1.2 
#'                      evaluated at LIN, BAR, BER, DAR at M44_BSRN_RsdpotExponent.r
#'                      
#' @version 2018-11-27 Update of the data cleaning and aggregation procedure (largely following Roesch 2011 M7) 
#'                      Fill IncomingShortwave SWD data with SumSW  
#' 
#' @version 2018-11-28 performed QR based on b = 1.2 for all sites IMPROVED comparison with LA2000
#' @version 2018-11-28 calc daily avg 
#' 
#' @version 2018-12-13 update the flagging procedure to account for excentricity of solar constant 
#' @version 2018-12-13 update negative fluxes set to 0 when SZA > 80 
#'  
#' @version 2019-08-28 update stuff for public access 

#' @TODO   
#' @version 2018-11-26 TODO METEO measurements are sometimes only reported every 10min 60 min and need a different nmin for aggregating 
#' 



# R

##' packages ----
library(data.table)

# install.packages("devtools")
# library(devtools)
# install_github("laubblatt/cleaRskyQunatileRegression")
if (inherits(try(library(cleaRskyQuantileRegression)),"try-error") ) { 
  print("library cleaRskyQuantileRegression not available, use source instead")
  source("~/bgc/gitbgc/clearskyquantileregression/R/PotentialRadiation.R")
  source("~/bgc/gitbgc/clearskyquantileregression/R/meann.R")
}



##' source functions ---- 
# pr = "~/Dissertation/maiphd/R/"
# source(paste(pr,"data.table.utils.R",sep=""))
# # source(paste(pr,"data.table.regression.fun.R",sep=""))
# source(paste(pr,"summary_fun.R",sep=""))
# source(paste(pr,"utils.ncdf.R",sep=""))
# source(paste(pr,"utils.rbash.R",sep=""))
# source(paste(pr,"axis_labels.R",sep=""))
# source(paste(pr,"fCalcPotRadiation_CosineResponsePower.R",sep=""))
# source(paste(pr,"utils_solarradiation_astro.R",sep=""))

#### I cannot install devtools on Rstudio server version 
# install.packages("devtools")
# library(devtools)
# install_github("laubblatt/phaselag")
# source("~/bgc/github/laubblatt/phaselag/R/data.table.regression.utils.R")


## path to the high resolution data 
pdat = "/Net/Groups/C-Side/BTM/mrenner/scratch/data/fielddata/BSRN/"

prdata = "~/bgc/gitbgc/clearsky_bsrn_esspaper/rdata/"
pfig = "~/bgc/gitbgc/clearsky_bsrn_esspaper/figures/"
pdataETHZ = "~/bgc/gitbgc/clearsky_bsrn_esspaper/dataETHZ/"

(BSRNStationMeta  = fread(paste0(prdata, "BSRN_Stations.csv")) )
BSRNStationMeta[ , URIofevent := NULL]
# tion no: 11; Surface type: tundra; Topography type: mountain valley, rural; Horizon: doi:10.1594/PANGAEA.669522; Station scientist: Marion.Maturilli@awi.de
BSRNStationMeta[ , Comment := NULL]
str(BSRNStationMeta)
setkey(BSRNStationMeta,SiteCode)
BSRNStationMeta[ , .(SiteCode,Elevation)]


(dtfiles = data.table(filename = list.files(paste0(pdat,"rdata/originaltimestep/") ) ))
dtfiles[ , SiteCode :=  tstrsplit(filename, split = "\\.")[1]]

dtfiles
sico = "SOV"

load(paste0(pdat,"rdata/originaltimestep/", sico, ".rdata"))
dtbsrn


#### loop through BSRN data, aggrgating to 30 min ####

# dtfiles[ , unique(SiteCode)] %in% dt60[ , unique(SiteCode)]
# dtfiles[ !SiteCode %in% dt60[ , unique(SiteCode)]  , unique(SiteCode)]

#### TODO Implement aggreagting strategy Roesch et al, 2011 AMT, method M7
# [x] get solar zenit angle at 1min  
# [x] set nighttime 93° to 0
# [x] calc sum of DIF and DIR https://bsrn.awi.de/data/calculation-of-global-radiation/
# [x] remove physical possible outliers 
# [x] get 15 min avg with nmin = 20%
# [x] get 30 min data with na.rm = FALSE
# [x] get 60 min data with nmin = 3


dt30 = data.table()
dt60 = data.table()
for (sico in dtfiles[ , unique(SiteCode)]) {
  # for (sico in dtfiles[ !SiteCode %in% dt60[ , unique(SiteCode)]  , unique(SiteCode)]  ) {
  load(paste0(pdat,"rdata/originaltimestep/",  sico, ".rdata"))
  if ("LWD" %in% colnames(dtbsrn) )   setnames(dtbsrn, "LWD", "IncomingLongwave")
  setnames(dtbsrn, "SWD", "IncomingShortwave")
  
  dtbsrn[ , SolElevRad := SolElev_rad(doy = yday(Date), hour = Time/3600, latDeg = BSRNStationMeta[ SiteCode == sico, Latitude ], longDeg = BSRNStationMeta[ SiteCode == sico, Longitude], timeZone = 0,  isCorrectSolartime = TRUE )]
  
  # ### check nighttime negative values after flagging 
  # # set negative IncomingSolar to 0 when nighttime, by Solar Zenit Angle > 93°
  # # dtbsrn[hour(Time) %in% 4,][1:60,(SolElevRad * 180 / pi), by = .(Time, IncomingShortwave)]
  # dtbsrn[(SolElevRad * 180 / pi) < 3  & IncomingShortwave < 0,]
  # 
  # dtbsrn[(SolElevRad * 180 / pi) < -3  & IncomingShortwave < 0,]
  # dtbsrn[(SolElevRad * 180 / pi) < -3  & IncomingShortwave < 0,IncomingShortwave := 0]
  # 
  # if ("DIF" %in% colnames(dtbsrn))  dtbsrn[(SolElevRad * 180 / pi) < -3  & DIF < 0, DIF := 0]
  # if ("DIR" %in% colnames(dtbsrn))  dtbsrn[(SolElevRad * 180 / pi) < -3  & DIR < 0, DIR := 0]
  # # dtbsrn
  
  #### flag and set to NA ####  
  #### EXTREMELY RARE LIMITS to minute data after Long and Dutton recommended in BSRNtoolbox
  #    and Roesch et al. 2011, Table 2
  # @update 20181213 Long and Shi 2008, excentricity must be accounted for in the limits  
  dtbsrn[ , S_a := 1368 * utils_eccentricitycorrectionfactor(yday(Date))]
  dtbsrn[ , S0mupower12 :=  S_a * ifelse(sin(SolElevRad) > 0,  abs(sin(SolElevRad))^1.2, -abs(sin(SolElevRad))^1.2) ]
  dtbsrn[ , S0mupower12 := ifelse(S0mupower12 < -4, 0, S0mupower12) ]
  
  ## extremely rare limits Roesch et al., 2011 Table 2 Long and Shi 2008
  # dtbsrn[IncomingShortwave < -2,]
  #' @update 20181213 thermal offset (Dutton et al., 2001), 
  #'   Extremely rare minimum is set to -2W m-2 , however,
  #'    some sites have consistently lower values resulting
  #'     in missing values during night time, especially (SBO, SOV, TAM) 
  #' To cope with these issues we set negative values < -10 to missing.
  
  dtbsrn[IncomingShortwave < -10, IncomingShortwave := NA]
  # dtbsrn[IncomingShortwave > 20 & IncomingShortwave >  1.2 * S0mupower12 + 50, ]
  dtbsrn[IncomingShortwave > 20 & IncomingShortwave >  1.2 * S0mupower12 + 50, 
         IncomingShortwave := NA]
  
  if ("DIF" %in% colnames(dtbsrn)) {
    dtbsrn[DIF < -10, DIF := NA]
    dtbsrn[DIF > 10 &  DIF >  0.75 * S0mupower12 + 30,      DIF := NA]
    # xyplot(DIF + I(0.75 * S0mupower12 + 30) ~ Time | Date, data = dtbsrn[Date == "2005-04-14", ])
    # xyplot(DIF + I(0.75 * S0mupower12 + 30) ~ Time | Date, data = dtbsrn[Date == "2014-06-17", ])
  }  
  
  ### LIMITS TO DIR are meant as cos(SZA) * DIR !!! ####
  ### this is more complicated since a nighttime cos(SZA) gets negative and then there will be too many outliers
  if ("DIR" %in% colnames(dtbsrn)) {
    dtbsrn[DIR < -10, DIR := NA]
    dtbsrn[DIR > S_a, DIR := NA]
    ### remove too large DIR when daytime  
    # dtbsrn[DIR > 30 &  (sin(SolElevRad) *  DIR) >  0.95 * S0mupower12 + 10 , ]
    dtbsrn[DIR > 30 &  (sin(SolElevRad) *  DIR) >  0.95 * S0mupower12 + 10 , DIR := NA]
    
    ### remove large values of DIR when nighttime 
    dtbsrn[DIR > 30 &  (SolElevRad * 180 / pi) < -3  ,  DIR := NA]
    ### remove very large values of DIR after sunset  
    dtbsrn[DIR > 150 &  (SolElevRad * 180 / pi) <= 0 ,   DIR := NA]
    ### remove very large values of DIR before sunset  
    dtbsrn[DIR > 250 &  (SolElevRad * 180 / pi) <= 1  , DIR := NA]
    
    # xyplot( DIR + I(0.95 * S0mupower12 + 10)   ~ Time | Date, data = dtbsrn[Date == "2006-06-13", ], type = "p")
    # xyplot( DIR + I(0.95 * S0mupower12 + 10)   ~ Time | Date, data = dtbsrn[Date == "2009-03-31", ], type = "p")
  }  
  
  ### check for negative nighttimers
  ### update 2018-12-13 set all remaining negative solar rad values to 0 when SZA > 80° to avoid low solar angle contamination with negative values 
  # dtbsrn[(SolElevRad * 180 / pi) < 10  & IncomingShortwave < 0,]
  # dtbsrn[(SolElevRad * 180 / pi) < 10  & IncomingShortwave < 0, min(IncomingShortwave)]
  # dtbsrn[(SolElevRad * 180 / pi) < 10  & IncomingShortwave < 0, sum(IncomingShortwave< -4)]
  # dtbsrn[(SolElevRad * 180 / pi) < 10  & IncomingShortwave < 0, median(IncomingShortwave), by = year(Date)]
  # dtbsrn[(SolElevRad * 180 / pi) < 10  & DIF < 0,]
  # dtbsrn[(SolElevRad * 180 / pi) < 10  & DIR < 0,]
  dtbsrn[(SolElevRad * 180 / pi) < 10  & IncomingShortwave < 0, IncomingShortwave := 0]
  if ("DIF" %in% colnames(dtbsrn))  dtbsrn[(SolElevRad * 180 / pi) < 10  & DIF < 0, DIF := 0]
  if ("DIR" %in% colnames(dtbsrn))  dtbsrn[(SolElevRad * 180 / pi) < 10  & DIR < 0, DIR := 0]
  
  
  # [x] calc sum of DIF and DIR https://bsrn.awi.de/data/calculation-of-global-radiation/
  if ( all( c("DIF","DIR") %in% colnames(dtbsrn)) ) 
    dtbsrn[  , SumSW :=  as.integer(DIF + DIR * sin(SolElevRad))]
  
  # dtbsrn[  , SumSWZenit :=  as.integer(DIF + DIR * cos(SolZenitRad))]
  # xyplot(IncomingShortwave + SumSW + SumSWZenit + DIR  ~ Time | Date, data = dtbsrn[Date == "1997-01-18", ], type = "l")
  # dtbsrn[is.na(IncomingShortwave) & (!is.na(SumSW)) , ]
  
  ## original version substitution with 
  ##    -	MW Generally the BSRN recommendation is to use the sum of diff + dir whenever possible as it is considered more accurate than the global radiation measurements. 
  ## MR: Since we here aim to use global radiation measurements to estimate clear sky conditions I think this choice is in agreement with the aims of the study
  dtbsrn[is.na(IncomingShortwave) & (!is.na(SumSW)) , IncomingShortwave := SumSW  ]
  
  
  # dtbsrn
  
  ### remove columns not required anymore 
  dtbsrn[ , S_a := NULL]
  dtbsrn[ , S0mupower12 := NULL]
  dtbsrn[ , SolElevRad := NULL]
  
  
  # [x] get 15 min avg with nmin = 20%, hence 3 obs , round to 1 digit 
  dt15 = dtbsrn[ , lapply(.SD, function(x) round(meann(x, nmin = 3, na.rm = TRUE),1) ),
                 by = list(SiteCode,Date, Time = as.ITime(floor(as.integer(Time) / (900)) * 900) ) ]
  
  # [x] get 30 min data with na.rm = FALSE
  # [ ] meteo data may have only 10min / 30 min resolution ! 
  dt30 = rbind(dt30, 
               dt15[ , lapply(.SD, function(x) round(mean(x, na.rm = FALSE),2) ),
                     by = list(SiteCode,Date, Time = as.ITime(floor(as.integer(Time) / (1800)) * 1800) ) ]
               , fill = TRUE)
  
  # [x] get 60 min data with nmin = 3
  dt60 = rbind(dt60, 
               dt15[ , lapply(.SD, function(x) round(meann(x, nmin = 3, na.rm = TRUE),2) ), 
                     by = list(SiteCode,Date,Time = as.ITime(floor(as.integer(Time) / (3600)) * 3600) ) ]
               , fill = TRUE)
  
  #[ ] save 15 min data per site 
  # save(dt15, file = paste0(prdata, sico,"_dt15.rdata" ))
  print(sico)
}

dt30
dt60

save(dt30, file = paste0(prdata,"BSRN_dt30.rdata"))
save(dt60, file = paste0(prdata,"BSRN_dt60.rdata"))