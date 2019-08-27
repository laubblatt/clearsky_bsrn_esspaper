#' Perform Quantile Regression runs on (half-)hourly BSRN data
#'
#' @filename clearsky_BSRN_quantreg.r
#' @version 2019-08-26 copied stuff from M44_BSRN_hourly_quantreg_85.r for public access
#'
#' @TODO which source scripts are really required?

##' packages ----
library(quantreg)
library(lubridate)
library(data.table)
library(parallel)

library(lattice)

#### I cannot install devtools on Rstudio server version
# install.packages("devtools")
# library(devtools)
source("~/bgc/gitbgc/clearskyquantileregression/R/ClearSkyQuantileRegression.R")
source("~/bgc/gitbgc/clearskyquantileregression/R/meann.R")
source("~/bgc/gitbgc/clearskyquantileregression/R/PotentialRadiation.R")
source("~/bgc/gitbgc/clearskyquantileregression/R/data_table_quantile_regression_utils.R")


si <- Sys.info()  # to get to know on which PC I am running
unlist(si)

##' PATH  ubuntu konya, birne oder cluster ----
if (is.element(unlist(si[[4]]), c("konya-ubuntu","birne","quillotavbox"))) {
  phom <- "~/Dissertation/daten/"
  pbtm = "mrenner@dialog.bgc-jena.mpg.de:/Net/Groups/C-Side/BTM/mrenner/scratch/data/"
  #   pdatnc = "~/Dissertation/daten/fielddata/ameriflux/"
  prdata <- "~/Dissertation/daten/fielddata/BSRN/"
  prdata <- "~/bgc/ownCloud/work/manuscripts_mr/M44_ClearSky/rdata/"
  pdat <- "~/bgc/ownCloud/work/manuscripts_mr/M44_ClearSky/"
}  else {

  pdat = "/Net/Groups/C-Side/BTM/mrenner/scratch/data/fielddata/BSRN/"
  prdata = "~/bgc/gitbgc/clearsky_bsrn_esspaper/rdata/"
  pfig = "~/bgc/gitbgc/clearsky_bsrn_esspaper/figures/"
}


#### LOAD META DATA ####
(BSRNmeta  = fread(paste0(pdat, "dataETHZ/BSRN_state.csv")) )
BSRNmeta[ , SiteCode := toupper(site)]
str(BSRNmeta)

(BSRNStationMeta  = fread(paste0(pdat, "dataBSRN/BSRN_Stations.csv")) )
BSRNStationMeta[ , URIofevent := NULL]
# tion no: 11; Surface type: tundra; Topography type: mountain valley, rural; Horizon: doi:10.1594/PANGAEA.669522; Station scientist: Marion.Maturilli@awi.de
BSRNStationMeta[ , Comment := NULL]
str(BSRNStationMeta)
setkey(BSRNStationMeta,SiteCode)

list.files(pdat)
(BSRNclearsky  = fread(paste0(pdat, "dataETHZ/CSW_mm_Swclidnz_nighttime0_201709.dat"), na.strings = "-999.9") )
setnames(BSRNclearsky, c("V1","V2", "V15"), c("site", "year", "annual") )
BSRNclearsky[ , SiteCode := toupper(site)]

dtclearsky =  melt.data.table(BSRNclearsky, id.vars = c("SiteCode", "year", "site") , value.name = "Swclidnz" )
dtclearsky[ , month := as.numeric(substr(variable,2,4)) - 2]
dtclearsky[, variable := NULL]
dtclearsky[, site := NULL]
setkey(dtclearsky,SiteCode, year, month)
dtclearsky

#### LOAD data ####
# save(dt30, file = paste0(prdata,"BSRN_dt30.rdata"))
# save(dt60, file = paste0(prdata,"BSRN_dt60.rdata"))
list.files(prdata)
load(paste0(prdata,"BSRN_dt30.rdata"))
load(file = paste0(prdata,"BSRN_dt60.rdata"))
dt30 <- dt30[!is.na(Date), ]


dt30 = merge(dt30, BSRNmeta[ , list(SiteCode, lat, lon) ], by = "SiteCode", all.x = TRUE)
dt30[ , Rsdpot_12 := calc_PotRadiation_CosineResponsePower(doy = yday(Date), hour = Time/3600 + 0.25,
                                                           latDeg = lat ,
                                                           longDeg = lon,
                                                           timeZone = 0, isCorrectSolartime = TRUE,
                                                           cosineResponsePower = 1.2 )]

dt30[ , Rsdpot_tzlon := calc_PotRadiation_CosineResponsePower(doy = yday(Date), hour = Time/3600 + 0.25,
                                                           latDeg = lat ,
                                                           longDeg = lon,
                                                           timeZone = 0, isCorrectSolartime = TRUE,
                                                           cosineResponsePower = 1 )]


#### DO Quantile Regression ####
## qr takes some time to finish
dt30[!is.na(IncomingShortwave) , nyrmon := .N , by = list(SiteCode, year(Date), month(Date))]
dt30[nyrmon > 2*24*26,]

## @update 20181204 better treshold for winter Rsdpot
dt30[  , nyrmon_Rp10 := sum(!is.na(IncomingShortwave) & Rsdpot_12 > 10) , by = list(SiteCode, year(Date), month(Date))]
# xyplot(nyrmon_Rp10 ~ I(month(Date)) | SiteCode, data = dt30[SiteCode %in% c("SPO","LIN","GVN"), ])
dt30[nyrmon_Rp10 > 100,]
dt30[nyrmon_Rp10 <= 100,]

#### AGGREGATE TO MONTHLYdiurnal cylce with 26 days required and then average the monthly mean diurnal cycle ####
dtyrmondiurnal = dt30[ , lapply(.SD, meann, nmin = 26, na.rm = TRUE),  by = list(SiteCode, year(Date), month(Date), Time)]
dtyrmon = dtyrmondiurnal[ , lapply(.SD, mean, na.rm = FALSE),  by = list(SiteCode, year, month)]
key(dtyrmon)
setkey(dtyrmon,SiteCode,year,month)
# merge with dtyrmon
dtyrmon =  merge( dtyrmon, dtclearsky)
dtyrmon =  merge( dtyrmon, BSRNStationMeta)

detectCores()
setkey(dt30, "SiteCode")
system.time(
  resultstau <- mclapply( unique(dt30[[ "SiteCode"]]),
                          function(x) {
                            dt30[ .(x), ][ nyrmon_Rp10 > 100,
                                           mlm.output.statlong.call.rq("IncomingShortwave ~ Rsdpot_12",
                                                                  data = .SD, tau = seq(0.7,0.99,0.01)),
                                           by = list(SiteCode, year(Date), month(Date))]
                          }, mc.cores = 10)
)
# user   system  elapsed
# 1902.499    5.448  237.105
resultstau[1]
resultstau[46]
#sapply(resultssite, is.list)
dt30qrtaus = rbindlist(resultstau,fill = TRUE)


dtyrmontau =  merge( dtyrmon, dcast(dt30qrtaus, ... ~ paste0("dt30rq","_",statistic)), by = c("SiteCode", "year", "month"))
dtyrmontau
dtyrmontau[ , IncomingShortwaveClearSky :=   dt30rq_slope1 * Rsdpot_12]

save(dtyrmontau, file = paste0(prdata,"dtyrmontau.rdata"))

### run per site with windowing
system.time(
qr30yrmon <- dt30[nyrmon_Rp10 > 100 , calc_ClearSky_QuantileRegression_MonthlyTimeWindow(Date = Date,
                           Time = Time, IncomingShortwave = IncomingShortwave,
                           IncomingShortwavePotential = Rsdpot_12,
                           mc.cores = 10,tau = 0.85, pdev = 0.25), by = SiteCode]
)

qr30yrmon

# Error in `[.data.table`(dth, , calc_ClearSky_QuantileRegression(IncomingShortwave,  :
#                                                                   j doesn't evaluate to the same number of columns for each group
#
# [.data.table`(dth, , calc_ClearSky_QuantileRegression(IncomingShortwave,
#     IncomingShortwavePotential, tau), by = list(year(Date), month(Date))) at
#  ClearSkyQuantileRegression.R#155



load(file = paste0(prdata,"BSRN_dt60.rdata"))
dt60 = merge(dt60, BSRNmeta[ , list(SiteCode, lat, lon) ], by = "SiteCode", all.x = TRUE)

dt60[ , Rsdpot_12 := calc_PotRadiation_CosineResponsePower(doy = yday(Date), hour = Time/3600 + 0.5,
                                                           latDeg = lat ,
                                                           longDeg = lon,
                                                           timeZone = 0, isCorrectSolartime = TRUE,
                                                           cosineResponsePower = 1.2 )]

dt60[  , nyrmon_Rp10 := sum(!is.na(IncomingShortwave) & Rsdpot_12 > 10) , by = list(SiteCode, year(Date), month(Date))]

system.time(
  qr60yrmon <- dt60[nyrmon_Rp10 > 50 , calc_ClearSky_QuantileRegression_MonthlyTimeWindow(Date = Date,
                                                                        Time = Time, IncomingShortwave = IncomingShortwave,
                                                                        IncomingShortwavePotential = Rsdpot_12,
                                                                        mc.cores = 10,tau = 0.85, pdev = 0.25), by = SiteCode]
)
# User      System verstrichen
# 1110.047     203.410     652.234
qr60yrmon
warnings()

library(latticeExtra)
xyplot(ftau ~ I(year + month/12) | SiteCode, data = qr60yrmon[!is.na(Window) ,], type = "l") +
  xyplot(ftau ~ I(year + month/12) | SiteCode, data = qr30yrmon[!is.na(Window) ,], type = "l", col = 2)






