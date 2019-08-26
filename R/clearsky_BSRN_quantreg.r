#' Perform Quantile Regression runs on (half-)hourly BSRN data 
#' 
#' @filename clearsky_BSRN_quantreg.r
#' @version 2019-08-26 copied stuff from M44_BSRN_hourly_quantreg_85.r for public access 
#' 
#' 

##' packages ----
library(quantreg)
library(lubridate)
library(data.table)
library(parallel)


#### I cannot install devtools on Rstudio server version
# install.packages("devtools")
# library(devtools)
source("~/bgc/gitbgc/clearskyquantileregression/R/ClearSkyQuantileRegression.R")




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
  prdata = "/Net/Groups/C-Side/BTM/mrenner/scratch/data/fielddata/BSRN/rdata/"
  pfig = "/Net/Groups/C-Side/BTM/mrenner/scratch/data/fielddata/BSRN/figures/"
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
load(paste0(prdata,"BSRN_dt30.rdata"))
load(file = paste0(prdata,"BSRN_dt60.rdata"))
dt30 <- dt30[!is.na(Date), ]

