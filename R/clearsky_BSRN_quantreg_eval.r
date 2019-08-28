#' Evaluation of the Quantile regression output with the Long and Ackerman 2000 estimates 
#' 
#' @filename clearsky_BSRN_quantreg_eval.r
#' @depends clearsky_BSRN_quantreg.r
#' 
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
#' @references  Renner, M., M. Wild, M. Schwarz, and A. Kleidon.
#'   "Estimating Shortwave Clear-Sky Fluxes from Hourly Global
#'     Radiation Records by Quantile Regression."
#'       Earth and Space Science, 2019.
#'       \url{https://doi.org/10.1029/2019EA000686}

#' 

library(data.table)
library(knitr)
library(phaselag)
library(latticeExtra)


#### I cannot install devtools on Rstudio server version
# install.packages("devtools")
# library(devtools)
# install_github("laubblatt/phaselag")
source("~/bgc/github/laubblatt/phaselag/R/data.table.regression.utils.R")


pdat = "/Net/Groups/C-Side/BTM/mrenner/scratch/data/fielddata/BSRN/"
prdata = "~/bgc/gitbgc/clearsky_bsrn_esspaper/rdata/"
pfig = "~/bgc/gitbgc/clearsky_bsrn_esspaper/figures/"


load(file = paste0(prdata,"qr30yrmon.rdata"))
load(file = paste0(pdat,"rdata/dtyrmon.rdata"))

setkey(dtyrmon,SiteCode,year,month)
setkey(qr30yrmon,SiteCode,year,month)

compare4package = merge(qr30yrmon, dtyrmon[ , .(SiteCode,year,month, ftau85dt30)])
xyplot(ftau ~ ftau85dt30 | SiteCode, data  = compare4package)
compare4package[abs(ftau - ftau85dt30) > 0.01,   ]
xyplot(ftau ~ ftau85dt30, data  = compare4package[abs(ftau - ftau85dt30) > 0.01,   ])

# SiteCode year month IncomingShortwave IncomingShortwavePotential      ftau IncomingShortwaveClearSky Window  tau ftau85dt30
# 1:      ALE 2008     3         29.470228                  26.684923 0.9912509                 26.451453   3mon 0.85  1.0141337
# 2:      ALE 2010    10          2.350706                   2.733753 0.8630349                  2.359324   3mon 0.85  0.8769211
# 3:      ALE 2011    10          1.906384                   2.733753 0.8531451                  2.332288   3mon 0.85  0.8652986
# 4:      BAR 1994     9                NA                 138.530830 0.7904661                109.503919   7mon 0.85  0.7625009
# 5:      BAR 2002    10         16.393817                  40.830352 0.8194527                 33.458543   7mon 0.85  0.6878413
# 6:      BAR 2005    10         16.647989                  40.830352 0.8283206                 33.820622   7mon 0.85  0.6694574
# 7:      BAR 2007    10         17.317436                  40.830352 0.8129472                 33.192921   7mon 0.85  0.6917041
# 8:      DOM 2006     8          4.450501                   3.322066 1.2837302                  4.264637   3mon 0.85  1.3053122
# 9:      DOM 2007     8          5.655687                   3.322066 1.2261801                  4.073451   3mon 0.85  1.2486729
# 10:      DOM 2008     8          5.614348                   3.875492 1.2473199                  4.833979   3mon 0.85  1.2616005
# 11:      DOM 2009     8          6.070957                   3.455931 1.1882435                  4.106487   5mon 0.85  1.2461970
# 12:      EUR 2008     3         43.587140                  38.688487 1.1132326                 43.069283   3mon 0.85  1.1380875
# 13:      EUR 2009     3         39.022324                  35.876915 1.0339574                 37.095203   3mon 0.85  1.0675016
# 14:      EUR 2010     3         38.920430                  35.876915 1.0492918                 37.645354   3mon 0.85  1.0786279
# 15:      EUR 2011     3         34.925941                  35.876915 1.0721293                 38.464691   3mon 0.85  1.0972845
# 16:      EUR 2011    10          5.420766                   7.105435 0.8714444                  6.191992   5mon 0.85  0.8983007
# 17:      GVN 2005     8         15.165874                  13.821590 1.1204769                 15.486772   3mon 0.85  1.1362084
# 18:      GVN 2009     8         14.762493                  13.821590 1.1049603                 15.272308   3mon 0.85  1.1213720
# 19:      SPO 1993     9          5.545826                   3.004640 1.1990445                  3.602697   3mon 0.85  1.2590854
# 20:      SPO 1995     9          6.452917                   3.004640 1.1166532                  3.355141   3mon 0.85  1.1524446
# 21:      SPO 1996     9          8.080520                   4.150829 1.0600444                  4.400063   3mon 0.85  1.1104346
# 22:      SPO 1997     9          5.000876                   3.004640 1.1352014                  3.410871   3mon 0.85  1.1508503

########################
#### data submission 
ClearSky_BSRN_monthly = dtyrmon[! is.na(Rsdpot_12) , .(SiteCode, year, month, IncomingShortwave, 
             IncomingShortwavePotential = Rsdpot_12, beta_qr = ftau85dt30, IncomingShortwaveClearSky, IncomingShortwaveClearSky_LA2000 = Swclidnz)]

write.table(format(ClearSky_BSRN_monthly, digits=2), file = paste0(prdata,"ClearSky_BSRN_monthly.tab"), sep = "\t", row.names = FALSE, quote = FALSE)

##################################################
#### START EVALUATION CODE #####
#### this will merge with the dtclearsky data to get the Long and Ackerman data 

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

dt30yrmon =  merge( qr30yrmon, dtclearsky)

dt30yrmon[ ,  ftau_LA2000 := Swclidnz/IncomingShortwavePotential]
dt30yrmon[ , SWCRE_QR := IncomingShortwave - IncomingShortwaveClearSky]
save(dt30yrmon, file = paste0(prdata,"dt30yrmon.rdata"))


dt60yrmon =  merge( qr60yrmon, dtclearsky, by = c("SiteCode","year", "month"))
dt60yrmon[ ,  ftau_LA2000 := Swclidnz/IncomingShortwavePotential]
dt60yrmon[ , SWCRE_QR := IncomingShortwave - IncomingShortwaveClearSky]
save(dt60yrmon, file = paste0(prdata,"dt60yrmon.rdata"))


#### calc a skill score ####
dt30yrmon[ , MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = IncomingShortwavePotential  * 0.81 )]
# [1] 0.718883
# xyplot(Swclidnz + IncomingShortwaveClearSky + IncomingShortwave ~ I(year+month/12), data = dt30yrmon[SiteCode == "DAA", ], type = "b")
# xyplot(IncomingShortwavePotential + IncomingShortwave ~ Time | Date, data = dt30[SiteCode == "DAA" & year(Date) == 2005, ])

(dt30yrmon_skill =  dt30yrmon[ , .(SSMSE =  MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = IncomingShortwavePotential* 0.81 ),
                               .N, RMSE = rmse(o = Swclidnz, p = IncomingShortwaveClearSky), r = cor(Swclidnz, IncomingShortwaveClearSky, use = "pair")  ), by = SiteCode][order(SSMSE),])

#### Monthly mean statistics reported in text ####
#' @update 20181213 remove condition of R1 in stats
(ymon_Rsdcs_lmstats = dcast(
  dt30yrmon[ , mlm.output.statlong.call(IncomingShortwaveClearSky~Swclidnz, .SD), by = SiteCode]
  , ... ~ statistic))

## NAIVE estimate using a fixed fration of Rsdpot (may need a skill score)
(ymon_Rsdpotcs_lmstats = dcast(
  dt30yrmon[ , mlm.output.statlong.call(I(IncomingShortwavePotential * 0.81) ~Swclidnz, .SD), by = SiteCode]
  , ... ~ paste0("Rsdpot_",statistic)) )
ymon_Rsdpotcs_rmse = dt30yrmon[ , .(Rsdpot_RMSE = rmse(I(IncomingShortwavePotential * 0.81), Swclidnz)) , by = SiteCode]

ymon_Rsdcs_rmse = dt30yrmon[ , .(RMSE = rmse(IncomingShortwaveClearSky, Swclidnz)) , by = SiteCode]
ymon_Rsdcs_lmstats =  merge(ymon_Rsdcs_lmstats,ymon_Rsdcs_rmse)

ymon_Rsdcs_lmstats[abs(slope1 - 1) < 0.1 & R2 > 0.95  , .(.N, range(slope1), range(R2), range(RMSE))]
ymon_Rsdcs_lmstats[RMSE < 10 & R2 > 0.95  , .(.N, range(slope1), range(R2), range(RMSE))]
ymon_Rsdcs_lmstats[RMSE < 20 & R2 > 0.95  , .(.N, range(slope1), range(R2), range(RMSE))]
# N        V2        V3        V4
# 1: 47 0.9209697 0.9685730  3.141172
# 2: 47 1.0886030 0.9990627 17.220763

ymon_Rsdcs_lmstats[ !(RMSE < 20 & R2 > 0.95)  , ]
ymon_Rsdcs_lmstats[ !(RMSE < 20 & R2 > 0.95)  , .(.N, range(slope1), range(R2), range(RMSE))]
ymon_Rsdcs_lmstats[order(R2), .(SiteCode,  RMSE, R2)]


### annual series and potential trends #####
## for high latitude sites with low winter values no missing values should be allowed !
## needs a grouping of sites to plot this
# dtyr =  dt30yrmon[ , lapply(.SD, meann, nmin = 11, na.rm = TRUE),  by = list(SiteCode, year)]
# dtyr =  dt30yrmon[ , lapply(.SD, base::mean),  by = list(SiteCode, year)]
# dtsite =  dtyr[ , lapply(.SD, meann, nmin = 2, na.rm = TRUE) ,  by = list(SiteCode)]

## IMPORTANT only use pairs of data of both data sets and remove strange values of QR

dt30yrmon[   !is.na(IncomingShortwaveClearSky - Swclidnz), ]
dt30yrmon[!is.na(IncomingShortwave)  &  !is.na(IncomingShortwaveClearSky - Swclidnz), ]
dtyr =  dt30yrmon[!is.na(IncomingShortwave)  &  !is.na(IncomingShortwaveClearSky - Swclidnz), lapply(.SD, base::mean),  by = list(SiteCode,year)]
dtyr
dtsite =  dtyr[ , lapply(.SD, meann, nmin = 2, na.rm = TRUE) ,  by = list(SiteCode)]
dtsite
dtsite[ , .(  rmse(IncomingShortwaveClearSky, Swclidnz), cor(IncomingShortwaveClearSky, Swclidnz)^2 )]
setkey(dtsite,SiteCode)

(dtsiteregSWcs = dcast(dt30yrmon[ , mlm.output.statlong.call("IncomingShortwaveClearSky ~ Swclidnz", data = .SD), by = list(SiteCode)], SiteCode ~ paste0("SWcsreg_",statistic)) )

(dtsiteRMSESWcs = dt30yrmon[ , .(RMSE_Rsdcs = rmse(IncomingShortwaveClearSky,Swclidnz) ), by = list(SiteCode)])

(dtsiteregtault1 = dcast(dt30yrmon[ftau < 1 & ftau_LA2000 < 1,  mlm.output.statlong.call("ftau ~ ftau_LA2000", data = .SD), by = list(SiteCode)], SiteCode ~ paste0("tault1reg_",statistic)) )


dtsitestats =  merge(dtsite,dtsiteregSWcs)
dtsitestats =  merge(dtsitestats,dtsiteRMSESWcs)
dtsitestats =  merge(dtsitestats,dtsiteregtault1)

### Order sites by their cloud effects
dtsite[ , .(SiteCode,ftau, SWCRE_QR, fCRE = SWCRE_QR/IncomingShortwave)][order(fCRE),]
dtsite[ , .(SiteCode,ftau, SWCRE_QR, fCRE = SWCRE_QR/IncomingShortwaveClearSky  )][order(SWCRE_QR),]

