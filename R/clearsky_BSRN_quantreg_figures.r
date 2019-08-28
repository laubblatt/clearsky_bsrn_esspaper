#' Calc quantile regression based on half-hourly BSRN data and compare with LA2000
#'

#' @filename clearsky_BSRN_quantreg_figures.r
#' @depends clearsky_BSRN_quantreg.r
#' @depends read_BSRN_aggregate.r
#' @depends read_BSRN_2datatable.r



#' 
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
#' @references  Renner, M., M. Wild, M. Schwarz, and A. Kleidon.
#'   "Estimating Shortwave Clear-Sky Fluxes from Hourly Global
#'     Radiation Records by Quantile Regression."
#'       Earth and Space Science, 2019.
#'       \url{https://doi.org/10.1029/2019EA000686}

#' @version 2019-08-28 This file is adapted from M44_ClearSky_figures.r  
#' @version 2019-08-28 code produces now figures and tables using the R package  # library(cleaRskyQuantileRegression)

#' @filename M44_ClearSky_figures.r
#' @depends M44_BSRN_hourly_quantreg_85.r
#' @depends read.BSRN.aggregate.r ## providing dt30 dt60 and dt15 aggregates
#' @depends read.bsrn.pangaea.snapshot.sh
#' 
#' @TODO DATA Dependencies 
#' load(file = paste0(prdata,"LIN_200308.rdata"))
#'  
#' 



##' packages ----
# library(tidyverse)
library(knitr)
# library(kableExtra)
library(data.table)
library(latticeExtra)
library(quantreg)
library(hexbin)
library(directlabels)
library(sp)



# ##' source functions ----
# pr = "~/Dissertation/maiphd/R/"
# source(paste(pr,"data.table.utils.R",sep=""))
# # source(paste(pr,"data.table.regression.fun.R",sep=""))
# source(paste(pr,"summary_fun.R",sep=""))
# source(paste(pr,"utils.ncdf.R",sep=""))
# source(paste(pr,"utils.rbash.R",sep=""))
# source(paste(pr,"axis_labels.R",sep=""))
# source(paste(pr,"fCalcPotRadiation_CosineResponsePower.R",sep=""))

#### I cannot install devtools on Rstudio server version
# install.packages("devtools")
# library(devtools)
# install_github("laubblatt/phaselag")

if (inherits(try(library(phaselag)),"try-error") ) { 
  print("library phaselag not available, use source instead")
  source("~/bgc/github/laubblatt/phaselag/R/data.table.regression.utils.R")
}
    
# install.packages("devtools")
# library(devtools)
# install_github("laubblatt/cleaRskyQunatileRegression")
library(cleaRskyQuantileRegression)
if (inherits(try(library(cleaRskyQuantileRegression)),"try-error") ) { 
  print("library cleaRskyQuantileRegression not available, use source instead")
  source("~/bgc/gitbgc/clearskyquantileregression/R/PotentialRadiation.R")
}



coordinates2degreeminute = function(x, ..., sec = FALSE) { 
  (zz = dd2dms(x,...))
  # str(zz)
  if (sec == FALSE)
    zz@sec <- 0
  sub("d","°",as.character(zz))
}


prdata = "~/bgc/gitbgc/clearsky_bsrn_esspaper/rdata/"
pfig = "~/bgc/gitbgc/clearsky_bsrn_esspaper/figures/"
pdataETHZ = "~/bgc/gitbgc/clearsky_bsrn_esspaper/dataETHZ/"


#### LOAD META DATA ####
### BRSN_state is the summary file of the Long and Ackerman 2000 method done by Maria Hakuba
# (BSRNmeta  = fread(paste0(pdat, "dataETHZ/BSRN_state.csv")) )
# BSRNmeta[ , SiteCode := toupper(site)]
# str(BSRNmeta)

# list.files(prdata)
# load(file = paste0(prdata,"dtconf.rdata"))
# dtconf
## get doubtful from comments in lost column of BSRN_state-txt
doubtful = c("ALE","ILO","SPO","SOV")
# dtconf[ , LA2000ok := TRUE]
# dtconf[daily == 0 , LA2000ok := FALSE]
# #dtconf[SiteCode %in% doubtful , doubtful := TRUE ]
# #dtconf[, doubtful := NULL ]
# dtconf[SiteCode %in% doubtful , LA2000ok := FALSE]
# 
# dtLA2000conf =  dtconf[, .(SiteCode,LA2000ok,daily)]
# fwrite(dtLA2000conf, file = paste0(prdata,"dtLA2000conf.csv" ))
dtconf =  fread(file = paste0(prdata,"dtLA2000conf.csv" ))

# merge(dtconf[, .(SiteCode,LA2000ok,daily)],BSRNmeta , by = "SiteCode")
# merge(dtconf[, .(SiteCode,LA2000ok,daily)],BSRNmeta , by = "SiteCode")[LA2000ok == FALSE, ]
# merge(dtconf[, .(SiteCode,LA2000ok,daily)],BSRNmeta , by = "SiteCode")[LA2000ok == TRUE, ]
# merge(dtconf[, .(SiteCode,LA2000ok,daily)],BSRNmeta[ , .(SiteCode,QCconf,QC)] , by = "SiteCode")

(BSRNStationMeta  = fread(paste0(prdata, "BSRN_Stations.csv")) )
BSRNStationMeta[ , URIofevent := NULL]
# tion no: 11; Surface type: tundra; Topography type: mountain valley, rural; Horizon: doi:10.1594/PANGAEA.669522; Station scientist: Marion.Maturilli@awi.de
BSRNStationMeta[ , Comment := NULL]
str(BSRNStationMeta)
setkey(BSRNStationMeta,SiteCode)
BSRNStationMeta[ , .(SiteCode,Elevation)]

(BSRNclearsky  = fread(paste0(pdataETHZ, "CSW_mm_Swclidnz_nighttime0_201709.dat"), na.strings = "-999.9") )
setnames(BSRNclearsky, c("V1","V2", "V15"), c("site", "year", "annual") )
BSRNclearsky[ , SiteCode := toupper(site)]
dtclearsky =  melt.data.table(BSRNclearsky, id.vars = c("SiteCode", "year", "site") , value.name = "Swclidnz" )
dtclearsky[ , month := as.numeric(substr(variable,2,4)) - 2]
dtclearsky[, variable := NULL]
dtclearsky[, site := NULL]
setkey(dtclearsky,SiteCode, year, month)
dtclearsky


#### load data 
# load(paste0(prdata,"BSRN_dt30.rdata"))

load(file = paste0(prdata,"qr30yrmon.rdata"))
load(file = paste0(prdata,"qr60yrmon.rdata"))
setkey(qr30yrmon,SiteCode,year,month)
setkey(qr60yrmon,SiteCode,year,month)

##################################################
#### this will merge with the dtclearsky data to get the Long and Ackerman data 
dt30yrmon =  merge( qr30yrmon, dtclearsky)
dt30yrmon[ ,  ftau_LA2000 := Swclidnz/IncomingShortwavePotential]
dt30yrmon[ , SWCRE_QR := IncomingShortwave - IncomingShortwaveClearSky]
save(dt30yrmon, file = paste0(prdata,"dt30yrmon.rdata"))


dt60yrmon =  merge( qr60yrmon, dtclearsky, by = c("SiteCode","year", "month"))
dt60yrmon[ ,  ftau_LA2000 := Swclidnz/IncomingShortwavePotential]
dt60yrmon[ , SWCRE_QR := IncomingShortwave - IncomingShortwaveClearSky]
save(dt60yrmon, file = paste0(prdata,"dt60yrmon.rdata"))

dt30yrmon = merge(dt30yrmon, dtconf[, .(SiteCode,LA2000ok,daily)], by = "SiteCode") 
dt60yrmon = merge(dt60yrmon, dtconf[, .(SiteCode,LA2000ok,daily)], by = "SiteCode") 

dtyr =  dt30yrmon[!is.na(IncomingShortwave)  &  !is.na(IncomingShortwaveClearSky - Swclidnz), lapply(.SD, base::mean),  by = list(SiteCode,year)]
dtyr
dtsite =  dtyr[ , lapply(.SD, meann, nmin = 2, na.rm = TRUE) ,  by = list(SiteCode)]
dtsite
dtsite[ , .(  rmse(IncomingShortwaveClearSky, Swclidnz), cor(IncomingShortwaveClearSky, Swclidnz)^2 )]
setkey(dtsite,SiteCode)

(dtsiteregSWcs = dcast(dt30yrmon[ , mlm.output.statlong.call("IncomingShortwaveClearSky ~ Swclidnz", data = .SD), by = list(SiteCode)], SiteCode ~ paste0("SWcsreg_",statistic)) )
(dtsiteRMSESWcs = dt30yrmon[ , .(RMSE_Rsdcs = rmse(IncomingShortwaveClearSky,Swclidnz) ), by = list(SiteCode)])

#' @version 2019-07-03 change site scale tau regression to include the above 1 estimates 
#' @TODO this should be done with the scatterplot as well 
(dtsiteregtau = dcast(dt30yrmon[,  mlm.output.statlong.call("ftau ~ ftau_LA2000", data = .SD), by = list(SiteCode)], SiteCode ~ paste0("taureg_",statistic)) )
dtsiteregtau[order(taureg_R2), .(SiteCode,taureg_R2,taureg_slope1,taureg_n)]
dtsitestats =  merge(dtsite,dtsiteregSWcs)
dtsitestats =  merge(dtsitestats,dtsiteRMSESWcs)
dtsitestats =  merge(dtsitestats,dtsiteregtau)

#### calc the number of sites with good stats and those not discussed in the text in sect 4.4 ####
(sitespoorLAok = dtconf[LA2000ok == FALSE, SiteCode])
sitespoordiscuss = c("TAM","PTR", "XIA")
sitespoor = c(sitespoorLAok, sitespoordiscuss)
length(sitespoor)
dtsitestats[ !(SiteCode %in%  sitespoor),  ][, .(r2min = min(taureg_R2), RMSE_Rsdcs_max = max(RMSE_Rsdcs)) ]
dtsitestats[ !(SiteCode %in%  sitespoor),  ][, .(r2mean = mean(taureg_R2), RMSE_Rsdcs_mean = mean(RMSE_Rsdcs)) ]


# load(file = paste0(prdata,"dtyrmon.rdata"))
## tau sensitivity output 
load(file = paste0(prdata,"dtyrmontau.rdata"))

### prepare the Table with all sites and the statistics 
colnames(dtsitestats)
#dtconf[ , .(SiteCode, daily, b, DiffRationNormSDlimit, Limit_difmax, Limit_NSW_max, LA2000ok)]
dtconf[ , table(LA2000ok)]
dtconf[LA2000ok == FALSE,  ]
#(dtsitestats = merge(dtsitestats,  dtconf[ , .(SiteCode, daily, b, DiffRationNormSDlimit, Limit_difmax, Limit_NSW_max, LA2000ok)] , by = "SiteCode"))
#dtsitestats[ , Location := NULL]
dtsitestats = merge(dtsitestats,BSRNStationMeta[, .(SiteCode,Location,SiteName,Longitude,Latitude)])

(dtsitestats_print =  dtsitestats[ , .(SiteCode,"Site Name" = substr(SiteName,1,12),Location, Lon = round(Longitude,1), Lat = round(Latitude,1),
                                       Rsd = round(IncomingShortwave), Rsdcs_QR = round(IncomingShortwaveClearSky), Rsdcs_LA = round(Swclidnz),
                                       RMSE = round(RMSE_Rsdcs,1),  
                                       beta_QR = round(ftau,2), beta_LA = round(ftau_LA2000,2),
                                       R2    = ifelse(LA2000ok == FALSE, NA, round(taureg_R2,2)), 
                                       slope = ifelse(LA2000ok == FALSE, NA,round(taureg_slope1,2)),
                                       n     = ifelse(LA2000ok == FALSE, NA,taureg_n),
                                       fCRE = round(-SWCRE_QR/IncomingShortwaveClearSky,2)  )][order(SiteCode),])

cat(kable(dtsitestats_print, format = "html"), sep = "\n", file = paste0(pfig,"dtsitestats_print.html"))

### requires kableExtra start <<<
dtsitestats_print[ , SiteCode := ifelse(dtconf[order(SiteCode), daily == 1], SiteCode, paste0(SiteCode, footnote_marker_alphabet(2, "html"))) ]
dtsitestats_print[SiteCode %in% doubtful , SiteCode := paste0(SiteCode, footnote_marker_alphabet(1, "html")) ]

### bold cells of slope when significant @version 2019-07-01 not used anymore ###
options(knitr.kable.na = '')
(dtsitestats_print_bold = dtsitestats_print%>%
    mutate(slope = ifelse(dtsitestats[order(SiteCode) , taureg_slope1_pvalue < 0.05], cell_spec(slope, bold = T), cell_spec(slope, bold = F) ) ) )

cat(kable(dtsitestats_print_bold, escape = FALSE, format = "html"), sep = "\n", file = paste0(pfig,"dtsitestats_print_bold.html"))

### add avg column 
(dtsitestats_avg =  dtsitestats_print[ , lapply(.SD, function(x) round(base::mean(x,na.rm = TRUE),2))])
(dtsitestats_avg[ , SiteCode := "Mean"] )
dtsitestats_avg[ , Lon := NA] 
dtsitestats_avg[ , Lat := NA] 
dtsitestats_avg[ , Rsd := round(Rsd) ]
dtsitestats_avg[ , Rsdcs_QR := round(Rsdcs_QR) ]
dtsitestats_avg[ , Rsdcs_LA := round(Rsdcs_LA) ]
dtsitestats_avg[ , RMSE := round(RMSE,1) ]
dtsitestats_avg[ , n := round(n) ]
dtsitestats_avg

dtsitestats_print_avg =  rbind(dtsitestats_print, dtsitestats_avg)

dtsitestats_print_avg_kable = 
  dtsitestats_print_avg %>% 
  kable(format = "html",escape = FALSE) %>%
  row_spec(55,bold = TRUE) %>%
  column_spec(2, width = "2cm") %>%
  column_spec(3, width = "3cm") %>%
  column_spec(4, width = "1.3cm") %>%
  column_spec(5, width = "1cm") %>%
  column_spec(6:9, width = "1cm") %>%
  add_header_above(c(" ", " "," "," "," "," ", "Shortwave fluxes" = 4, "Fractional transmission" = 5)," ")
dtsitestats_print_avg_kable
cat(dtsitestats_print_avg_kable, sep = "\n", file = paste0(pfig,"dtsitestats_print_avg_kable.html"))

## add footnotes 
dtsitestats_print_avg
dtsitestats_print_avg_kable_foot = 
  dtsitestats_print_avg %>%
  kable(format = "html",escape = FALSE) %>%
  kable_styling(font_size = 10) %>%
  row_spec(55,bold = TRUE) %>%
  column_spec(1, width = "1cm") %>%
  column_spec(2, width = "2cm") %>%
  column_spec(3, width = "3cm") %>%
  column_spec(4, width = "1.5cm") %>%
  column_spec(5, width = "1.2cm") %>%
  column_spec(6:9, width = "1cm") %>%
  column_spec(10:15, width = "1cm") %>%
  add_header_above(c(" ", " "," "," "," ", "Shortwave fluxes [W m^-2]" = 4, "Fractional transmission" = 6)," ") %>%
  footnote(alphabet = c("Long and Ackerman flagged doubtful", "Long and Ackerman approach in climatological mode"))

dtsitestats_print_avg_kable_foot
cat(dtsitestats_print_avg_kable_foot, sep = "\n", file = paste0(pfig,"dtsitestats_print_avg_kable_foot.html"))
### requires kableExtra end >>>



########################################
#### FIGURES ####
#### Axis labels #####
lab_Rsdpot12_name = expression(bold("Potential Shortwave Radiation "*R['sd,pot']==S[0]*cos(SZA)^1.2*" ")(W*m^-2))
lab_Rsdpot12_name = expression(bold("Potential Shortwave Radiation ")(W*m^-2))
lab_Rsdpot12 = expression(bold( R['sd,pot']==S[0]*cdot*cos(SZA)^1.2*" " )(W*m^-2))
lab_Rsd      = expression(bold("Observed Shortwave Radiation ")(W*m^-2))
lab_ftau = "Fractional Solar Transmission"
lab_ftau = "Fraction of potential Radiation"
lab_Rsdcs      = expression(bold("Clear-sky shortwave radiation "*R['sd,cs']*" ")(W*m^-2))
labhourday = expression(bold("Hour of day (UTC)"))
#### ILLUS QR ##### 
#### Figure 1 Illus QR #####
sico = "LIN"
dt30 = merge(dt30, BSRNStationMeta[ , list(SiteCode, lat = Latitude, lon = Longitude) ], by = "SiteCode", all.x = TRUE)

(dat = dt30[SiteCode == sico & year(Date) == 2003 & month(Date) == 8, ])

dat[ , Rsdpot_12 := calc_PotRadiation_CosineResponsePower(doy = yday(Date), hour = Time/3600 + 0.25,
                                                           latDeg = lat ,
                                                           longDeg = lon,
                                                           timeZone = 0, isCorrectSolartime = TRUE,
                                                           cosineResponsePower = 1.2 )]

dat[ , Rsdpot_tzlon := calc_PotRadiation_CosineResponsePower(doy = yday(Date), hour = Time/3600 + 0.25,
                                                              latDeg = lat ,
                                                              longDeg = lon,
                                                              timeZone = 0, isCorrectSolartime = TRUE,
                                                              cosineResponsePower = 1 )]



(gname = paste(pfig,"BSRN_Rsd2Rsdpot12_P85fit30min_",sico,"_200308.pdf",sep=""))
pdf(gname,6,5)
xyplot(IncomingShortwave ~  Rsdpot_12, data = dat,
       xlab = list(lab_Rsdpot12_name,cex = 1.3), ylab = list(label=lab_Rsd, cex=1.3),
       type = c("p","g"), pch = ".", cex = 3,
       # main = paste0(dtyrmon[SiteCode == sico, unique(SiteName)]," ,Quantile Regression, tau = 0.90" ),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.abline(rq(y~x, tau = 0.85), col=1, lwd = 2)
         # panel.abline(0,1, col=1, lty = 2, lwd = 2)
         panel.ablineq(0,1, col=1, lty = 2, lwd = 2, at = 0.89, label = "1:1", rotate = TRUE, fontfamily = "sans", cex = 1.5,pos = 3)
         panel.ablineq(rq(y~x, tau = 0.85), col=1, at = 0.7, rotate =TRUE,
                       pos = 3, label = "Slope of 85% Quantile", cex = 1.5, font = "Helvetica", fontface = 2 )
         panel.text(700,-20, "half hourly data of one month", fontfamily = "Helvetica")
       },
       grid = TRUE)
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))


#### FIGURE 7 all monthly fluxes, with LA2000 ####
#' @version 2019-07-03 update by removing poor sites here 
(gname = paste(pfig,"BSRN_RClearSky85fit30min_vs_LongAckerman_allmonth_hexbin.pdf",sep=""))
pdf(gname,6,6)
hexbinplot(IncomingShortwaveClearSky ~ Swclidnz, data = dt30yrmon[LA2000ok == TRUE, ], type =  c("p", "r","g"),
           xlab = expression(bold("Standard method (Long and Ackerman 2000) ")(W*m^-2)),
           ylab = expression(bold("Quantile Regression method, "*omega==0.85*"  ")(W*m^-2)),
           main = "Monthly mean Clear Sky Solar Radiation") +
  layer(panel.abline(0,1,lty = 2, col = "darkgrey")) +
  layer(panel.ablineq(lm(y ~ x), rotate = FALSE, x = 110, y = 30, pos = 4, cex = 1, r.sq = FALSE, digits = 2, col = 4) )
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))

dt30yrmon[LA2000ok == TRUE , rmse(Swclidnz, IncomingShortwaveClearSky)]
# [1] 8.600331 ## reported in text ###

MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = Rsdpot_12* 0.81 )
dt30yrmon[LA2000ok == TRUE , MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = IncomingShortwavePotential* 0.81 )]
# 0.7396669  reduction in residual variance by 74% reportent in text 

(gname = paste(pfig,"BSRN_RClearSky85fit60min_vs_LongAckerman_allmonth_hexbin.pdf",sep=""))
pdf(gname,6,6)
hexbinplot(IncomingShortwaveClearSky ~ Swclidnz, data = dt60yrmon[LA2000ok == TRUE, ], type =  c("p", "r","g"),
           xlab = expression(bold("Standard method (Long and Ackerman 2000) ")(W*m^-2)),
           ylab = expression(bold("Quantile Regression method, "*omega==0.85*"  ")(W*m^-2)),
           main = "Monthly mean Clear Sky Solar Radiation") +
  layer(panel.abline(0,1,lty = 2, col = "darkgrey")) +
  layer(panel.ablineq(lm(y ~ x), rotate = FALSE, x = 110, y = 30, pos = 4, cex = 1, r.sq = FALSE, digits = 2, col = 4) )
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))



### show the improvement against a fixed fraction estimate
(gname = paste(pfig,"BSRN_Rsdcs_diff_LA2000_ref_monthly.pdf",sep=""))
pdf(gname,6,6)
xyplot( I(IncomingShortwavePotential * 0.81 - Swclidnz)  + I(IncomingShortwaveClearSky - Swclidnz) ~ Swclidnz, data = dtyrmon,
        pch = ".", cex = 2,
        xlab = expression(bold("Clear-sky global radiation, Standard method ")(W*m^-2) ),
        ylab = expression(bold("Difference to Standard method ")(W*m^-2) ),
        auto.key = list(corner = c(0,0), text = c("Constant fractional transmission","Quantile regression approach"))
) + layer(panel.abline(h = 0))
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))

#### Figure 8 updated with a marginal distribution of the residuals 
(gname = paste(pfig,"BSRN_Rsdcs_diff_LA2000_ref_monthly_margdens.pdf",sep=""))
pdf(gname,7.8,5)
densref = dtyrmon[LA2000ok == TRUE ,  density(IncomingShortwavePotential * 0.81 - Swclidnz, na.rm = TRUE)]
densqr = dtyrmon[LA2000ok == TRUE ,  density(IncomingShortwaveClearSky - Swclidnz, na.rm = TRUE)]

ylims = c(-130,100)
ylims = c(-70,70)
op = par(mar=c(4.1, 4.1, 2, 0), las = 1, mgp = c(2,0.15,0), tck = 0.01)
par(fig=c(0,0.8,0,1), new=FALSE)
plot( I(IncomingShortwavePotential * 0.81 - Swclidnz)  ~ Swclidnz, data = dtyrmon[LA2000ok == TRUE,],
      pch = ".", cex = 2, ylim = ylims, col = "grey",
      xlab = expression(bold("Clear-sky shortwave radiation, Long and Ackerman method ")(W*m^-2) ),
      # xlab = "",
      ylab = expression(bold("Difference to Long and Ackerman method ")(W*m^-2) )
)
# opx = par(xpd = TRUE)
# mtext(text = expression(bold("Clear-sky shortwave radiation, Long and Ackerman method ")(W*m^-2) ),side = 1,)
# par(opx)
abline(h = 0, lty = 2)
grid()
points( I(IncomingShortwaveClearSky - Swclidnz)  ~ Swclidnz, data = dtyrmon[LA2000ok == TRUE,],pch = ".", cex = 2, ylim = ylims, col = 2 )
legend("bottomleft", c("Reference with fractional transmission","Quantile regression approach"), col = c("grey","red"), pch = ".", pt.cex = 4, bty = "n")

# auto.key = list(corner = c(0,0), text = c("Constant fractional transmission","Quantile regression approach")
axis(4, labels = NA, las = 2, tck = 0.01, mgp = c(1,0.25,0))
# par(fig=c(0.8,1,0,1), new=TRUE, mar=c(0, 6.1, 2.1, 8.1), xpd = TRUE)
# par(fig=c(0.5,1,0,1), new=TRUE)
# par(fig=c(0.8,1,0,1), new=TRUE, xpd = TRUE, mar = c(4.1, 0.1, 2, 0))
par(fig=c(0.8,1,0,1), new=TRUE, xpd = FALSE, mar = c(4.1, 0.1, 2, 0))
plot(densqr$y , densqr$x, ylab = "", xlab = "", ylim = ylims, type = "l", col = 2, lwd = 2, axes = F)
abline(h = 0, lty = 2)
lines(densref$y , densref$x, col = "grey", lwd = 2)
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))




#### Figure 2 effect of exponent b ####
load(file = paste0(prdata,"LIN_200308.rdata"))
dat1 = dtbsrn[year(Date) == 2003 & month(Date) == 8 & mday(Date) == 11, ]

# lab_Rsdpot12_name = expression(bold("Potential Shortwave Radiation ")(W*m^-2))
# lab_Rsdpot12 = expression(bold( R['sd,pot']==S[0]*cdot*cos(SZA)^1.2*" " )(W*m^-2))
# lab_Rsd      = expression(bold("Observed Shortwave Radiation ")(W*m^-2))
#### Figure 2a
(gname = paste(pfig,"BSRN_Rsd2RsdpotExpo_clearskyday_",sico,".pdf",sep=""))
pdf(gname,5,5)
par(mar = c(4,5,3,2), las = 1, mgp = c(1.5,0.15,0), tck = 0.01,pty = "sq")
lim = c(0,1100)
plot(IncomingShortwave  ~ Rsdpot_tzlon , data = dat1, type = "p", pch = ".", col = "lightblue", cex = 2, ylim = lim, xlim = lim, xlab = lab_Rsdpot12_name, ylab = lab_Rsd) 
points(IncomingShortwave  ~ Rsdpot_12, data = dat1, type = "p", pch = ".", col = "#fdcc8a"  , cex = 4) 
fit12 = lm(IncomingShortwave  ~ Rsdpot_12, data = dat1)
summary(fit12)
abline(fit12, col = "#b30000", lwd = 2, lty = 2)
grid()
fit = lm(IncomingShortwave  ~ Rsdpot_tzlon, data = dat1)
summary(fit)
abline(fit, col = 4, lwd = 2, lty = 4)
anova(fit,fit12)
sum(resid(fit)^2)
legend("topleft", c(expression(R['sd,pot']==S[0]*cos(theta)^bold(1.2)), expression(R['sd,pot']==S[0]*cos(theta)^1.0) ), lty = c(2,4), lwd = 2, col =  c(2,4), bty = "n", inset = c(0,0.1))
# legend("topleft", c(expression(R['sd,pot']==S[0]*cos(SZA)^bold(1.2)), expression(R['sd,pot']==S[0]*cos(SZA)) ), lty = 2, col =  c(2,4), bty = "n", inset = c(0,0.1))
abline(0,1, lty = 3, col = "grey")
legend("topleft", "a)" , bty = "n", cex = 1.5, text.font = 2, inset = c(-0.07,-0.02))
gg = dev.off()
system(paste("pdfcrop ",gname,gname))


### plot time course 
(gname = paste(pfig,"BSRN_Rsd_RsdpotExpo_Time_clearskyday_",sico,".pdf",sep=""))
pdf(gname,5,5)
par(mar = c(4,5,3,2), las = 1, mgp = c(1.5,0.15,0), tck = 0.01,pty = "sq")
plot(IncomingShortwave  ~ I(Time/3600) , data = dat1, type = "p", pch = ".", xlab = labhourday, col = "grey", cex = 5, ylab = lab_Rsd, xaxt = "n") 
axis(1,at = seq(0,24,6))
# add the clear-sky estimate without exponent 
lines( I( coef(fit)[2]*Rsdpot_tzlon) ~ I(Time/3600) , data = dat1, col = 4, lty = 4, lwd = 2) 
lines( I( coef(fit12)[2]*Rsdpot_12) ~ I(Time/3600) , data = dat1, col = 2, lty = 2, lwd = 2) 

lines( I( coef(fit)[2]*Rsdpot_tzlon - IncomingShortwave) ~ I(Time/3600) , data = dat1, col = 3, lty = 1, lwd = 2) 
lines( I( coef(fit12)[2]*Rsdpot_12- IncomingShortwave) ~ I(Time/3600) , data = dat1, col = 3, lty = 1, lwd = 2) 
lines( I( coef(fit)[2]*Rsdpot_tzlon - IncomingShortwave) ~ I(Time/3600) , data = dat1, col = 4, lty = 4, lwd = 2) 
lines( I( coef(fit12)[2]*Rsdpot_12- IncomingShortwave) ~ I(Time/3600) , data = dat1, col = 2, lty = 2, lwd = 2) 

# legend("bottom", c("Clear Sky fit, b = 1.2", "Clear Sky fit, b = 1.0", "Observation"), lty = c(2,2,NA), pch = c(NA,NA, 15), col =  c(2,4,"grey"), bty = "n")
legend("topright", c(expression(R['sd']*"~"*cos(theta)^bold(1.2)), expression(R['sd']*"~"*cos(theta)^bold(1.0)), "Observation"), lty = c(2,4,NA), lwd = c(2,2,NA), pch = c(NA,NA, 15), col =  c(2,4,"grey"), bty = "n")
# legend("topright", c("fit, b = 1.2", "Clear Sky fit, b = 1.0", "Observation"), lty = c(2,2,NA), pch = c(NA,NA, 15), col =  c(2,4,"grey"), bty = "n")
legend("topleft", "b)" , bty = "n", cex = 1.5, text.font = 2, inset = c(-0.07,-0.02))
# title("Modelled Clear Sky Flux")
abline(h = seq(0,1000,200), lty = 3, col = "grey")
abline(v = seq(0,24,6), lty = 3, col = "grey")
gg = dev.off()
system(paste("pdfcrop ",gname,gname))

tz(dat1$Time)



##### Figure 5 Scatterplot frac ######
#' Figure 5 Comparison of monthly fractional clear-sky solar transmission obtained with the proposed QR method and the Standard method.
#' Both methods reveal estimates larger than 1 which are deemed unreliable and colored as grey dots. 
#' Blue dots and the blue linear best fit line are derived from values lower than 1. The dashed line shows the 1-1 perfect fit. 
#'
#' @version 2019-07-02 remove points from comparison when LA2000ok 
#' @version 2019-07-03 do not exclude ftau > 1 since this corrupts site scale ostats in table 1 

(gname = paste(pfig,"BSRN_fractrans85fit30min_vs_LongAckerman_allmonth.pdf",sep=""))
pdf(gname,5,5)
lims = c(0.4,1.3)
op = par(mar = c(4,5,3,2), las = 1, mgp = c(1.5,0.15,0), tck = 0.01,pty = "sq")
plot(ftau ~ ftau_LA2000, data = dt30yrmon[SiteCode %in% dtconf[LA2000ok==TRUE,SiteCode], ], type =  c("p"),
     xlab = expression(bold("Standard method (Long and Ackerman 2000) ")),
     ylab = expression(bold("Quantile Regression method  ")),
     xlim = lims, ylim = lims, col = 4, pch = "." )
# xlim = lims, ylim = lims, col = rgb(.5,.5,.5), pch = "." )
# points(ftau85dt30 ~ ftau_LA2000, data = dtyrmon[ftau85dt30 < 1 & ftau_LA2000 < 1, ][SiteCode %in% dtconf[LA2000ok==TRUE,SiteCode], ], type =  c("p"),
#        xlim = lims, ylim = lims, col = 4, pch = "." )
title("Fractional clear-sky solar transmission", line = 0.3)
grid()
abline(0,1,lty = 2, col = "darkgrey")
fitall = lm(ftau ~ ftau_LA2000, data = dt30yrmon[SiteCode %in% dtconf[LA2000ok==TRUE,SiteCode], ])
summary(fitall)
abline(fitall, col = 4)
sufiall = summary(fitall)
(textfitall =  bquote(y == .(sprintf('%.2f', sufiall$coef[1,1])) + .(sprintf('%.2f', sufiall$coef[2,1]))*x*", "* italic(r)^2*"="*.(sprintf('%.2f', sufiall$adj.r.squared) )*", n = "*.(length(sufiall$residuals)) ))
legend("bottom", as.expression(textfitall), col = 4, lty = 1, bty = "n")
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))

(gname = paste(pfig,"BSRN_fractrans85fit60min_vs_LongAckerman_allmonth.pdf",sep=""))
pdf(gname,5,5)
lims = c(0.4,1.3)
op = par(mar = c(4,5,3,2), las = 1, mgp = c(1.5,0.15,0), tck = 0.01,pty = "sq")
plot(ftau ~ ftau_LA2000, data = dt60yrmon[SiteCode %in% dtconf[LA2000ok==TRUE,SiteCode], ], type =  c("p"),
     xlab = expression(bold("Standard method (Long and Ackerman 2000) ")),
     ylab = expression(bold("Quantile Regression method using hourly data ")),
     xlim = lims, ylim = lims, col = 4, pch = "." )
title("Fractional clear-sky solar transmission", line = 0.3)
grid()
abline(0,1,lty = 2, col = "darkgrey")
fitall = lm(ftau ~ ftau_LA2000, data = dt60yrmon[SiteCode %in% dtconf[LA2000ok==TRUE,SiteCode], ])
summary(fitall)
abline(fitall, col = 4)
sufiall = summary(fitall)
(textfitall =  bquote(y == .(sprintf('%.2f', sufiall$coef[1,1])) + .(sprintf('%.2f', sufiall$coef[2,1]))*x*", "* italic(r)^2*"="*.(sprintf('%.2f', sufiall$adj.r.squared) )*", n = "*.(length(sufiall$residuals)) ))
legend("bottom", as.expression(textfitall), col = 4, lty = 1, bty = "n")
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))


##### map of the sites #####
library(maptools) # pointLabel positioning
### thats it! BEST land sea mask
library(rgdal)
# install.packages("rworldmap")
library(rworldmap)
data(coastsCoarse)

dtmetamap = copy(dtsitestats)
coordinates(dtmetamap) <- c("Longitude", "Latitude")
proj4string(dtmetamap) <- proj4string(coastsCoarse)
# proj4string(dtmetamap) <- proj4string(tz_world)
str(dtsite)
dtsite[ , summary(IncomingShortwaveClearSky)]
dtsite[ , IncomingShortwaveClearSky_cut := cut(IncomingShortwaveClearSky, seq(0,350,50))]
dtsite[ , levels(IncomingShortwaveClearSky_cut)]
nl = dtsite[ , nlevels(IncomingShortwaveClearSky_cut)]
RsdcsCols = colorRampPalette(brewer.pal(7,"Reds"))(nl)
RsdcsCols = colorRampPalette(brewer.pal(7,"PuRd"))(nl)
dtsite[ , RsdcsCols[ as.numeric(IncomingShortwaveClearSky_cut)]]

dtsite[ , summary(ftau)]
dtsite[ , ftau_cut := cut(ftau, c(0.6,seq(0.7,0.95,0.05),1.1))]
dtsite[ , levels(ftau_cut)]
nl = dtsite[ , nlevels(ftau_cut)]
RsdcsCols = rev(colorRampPalette(brewer.pal(7,"YlOrBr"))(nl))
dtsite[ , RsdcsCols[ as.numeric(ftau85dt30_cut)]]

#### Figure Map of the sites with ftau color ####
(gname = paste(pfig,"BSRN_worldmap_ftau.pdf",sep=""))
pdf(gname,10,7)
plot(coastsCoarse, ylim = c(-65,90), col = "grey")
points(dtmetamap, col = 1 , pch = 0)
points(dtmetamap, col = dtsite[ , RsdcsCols[ as.numeric(ftau_cut)]], pch = 15)
pointLabel(dtmetamap$Longitude, dtmetamap$Latitude, dtmetamap$SiteCode, cex = 0.8)
(leg = dtsitestats[order(Latitude, decreasing = TRUE) , paste(SiteCode, SiteName, sep = ", ")])
# (leg = dtsite[order(lat, decreasing = TRUE) , paste(SiteCode, SiteName, round(IncomingShortwaveClearSky), round(Swclidnz), sep = ", ")])
# legend("bottomleft", leg , col = 2 , pch = 0, cex = 0.5, bty = "n", inset = c(0,0.021))
legend("bottomleft", dtsite[ , levels(ftau_cut)] , col = RsdcsCols , pch = 15, cex = 1,
       bty = "n", inset = c(0.15,0.1), title = expression(bold(R['sd']/R['sd,pot'])) )

gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))



(gname = paste(pfig,"BSRN_worldmap.pdf",sep=""))
pdf(gname,16,10)
plot(coastsCoarse, ylim = c(-65,90), col = "grey")
points(dtmetamap, col = dtsite[ , RsdcsCols[ as.numeric(IncomingShortwaveClearSky_cut)]], pch = 15)
pointLabel(dtmetamap$Longitude, dtmetamap$Latitude, dtmetamap$SiteCode, cex = 0.8)
(leg = dtsitestats[order(Latitude, decreasing = TRUE) , paste(SiteCode, SiteName, sep = ", ")])
# (leg = dtsite[order(lat, decreasing = TRUE) , paste(SiteCode, SiteName, round(IncomingShortwaveClearSky), round(Swclidnz), sep = ", ")])
legend("bottomleft", leg , col = 2 , pch = 15, cex = 0.5, bty = "n", inset = c(0,0.021))
legend("bottomleft", dtsite[ , levels(IncomingShortwaveClearSky_cut)] , col = RsdcsCols , pch = 15, cex = 1,
       bty = "n", inset = c(0.15,0.1), title = expression(bold(R['sd,cs']*" ")(W*m^-2)) )

gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))


(gname = paste(pfig,"BSRN_ClearSky12_vs_LongAckerman_AVG.pdf",sep=""))
pdf(gname,6,6)
lims = c(150,350)
op = par(mar = c(4,5,3,2), las = 1, mgp = c(2,0.15,0), tck = 0.01,pty = "sq")
plot(IncomingShortwaveClearSky ~ Swclidnz, data = dtsite, xlim = lims, ylim =lims,
     xlab = expression(bold("Standard method based on Long and Ackerman (2000)")),
     ylab = expression(bold("Quantile Regression method, "*omega==0.85)),
     main = "Site average Clear Sky Solar Radiation")
abline(0,1, lty = 2)
grid()
# use the full site data to estimate the long term mean
# text(IncomingShortwaveClearSky ~ Swclidnz, data = dtsite, label = SiteCode)
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))


#### fractional solar transmission time series ####
dt30yrmon[ , cor(ftau,ftau_LA2000, use = "pair")^2] 

### FIGURE S3 ####
#' @version 2019-07-02 update and annotate problem sites in panel-label 
#' 
# dtyrmon[ , doubtful := NULL]

dtyrmon4S3 = copy(dt30yrmon)
dtyrmon4S3[, SiteCodeOnPanel := SiteCode]
dtyrmon4S3[daily == 0 , SiteCodeOnPanel :=  paste(SiteCode, '(b)')  ]
dtyrmon4S3[SiteCode %in% doubtful, SiteCodeOnPanel :=  paste(SiteCode, '(a)')  ]

(gname = paste(pfig,"BSRN_FractionalSolarTransmission12_ftau85dt30_vs_LongAckerman_timeseries.pdf",sep=""))
pdf(gname,8,12)
xyplot(ftau + ftau_LA2000 ~ I(year + month/12) | paste(SiteCodeOnPanel),
       data = dtyrmon4S3, ylim = c(0.6,1.3), lwd = c(2,1),
       type =  c("l","g") , xlab = list(label = "Time", cex = 1.3),
       ylab = list(label = "Fractional Solar Transmission", cex = 1.3),
       scales=list(tck=c(1,0), x=list(cex=1), y=list(cex=1)),
       as.table = TRUE,
       auto.key = list("top", columns = 2, points=FALSE, lines=TRUE, text = c("Quantile regression slope",expression("Long and Ackermann (2000) "*R['sd,cs']/R['sd,pot'])), cex = 1)
)
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))

#### Figure 6 beta series for selected sites ####
sics = c("PAY",  "E13", "BER", "MAN", "GVN")
sics = c("PAY",  "E13", "TAM", "DAR", "GVN")
ii = 1
colfilled = "magenta"
colfilled = "purple"

(gname = paste(pfig,"BSRN_ftau_timeseries_",paste(sics, collapse = "_"),".pdf",sep=""))
pdf(gname,5,7)
op = par(mar = c(1,3,1.2,1), las = 1, mgp = c(2,0.15,0), tck = 0.01, mfrow = c(5,1))
ylims = c(0.72,0.95)
ii = 1
for (sico in sics) {
  dat = BSRNStationMeta[dt30yrmon][SiteCode == sico , ]
  if (dat[ , unique(Latitude) < -50 ]) 
    (ylims = dat[ , range(c(ftau,ftau_LA2000),na.rm=TRUE)] + c(-0.02, 0.02))
  
  if (ii == 3) { ylab = expression(bold("Fractional transmission")) 
  } else ylab = ""
  plot(ftau_LA2000 ~ I(year + month/12),  data = dat , type = "l" , lwd = 2.5 , xlab = "", ylab = ylab, ylim = ylims, xlim = c(1992,2015))
  # plot(ftau_LA2000 ~ I(year + month/12),  data = dat , type = "l" , lwd = 2.5 , xlab = "", ylab = lab_ftau, ylim = ylims)
  grid()
  lines(ftau ~ I(year + month/12),  data = dat , col = 2, lwd = 1.5  )
  points(ftau ~ I(year + month/12),  data = dat[Window != "1mon" , ] , col = colfilled, pch = 3 , lwd = 2 )
  # title(paste(sico, dat[, unique(SiteName)], dat[!is.na(lon), unique(lon)], dat[!is.na(lat), unique(lat)]), line = 0.3)
  title(paste(sico, dat[, unique(SiteName)], dat[!is.na(Latitude), coordinates2degreeminute(unique(Latitude),NS =TRUE)], dat[!is.na(Longitude), coordinates2degreeminute(unique(Longitude))], sep = ", "), line = 0.3)
  if (ii == 4)
    legend("bottomleft", c("Long and Ackerman method" , "Quantile Regression approach", "Sampling window > 1 month"), col = c(1,2,colfilled), lty = c(1,1,NA), pch = c(NA,NA,3), bty = "n")
  # legend("topleft", c("Standard Method" , "Quantile Regression approach"), col = c(1,2), lty = 1, bty = "n")
  legend("topleft", paste0(letters[ii],")") , bty = "n", cex = 1.2, text.font = 2, inset = c(-0.02,-0.02))
  ii = ii + 1
}
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))


(gname = paste(pfig,"BSRN_FractionalSolarTransmission12_ftau85dt30_vs_LongAckerman_timeseries.pdf",sep=""))
pdf(gname,20,12)
xyplot(ftau + ftau_LA2000   ~ I(year + month/12) | paste(SiteCode, substr(SiteName,1,12)),
#       data = dtyrmon[SiteCode %in% sics, ], ylim = c(0.45,1.3),
       data = BSRNStationMeta[ dt30yrmon ], ylim = c(0.45,1.3),
       type =  c("l","g") , xlab = list(label = "year and month", cex = 1.5),
       ylab = list(label = "Fractional Solar Transmission", cex = 1.3),
       scales=list(tck=c(1,0), x=list(cex=1.2), y=list(cex=1.2)),
       par.strip.text=list(cex=1.3),
       auto.key = list("top", columns = 2, text = c("Quantile regression slope","SWclear/SWpot"), cex = 1.3)
)
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))


(gname = paste(pfig,"BSRN_SWClearSky12_P85fit30min_vs_LongAckerman_timeseries.pdf",sep=""))
pdf(gname,20,12)
xyplot(IncomingShortwaveClearSky + I(Swclidnz)  ~ I(year + month/12) | paste(SiteCode, substr(SiteName,1,12)),
       data = BSRNStationMeta[ dt30yrmon ], lwd = c(1.5,0.8),
       type =  c("l","g") , xlab = list(label = "year and month", cex = 1.5),
       ylab = list(label =lab_Rsdcs , cex = 1.3),
       scales=list(tck=c(1,0), x=list(cex=1.2), y=list(cex=1.2)),
       par.strip.text=list(cex=1.3),
       auto.key = list("top", columns = 2, text = c("Quantile regression Approach","Long and Ackerman 2000"), cex = 1.3)  )
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))



#### Sensitivity study on tau ######
#' @update 2019-03-06 tropical sites without seasonal variation in tau in standrad methdo and high latitude sites with outliers excluded 
#' @update 2019-03-07 these tropical sites still show variation, but really small, so I again include all 
#' @version 20190702 update only LA2000ok
#' 
dtyrmontau[dt30rq_R1 > 0.75 , .N, by = list(SiteCode,year,month)]
dtyrmontau_R1set = dtyrmontau[dt30rq_R1 > 0.75 , sum(dt30rq_R1 > 0.75) == 30 , by = list(SiteCode,year,month)][V1 == TRUE, ][ , V1:=NULL]
#merge
dtyrmontau[dtyrmontau_R1set]


# (dtyrmontau_skilltau =  dtyrmontau[dt30rq_R1 > 0.75 , MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = Rsdpot_12* 0.81 ), by = tau])
# (dtyrmontau_skilltau =  dtyrmontau[dtyrmontau_R1set][dt30rq_R1 > 0.75 , .(SWcs_MSES =  MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = Rsdpot_12* 0.81),
#                                                                           FST_r2 = cor(dt30rq_slope1, Swclidnz/Rsdpot_12, use = "pair")^2
# ), by = tau])

dtconf[LA2000ok==TRUE,SiteCode]

(dtyrmontau_skilltau =  dtyrmontau[dtyrmontau_R1set][SiteCode %in% dtconf[LA2000ok==TRUE,SiteCode], ] [dt30rq_R1 > 0.75 , .(SWcs_MSES =  MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = Rsdpot_12* 0.81),
                                                                          FST_r2 = cor(dt30rq_slope1, Swclidnz/Rsdpot_12, use = "pair")^2
), by = tau])

dtyrmontau_skilltau[ , skillsum := SWcs_MSES + FST_r2]
dtyrmontau_skilltau
dtyrmontau_skilltau[tau == 0.85, ]

xyplot(SWcs_MSES ~ tau, data = dtyrmontau_skilltau)

(gname = paste(pfig,"BSRN_skilltau.pdf",sep=""))
pdf(gname,6,6)
op = par(mar = c(4,5,3,2), las = 1, mgp = c(1.5,0.15,0), tck = 0.01,pty = "sq")
plot(SWcs_MSES ~ tau, data = dtyrmontau_skilltau, xlab = expression(bold("Quantile "*omega*" used in Quantile Regression" )), ylab = "Skill score", ylim = c(-0.2,1), pch = 0 )
legend("bottomleft", c("Explained Variance fractional transmission", "Mean Squared Error Skill Score"), col = c(2,1), pch = c(1,0), bty = "n")
points(FST_r2 ~ tau, data = dtyrmontau_skilltau, col = 2)
grid()
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))

dtyrmontauwide = dcast(dtyrmontau[dt30rq_R1 > 0.75 , ], SiteCode + year + month + Rsdpot_12 ~ tau, value.var = c("dt30rq_slope1", "IncomingShortwaveClearSky"))
dtyrmontauwide[ , SiteCodeRev := factor(SiteCode, levels = sort(unique(SiteCode), decreasing = TRUE))]
(gname = paste(pfig,"BSRN_tausensi_onRsdcs_boxplot.pdf",sep=""))
pdf(gname,5,8)
par(mar = c(4,5,3,2), las = 1, mgp = c(1.5,0.15,0), tck = 0.01)
bp = boxplot( I(IncomingShortwaveClearSky_0.9 - IncomingShortwaveClearSky_0.85 )  ~ SiteCodeRev , data = dtyrmontauwide, horizontal = TRUE,yaxt="n", xlab = expression(bold("Change in monthly mean clear-sky flux ")(W*m^-2)))
axis(2, at = 1:54, labels =  bp$names , cex.axis=0.7)
grid()
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))

#' @version 20190624 combine the two plots in a multi panel
#' Figure 4 Sensitivity of the quantile regression method to the choice of the quantile. 
#' Panel a) Skill scores across all sites and months using the estimates of the standard method for evaluation. 
#' Red shows the squared correlation of the fractional solar transmission and black shows the skill score with a constant fractional solar transmission as reference. 
#' Panel b) shows the difference in the monthly estimate of the clear-sky shortwave flux when the quantile ??? = 90% is compared to  ??? = 85%.
#'   The sensitivities only cover months where the quantile regressions yielded a R1>0.75.  
#' 
(gname = paste(pfig,"BSRN_skilltau_tausensi_panels.pdf",sep=""))
pdf(gname,8,6)
nf <- layout(matrix(1:2,1,2), widths = c(5,3.5), heights = c(5))
# layout.show(nf)
# op = par(mar = c(4,5,3,2), las = 1, mgp = c(1.5,0.15,0), tck = 0.01)
op = par(mar = c(4,4,3,0.4), las = 1, mgp = c(1.5,0.15,0), tck = 0.01)
plot(SWcs_MSES ~ tau, data = dtyrmontau_skilltau, xlab = expression(bold("Quantile "*omega*" used in Quantile Regression" )), ylab = "Skill score", ylim = c(-0.2,1), pch = 0 )
legend("bottomleft", c("Explained Variance fractional transmission", "Mean Squared Error Skill Score"), col = c(2,1), pch = c(1,0), bty = "n")
points(FST_r2 ~ tau, data = dtyrmontau_skilltau, col = 2)
grid()
legend("topleft",expression(bold("a)")), bty = "n", cex = 1.5, inset = c(-0.04,-0,01))

bp = boxplot( I(IncomingShortwaveClearSky_0.9 - IncomingShortwaveClearSky_0.85 )  ~ SiteCodeRev , data = dtyrmontauwide, horizontal = TRUE,yaxt="n",range = 0,
              xlab = expression(bold("Clear-sky flux difference ")(W*m^-2)))
axis(2, at = 1:54, labels =  bp$names , cex.axis=0.4)
grid()
legend("topright",expression(bold("b)")), bty = "n", cex = 1.5, inset = c(0,0))
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))



(dtyrmontau_yrmonreg = dcast(dtyrmontau[dtyrmontau_R1set][ , mlm.output.statlong.call("IncomingShortwaveClearSky ~ Swclidnz", .SD), by = tau ], tau ~ statistic))

(gname = paste(pfig,"BSRN_tau_SWcs_xy.pdf",sep=""))
pdf(gname,7,7)
op = par(mar = c(4,5,3,2), las = 1, mgp = c(1.5,0.15,0), tck = 0.01,pty = "sq")
plot(IncomingShortwaveClearSky ~ Swclidnz, data = dtyrmontau[dtyrmontau_R1set][round(tau,2) == 0.90, ], pch = ".")
abline(0,1,lty = 2)
abline(lm(IncomingShortwaveClearSky ~ Swclidnz, data = dtyrmontau[dtyrmontau_R1set][round(tau,2) == 0.90, ]), col = 1)
abline(lm(IncomingShortwaveClearSky ~ Swclidnz, data = dtyrmontau[dtyrmontau_R1set][round(tau,2) == 0.80, ]), col = 2)
abline(lm(IncomingShortwaveClearSky ~ Swclidnz, data = dtyrmontau[dtyrmontau_R1set][round(tau,2) == 0.95, ]), col = 4)
#points(FST_r2 ~ tau, data = dtyrmontau_skilltau, col = 2)
grid()
# xyplot(slope1 + R2adj ~ tau, data = dtyrmontau_yrmonreg)
gg = dev.off()
system(paste("pdfcrop ",gname,gname))


### the skill score is quite sensitive to outliers in the slopes, so I use the R1
# (dtyrmontau_skilltausite =  dtyrmontau[dt30rq_R1 > 0.75 , MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = Rsdpot_12* 0.81 ), by = list(SiteCode,tau)])

(dtyrmontau_skilltausite =  dtyrmontau[dtyrmontau_R1set][dt30rq_R1 > 0.75 , .(MSES = MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = Rsdpot_12* 0.81 ),
                                                                              SWcs_cor = cor(Swclidnz,IncomingShortwaveClearSky, use = "pair")^2 ,
                                                                              FST_r2 = cor(dt30rq_slope1, Swclidnz/Rsdpot_12, use = "pair")^2)   , by = list(SiteCode,tau)])

(dtyrmontau_skilltausite =  dtyrmontau[dtyrmontau_R1set][SiteCode %in% dtconf[LA2000ok==TRUE,SiteCode], ][dt30rq_R1 > 0.75 , .(MSES = MSE_SkillScore(o = Swclidnz, p = IncomingShortwaveClearSky, ref = Rsdpot_12* 0.81 ),
                                                                              SWcs_cor = cor(Swclidnz,IncomingShortwaveClearSky, use = "pair")^2 ,
                                                                              FST_r2 = cor(dt30rq_slope1, Swclidnz/Rsdpot_12, use = "pair")^2)   , by = list(SiteCode,tau)])


# dtyrmontau_skilltausite[ , which.max(V1), by = SiteCode]

## LER and DAA are very negative 
dtyrmontau_skilltausite[SiteCode == "LER", ]
dtyrmontau_skilltausite[SiteCode == "DAA", ]

## Supplement Figure S1 ####
(gname = paste(pfig,"BSRN_skilltau_site.pdf",sep=""))
pdf(gname,8,10)
xyplot(MSES ~ tau | paste(SiteCode), data = merge(dtyrmontau_skilltausite, BSRNStationMeta, by = "SiteCode"),  as.table = TRUE,
       xlab = expression(bold("Quantile "*omega*" used in Quantile Regression")), ylab =  "Mean Squared Error Skill Score",
       ylim = c(-0.5,1.02), axis = axis.grid) + layer(panel.abline(v=0.85, col = 2))
gg = dev.off()
system(paste("pdfcrop ",gname,gname))
system(paste("pdftoppm -singlefile -rx 300 -ry 300 -png ",gname, gname))

(gname = paste(pfig,"BSRN_SWcscor_tau_site.pdf",sep=""))
pdf(gname,8,12)
xyplot(SWcs_cor ~ tau | paste(SiteCode, substr(SiteName,1,12)), data = merge(dtyrmontau_skilltausite, BSRNStationMeta, by = "SiteCode"), xlab = expression(bold("Quantile used in Quantile Regression "*tau)), ylab =  "Explained Variance ", ylim = c(0.5,1.02), axis = axis.grid) + layer(panel.abline(v=0.86, col = 2))
gg = dev.off()
system(paste("pdfcrop ",gname,gname))

(gname = paste(pfig,"BSRN_FST_r2tau_site.pdf",sep=""))
pdf(gname,15,10)
xyplot(FST_r2 ~ tau | paste(SiteCode, substr(SiteName,1,12)), data = merge(dtyrmontau_skilltausite, BSRNStationMeta, by = "SiteCode"), xlab = expression(bold("Quantile used in Quantile Regression "*tau)), ylab =  "Explained Variance of fractional transmission ", ylim = c(0.5,1.02), axis = axis.grid) + layer(panel.abline(v=0.85, col = 2))
gg = dev.off()
system(paste("pdfcrop ",gname,gname))



