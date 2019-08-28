#' read BSRN snapshot 2015-09 into R data.table format 
#'
#' This reads in data from a snapshot provided on pangaea. This snapshot has a different format than the station files. 
#' Here only the radiation data is being used.
#' The script will merge data per site into one file. 
#' @references \url{http://hs.pangaea.de/Projects/BSRN/2015-09_snapshot/radiation.zip}
#' König-Langlo, Gert; Driemel, Amelie; Raffel, Bonnie; Sieger, Rainer (2015): BSRN snapshot 2015-09, links to zip archives. PANGAEA, https://doi.org/10.1594/PANGAEA.852720
#' @filename read_BSRN_2datatable.r
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de

#### ---------- download and extract data 
# cd /Net/Groups/C-Side/BTM/mrenner/scratch/data/fielddata/BSRN/
# ls
# wget http://hs.pangaea.de/Projects/BSRN/2015-09_snapshot/radiation.zip
# unzip radiation.zip
# filestructure: radiation/
#   radiation/XIA/XIA_radiation_2015-03.txt
# less radiation/XIA/XIA_radiation_2015-03.txt
#   
# grep -n "Height [m]"  radiation/LIN/LIN_radiation_1999-10.txt
# grep -n ^Date  radiation/LIN/LIN_radiation_1999-10.txt
# #### identify line number when data starts
# grep -n ^Date  radiation/LIN/LIN_radiation_1999-10.txt |cut -f1 -d:
  

#### ---------- merge data per site into one file 
# R
library(data.table)

#### SET PATH WHERE DATA is being downloaded 
pdat = "/Net/Groups/C-Side/BTM/mrenner/scratch/data/fielddata/BSRN/"

dtfiles =  data.table(pathfile = list.files(pdat,pattern = glob2rx("*_radiation*.txt"), recursive = TRUE))
dtfiles
dtfiles[ , c("type", "SiteCode", "filename") := tstrsplit(pathfile, split = "/")]
sicos = dtfiles[ , unique(SiteCode)]
sicos

## change colnames
# cona = colnames(dat)
# substr(cona,1,3)
# coneu = c("DateTime", "Height", "IncomingShortwave", "SWD_stddev", "SWD_min","SWD_max",
#  "DIR" ,        "DIR_stddev",  "DIR_min", "DIR_max" ,    "DIF", "DIF_stddev",
# "DIF_min",  "DIF_max",      "IncomingLongwave", "LWD_stddev",  "LWD_min" ,   "LWD_max",
#  "T2"    ,          "RelativeHumidity"  ,             "StationPressure" )

sico = "LIN"
setwd(pdat)
for (sico in sicos ) {
  dtbsrn = data.table()
  # datei = "radiation/LIN/LIN_radiation_2006-03.txt"
  for( datei in dtfiles[SiteCode == sico, pathfile]) {
    # identity how many lines I need to skip
    (linestart = as.numeric(system(paste0("grep -n ^Date " ,  datei , " |cut -f1 -d:"), intern = TRUE)))
    dat = fread(paste0(pdat,datei), skip = linestart-1)
    str(dat)
    cona = colnames(dat)
    cona = gsub("\\s*\\[[^\\)]+\\]", "" , cona)
    cona = gsub("/" , "" , cona)
    (cona = gsub(" " , "" , cona))
    
    setnames(dat,cona)
    
    dat[ , c("Date" , "Time") := IDateTime(as.POSIXct(DateTime, format = "%Y-%m-%dT%H:%M",tz = "GMT")) ]
    dat[ , SiteCode := sico ]
    ### write output
    wfields = c( "SWD", "LWD",
                 "DIR", "DIF", "T2", "RH", "PoPoPoPo" )
    # (cout = c(  wfields[wfields %in% cona])
    (cout = c("SiteCode", "Date", "Time",  wfields[wfields %in% cona]) )
    dtbsrn = rbind(dtbsrn, dat[ , cout, with = FALSE], fill = TRUE)
    print(paste0("Done ", datei, ", ", nrow(dtbsrn)))
    rm(dat)
  }
  print(paste0("############Done ", sico))
  save(dtbsrn, file = paste0(pdat,"rdata/", sico,".rdata" ) )
  
} ## done read all raw data into sitecode files

