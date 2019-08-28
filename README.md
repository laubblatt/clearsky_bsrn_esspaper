# R script to calculate and validate monthly clear sky fluxes at BSRN stations

**_By Maik Renner, Max-Planck Institute for Biogeochemistry, Jena, Germany_**
> 2019-08-26

This script contains all code to reproduce the analysis which is described in 
    
Renner, M., M. Wild, M. Schwarz, and A. Kleidon.
    **Estimating Shortwave Clear-Sky Fluxes from Hourly Global
    Radiation Records by Quantile Regression.**
    Earth and Space Science, 2019.
    https://doi.org/10.1029/2019EA000686. 
    
The script will merge the BSRN radiation data, perform aggregation to half-hourly values and apply the novel quantile regression approach to estimate the fractional clear-sky solar radiation from sub-daily records of global radiation.
The script also allows to reproduce the figures in the above mentioned manuscript. 

## Instructions    
### Install dependencies 
To use the script you need to install the related R package **cleaRskyQuantileRegression**, and the **phaselag** package from github:
```R
library(devtools)
install_github("laubblatt/cleaRskyQuantileRegression")
install_github("laubblatt/phaselag")
 ```

## Download the data from BSRN
Here I used a snapshot which contains all data in a different format. Note, that you need about 3.7GB space.

### Source Reference:
König-Langlo, Gert; Driemel, Amelie; Raffel, Bonnie; Sieger, Rainer (2015): BSRN snapshot 2015-09, links to zip archives. PANGAEA, https://doi.org/10.1594/PANGAEA.852720

```bash
cd YOUR_DATA_FOLDER 
wget http://hs.pangaea.de/Projects/BSRN/2015-09_snapshot/radiation.zip
unzip radiation.zip
```


## Merge data 
into a R data.table format using (R/read_BSRN_2datatable.r) 
You will need to set the variable *pdat* in this file to let R know where the data is stored in your system. 
```R
source("R/read_BSRN_2datatable.r") 
 ```

## Aggregate data 
BSRN data is stored in minutes. This script will process these and aggregate data into 30min and 60min values, using the aggregation scheme of Driemel et al., 2011. 
```R
source("R/read_BSRN_aggregate.r") 
 ```

## Perfrom the Quantile regression 
This will use the 30min and 60min data and derived the clear sky estimates for each month:

```R
source("R/clearsky_BSRN_quantreg.r") 
 ```

## Analyse output and create figures
This will compare with estimates of the quasi-standard *Long and Ackerman 2000 JGR* method. 
The data for comparison can be obtained from Prof. Martin Wild, ETH. 
```R
source("R/clearsky_BSRN_quantreg_figures.r") 
 ```


## for questions please open an issue in the project page 

