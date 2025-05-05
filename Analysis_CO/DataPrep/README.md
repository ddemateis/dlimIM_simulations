# Data Sharing for Colorado Data

All data for the Colorado analysis is publicly available. Some data requires a data use agreement and IRB approval. Below is the infromation on each data soure and instructions for data processing.

All data processing scripts are in the folder Analysis_CO/DataPrep.

## Data Sources and Pre-Processing

### Health data
- Colorado vital statistics data was obtained from the Colorado Department of Public Health and Environment (CDPHE) Vital Statistics Program. Data can be requested for qualifying purposes through a data use agreement. Additional informtion is available at: https://cdphe.colorado.gov/request-data-for-research.

### Colorado EnviroScreen data
- Form: Colorado_EnviroScreen_v1_CensusTract.csv 
- Downloaded from https://cdphe.colorado.gov/enviroscreen. 
- This is publicly available without restriction.
- No processing was done.

### PM2.5 data
- Form: YYYY_pm25_daily_average.txt
- Downloaded from https://www.epa.gov/hesc/rsig-related-downloadable-data-files. 
- This is publicly available without restriction.
- Processed with scripts: 1) CMAQdownscalar_CO.R, 2) CMAQdownscalar_CO_lag_tract.R

### Temperature data
- Form: tmmx_YYYY.nc
- Downloaded from https://www.climatologylab.org. 
- This is publicly available without restriction.
- Processed with scripts: 1) gridmet_temp_lag_tract.R, 2) gridmet_temp.R

## Constrution of Data Analysis File

- Processed with the script: data_prep_Demateis_Enviroscreen.R
- The script takes the health data and merges it with the PM2.5 and temperature data, as well as the Colorado EnviroScreen data.






