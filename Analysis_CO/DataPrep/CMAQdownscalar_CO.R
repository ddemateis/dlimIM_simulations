# Create clean file of CMAQ downscalar data at the census tract level
# also create population weighted census tract level files.
# Written by Ander Wilson
# Date created: 11/17/2022
# Date last edited: 8/9/2024

library(tidyverse)

year_range <- 2002:2019


# load and combine data for all years
# these files are saved as separate .txt files

# initiate combined file 
cmaqdownscaler_colorado_tract <- tribble()

for(year in year_range){
  print(year)
  
  #sample file to know layout
  dta <- read_csv(paste0("ExposureData/CMAQ_ozone_pm/",year,"_pm25_daily_average.txt"),n_max = 1)
  
  # load PM data
  if(class(attr(dta,"spec")$cols[[1]])[1]=="collector_date"){
    pm_in <- read_csv(paste0("ExposureData/CMAQ_ozone_pm/",year,"_pm25_daily_average.txt"),
                    col_select = c(1,2,5),
                    col_types=list(col_date(),
                                   col_character(),
                                   col_double(),
                                   col_double(),
                                   col_double(),
                                   col_double()))
  }else{
    pm_in <- read_csv(paste0("ExposureData/CMAQ_ozone_pm/",year,"_pm25_daily_average.txt"),
                    col_select = c(1,2,5),
                    col_types=list(col_date("%b-%d-%Y"),
                                   col_character(),
                                   col_double(),
                                   col_double(),
                                   col_double(),
                                   col_double()))
  }
  
  # some years have different names
  colnames(pm_in) <- c("Date","FIPS","cmaq_pm25")
  
  pm_in <- pm_in %>% 
    mutate(FIPS = if_else(nchar(FIPS)==10, paste0("0",FIPS), FIPS)) %>% # add leading 0s that were dropped
    filter(substr(FIPS,1,2)=="08") %>% # keep colorado only
    mutate(fipscoor = substr(FIPS,3,5), # make identifiers match birth data
           tract = substr(FIPS,6,11)) %>%
    select(-FIPS) # drop full FIPS
  
  #sample file to know layout
  dta <- read_csv(paste0("ExposureData/CMAQ_ozone_pm/",year,"_ozone_daily_8hour_maximum.txt"),n_max = 1)
  
  # load ozone data
  if(class(attr(dta,"spec")$cols[[1]])[1]=="collector_date"){
    ozone_in <- read_csv(paste0("ExposureData/CMAQ_ozone_pm/",year,"_ozone_daily_8hour_maximum.txt"),
                      col_select = c(1,2,5),
                      col_types=list(col_date(),
                                     col_character(),
                                     col_double(),
                                     col_double(),
                                     col_double(),
                                     col_double()))
  }else{
    ozone_in <- read_csv(paste0("ExposureData/CMAQ_ozone_pm/",year,"_ozone_daily_8hour_maximum.txt"),
                      col_select = c(1,2,5),
                      col_types=list(col_date("%b-%d-%Y"),
                                     col_character(),
                                     col_double(),
                                     col_double(),
                                     col_double(),
                                     col_double()))
  }
  
  
  # some years have different names
  # this fixes that
  colnames(ozone_in) <- c("Date","FIPS","cmaq_o3")
  
  ozone_in <- ozone_in %>% 
    mutate(FIPS = if_else(nchar(FIPS)==10, paste0("0",FIPS), FIPS)) %>% # add leading 0s that were dropped
    filter(substr(FIPS,1,2)=="08") %>% # keep colorado only
    mutate(fipscoor = substr(FIPS,3,5), # make identifiers match birth data
           tract = substr(FIPS,6,11)) %>%
    select(-FIPS) # drop full FIPS
  
  # merge ozone and pm data
  one_year <- 
    full_join(pm_in,
              ozone_in, 
              by=c("Date","fipscoor","tract")) %>% 
    relocate(cmaq_pm25, .after = last_col())
  
  # save to all years
  cmaqdownscaler_colorado_tract <- bind_rows(cmaqdownscaler_colorado_tract,
                                       one_year)
}


table(lubridate::year(cmaqdownscaler_colorado_tract$Date))


# this creates population weighted averages at the county level

# get county population data
library(tidycensus)
sf1 <- load_variables(2010, "sf1", cache = TRUE)

tractpop <- get_decennial(geography = "tract",
                          state = "CO",
                          variables = "P001001",
                          year = 2010) %>%
  rename("pop2010"="value") %>%
  mutate(fipscoor = substr(GEOID,3,5), # make identifiers match birth data
         tract = substr(GEOID,6,11)) %>%
  select(fipscoor, tract, pop2010) %>%
  arrange(fipscoor, tract)

rm("sf1")

# merge pop and cmaq data
# summarize to county level
cmaqdownscaler_colorado_county <- left_join(cmaqdownscaler_colorado_tract,
                                            tractpop,
                                            by = c("fipscoor", "tract")) %>%
  group_by(fipscoor, Date) %>%
  summarise(cmaq_pm25 = weighted.mean(cmaq_pm25,pop2010),
            cmaq_o3 = weighted.mean(cmaq_o3,pop2010))



# save files
save(cmaqdownscaler_colorado_tract,file="ExposureData/CMAQ_ozone_pm/cmaqdownscaler_colorado_tract.rda")
save(cmaqdownscaler_colorado_county,file="ExposureData/CMAQ_ozone_pm/cmaqdownscaler_colorado_county.rda")

