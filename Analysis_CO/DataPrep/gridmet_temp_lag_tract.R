# Create lagged variables for temperature and trimester and pregnancy average exposure for temperature
# this is done at the census tract level
# Written by Ander Wilson
# Date created: 8/1/2024
# Date last edited: 8/9/2024


library(tidyverse)
library(data.table)

load(file="ExposureData/gridmet_temp/gridmet_tmmx_tract.rda")
ls()


# function to make weekly lags
wklead <- function(x,wk){
  a <- (wk-1)*7
  b <- rowMeans(cbind(lead(x, n=0+a),
                      lead(x, n=1+a),
                      lead(x, n=2+a),
                      lead(x, n=3+a),
                      lead(x, n=4+a),
                      lead(x, n=5+a),
                      lead(x, n=6+a)), na.rm=FALSE)
  return(b)
}



# make weekly lags for PM
gridmet_tmmx_tract_weeklag <- gridmet_tmmx_tract %>% 
  data.table()

setorder(gridmet_tmmx_tract_weeklag, fipscoor, tract, date)

for(l in 1:42){
  gridmet_tmmx_tract_weeklag[, (paste0("tmmx_lag",l)) := wklead(tmmx, l), by=c("fipscoor","tract")]  
}



# save files
gridmet_tmmx_tract_weeklag <- gridmet_tmmx_tract_weeklag %>% 
  tibble() %>% select(-c("tmmx"))


save(gridmet_tmmx_tract_weeklag,file="ExposureData/gridmet_temp/gridmet_tmmx_tract_weeklag.rda")



# make trimester average exposures files and pregnancy average
gridmet_tmmx_tract_trimester <- gridmet_tmmx_tract_weeklag %>%
  mutate(tmmx_tri1 = (  tmmx_lag1
    +tmmx_lag2
    +tmmx_lag3
    +tmmx_lag4
    +tmmx_lag5
    +tmmx_lag6
    +tmmx_lag7
    +tmmx_lag8
    +tmmx_lag9
    +tmmx_lag10
    +tmmx_lag11
    +tmmx_lag12)/12,
    tmmx_tri2 = ( tmmx_lag13
      +tmmx_lag14
      +tmmx_lag15
      +tmmx_lag16
      +tmmx_lag17
      +tmmx_lag18
      +tmmx_lag19
      +tmmx_lag20
      +tmmx_lag21
      +tmmx_lag22
      +tmmx_lag23
      +tmmx_lag24
      +tmmx_lag25)/13,
    tmmx_tri3 = ( tmmx_lag26
      +tmmx_lag27
      +tmmx_lag28
      +tmmx_lag29
      +tmmx_lag30
      +tmmx_lag31
      +tmmx_lag32
      +tmmx_lag33
      +tmmx_lag34
      +tmmx_lag35
      +tmmx_lag36
      +tmmx_lag37)/12,
    tmmx_preg = (  tmmx_lag1
                      +tmmx_lag2
                      +tmmx_lag3
                      +tmmx_lag4
                      +tmmx_lag5
                      +tmmx_lag6
                      +tmmx_lag7
                      +tmmx_lag8
                      +tmmx_lag9
                      +tmmx_lag10
                      +tmmx_lag11
                      +tmmx_lag12
                      +tmmx_lag13
                      +tmmx_lag14
                      +tmmx_lag15
                      +tmmx_lag16
                      +tmmx_lag17
                      +tmmx_lag18
                      +tmmx_lag19
                      +tmmx_lag20
                      +tmmx_lag21
                      +tmmx_lag22
                      +tmmx_lag23
                      +tmmx_lag24
                      +tmmx_lag25
                      +tmmx_lag26
                      +tmmx_lag27
                      +tmmx_lag28
                      +tmmx_lag29
                      +tmmx_lag30
                      +tmmx_lag31
                      +tmmx_lag32
                      +tmmx_lag33
                      +tmmx_lag34
                      +tmmx_lag35
                      +tmmx_lag36
                      +tmmx_lag37)/37) %>%
  select(-starts_with("tmmx_lag"))



# save files
save(gridmet_tmmx_tract_trimester,file="ExposureData/gridmet_temp/gridmet_tmmx_tract_trimester.rda")

