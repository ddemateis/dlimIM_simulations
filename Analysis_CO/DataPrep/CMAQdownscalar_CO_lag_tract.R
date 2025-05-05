# Create lagged variables for pm and ozone from CMAQ and trimester and pregnancy average exposures
# this is done at the census tract level
# Written by Ander Wilson
# Date created: 8/1/2024
# Date last edited: 8/9/2024

library(tidyverse)
library(data.table)

load(file="ExposureData/CMAQ_ozone_pm/cmaqdownscaler_colorado_tract.rda")
ls()


# load tract level data
cmaqdownscaler_colorado_tract <- cmaqdownscaler_colorado_tract %>% 
  drop_na(Date)  

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
cmaqdownscaler_colorado_tract_pm25_weeklag <- cmaqdownscaler_colorado_tract %>% 
  data.table()

setorder(cmaqdownscaler_colorado_tract_pm25_weeklag, fipscoor, tract, Date)

for(l in 1:42){
  cmaqdownscaler_colorado_tract_pm25_weeklag[, (paste0("cmaq_pm25_lag",l)) := wklead(cmaq_pm25, l), by=c("fipscoor","tract")]  
}


# make weekly lags for Ozone
cmaqdownscaler_colorado_tract_o3_weeklag <- cmaqdownscaler_colorado_tract %>% 
  data.table()

setorder(cmaqdownscaler_colorado_tract_o3_weeklag, fipscoor, tract, Date)

for(l in 1:42){
  cmaqdownscaler_colorado_tract_o3_weeklag[, (paste0("cmaq_o3_lag",l)) := wklead(cmaq_o3, l), by=c("fipscoor","tract")]  
}


# save files
cmaqdownscaler_colorado_tract_o3_weeklag <- cmaqdownscaler_colorado_tract_o3_weeklag %>% 
  tibble() %>% select(-c("cmaq_o3", "cmaq_pm25"))
cmaqdownscaler_colorado_tract_pm25_weeklag <- cmaqdownscaler_colorado_tract_pm25_weeklag %>% 
  tibble() %>% select(-c("cmaq_o3", "cmaq_pm25"))


save(cmaqdownscaler_colorado_tract_o3_weeklag,file="ExposureData/CMAQ_ozone_pm/cmaqdownscaler_colorado_tract_o3_weeklag.rda")
save(cmaqdownscaler_colorado_tract_pm25_weeklag,file="ExposureData/CMAQ_ozone_pm/cmaqdownscaler_colorado_tract_pm25_weeklag.rda")



# make trimester average exposures files and pregnancy average
cmaqdownscaler_colorado_tract_o3_trimester <- cmaqdownscaler_colorado_tract_o3_weeklag %>%
  mutate(cmaq_o3_tri1 = (  cmaq_o3_lag1
    +cmaq_o3_lag2
    +cmaq_o3_lag3
    +cmaq_o3_lag4
    +cmaq_o3_lag5
    +cmaq_o3_lag6
    +cmaq_o3_lag7
    +cmaq_o3_lag8
    +cmaq_o3_lag9
    +cmaq_o3_lag10
    +cmaq_o3_lag11
    +cmaq_o3_lag12)/12,
    cmaq_o3_tri2 = ( cmaq_o3_lag13
      +cmaq_o3_lag14
      +cmaq_o3_lag15
      +cmaq_o3_lag16
      +cmaq_o3_lag17
      +cmaq_o3_lag18
      +cmaq_o3_lag19
      +cmaq_o3_lag20
      +cmaq_o3_lag21
      +cmaq_o3_lag22
      +cmaq_o3_lag23
      +cmaq_o3_lag24
      +cmaq_o3_lag25)/13,
    cmaq_o3_tri3 = ( cmaq_o3_lag26
      +cmaq_o3_lag27
      +cmaq_o3_lag28
      +cmaq_o3_lag29
      +cmaq_o3_lag30
      +cmaq_o3_lag31
      +cmaq_o3_lag32
      +cmaq_o3_lag33
      +cmaq_o3_lag34
      +cmaq_o3_lag35
      +cmaq_o3_lag36
      +cmaq_o3_lag37)/12,
    cmaq_o3_preg = (  cmaq_o3_lag1
                      +cmaq_o3_lag2
                      +cmaq_o3_lag3
                      +cmaq_o3_lag4
                      +cmaq_o3_lag5
                      +cmaq_o3_lag6
                      +cmaq_o3_lag7
                      +cmaq_o3_lag8
                      +cmaq_o3_lag9
                      +cmaq_o3_lag10
                      +cmaq_o3_lag11
                      +cmaq_o3_lag12
                      +cmaq_o3_lag13
                      +cmaq_o3_lag14
                      +cmaq_o3_lag15
                      +cmaq_o3_lag16
                      +cmaq_o3_lag17
                      +cmaq_o3_lag18
                      +cmaq_o3_lag19
                      +cmaq_o3_lag20
                      +cmaq_o3_lag21
                      +cmaq_o3_lag22
                      +cmaq_o3_lag23
                      +cmaq_o3_lag24
                      +cmaq_o3_lag25
                      +cmaq_o3_lag26
                      +cmaq_o3_lag27
                      +cmaq_o3_lag28
                      +cmaq_o3_lag29
                      +cmaq_o3_lag30
                      +cmaq_o3_lag31
                      +cmaq_o3_lag32
                      +cmaq_o3_lag33
                      +cmaq_o3_lag34
                      +cmaq_o3_lag35
                      +cmaq_o3_lag36
                      +cmaq_o3_lag37)/37) %>%
  select(-starts_with("cmaq_o3_lag"))




cmaqdownscaler_colorado_tract_pm25_trimester <- cmaqdownscaler_colorado_tract_pm25_weeklag %>% 
  mutate(cmaq_pm25_tri1 = ( cmaq_pm25_lag1
                   +cmaq_pm25_lag2
                   +cmaq_pm25_lag3
                   +cmaq_pm25_lag4
                   +cmaq_pm25_lag5
                   +cmaq_pm25_lag6
                   +cmaq_pm25_lag7
                   +cmaq_pm25_lag8
                   +cmaq_pm25_lag9
                   +cmaq_pm25_lag10
                   +cmaq_pm25_lag11
                   +cmaq_pm25_lag12)/12,
         cmaq_pm25_tri2 = ( cmaq_pm25_lag13
                   +cmaq_pm25_lag14
                   +cmaq_pm25_lag15
                   +cmaq_pm25_lag16
                   +cmaq_pm25_lag17
                   +cmaq_pm25_lag18
                   +cmaq_pm25_lag19
                   +cmaq_pm25_lag20
                   +cmaq_pm25_lag21
                   +cmaq_pm25_lag22
                   +cmaq_pm25_lag23
                   +cmaq_pm25_lag24
                   +cmaq_pm25_lag25)/13,
         cmaq_pm25_tri3 = ( cmaq_pm25_lag26
                   +cmaq_pm25_lag27
                   +cmaq_pm25_lag28
                   +cmaq_pm25_lag29
                   +cmaq_pm25_lag30
                   +cmaq_pm25_lag31
                   +cmaq_pm25_lag32
                   +cmaq_pm25_lag33
                   +cmaq_pm25_lag34
                   +cmaq_pm25_lag35
                   +cmaq_pm25_lag36
                   +cmaq_pm25_lag37)/12,
         cmaq_pm25_preg = (  cmaq_pm25_lag1
                           +cmaq_pm25_lag2
                           +cmaq_pm25_lag3
                           +cmaq_pm25_lag4
                           +cmaq_pm25_lag5
                           +cmaq_pm25_lag6
                           +cmaq_pm25_lag7
                           +cmaq_pm25_lag8
                           +cmaq_pm25_lag9
                           +cmaq_pm25_lag10
                           +cmaq_pm25_lag11
                           +cmaq_pm25_lag12
                           +cmaq_pm25_lag13
                           +cmaq_pm25_lag14
                           +cmaq_pm25_lag15
                           +cmaq_pm25_lag16
                           +cmaq_pm25_lag17
                           +cmaq_pm25_lag18
                           +cmaq_pm25_lag19
                           +cmaq_pm25_lag20
                           +cmaq_pm25_lag21
                           +cmaq_pm25_lag22
                           +cmaq_pm25_lag23
                           +cmaq_pm25_lag24
                           +cmaq_pm25_lag25
                           +cmaq_pm25_lag26
                           +cmaq_pm25_lag27
                           +cmaq_pm25_lag28
                           +cmaq_pm25_lag29
                           +cmaq_pm25_lag30
                           +cmaq_pm25_lag31
                           +cmaq_pm25_lag32
                           +cmaq_pm25_lag33
                           +cmaq_pm25_lag34
                           +cmaq_pm25_lag35
                           +cmaq_pm25_lag36
                           +cmaq_pm25_lag37)/37) %>%
  select(-starts_with("cmaq_pm25_lag"))





# save files
save(cmaqdownscaler_colorado_tract_o3_trimester,file="ExposureData/CMAQ_ozone_pm/cmaqdownscaler_colorado_tract_o3_trimester.rda")
save(cmaqdownscaler_colorado_tract_pm25_trimester,file="ExposureData/CMAQ_ozone_pm/cmaqdownscaler_colorado_tract_pm25_trimester.rda")

