rm(list=ls())
gc()
library(tidyverse)
library(ncdf4)
library(sf)
library(tigris)

# library(raster)
# library(sp)
# library(tidync)

# data source
# https://www.nature.com/articles/s41597-021-00973-0
# https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=2129
# downloaded 12/6/2022


#------------------
# FIND COLORADO BOUNDARY
#------------------

# state
co <- states() %>% filter(NAME=="Colorado")
boundary <- st_bbox(st_sf(co))
boundary

lat_range <- boundary[c(2,4)]
lon_range <- boundary[c(1,3)]

lat_range
lon_range


# get tracts to merge with
tracts <- tracts(state = "CO", year=2010) %>%
  select(GEOID=GEOID10, geometry)


#------------------
# Check that all grids are the same
#------------------


check_grid <- data.frame(year=2000:2020, equal=NA)
year <- check_grid$year[1]
print(year)
nc_data <- nc_open(paste0('ExposureData/gridmet_temp/tmmx_',year,'.nc'))
lon_grid <- ncvar_get(nc_data, "lon")
lat_grid <- ncvar_get(nc_data, "lat")
nc_close(nc_data)
check_grid[check_grid$year==year,"equal"] <-  TRUE


for(year in check_grid$year[-1]){
  
  print(year)
  nc_data <- nc_open(paste0('ExposureData/gridmet_temp/tmmx_',year,'.nc'))
  lon_grid1 <- ncvar_get(nc_data, "lon")
  lat_grid1 <- ncvar_get(nc_data, "lat")
  
  check_grid[check_grid$year==year,"equal"] <- all.equal(lon_grid1,lon_grid) & all.equal(lat_grid1,lat_grid)
  nc_close(nc_data)
}

all(check_grid$equal)


#------------------
# FIND OVERLAP BETWEEN GRID AND TRACTS
#------------------

keep_lon <- which(lon_grid>lon_range[1] & lon_grid < lon_range[2])
keep_lat <- which(lat_grid>lat_range[1] & lat_grid < lat_range[2])
lon <- lon_grid[which(lon_grid>lon_range[1] & lon_grid < lon_range[2])]
lat <- lat_grid[which(lat_grid>lat_range[1] & lat_grid < lat_range[2])]

# make to data.frame with lon and lats
df1 <- data.frame(id=1:(length(lon)*length(lat)), lon=rep(lon,length(lat)), lat=rep(lat,each=length(lon)))

grid_locations_sf <- st_as_sfc(paste0("POINT(",df1$lon," ",df1$lat,")")) %>%
  st_sf(id = df1$id,
        lat = df1$lat,
        lon = df1$lon) %>%
  st_set_crs("NAD83")

# find distance between tracts and each grid point
# this is done because some tracts are small and contain no grid points. 
# treat grid as circles
dist <- st_distance(grid_locations_sf,tracts)
# define radius to search within. this is the size of one grid cell.
a <- sqrt(2*2000^2)
units(a) <- "m"
link <- apply(dist<a, 2,which)
names(link) <- tracts$GEOID
linkvec <- unlist(link)
grid_fips <- data.frame(GEOID = substr(names(linkvec),1,11), id = linkvec)


#------------------
# Get data year by year
#------------------

for(year in 2000:2020){
  
  print(year)
  
  # open data
  nc_data <- nc_open(paste0('ExposureData/gridmet_temp/tmmx_',year,'.nc'))
  # print(nc_data)
  # lon_grid <- ncvar_get(nc_data, "lon")
  # lat_grid <- ncvar_get(nc_data, "lat")
  time_grid <- ncvar_get(nc_data, "day")
  time_grid_date <- lubridate::as_date(as.POSIXct('1900-01-01 00:00:00') + as.difftime(nc_data$dim$day$vals,units="days"))
  
  # get data for bounding box around Colorado
  # keep_lon <- which(lon_grid>lon_range[1] & lon_grid < lon_range[2])
  # keep_lat <- which(lat_grid>lat_range[1] & lat_grid < lat_range[2])
  # lon <- lon_grid[which(lon_grid>lon_range[1] & lon_grid < lon_range[2])]
  # lat <- lat_grid[which(lat_grid>lat_range[1] & lat_grid < lat_range[2])]

  # # make to data.frame with lon and lats
  # df <- data.frame(id=1:(length(lon)*length(lat)), lon=rep(lon,length(lat)), lat=rep(lat,each=length(lon)))
  # 
  # grid_locations_sf <- st_as_sfc(paste0("POINT(",df$lon," ",df$lat,")")) %>%
  #   st_sf(id = df$id,
  #         lat = df$lat,
  #         lon = df$lon) %>%
  #   st_set_crs("NAD83")

  # find overlap between grid and census tracts
  # grid_fips <- st_join(grid_locations_sf, t2, join = st_intersects) %>%
  #   data.frame() %>%
  #   select(id, GEOID) %>%
  #   drop_na()

  # # find distance between tracts and each grid point
  # # this is done because some tracts are small and contain no grid points. 
  # # treat grid as circles
  # dist <- st_distance(grid_locations_sf,tracts)
  # # define radius to search within. this is the size of one grid cell.
  # a <- sqrt(2*2000^2)
  # units(a) <- "m"
  # link <- apply(dist<a, 2,which)
  # names(link) <- tracts$GEOID
  # linkvec <- unlist(link)
  # grid_fips <- data.frame(GEOID = substr(names(linkvec),1,11), id = linkvec)

  
  # get temperature data for full bounding box
  tmmx <- ncvar_get(nc_data, "air_temperature",
                    start = c(min(keep_lon),min(keep_lat),1),
                    count = c(length(keep_lon),length(keep_lat),length(time_grid)))
  
  
  df <- data.frame(id=rep(df1$id,length(time_grid)),
                   tmmx=c(tmmx),
                   date=rep(time_grid_date, each=nrow(df1)))
  
  # merge data
  full_year <- left_join(df, grid_fips, by="id")
  
  tract_day_level <- full_year %>%
    drop_na() %>%
    group_by(GEOID, date) %>%
    summarise(tmmx = mean(tmmx) - 273.15)  # subtraction to go from Kelvin to Celsius
  
  nc_close( nc_data )
  
  save(tract_day_level, file=paste0('ExposureData/gridmet_temp/gridmet_tmmx_temp_',year,'.rda'))
  
}


#------------------
# LOAD & COMBINE FILES
#------------------

rm(list=ls())

gridmet_tmmx_tract <- NULL
for(year in 2000:2020){
  load(file=paste0('ExposureData/gridmet_temp/gridmet_tmmx_temp_',year,'.rda'))
  
  gridmet_tmmx_tract <- bind_rows(gridmet_tmmx_tract,
                                 tract_day_level)
  
}


#------------------
# TRACT LEVEL FILE
#------------------
gridmet_tmmx_tract <- gridmet_tmmx_tract %>% 
  ungroup() %>% 
  mutate(fipscoor = substr(GEOID,3,5), 
         tract = substr(GEOID,6,11)) %>%
  select(fipscoor, tract, date, tmmx)

# tract_day_level %>% filter(GEOID=="08001007900") 
# 
# test <- daymet_temp_tract %>% filter(GEOID=="08001007900") 
# 
# daymet_temp_tract %>% filter(GEOID=="08001007900") %>%
#   ggplot(aes(x=date, y=tmax)) + geom_point()


save(gridmet_tmmx_tract, file=paste0('ExposureData/gridmet_temp/gridmet_tmmx_tract.rda'))


#------------------
# COUNTY LEVEL FILE
#------------------

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


# merge pop and gridmet data
# summarize to county level
gridmet_tmmx_county <- left_join(gridmet_tmmx_tract,
                                tractpop,
                                by = c("fipscoor", "tract")) %>%
  group_by(fipscoor, date) %>%
  summarise(tmmx = weighted.mean(tmmx,pop2010))



save(gridmet_tmmx_county, file=paste0('ExposureData/gridmet_temp/gridmet_tmmx_county.rda'))
save(gridmet_tmmx_county, file=paste0('~/Dropbox/gridmet_temp/gridmet_tmmx_county.rda'))
save(gridmet_tmmx_tract, file=paste0('~/Dropbox/gridmet_temp/gridmet_tmmx_tract.rda'))

# 
gridmet_tmmx_county %>% filter(fipscoor=="069") %>%
  ggplot(aes(x=date, y=tmmx)) + geom_point()





# 
# 
# 
# ggplot(tract_day_level %>% filter(GEOID=="08001008353")) +
#   geom_point(aes(x=date, y=tmmx))
# 
# 
