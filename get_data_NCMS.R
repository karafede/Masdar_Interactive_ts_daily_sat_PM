
library(readr)
library(dplyr)
library(threadr)
library(tidyr)
library(dygraphs)

# function to generate time-series based on data for each year in the UAE

#setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/Interactive_plots_R")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV")


# station data

# NCMS_2013 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_NCMS_ 2013 _daily_filtered.csv")
# NCMS_2014 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_NCMS_ 2014 _daily_filtered.csv")
# NCMS_2015 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_NCMS_ 2015 _daily_filtered.csv")
# NCMS_2016 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_NCMS_ 2016 _daily_filtered.csv")

# load hourly data and filter only one specific time~~~~~~~~~~~~~~~~~~~
NCMS_2013_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_NCMS_ 2013 _hourly_filtered.csv")
NCMS_2014_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_NCMS_ 2014 _hourly_filtered.csv")
NCMS_2015_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_NCMS_ 2015 _hourly_filtered.csv")
NCMS_2016_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_NCMS_ 2016 _hourly_filtered.csv")


NCMS_2013_filtered_time <- filter(NCMS_2013_filtered, grepl('12:', DateTime))
NCMS_2014_filtered_time <- filter(NCMS_2014_filtered, grepl('12:', DateTime))
NCMS_2015_filtered_time <- filter(NCMS_2015_filtered, grepl('12:', DateTime))
NCMS_2016_filtered_time <- filter(NCMS_2016_filtered, grepl('11:', DateTime))




# DM_2013 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_DM_ 2013 _daily_filtered.csv")
# DM_2014 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_DM_ 2014 _daily_filtered.csv")
# DM_2015 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_DM_ 2015 _daily_filtered.csv")
# DM_2016 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_DM_ 2016 _daily_filtered.csv")

# load hourly data and filter only one specific time~~~~~~~~~~~~~~~~~~~
DM_2013_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_DM_ 2013 _hourly_filtered.csv")
DM_2014_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_DM_ 2014 _hourly_filtered.csv")
DM_2015_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_DM_ 2015 _hourly_filtered.csv")
DM_2016_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_DM_ 2016 _hourly_filtered.csv")

DM_2016_filtered$Site  <- ifelse(grepl("DUBAIAIRPORT", DM_2016_filtered$Site, ignore.case = TRUE), 
                    "DUBAI AIR PORT", DM_2016_filtered$Site)


DM_2013_filtered_time <- filter(DM_2013_filtered, grepl('12:', DateTime))
DM_2014_filtered_time <- filter(DM_2014_filtered, grepl('12:', DateTime))
DM_2015_filtered_time <- filter(DM_2015_filtered, grepl('12:', DateTime))
DM_2016_filtered_time <- filter(DM_2016_filtered, grepl('11:', DateTime))



# EAD_2013 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_EAD_ 2013 _daily_filtered.csv")
# EAD_2014 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_EAD_ 2014 _daily_filtered.csv")
# EAD_2015 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_EAD_ 2015 _daily_filtered.csv")
# EAD_2016 <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/daily data/daily moved/daily_filtered_4_box/database_EAD_ 2016 _daily_filtered.csv")

# load hourly data and filter only one specific time~~~~~~~~~~~~~~~~~~~
EAD_2013_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_EAD_ 2013 _hourly_filtered.csv")
EAD_2014_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_EAD_ 2014 _hourly_filtered.csv")
EAD_2015_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_EAD_ 2015 _hourly_filtered.csv")
EAD_2016_filtered <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box/database_EAD_ 2016 _hourly_filtered.csv")


EAD_2013_filtered_time <- filter(EAD_2013_filtered, grepl('12:', DateTime))
EAD_2014_filtered_time <- filter(EAD_2014_filtered, grepl('12:', DateTime))
EAD_2015_filtered_time <- filter(EAD_2015_filtered, grepl('12:', DateTime))
EAD_2016_filtered_time <- filter(EAD_2016_filtered, grepl('11:', DateTime))

# bind data together
# AQ_data <- rbind(EAD_2013, EAD_2014, EAD_2015, EAD_2016,
#                  DM_2013, DM_2014, DM_2015, DM_2016,
#                  NCMS_2013, NCMS_2014, NCMS_2015, NCMS_2016)


AQ_data <- rbind(EAD_2013_filtered_time, EAD_2014_filtered_time, EAD_2015_filtered_time, EAD_2016_filtered_time,
                 DM_2013_filtered_time, DM_2014_filtered_time, DM_2015_filtered_time, DM_2016_filtered_time,
                 NCMS_2013_filtered_time, NCMS_2014_filtered_time, NCMS_2015_filtered_time, NCMS_2016_filtered_time)


AQ_data_PM25 <- AQ_data %>%
  filter(Pollutant == "PM2.5") 

AQ_data_PM10 <- AQ_data %>%
  filter(Pollutant == "PM10")

AQ_data_PM <- rbind(AQ_data_PM25,
                    AQ_data_PM10)


dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/dawit Data/Hourly Database format CSV/Arranged dates/R files/filtered_4_box"
write_csv(AQ_data_PM, paste0(dir, "/","PM25_PM10_data_filtered_4_box.csv"))


# shift time back of three hours~~~~~~~~~~~~~~~~~~~
# AQ_data$DateTime <- (AQ_data$DateTime) - 4*60*60


# display only the date~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AQ_data <- AQ_data %>%
  mutate(DateTime = ymd_hms(DateTime, tz = "UTC"),
         Date = date(DateTime)) 


# satellite data 
Sat_data <- read_csv("Z:/_SHARED_FOLDERS/Air Quality/Phase 1/Pathflow of Phase I_DG/Sat_AOD_Correlation/PM25_from_AOD_MODIS.csv")



# location of the stations with the PM2.5 

get_sites <- function(var) {
  NCMS_PM25 <- AQ_data %>%
    filter(Pollutant == var) %>%
    distinct(Site, Latitude, Longitude)
  
# Return
NCMS_PM25
}

 
get_measurement_time_series <- function(station, pollutant) {

  # Import hourly data from several years
  
  # NCMS[sapply(NCMS,is.na)] = NA 
  

  # NCMS <- AQ_data %>%
  #  #  mutate(date = mdy_hms(DateTime, tz = "UTC")) %>%
  #    mutate(date = ymd(Date)) %>%
  #    dplyr:: select(date,
  #           Site,
  #           Pollutant,
  #           Daily_mean) %>%
  #   filter(Site == station)
  
  
  NCMS <- AQ_data %>%
    #  mutate(date = mdy_hms(DateTime, tz = "UTC")) %>%
    mutate(date = ymd(Date)) %>%
    dplyr:: select(date,
                   Site,
                   Pollutant,
                   Value) %>%
    filter(Site == station)
  
  
  
  # replace NaN (not a Number with NA that is a missing value)
  # NCMS[sapply(NCMS,is.na)] = NA 
  NCMS_filtered <- Sat_data %>%
    #  mutate(date = mdy_hms(DateTime, tz = "UTC")) %>%
    mutate(date = ymd(Date)) %>%
    dplyr:: select(date,
                   Site,
                   AOD_PM25) %>%
    filter(Site == station)
  
  
  # data_time <- NCMS %>%
  #   spread(Pollutant, Daily_mean)
  
  data_time <- NCMS %>%
    spread(Pollutant, Value)
  
  data_time_filtered <- NCMS_filtered %>%
    select(-Site)
  
  # data_time <- data_time %>%
  #   left_join(data_time_filtered, by = "date")
  
  # Build timeseries for plots
  time_series <- data_frame_to_timeseries(data_time,  tz = "UTC")
  time_series_filtered <- data_frame_to_timeseries(data_time_filtered,  tz = "UTC")
  
  # Return
  #time_series
  #time_series_filtered
  # bind two time series together
  All_data <- cbind(time_series,time_series_filtered)
  
  #return
  
  All_data
}

# All_data
# data_both<- All_data[pollu,]
# data_both_ts<-as.ts(data_both[1])
# data_both_ts_2<-as.ts(data_both[2])
# data_both_ts$time_series
# data_both_ts_2$time_series_filtered
############################################################################
############################################################################

# to create grouped interactive dygraphs
# pollutant <- "PM<sub>2.5</sub>"
# station<-"Deira"
# 
# ts<-All_data
# station = "Zabeel"
# group = pollutant
# pollu= "PM2.5"
# da<- is.na(data_time$AOD_PM25 )
# dada <-  which(da , arr.ind = T,useNames = TRUE)
# 

# station<- "Al Ain Street"
# pollutant <- "PM<sub>2.5</sub>"
# 
# All_data<-get_measurement_time_series(station, pollutant)
# ts<-All_data
# pollu<-"PM2.5"
# group =   pollutant
# 
# 
# station = "Al Ain Islamic Ins"
# group =   pollutant
# pollu="PM2.5"
# data_both<- All_data[pollu,]
# data_both_ts<-as.ts(data_both[1])
# data_both_ts_2<-as.ts(data_both[2])

#ts_xxx <- cbind(data_both_ts$time_series, data_both_ts_2$time_series_filtered)
#dawit<-interactive_plot(time_series_BainAljesrain, station, group, pollu)


interactive_plot <- function(ts, station, group, pollu) {
  check_<-row.names(ts)
  
  if (!is.null(ts) & is.element(pollu, check_) ) {
    
    # Get colour vector
    colour_vector <- threadr::ggplot2_colours(45)
    
    
    #PM10
    
    if (pollutant == "PM<sub>10</sub>") {
      data_both<- ts[pollu,]
      data_both_ts<-as.ts(data_both[1])
      data_both_ts_2<-as.ts(data_both[2])
      
      ts_xxx <- cbind(data_both_ts$time_series, data_both_ts_2$time_series_filtered)
      
      
      plot <- dygraph(ts_xxx, group = group, main = paste(station, " - ", pollutant)) %>% 
        dySeries("..1",label = "Station", color = "red") %>%
        dySeries("..2",label = "SAT.", color = "blue") %>%
        dyAxis("y", label = "Hourly PM<sub>10</sub> (&#956;g m<sup>-3</sup>)") %>% 
        dyRangeSelector()
    }
    #PM2.5
    
    if (pollutant == "PM<sub>2.5</sub>") {
      data_both<- ts[pollu,]
      data_both_ts<-as.ts(data_both[1])
      data_both_ts_2<-as.ts(data_both[2])
      
      ts_xxx <- cbind(data_both_ts$time_series, data_both_ts_2$time_series_filtered)
      
      
      plot <- dygraph(ts_xxx, group = group, main = paste(station, " - ", pollutant)) %>% 
        dySeries("..1",label = "Station", color = "red") %>%
        dySeries("..2",label = "SAT.", color = "blue") %>%
        dyAxis("y", label = "Hourly PM<sub>2.5</sub> (&#956;g m<sup>-3</sup>)") %>% 
        dyRangeSelector()
    }
    
    
    
    
    # Return
    plot
    
  }
  
}

        
# pollutant <- "NO2"
# plot <- interactive_plot(time_series_Ghalilah$NO2, station = "Ghalilah", group = pollutant)
# plot



interactive_map_index <- function(df) {

  # Map
  map <- leaflet() %>%
    addTiles(group = "OpenStreetMap") %>%
    addProviderTiles("Stamen.Toner", group = "Toner") %>%
    addProviderTiles("Esri.WorldImagery", group = "Images") %>%
    addMarkers(data = df, lng = ~ Longitude, lat = ~ Latitude,
               popup = ~ Site, group = "Sites") %>%
    addPolygons(stroke = TRUE, data = shp_UAE,
                weight = 1.5, color = ~ colorNumeric(c("#a56e6e", "#7a7acc", "#FFFF00", "#ff0000", "#be68be", "#7fbf7f", "#008000", "#0000ff"), shp_UAE$ID_1)(ID_1),
                fillOpacity = 0.5,
                group = "shape_UAE") %>%
    addLayersControl(baseGroups = c("OpenStreetMap", "Toner", "Images"),
                     overlayGroups = c("Sites"))
  
  # Return
  map

}

# 

interactive_map <- function(df) {
  
  # Map
  map <- leaflet() %>%
    setView(lng = 55.9971, lat = 25.3302, zoom = 9) %>%  
    addTiles(group = "OpenStreetMap") %>%
    addProviderTiles("Stamen.Toner", group = "Toner") %>%
    addProviderTiles("Esri.WorldImagery", group = "Images") %>%
    addMarkers(data = df, lng = ~ Longitude, lat = ~ Latitude,
               popup = ~ Site, group = "Sites") %>%
    addLayersControl(baseGroups = c("OpenStreetMap", "Toner", "Images"),
                     overlayGroups = c("Sites"))
  
  # Return
  map
  
}


