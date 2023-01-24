if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,dplyr,patchwork,rnoaa,dplyr,hstr)
system("git config --global user.email 'tylersmith22@vt.edu' ") 
system("git config --global user.name 'tyler-26-smith' ")

coeff <- 10

# A few constants
temperatureColor <- "#69b3a2"
  priceColor <- rgb(0.2, 0.6, 0.9, 1)
  
  ggplot(data, aes(x=day)) +
    
    geom_line( aes(y=temperature), size=2, color=temperatureColor) + 
    geom_line( aes(y=price / coeff), size=2, color=priceColor) +
    
    scale_y_continuous(
      
      # Features of the first axis
      name = "Temperature (Celsius °)",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*coeff, name="Price ($)")
    ) + 
    
    theme_ipsum() +
    
    theme(
      axis.title.y = element_text(color = temperatureColor, size=13),
      axis.title.y.right = element_text(color = priceColor, size=13)
    ) +
    
    ggtitle("Temperature down, price up")
  
  ########## USGS Stage Data
  
  source("https://goo.gl/Cb8zGn")
  myflowgage_id="01672500"
  myflowgage=get_usgs_gage(myflowgage_id,
                           begin_date="2017-02-01",end_date="2023-02-01")
  class(myflowgage)
  View(myflowgage$flowdata)
  
  plot(myflowgage$flowdata$mdate,myflowgage$flowdata$flow,
       main=myflowgage$gagename,xlab = "Date",
       ylab="Flow m^3/day",type="l")
  
  ########## NOAA
  station_data <- ghcnd_stations()
  meteo_distance(station_data, 37.8, -77.5, radius = 10, limit = 10)
  
  
  stns=meteo_distance(
    station_data=ghcnd_stations(),
    lat=myflowgage$declat,
    long=myflowgage$declon,
    units = "deg",
    radius = 20,
    limit = NULL
  )
  
 ##USC00440327
  WXData=meteo_pull_monitors(
    monitors=stns[25,1],    # replace the *** with index you find
    keep_flags = FALSE,
    date_min = "2016-01-01",
    date_max = NULL,
    var = c("TMAX","TMIN","PRCP")
    
  )
  
  plot(WXData$date,WXData$tmax,type="l", col="red")
  lines(WXData$date,WXData$tmin,type="l", col="blue")
  points(WXData$date,WXData$prcp, col="black")
  