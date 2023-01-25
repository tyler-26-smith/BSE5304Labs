if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,dplyr,patchwork,rnoaa,dplyr,hstr)
system("git config --global user.email 'tylersmith22@vt.edu' ") 
system("git config --global user.name 'tyler-26-smith' ")

########## USGS Stage Data
  
  source("https://goo.gl/Cb8zGn")
  myflowgage_id="01672500"
  myflowgage=get_usgs_gage(myflowgage_id,
                           begin_date="2022-01-01",end_date="2023-01-24")
  class(myflowgage)
  View(myflowgage$flowdata)

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
  
 #USC00440327
  WXData=meteo_pull_monitors(
    monitors=stns[25,1],    # replace the *** with index you find
    keep_flags = FALSE,
    date_min = "2022-01-01",
    date_max = NULL,
    var = c("TMAX","TMIN","PRCP")
    
  )
  
########## Graphing
    coeff <- 10
    # A few constants
    tminColor <- "#0000ff"
      tmaxColor <- "#ff0000"
        prcpColor <- "#000000"
        dir.create("pdfs")
        basestr=format(Sys.time(),"./pdfs/%Y%m%d%H%M")
        
      p1 <- ggplot(WXData, aes(x=date)) +
          geom_line( aes(y=tmin / coeff), linewidth=0.5, color=tminColor) + 
          geom_line( aes(y=tmax / coeff), linewidth=0.5, color=tmaxColor) + 
          geom_line( aes(y=prcp / coeff), linewidth=0.5, color=prcpColor) +
          scale_y_continuous(
            # Features of the first axis
            name = "Temperature (C)",
            
            # Add a second axis and specify its features
            sec.axis = sec_axis(~., name="Precipitation (mm)")
          ) + 
          theme(
            axis.title.y = element_text(color = "black", size=13),
            axis.title.y.right = element_text(color = prcpColor, size=13)
          ) +
          ggtitle("Ashland, VA Temperature and Precipitation")
      
      p2 <- ggplot(myflowgage$flowdata, aes(x=mdate)) +
          geom_line( aes(y=flow), linewidth=0.5, color="#000000") +
          labs(y = "Flow (m^3/day)", x = "Date") +
          theme(axis.title.y = element_text(color = "black", size=13)
          ) +
          ggtitle("South Anna River Near Ashalnd, VA Flow Data")
        
plot(p1)
plot(p2)
plot(p1+p2)
        