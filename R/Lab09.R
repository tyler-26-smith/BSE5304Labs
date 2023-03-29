 Sys.getenv('USER')
 LabNo="/Lab09"
#
# What needs to be loaded
#
 if (!require("pacman")) install.packages("pacman")
 myhomedir=Sys.getenv("HOME")
 datadir=paste0(myhomedir,"/data",LabNo)
 dir.create(datadir,recursive = T)
 srcdir=paste0(myhomedir,"/src")
 dir.create(srcdir,recursive = T)
# Setting the directory for where the GitHub project exists. 
# This depends on where you set up your git, and what you called it locally, 
# but when you start a new git project, it will be the first directory you 
# are placed in... or if later in the project:
# 
 mygitdir=rstudioapi::getActiveProject()
 mypdfdir=paste0(mygitdir,"/pdfs",LabNo)
 dir.create(mypdfdir)
# 
 setwd(mygitdir)
 system("git config --global user.email 'tylersmith22@vt.edu' ") 
 system("git config --global user.name 'tyler-26-smith' ")
 system("git config pull.rebase false")

 if (!require("pacman")) install.packages("pacman")
 pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
                 data.table,foreign,maptools,dataRetrieval,gdistance)
 setwd(datadir)
#
# Note we have a new library to access USGS Waterdata
# https://owi.usgs.gov/R/dataRetrieval.html
# https://owi.usgs.gov/R/training-curriculum/usgs-packages/dataRetrieval-readNWIS/
#

  ?dataRetrieval  # Review the man page for this package
  ?readNWISuv
  ?readNWISdv
  ?readNWISdata
 
  url="https://help.waterdata.usgs.gov/parameter_cd?group_cd=%"
  browseURL(url)
  View(parameterCdFile)
 
 ##############################################
 # 0205551460 LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
 ##############################################
  siteNo = "0205551460"
  parameterCd = c("00060","00065")
  start.date = "2017-05-01"  # Not frozen to not frozen
  end.date = "2017-11-01"    # to still not frozen
 #
 # For each gage location, let's keep the data organized as a 
 # list.
  USGS0205551460=list()   # Organize the data in a nice list as in previous labs
   USGS0205551460[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
   head(USGS0205551460$flowdata)  # Note that we have 00060 and 00065...
  #  agency_cd	site_no        	dateTime X_00060_00000 X_00060_00000_cd
  #1  	USGS 0205551460 2017-05-01 04:00:00      	6.38            	A
  #2  	USGS 0205551460 2017-05-01 04:05:00      	6.38            	A
  #  X_00065_00000 X_00065_00000_cd tz_cd
  #1      	2.74            	A   UTC
  #2      	2.74            	A   UTC
  #
  # And of course we want to work in SI units so:
   USGS0205551460$flowdata$depth_m=USGS0205551460$flowdata$X_00065_00000*0.3048
  # m/ft depth
   USGS0205551460$flowdata$cms=USGS0205551460$flowdata$X_00060_00000*.02832
  # m3/ft3 flow
  #
  # Let's add in the USGS gage site information to the list and inspect
   USGS0205551460[["site"]]=readNWISsite(siteNo)
   head(USGS0205551460$site)
   class(USGS0205551460$site$dec_lat_va)
  #
  # Set the Manning Coefficient in the USGS Gage's Site Table
  #
   url="https://www.google.com/search?q=manning%27s+n+for+stream"
   browseURL(url)
   url="https://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm"
   browseURL(url)
   USGS0205551460$site$man_n=.035/1.49
  #
  # Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
   coordinates(USGS0205551460$site)=~dec_long_va+dec_lat_va
  
