# make contextual data for ncov19 outbreak
library(sp)
library(rgdal)
library(raster)
library(dplyr)

#import temp monitor data (a sample for now, need to read more in)
temp3 <- read.csv(file = "../data/context/2013051.csv")[, c(3, 4, 6, 11:14)]
temp4 <- read.csv(file = "../data/context/2013054.csv")[, c(3, 4, 6, 11:14)]
temp5 <- read.csv(file = "../data/context/2013024.csv")[, c(3, 4, 6, 27, 29, 31, 33)]
temp6 <- read.csv(file = "../data/context/2016926.csv")[, c(3, 4, 6, 8:11)]
#limit to february temps
temps <- rbind(temp3, temp4, temp5, temp6)
temps <- temps[substr(temps$DATE, 6, 7) == "02",]

#get administration boundaries
admU <-getData('GADM', country='CHINA', level=1)

#get average temps in each province
tempsdat2 <- data.frame(temps[, c(1:2,5)])
tempsdat <- data.frame(temps[, c(1:2,5)])
coordinates(tempsdat) <- c("LONGITUDE", "LATITUDE")
proj4string(tempsdat) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
tempsdat2$prov <- over(tempsdat, admU)$NAME_1
tempcsv <- aggregate(TAVG ~ prov, data = tempsdat2, FUN = mean)


#get average precipitation in each province
precipdat2 <- data.frame(temps[, c(1:2,4)])
precipdat <- data.frame(temps[, c(1:2,4)])
coordinates(precipdat) <- c("LONGITUDE", "LATITUDE")
proj4string(precipdat) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
precipdat2$prov <- over(precipdat, admU)$NAME_1
precipcsv <- aggregate(PRCP ~ prov, data = precipdat2, FUN = mean)

#estimate province areas in km^2
admU$area_sqkm <- area(admU) / 1000000
areacsv <- data.frame(prov=admU$NAME_1, area = admU$area_sqkm)
tacsv <- merge(tempcsv, areacsv, by = "prov", all = TRUE)

#estimate populations per province from world pop
#for now, use china statistics bureau data
pop <- read.csv("../data/context/pop2.csv", header = TRUE, 
                colClasses = c(popper10k = "numeric", PerCapitaDisposibleIncome = "numeric",
                               PassengerTraffic.in10000. = "numeric", BedsInHealthcareInst = "numeric"))
datcsv <- merge(tacsv, pop, by = "prov", all = TRUE)
datcsv$popdensity <- datcsv$popper10k*10000/datcsv$area
datcsv2 <- merge(datcsv, precipcsv, by = "prov", all = TRUE)

write.csv(datcsv2, file = "../data/provincedat.csv")

