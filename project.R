#------------------------------------------------------------------------------
#   ENSO AND AUSTRALIAN PRECIPITATION VARIABILITY (before and after 1880)
library(ncdf4)
library(maps)
library(RColorBrewer)
library(paletteer)

#------------------------------------------------------------------------------


# read data -----------------------------------
f1 <- nc_open("../data/ModE-RA_ensmean_temp2_abs_1420-2009nc.sec")
f2 <- nc_open("../data/ModE-RA_ensmean_totprec_abs_1420-2009nc.sec")
LaNina <- unlist(read.table("../exercises/LaNinaYears.txt", header = T))
ElNino <- unlist(read.table("../exercises/ElNinoYears.txt", header = T)) 

time <- f1$dim[[1]]$vals
lon <- f1$dim[[2]]$vals
lat <- f1$dim[[3]]$vals


# temp composites -------------------------
temp <- ncvar_get(f1,varid="temp2")
temp2 <- array(temp,dim=c(length(lon),length(lat),12,length(time)/12))
temp3 <- aperm(temp2,c(1,2,4,3))
temp4 <- apply(temp3,c(1:3),mean)
temp_mean <- apply(temp3, c(1,2), mean)

# prec composites --------------------------
prec <- ncvar_get(f2,varid="totprec")
prec2 <- array(prec,dim=c(length(lon),length(lat),12,length(time)/12))
prec3 <- aperm(prec2,c(1,2,4,3)) # EINHEIT kg/m^2s
# wir brauchen kg/m^2 für monate
# 60*60*24*365 pro jahr
prec4 <- apply(prec3,c(1:3),mean)
prec4 = prec4*60*60*24*365
prec_mean <- apply(prec4, c(1,2), mean)
prec_sd <- apply(prec4, c(1,2), sd)





# PLOTS Composites
x <- lon
y <- rev(lat)

## average temperature to see ENSO
z <- temp_mean[,rev(1:length(lat))]

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- seq(210, 310, 10)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Mean T all years")## average prec


## average precipitation
z <- prec_mean[,rev(1:length(lat))]

mycol <- c("#c6dbef", "#9ecae1", "#4292c6", "#2171b5", "#08519c", "#08306b")
mylevs <- seq(0,6000,1000)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Mean yearly prec [mm]")


## sd precipitation
z <- prec_sd[,rev(1:length(lat))]

mylevs <- seq(0,600,100)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Sd of yearly prec [mm]")



# PLOT correlation

## correlation of global temperature with precipitation in australia 
## -26.853388 S, 133.275154 E
clon <- 133.275154
clat <- -26.853388
dlon <- lon[2]-lon[1]
dlat <- lat[1]-lat[2]

## now selcet that gridcell: 
sellon <- (lon<(clon+(dlon/2)))&(lon>=(clon-(dlon/2)))
sellat <- (lat<(clat+(dlat/2)))&(lat>=(clat-(dlat/2)))
ausie <- prec4[sellon,sellat,]

par(mfrow=c(1,2))

plotfield <- apply(temp4,c(1,2),cor,ausie)
z <- plotfield[,rev(1:length(lat))]

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor of prec in Australia with world temperature")




# correlation of average Elnino temperature vs
# precipitation in Australia
yr.jfm <- c(1421:2008)
selnino <- match(ElNino,yr.jfm)
ausie <- prec4[sellon,sellat,selnino]
par(mfrow=c(1,2))

plotfield <- apply(temp4[,,selnino],c(1,2),cor,ausie)
z <- plotfield[,rev(1:length(lat))]

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor of prec in Australia with world temperature for ElNino years")



# correlation of average LaNina temperature 
# vs precipitation in Australia

selnina <- match(LaNina,yr.jfm)                     
ausie <- prec4[sellon,sellat,selnina]
plotfield <- apply(temp4[,,selnina],c(1,2),cor,ausie)
# y <- rev(sel.lat.num)
# z <- plotfield[,rev(1:length(lat[sellat]))]

y <- rev(lat)
z <- plotfield[,rev(1:length(lat))]

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor of prec in Australia with world temperature for LaNina years")




## Cropping of Map function
sellat <- which(lat<20)
sel.lat.num <- lat[sellat]




# Erstellen der Timeseries aus Mittelwerten für australien 
# Range der Koordinaten:
# S -30.52     -     -19.109
# E 123.11     -      142.34
l = lat >= -30.52 & lat <= -19.109
b = lon >= 123.11 & lon <= 142.34

prec_time = data.frame(prec = apply(prec4[b,l,], 3, mean),
                       time = 1421:2008)
plot(prec_time$time,prec_time$prec)

ElNino <- sort(ElNino)
time_cuts = data.frame("1425" = rep(NA,3),
                       "1485" = rep(NA,3),
                       "1610" = rep(NA,3),
                       "1647" = rep(NA,3),
                       "1735" = rep(NA,3),
                       "1747" = rep(NA,3),
                       "1759" = rep(NA,3),
                       "1783" = rep(NA,3),
                       "1804" = rep(NA,3),
                       "1869" = rep(NA,3),
                       "1878" = rep(NA,3),
                       "1889" = rep(NA,3),
                       "1914" = rep(NA,3),
                       "1920" = rep(NA,3),
                       "1931" = rep(NA,3),
                       "1973" = rep(NA,3),
                       "1987" = rep(NA,3),
                       "1988" = rep(NA,3),
                       "1992" = rep(NA,3),
                       "1993" = rep(NA,3),
                       "1998" = rep(NA,3))

for (i in 1:length(ElNino)){
  ind = which(prec_time$time == ElNino[i])
  time_cuts[1,i] = prec_time$prec[ind]
  time_cuts[2,i] = prec_time$prec[ind+1]
  time_cuts[3,i] = prec_time$prec[ind+2]
  # time_cuts[4,i] = prec_time$prec[ind+3]
  # time_cuts[5,i] = prec_time$prec[ind+4]
}

par(mfrow=c(1,2))
plot(time_cuts[,1], type = "l", ylim = c(50,200), main = "bis 1800")
for (i in 2:8) {
  lines(time_cuts[,i])
  
}
plot(time_cuts[,9], type = "l", ylim = c(50,200), main = "ab 1800")
for (i in 10:length(ElNino)) {
  lines(time_cuts[,i])
  
}

dev.off()

plot(prec_time$prec, type="l")
