
#   ENSO AND AUSTRALIAN PRECIPITATION VARIABILITY (before and after 1880)

# libraries ------------------------
library(ncdf4)
library(maps)
library(RColorBrewer)
library(paletteer)

# To Dos ----------------------------

# 1) crop composite map to australia to visualize local precipitation effects
# 2) do a bit prettier timeseries composites
# 3) calculate own el nino index


# read data -----------------------------------
f1 <- nc_open("../data_raw/ModE-RA_ensmean_temp2_abs_1420-2009.nc")
f2 <- nc_open("../data_raw/ModE-RA_ensmean_totprec_abs_1420-2009.nc")
LaNina <- read.table("../data_raw/LaNinaYears.txt", header = T)
ElNino <- read.table("../data_raw/ElNinoYears.txt", header = T)

# read dimensions
time <- f1$dim[[1]]$vals
lon <- f1$dim[[2]]$vals
lat <- f1$dim[[3]]$vals

# change colnames to something writeable
colnames(ElNino) <- "yr"
colnames(LaNina) <- "yr"

# select years 1800-2000
selyrs <- c(380:580) # indices
yrs <- c(1421:2008)

# select el nino/la nina indeces
selnino <- match(ElNino$yr,yrs)
selnina <- match(LaNina$yr,yrs)


# temp composites -------------------------

# get variables lon, lat, year, mon
temp <- ncvar_get(f1, varid="temp2")
temp2 <- array(temp, dim=c(length(lon), length(lat), 12, length(time)/12))
temp3 <- aperm(temp2, c(1,2,4,3))

# compute annual mean 
# for months: use temp3[,,,1:3]
temp4 <- apply(temp3, c(1:3), mean) # annual mean for each cell

# reference period
temp_mean <- apply(temp4, c(1,2), mean) # global mean for entire period

# composites and difference to reference period
temp_comp_nino <- apply(temp4[,,selnino], c(1,2), mean)
diff_temp_nino <- temp_comp_nino - temp_mean

# select years of choice
# temp4 <- temp4[,,selyrs]


# prec composites --------------------------
prec <- ncvar_get(f2, varid="totprec")
prec2 <- array(prec, dim=c(length(lon), length(lat), 12, length(time)/12))
prec3 <- aperm(prec2, c(1,2,4,3)) # UNIT kg/m^2 s

# compute precip sum for each year
# convert to kg/m^2 per year --> *60*60*24*365
prec4 <- apply(prec3,c(1:3), mean)
prec4 = prec4*60*60*24*365

# compute mean sum per grid cell
prec_mean <- apply(prec4, c(1,2), mean)
prec_sd <- apply(prec4, c(1,2), sd)

# composites and difference to reference period
prec_comp_nino <- apply(prec4[,,selnino], c(1,2), mean)
diff_prec_nino <- prec_comp_nino - prec_mean


# PLOTS --------------------------------------

# Cropping of Map function
# careful - if done before analysis, analysis only looks at chosen area, 
# otherwise only the map is cropped
sellat <- which(lat<20)
sel.lat.num <- lat[sellat]

x <- lon
y <- rev(lat)

# composites for el nino years ---------------------------

z <- diff_prec_nino[,rev(1:length(lat))]

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Composites Precipitation El Nino years")


# annual means --------------------------

# average temperature whole time period
z <- temp_mean[,rev(1:length(lat))]

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- seq(210, 310, 10)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Mean T all years")## average prec


# average precipitation whole time period
z <- prec_mean[,rev(1:length(lat))]

mycol <- c("#c6dbef", "#9ecae1", "#4292c6", "#2171b5", "#08519c", "#08306b")
mylevs <- seq(0,6000,1000)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Mean yearly prec [mm]")


# sd precipitation
z <- prec_sd[,rev(1:length(lat))]

mylevs <- seq(0,600,100)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Sd of yearly prec [mm]")


# correlations --------------------------------

# How does every temperature gridcell correlate with the australian precip timeseries?

# correlation of global temperature with precipitation in australia 
# random gridcell in middle australia:  -26.853388 S, 133.275154 E
clon <- 133.275154
clat <- -26.853388
dlon <- lon[2]-lon[1]
dlat <- lat[1]-lat[2]

# now selcet that gridcell: 
sellon <- (lon<(clon+(dlon/2)))&(lon>=(clon-(dlon/2)))
sellat <- (lat<(clat+(dlat/2)))&(lat>=(clat-(dlat/2)))

# timeseries precip at that gridcell
ausie <- prec4[sellon,sellat,]

# ausie has negative precip values, wtf?

# format 2 plots
par(mfrow=c(1,2))

plotfield <- apply(temp4,c(1,2),cor,ausie)
z <- plotfield[,rev(1:length(lat))]

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor of prec in Australia with world temperature")





# correlation of average Elnino temperature with precipitation in Australia
# only comparison of intensities of el nino
yr.jfm <- c(1421:2008)
selnino <- match(ElNino$yr,yr.jfm) # indices of the years
ausie <- prec4[sellon,sellat,selnino]

par(mfrow=c(1,2))

plotfield <- apply(temp4[,,selnino],c(1,2),cor,ausie) 
z <- plotfield[,rev(1:length(lat))]

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor of prec in Australia with world temperature for ElNino years")



# correlation of average LaNina temperature vs precipitation in Australia
# only comparison of intensities of la nina

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


# COMPOSITES TIMESERIES --------------------------------

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
