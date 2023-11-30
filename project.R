
#   ENSO AND AUSTRALIAN PRECIPITATION VARIABILITY (before and after 1880)

# libraries ------------------------
library(ncdf4)
library(maps)
library(RColorBrewer)
library(paletteer)
library(tidyverse)
library(ggplot2)

# To Dos ----------------------------

# 1) crop composite map to australia to visualize local precipitation effects
# 2) calculate own el nino index
# 3) use other months for analysis

# 4) figure out why there are 
# --> negative precip values
# --> mean temp values >300


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

# sort El Nino
ElNino <- ElNino |> arrange(yr)
LaNina <- LaNina |> arrange(yr)

# select years 1800-2000
yrs <- c(1421:2008)
selyrs_1800 <- which(yrs < 1800)
selyrs_2008 <- which(yrs >= 1800) # indices

# select el nino/la nina indices
selnino <- match(ElNino$yr,yrs)
selnina <- match(LaNina$yr,yrs)
selnino1800 <- selnino[selnino %in% selyrs_1800]
selnino2008 <- selnino[selnino %in% selyrs_2008]



# temp composites -------------------------

# get variables lon, lat, year, mon
temp <- ncvar_get(f1, varid="temp2")
temp2 <- array(temp, dim=c(length(lon), length(lat), 12, length(time)/12))
temp3 <- aperm(temp2, c(1,2,4,3))

# compute annual mean 
# for months: use temp3[,,,1:3]
temp4 <- apply(temp3[,,,1:12], c(1:3), mean) # annual mean for each cell

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


# El nino index ----------------------

# lon lat indices
nino_34_lon <- which(lon > -170 & lon < -120)
nino_34_lat <- which(lat > -5 & lat < 5)

# annual mean temp value for nino 3.4 region
# maybe use different months?
temp_nino_34 <- apply(temp3[nino_34_lon,nino_34_lat,,1:12], c(1:3), mean)


temp_nino_34


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

# version 1 -----------------------

# Erstellen der Timeseries aus Mittelwerten für australien 
# Range der Koordinaten:
# S -30.52     -     -19.109
# E 123.11     -      142.34
l <-  lat >= -30.52 & lat <= -19.109
b <-  lon >= 123.11 & lon <= 142.34


# precip and precip difference (split in 1800)
df_prec_time <-  tibble(yr = 1421:2008,
                     prec = apply(prec4[b,l,], 3, mean),
                     prec_diff_1800 = prec - mean(prec[selyrs_1800]),
                     prec_diff_2008 = prec - mean(prec[selyrs_2008])
                     )


# make df of composites
df_prec_comp <- tibble()
compyrs <- c(-2:5) # years for the analysis

# add the years before 1800
for (i in selnino1800) {
  comp_indices <- compyrs + i # indices for the years before and after the event
  
  # calculate composites
  composite_df <- tibble(yr = c(-2:5), # comp years
                         nino_yr = yrs[i], # el nino event
                         comps = df_prec_time$prec_diff_1800[comp_indices],  # precip values
                         period = "before 1800")
  
  # add to complete dataframe
  df_prec_comp <- rbind(df_prec_comp, composite_df)
}

# add the years after 1800
for (i in selnino2008) {
  comp_indices <- compyrs + i # indices for the years before and after the event
  
  # calculate composites
  composite_df <- tibble(yr = c(-2:5), # comp years
                         nino_yr = yrs[i], # el nino event
                         comps = df_prec_time$prec_diff_2008[comp_indices],  # precip values
                         period = "after 1800")
  
  # add to complete dataframe
  df_prec_comp <- rbind(df_prec_comp, composite_df)
}

# if we want the color in the plot to be discrete
# df_prec_comp$nino_yr <- factor(df_prec_comp$nino_yr, levels = ElNino$yr)
df_prec_comp$period <- factor(df_prec_comp$period, levels = c("before 1800", "after 1800")) # sorts the facet wrap

# plot
ggplot(df_prec_comp, aes(x=yr, y=comps, group=nino_yr, color=nino_yr)) + # line plot by nino event
  geom_line() +
  facet_wrap(~period) # distinguish period

# version 2 ------------------------

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
  ind = which(df_prec_time$yr == ElNino[i])
  time_cuts[1,i] = df_prec_time$prec[ind]
  time_cuts[2,i] = df_prec_time$prec[ind+1]
  time_cuts[3,i] = df_prec_time$prec[ind+2]
  # time_cuts[4,i] = df_prec_time$prec[ind+3]
  # time_cuts[5,i] = df_prec_time$prec[ind+4]
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
