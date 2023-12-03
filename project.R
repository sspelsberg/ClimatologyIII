
#   ENSO AND AUSTRALIAN PRECIPITATION VARIABILITY (before and after 1880)

# libraries ------------------------
library(ncdf4)
library(maps)
library(RColorBrewer)
library(paletteer)
library(tidyverse)
library(ggplot2)
library(stats)

# To Dos ----------------------------

# 1) crop composite map to australia to visualize local precipitation effects 
# 2) calculate own el nino index (done?)
# 3) use other months for analysis 

# 4) figure out why there are 
# --> negative precip values -> because of the moving average
# --> mean temp values >300  -> because unit is Kelvin


# read data -----------------------------------
f1 <- nc_open("../data_raw/ModE-RA_ensmean_temp2_abs_1420-2009.nc")
f2 <- nc_open("../data_raw/ModE-RA_ensmean_totprec_abs_1420-2009.nc")
LaNina <- read.table("../data_raw/LaNinaYears.txt", header = T)
ElNino <- read.csv("./data/ElNino_years_sophie.csv", header = T)

# compute selection indices etc. ----------------------

# read dimensions
time <- f1$dim[[1]]$vals
lon <- f1$dim[[2]]$vals
lat <- f1$dim[[3]]$vals

# change colnames to something writeable
# colnames(ElNino) <- "yr" ;we have our own series now
colnames(LaNina) <- "yr"

# sort El Nino
# ElNino <- ElNino$year |> arrange(yr); we have our own series now
ElNino <- data.frame(yr = as.numeric(ElNino$year))
LaNina <- LaNina |> arrange(yr)

# select years 1800-2000
yrs <- c(1421:2008)
selyrs_1800 <- which(yrs < 1800)
selyrs_2008 <- which(yrs >= 1800) # indices
selyrs_1850 <- which(yrs >= 1850) # skip period with few measurements
selyrs_1800_1850 <- which(yrs >= 1800 & yrs < 1850)

# select el nino/la nina indices
selnino <- match(ElNino$yr,yrs)
selnina <- match(LaNina$yr,yrs)
selnino1800 <- selnino[selnino %in% selyrs_1800]
selnino2008 <- selnino[selnino %in% selyrs_2008]

append(selyrs_1800, c(380:395))

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

min(prec2) # ignore negative values when looking at absolutes, they make sense when looking at deviation

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

# 1) standardize after 1850
# 2) look at what sd corresponds to modern elnino threshhold 
# 3) apply that sd threshhold to time before 1800 to define el nino
# official NOAA index: 30 years base period (implemented), 
# 3 months running mean +/- 0.5 °C for at least 5 consecutive months = el Nino

# lon lat indices
nino_34_lon <- which(lon > -170 & lon < -120)
nino_34_lat <- which(lat > -5 & lat < 5)

# annual mean temp value for nino 3.4 region
# maybe use different months?
temp_nino_34 <- apply(temp3[nino_34_lon,nino_34_lat,,1:12], 3, mean) # in K
temp_nino_34 <- tibble(yr = yrs,
                       nino34 = temp_nino_34)

# overall mean temp
mean_nino_34_1800 <- mean(temp_nino_34$nino34[selyrs_1800])
mean_nino_34_2008 <- mean(temp_nino_34$nino34[selyrs_2008])
mean_nino_34_1850 <- mean(temp_nino_34$nino34[selyrs_1850])
mean_nino_34_1800_1850 <- mean(temp_nino_34$nino34[selyrs_1800_1850])

# standardized data for 1421-1800, 1800-1850, 1850-2008
std_nino_34_1800 <- (temp_nino_34$nino34[selyrs_1800] - mean_nino_34_1800) / sd(temp_nino_34$nino34[selyrs_1800])
std_nino_34_1800_1850 <- (temp_nino_34$nino34[selyrs_1800_1850] - mean_nino_34_1800_1850) / sd(temp_nino_34$nino34[selyrs_1800_1850])
std_nino_34_1850 <- (temp_nino_34$nino34[selyrs_1850] - mean_nino_34_1850) / sd(temp_nino_34$nino34[selyrs_1850])

# compute threshhold for el nino events after 1850 --> el nino event = 0.905 sigma
nino_threshhold <- 0.5 / sd(temp_nino_34$nino34[selyrs_1850])

# add standardized column to dataframe
temp_nino_34$std = append(std_nino_34_1800, append(std_nino_34_1800_1850, std_nino_34_1850))

# calculate a moving average for a given number of years
moving_average <- function(vals, years = 30){
  return(stats::filter(vals, rep(1/years, years), sides = 2))
}

# compute running means
temp_nino_34$moving_average = moving_average(temp_nino_34$nino34)
temp_nino_34$moving_average_std = moving_average(temp_nino_34$std)

# compute residuals of std data
temp_nino_34$res_std = temp_nino_34$std - temp_nino_34$moving_average_std 


# plot original data with moving average
png("./plots/temp_nino34_movingaverage.png", width = 1000, height=700)
ggplot(temp_nino_34, aes(x=yr)) +
  geom_line(aes(y=nino34)) +
  geom_line(aes(y=moving_average)) +
  geom_segment(aes(x = 1421, xend = 1800, y = mean_nino_34_1800, yend = mean_nino_34_1800), color="red") +
  geom_segment(aes(x = 1800, xend = 2008, y = mean_nino_34_2008, yend = mean_nino_34_2008), color="blue")
dev.off()

# plot standardized data
png("./plots/temp_nino34_stand_data.png", width = 1000, height=700)
ggplot(temp_nino_34, aes(x=yr)) +
  geom_line(aes(y=std)) +
  geom_line(aes(y=moving_average_std))
dev.off()

# plot residuals (data - moving average)
png("./plots/temp_nino34_residuals.png", width = 1000, height=700)
ggplot(temp_nino_34, aes(x=yr, y=res_std)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = nino_threshhold, linetype="dashed", color = "red") +
  geom_hline(yintercept = -nino_threshhold, linetype="dashed", color = "blue") +
  geom_line()
dev.off()

# PLOTS --------------------------------------

# Cropping of Map function
# careful - if done before analysis, analysis only looks at chosen area, 
# otherwise only the map is cropped

# latitude < -8 - -43 for only Australia
# longitude          
sellat <- which(lat < -8 & lat > -43 )
sellon <- which(lon < 157 & lon > 109 )

x <- lon[sellon]
y <- rev(lat[sellat])

# composites for el nino years ---------------------------

z <- diff_prec_nino[1:length(lon[sellon]),rev(1:length(lat[sellat]))]

png("./plots/prec_composites_elninoyrs.png", width = 1000, height=700)
mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Composites Precipitation El Nino years")
dev.off()

# annual means --------------------------

# average temperature whole time period

z <- temp_mean[1:length(lon[sellon]),rev(1:length(lat[sellat]))]

png("./plots/temp_average_allyr.png", width = 1000, height=700)
mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- seq(210, 310, 10)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Mean T all years")
dev.off()

# average precipitation whole time period

z <- prec_mean[1:length(lon[sellon]),rev(1:length(lat[sellat]))]

png("./plots/prec_average_allyr.png", width = 1000, height=700)
mycol <- c("#c6dbef", "#9ecae1", "#4292c6", "#2171b5", "#08519c", "#08306b")
mylevs <- seq(0,6000,1000)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Mean yearly prec [mm]")
dev.off()

# sd precipitation

z <- prec_sd[1:length(lon[sellon]),rev(1:length(lat[sellat]))]

png("./plots/prec_sd_allyr.png", width = 1000, height=700)
mylevs <- seq(0,600,100)
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main="Sd of yearly prec [mm]")
dev.off()

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



# timeseries average precip at those gridcells BEFORE 1800
ausie <- prec4[sellon,sellat,selyrs_1800]

par(mfrow=c(1,2))

plotfield <- apply(temp4[,,selyrs_1800],c(1,2),cor,ausie)

x <- lon
y <- rev(lat)
z <- plotfield[,rev(1:length(lat))]

png("./plots/correlation_prec_australia_global_temp_1800.png", width = 1000, height=700)

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor of prec in Australia with world temperature (1800)")

dev.off()

# timeseries average precip at those gridcells AFTER 1800
ausie <- prec4[sellon,sellat,selyrs_2008]

par(mfrow=c(1,2))

plotfield <- apply(temp4[,,selyrs_2008],c(1,2),cor,ausie)

z <- plotfield[,rev(1:length(lat))]

png("./plots/correlation_prec_australia_global_temp_2008.png", width = 1000, height=700)

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor of prec in Australia with world temperature (2008)")

dev.off()



# correlation of Elnino intensity with average precipitation in Australia BEFORE 1800
# get prec in el nino years for that grid cell
ausie <- prec4[sellon,sellat,selnino[selnino<380]]

par(mfrow=c(1,2))

plotfield <- apply(temp4[,,selnino[selnino<380]],c(1,2),cor,ausie) 

z <- plotfield[,rev(1:length(lat))]

png("./plots/correlation_prec_australia_ElNinointensity_1800.png", width = 1000, height=700)

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor ElNino intensity with prec in Australia")

dev.off()

# correlation of Elnino intensity with average precipitation in Australia AFTER 1800
# get prec in el nino years for that grid cell
ausie <- prec4[sellon,sellat,selnino[selnino>=380]]

par(mfrow=c(1,2))

plotfield <- apply(temp4[,,selnino[selnino>=380]],c(1,2),cor,ausie) 

z <- plotfield[,rev(1:length(lat))]

png("./plots/correlation_prec_australia_ElNinointensity_2008.png", width = 1000, height=700)

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor ElNino intensity with prec in Australia")

dev.off()



# correlation of LaNina intensity vs precipitation in Australia BFORE 1800
ausie <- prec4[sellon,sellat,selnina[selnina<380]]

par(mfrow=c(1,2))

plotfield <- apply(temp4[,,selnina[selnina<380]],c(1,2),cor,ausie)

z <- plotfield[,rev(1:length(lat))]

png("./plots/correlation_prec_australia_LaNinaintensity_1800.png", width = 1000, height=700)

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor LaNina intensity with prec in Australia")

dev.off()



# correlation of LaNina intensity vs precipitation in Australia AFTER 1800
ausie <- prec4[sellon,sellat,selnina[selnina>=380]]

plotfield <- apply(temp4[,,selnina[selnina>=380]],c(1,2),cor,ausie)

z <- plotfield[,rev(1:length(lat))]

png("./plots/correlation_prec_australia_LaNinaintensity_2008.png", width = 1000, height=700)

mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <- max(max(z),abs(min(z)))*(c(0:10)-5)/5
filled.contour(x,y,z,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)}, main = "Cor LaNina intensity with prec in Australia")

dev.off()

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
png("./plots/timeseries_prec_elninoyrs.png", width = 1000, height=700)
ggplot(df_prec_comp, aes(x=yr, y=comps, group=nino_yr, color=nino_yr)) + # line plot by nino event
  geom_line() +
  facet_wrap(~period) # distinguish period
dev.off()


# save Sophies ElNino yrs in csv for above analysis
timeseries_elnino <- data.frame(index = 1:length(unique(df_prec_comp$nino_yr)),
                                year = unique(df_prec_comp$nino_yr))
write.csv(timeseries_elnino, "./data/ElNino_years_sophie.csv")





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
