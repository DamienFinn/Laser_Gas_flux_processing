#Laser gas flux linearity script
#From input laser data and appropriate metadata, calculates gas fluxes from raw data and
#data cropped to remove erroneous measurements from the beginning and end of field experiments,
#where disturbances to the chamber can affect linear model fits

#The code will work with any gas, however the examples here are for methane

#Written by Damien Finn
#email: damien.finn@uqconnect.edu.au

#Name of your gas of interest and its units 
mygas <- "Methane (ppb)"

#Molecular weight of your gas
MW <- 16.04 #g mol-1

#Set working directory
dir <- ".../mydir"
setwd(dir)

#Write lazer file name
lazer_file <- "CH4_Pico.txt"

#Write mapping file name
mapping_file <- "MapFile_CH4__testing.txt"

#Set a general file name
file_name <- as.character(gsub(".txt", "", lazer_file))

#Set up outputs
lazer_data <- read.table(lazer_file, sep = ",", header = T)
mapping_data <- read.table(mapping_file, sep = "\t", header = T)
#
#print header in output tables
Z <- data.frame("Name","Start time", "End time", "Cropped start time", "Cropped end time", "Slope (ppb sec-1)", "Cropped slope (ppb sec-1)", "Delta slope", "R2", "Cropped R2", "Delta R2", "Flux (mg m2 h-1)", "Denoised flux (mg m2 h-1)", "Field Corr Flux (mg m2 h-1)", "Denoised Field Corr Flux (mg m2 h-1)", "Start gas", "End gas", "Cropped start gas", "Cropped end gas", "P value", "Cropped P value", "For window length", "Rev window lenth", "Comments")
Y <- data.frame("plot_code", "sub_plot", "plot_corner", "collar_number", "replica", "measurement_code", "treatment_code_partitioning", "disturbance_code_control", "litter_code", "year", "month", "day", "hour", "min", "instr", "chmb_type", "d_chamber","smpl_date", "raw_start_time", "raw_end_time", "Time", "CH4ppb", "total_vol_L", "collar_area_m2", "collar_hght_cm", "air_temp_c", "soil_temp_c", "atmp_mb", "atmp_Kpa","RawReads", "Flux_name", "Comments")
#set name of outputs
output_name <- paste(file_name, "DX_calc.txt", sep = ".")
second_output <- paste(file_name, "RData4Flux.txt", sep = ".")
write.table(Z, output_name, sep = "\t", row.names = F, col.names = F)
write.table(Y, second_output, sep = "\t", row.names = F, col.names = F)



#Define linearity function for methane analyses
check_linearity <- function(lazer_data){
  par(mfrow = c(2,2))
  tmp <- lazer_data[begin:end,]
  len <- (end - begin) + 1
  time <- tmp[, 1]
  time <- gsub("^.*\\s", "", time)
  time <- gsub("[.^].*", "", time)
  ticknames <- time[seq(1, len, 50)]
  x <- c(1:len)
  methane <- as.numeric(tmp$CH4..ppb.)
  Uname <- paste("Flux", name, sep = " ")
  plot(methane ~ x, xaxt = "n", ylab = mygas, xlab = "", main = Uname)
  axis(1, at = seq(1, len, 50), labels = ticknames, las =2, cex = 0.2)
  #
  lm1 <- lm(methane ~ x) #This calculates a linear model on the uncorrected data
  #
  abline(lm1, col = "red", lty = 2, lwd = 3)
  print(name)
  print(summary(lm1))
  lm1.p <- summary(lm1)$coefficients[2, 4]
  Residuals <- resid(lm1)
  plot(x, Residuals, ylab = "Residuals", xaxt = "n", xlab = "", main = "Residual deviation from model")
  axis(1, at = seq(1, len, 50), labels = ticknames, las =2, cex = 0.2)
  abline(0, 0, col = "red", lty = 2, lwd = 3)
  resid.df <- data.frame(Residuals)
  resid.start <- head(resid.df, forward.slice)
  resid.end <- tail(resid.df, reverse.slice)
  max.resid <- max(Residuals)
  min.resid <- min(Residuals)
  low.start <- subset(resid.start, Residuals < 0.1*min.resid)
  high.start <- subset(resid.start, Residuals > 0.1*max.resid)
  low.end <- subset(resid.end, Residuals < 0.1*min.resid)
  high.end <- subset(resid.end, Residuals > 0.1*max.resid)
  low.start.names <- as.numeric(rownames(low.start))
  low.end.names <- as.numeric(rownames(low.end))
  high.start.names <- as.numeric(rownames(high.start))
  high.end.names <- as.numeric(rownames(high.end))
  slice.obj <- c(low.start.names, low.end.names, high.start.names, high.end.names)
  sliced.df <- resid.df[-slice.obj,]
  sliced.time <- time[-slice.obj]
  sliced.len <- length(sliced.df)
  sliced.ticknames <- sliced.time[seq(1, sliced.len, 50)]
  w <- c(1:length(sliced.df))
  methane2 <- as.numeric(tmp$CH4..ppb.)
  sliced.methane <- methane2[-slice.obj]
  mylen <- length(sliced.time)
  #
  #Here is where the values in the listed order will be printed into RDATA4Flux table
  tmpFluxdf <- data.frame("plot_code"=rep(plotcode, mylen), "sub_plot"=rep(plotno, mylen), "plot_corner"=rep(plotcorn, mylen), "collar_number"=rep(collarno, mylen), "replica"=rep(replica, mylen), "measurement_code"=rep(meascode, mylen), "treatment_code_partitioning"=rep(treatcodepart, mylen), "disturbance_code_control"=rep(distcodepart, mylen), "litter_code"=rep(littercode, mylen), "year"=rep(year, mylen), "month"=rep(month, mylen), "day"=rep(day, mylen), "hour"=rep(hour, mylen), "min"=rep(minute, mylen), "instr"=rep(instrument, mylen), "chmb_type"=rep(chambertype, mylen), "d_chamber"=rep(name, mylen), "smpl_date"=rep(Date, mylen), "raw_start_time"=rep(Commencer, mylen), "raw_end_time"=rep(Finir, mylen), "Time"=rep(NA, mylen), "CH4ppb"=rep(NA, mylen), "total_vol_L"=rep(volume, mylen), "collar_area_m2"=rep(area,mylen), "collar_hght_cm"=rep(height, mylen), "air_temp_c"=rep(temp, mylen), "soil_temp_c"=rep(soiltemp, mylen),  "atmp_mb"=rep(pressmb, mylen), "atmp_Kpa"=rep(press, mylen), "CH4_RawReads"=rep(CH4readstxt, mylen), "CH4_flux_name"=rep(CH4fluxtxt, mylen), "Comments"=rep(Comments, mylen))
  tmpFluxdf$Time <- sliced.time
  tmpFluxdf$CH4ppb <- sliced.methane
  write.table(tmpFluxdf, second_output, sep = "\t", row.names = F, quote = F, append = T, col.names = F)
  if (length(sliced.df) < 1){
    print("Non-linear process detected.")
  }else{
    lm2 <- lm(sliced.methane ~ w) #This calculates a linear model on the corrected data
    print(name)
    print(summary(lm2))
    lm2.p <- summary(lm2)$coefficients[2,4]
    #if (lm1.p < lm2.p){
    #  print("Uncorrected data has a better fit.")
    #}else{
    Cname <- paste("Denoised Flux", name, sep = " ")
    plot(sliced.methane ~ w, xaxt = "n", ylab = mygas, xlab = "", main = Cname)
    axis(1, at = seq(1, sliced.len, 50), labels = sliced.ticknames, las =2, cex = 0.2)
    abline(lm2, col = "red", lty = 2, lwd = 3)
    plot(w, sliced.df, ylab = "Residuals", xaxt = "n", xlab = "", main = "Residual deviation from model")
    axis(1, at = seq(1, sliced.len, 50), labels = sliced.ticknames, las =2, cex = 0.2)
    abline(0, 0, col = "red", lty = 2, lwd = 3)
    start.time <- sliced.time[1]
    s <- as.POSIXlt(paste(Sys.Date(), start.time))
    Ustart.time <- time[1]
    Us <- as.POSIXlt(paste(Sys.Date(), Ustart.time))
    end.time.point <- (length(sliced.time)-1)
    end.time <- sliced.time[end.time.point]
    e <- as.POSIXlt(paste(Sys.Date(), end.time))
    Uend.time.point <- (length(time)-1)
    Uend.time <- time[Uend.time.point]
    Ue <- as.POSIXlt(paste(Sys.Date(), Uend.time))
    start.methane <- sliced.methane[1]
    end.methane.point <- (length(sliced.methane)-1)
    end.methane <- sliced.methane[end.methane.point]
    V <- volume #Volume of chamber (L)
    UCH4ppbsec <- lm1$coefficients[2]
    CCH4ppbsec <- lm2$coefficients[2]
    lm1R2 <- summary(lm1)$adj.r.squared
    lm2R2 <- summary(lm2)$adj.r.squared
    t <- temp #temperature (Celsius)
    R <- 8.3144598
    A <- area #m2
    Uflux <- (UCH4ppbsec*MW*V*press*3.6)/(R*(t + 273.13)*1000*A)
    Cflux <- (CCH4ppbsec*MW*V*press*3.6)/(R*(t + 273.13)*1000*A) #Corrected flux in mg m2 h-1
    blank <- c("")
    FieldUflux <- ((V/((101.325/press)*22.4*((t + 273.13)/298.15)))/A)*(UCH4ppbsec/1000)*MW*3.6 
    FieldCflux <- ((V/((101.325/press)*22.4*((t + 273.13)/298.15)))/A)*(CCH4ppbsec/1000)*MW*3.6 #Field corrected flux in mg m2 h-1
    tmp.df <- data.frame(name, as.character(time[1]), as.character(time[len]),  as.character(sliced.time[1]), as.character(sliced.time[sliced.len]), round(UCH4ppbsec, digits = 4), round(CCH4ppbsec, digits = 4), round(abs(UCH4ppbsec - CCH4ppbsec), digits = 4), round(lm1R2, digits = 3), round(lm2R2, digits = 3), round(abs(lm1R2 - lm2R2), digits = 4),  round(Uflux, digits = 4),  round(Cflux, digits = 4), round(FieldUflux, digits = 4), round(FieldCflux, digits = 4), as.character(round(methane[1], digits = 1)), as.character(round(methane[len], digits = 1)), round(start.methane, digits = 1), round(end.methane, digits = 1), summary(lm1)$coefficients[2,4], summary(lm2)$coefficients[2,4], as.character(forward.slice), as.character(reverse.slice), blank)
    write.table(tmp.df, output_name, sep = "\t", row.names = F, quote = F, append = T, col.names = F)
  }
}

#Choose a window at which to search for linearity in the methane curve
forward.slice <- 50 #Units in data points of lazer measurements
reverse.slice <- 50

#Apply the methane function to all the sampling points in the mapping file
for(i in 1:nrow(mapping_data)){
  name <- as.character(mapping_data$d_chamber[i])
  statement <- ""
  Date <- as.character(mapping_data$smpl_date[i])
  Commencer <- as.character(mapping_data$raw_start_time[i])
  Finir <- as.character(mapping_data$raw_end_time[i])
  Startdatetime <- paste(Date, Commencer, sep = " ")
  print(Startdatetime)
  Enddatetime <- paste(Date, Finir, sep = " ")
  print(Enddatetime)
  begin <- grep(pattern = Startdatetime, lazer_data$Time.Stamp)
  end <- grep(pattern = Enddatetime, lazer_data$Time.Stamp)
  if(length(begin) < 1){
    begin <- end - 350
  }
  if(length(end) < 1){
    end <- begin + 350
  } 
  if((length(begin) + length(end)) < 1){
    statement <- paste("***Cannot find", name, "***", sep = " ")
    print(statement)
  }else{
    temp <- mapping_data$air_temp_c[i] #Celsius
    volume <- mapping_data$total_vol_L[i] #L
    press <- mapping_data$atmp_Kpa[i] #kPa
    area <- mapping_data$collar_area_m2[i] #m2
    plotcode <- mapping_data$plot_code[i]
    plotno <- mapping_data$sub_plot[i]
    plotcorn <- mapping_data$plot_corner_code[i]
    collarno <- mapping_data$collar_number[i]
    replica <- mapping_data$replica[i]
    meascode <- mapping_data$measurement_code[i]
    treatcodepart <- mapping_data$treatment_code_partitioning[i]
    distcodepart <- mapping_data$disturbance_code_control[i]
    littercode <- mapping_data$litter_code[i]
    year <- mapping_data$year[i]
    month <- mapping_data$month[i]
    day <- mapping_data$day[i]
    hour <- mapping_data$hour[i]
    minute <- mapping_data$min[i]
    pressmb <- mapping_data$atmp_mb[i]
    instrument <- mapping_data$instr[i]
    chambertype <- mapping_data$chmb_type[i]
    height <- mapping_data$collar_hght_cm[i]
    soiltemp <- mapping_data$soil_temp_c[i]
    CH4readstxt <- mapping_data$CH4_RawReads[i]
    CH4fluxtxt <- mapping_data$CH4_flux_name[i]
    CO2fluxtxt <- mapping_data$C02_flux_Filename[i]
    Dissgas <- mapping_data$Dissolve_gas[i]
    Dissgasname <- mapping_data$Dissolve_gas_name[i]
    Comments <- mapping_data$comments[i]
    res <- check_linearity(lazer_data)
  }
}

##
##
##
##
##

#Import corrected data for R Flux package analysis and the libraries
library(flux)

#The following is for CH4 fluxes
Fluxdf <- read.table(second_output, sep = "\t", header = T)

#Need a vector of pressure in Pascals
atmp_Pa <- (Fluxdf$atmp_Kpa)*1000
Vol_cbm <- (Fluxdf$total_vol_L/1000)
#Need a vector of CH4codes 
CH4code <- data.frame("CH4code"=rep(0, length(atmp_Pa)))

#Need unique integers for x axis, so must use seconds instead of minutes
Time_secs <- data.frame("Timestamp"=rep(NA, length(atmp_Pa)))
for(i in 1:nrow(Fluxdf)){
  timestamp <- as.character(paste(Fluxdf$year[i], Fluxdf$month[i], Fluxdf$day[i], sep = "-"))
  timestamp <- as.character(paste(timestamp, Fluxdf$Time[i], sep = " "))
  Time_secs$Timestamp[i] <- timestamp
}
Fluxdf <- cbind(Fluxdf, atmp_Pa, CH4code, Time_secs, Vol_cbm)
Time_from_0 <- data.frame("Time_from_0.mins"=rep(NA, nrow(Fluxdf))) #Note, this is transformed to Time from 0 in minutes, as Flux package uses these units for its calculations
for(j in 1:nrow(Fluxdf)){
  lev <- Fluxdf$d_chamber[j]
  startlev <- match(lev, Fluxdf$d_chamber)
  n1 <- as.numeric(as.POSIXct(Fluxdf$Timestamp[startlev]))
  ni <- as.numeric(as.POSIXct(Fluxdf$Timestamp[j]))
  res <- (ni - n1)/60
  Time_from_0$Time_from_0.mins[j] <- res
}
#Create Flux dataframe with the Time_from_0 vector required for flux calculations
Fluxdf <- cbind(Fluxdf, Time_from_0)
#Issue! There is too much data per chamber. Roughly 250 points per chamber
#must reduce to make model comparisons with Flux package realistic
##  this subsets data at a fix frequency of seconds in time, so if number below says 20 is every 20 seconds and so on
subset.Fluxdf <- Fluxdf[seq(1, nrow(Fluxdf), 60),]

#Separate chambers based on names
partitioneddf <- chop(subset.Fluxdf, factors = c("plot_code", "d_chamber", "smpl_date"), nmes = c("plot_code", "d_chamber", "smpl_date")) 
#Create list of parameters for flux calculation
vp.CH4 <- list(CH4 = "CH4ppb", time = "Time_from_0.mins", volume = "Vol_cbm", t.air = "air_temp_c", area = "collar_area_m2", p.air = "atmp_Pa", CH4.gcq = "CH4code")



#Perform flux calculations
## fixing units to mg: out.unit = "mg"
## identified GHG as CH4 (default): ghg = "CH4"
flux.CH4 <- flux(partitioneddf, var.par = vp.CH4, out.unit = "mg", elementar = T)
#Plot
plot(flux.CH4, dims = c(4,8))

#Write the raw Flux table out
fluxR_Output_calc <- paste(file_name, "fluxR_Program_Output.txt", sep = "_")
write.table(flux.CH4$flux.table, fluxR_Output_calc, sep = "\t", row.names = F, col.names = T)



