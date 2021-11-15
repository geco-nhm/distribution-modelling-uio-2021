#------------------------------------------------------------------------------------------------------------#
############################## Create LPJ-GUESS climate input diagrams #######################################
#------------------------------------------------------------------------------------------------------------#
##### Author: Lasse T. Keetz
##### Date:   2021-10-12 [YYYY-MM-DD]
##### Description: Reads in LPJ-GUESS Education climate input .txt files created by GetClim
#####              and produces basic figures to visualize the data.
#####              Obs! Requires R version >= 4.* to run! Update if necessary!!
#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

### Set file paths and names ### 

### Set working directory, must be the folder where you store the .txt climate file you want to analyze ###
# OBS! Do not remove the r"()" ('r' and parentheses)! #
setwd(r"(C:\Users\lassetk\OneDrive - Universitetet i Oslo\Teaching\GEO9915_DM_DGVM\DGVM\Lab\Guesswork\Testrun_tropics\)")

### Set name of the climate input file you want to analyze, see example ###
climfile_name <-
  "climate_tropics_rcp8-5.txt"

### For convenience, add ".txt" ending in case it's not included above
if(!endsWith(climfile_name, '.txt')) {
  climfile_name = paste0(climfile_name,".txt")
}

### Number of spinup years you used in your simulation ###
yrs_spinup <- 1000

#------------------------------------------------------------------------------------------------------------#
####################################### Load and prepare data ################################################
#------------------------------------------------------------------------------------------------------------#

### Choose one or more types of output from one of your simulations. Read in output file. ### 
# Annual net primary production #
climdata <- read.table(climfile_name, header=T, sep="", dec=".", skip = 8) 

# Create plot for transient phase and 200yrs of spinup
yrs_plot <- climdata$Year[climdata$Year >= (yrs_spinup - 100)]
yrs_ind <- which(climdata$Year %in% yrs_plot)


#------------------------------------------------------------------------------------------------------------#
############################################# Plot MAT #######################################################
#------------------------------------------------------------------------------------------------------------#

### Create data frame to store values ###
avg_tmp_df <- data.frame(year = yrs_plot, meanT = 0)

### Loop through defined period ###
for(i in 1:length(yrs_ind)){
  
  ### Calculate annual mean temperature in each year
  avg_tmp_df$meanT[i] <-
    mean(c(
      climdata$Jan[yrs_ind[i]],
      climdata$Feb[yrs_ind[i]],
      climdata$Mar[yrs_ind[i]],
      climdata$Apr[yrs_ind[i]],
      climdata$May[yrs_ind[i]],
      climdata$Jun[yrs_ind[i]],
      climdata$Jul[yrs_ind[i]],
      climdata$Aug[yrs_ind[i]],
      climdata$Sep[yrs_ind[i]],
      climdata$Oct[yrs_ind[i]],
      climdata$Nov[yrs_ind[i]],
      climdata$Dec[yrs_ind[i]]
    ))
  
}

## Plot time series of mean annual temperature ##
# Note: gsub function removes the ".txt" ending
png(filename = paste0(gsub('.{0,4}$', '', climfile_name),
                      "_MAT.png"),
    width = 1100,
    height = 600)

par(mar = c(5, 6, 2, 2))
# Plot #
plot(avg_tmp_df,
     cex=2, cex.axis= 1.6, cex.lab = 2.5,
     type="o", pch=21, col="black", bg="#DC143C",
     xlab="Simulation year",
     ylab="Mean annual temperature [deg. C]")

# Add vertical line to seperate scenario and spinup phases
abline(v = yrs_spinup, col = "black")
mtext("Spin-up phase", side=3, adj=0, cex=1.5)#x = start_year_plot_spinup+10, y = 0+max(anpp_transient$Total)*0.05,
#"Spin-up phase")
mtext("Scenario phase", side=3, adj=1, cex=1.5)

# close png
dev.off()



#------------------------------------------------------------------------------------------------------------#
############################################## Plot MAP ######################################################
#------------------------------------------------------------------------------------------------------------#

### Create data frame to store values ###
tot_prec_df <- data.frame(year = yrs_plot, MAP = 0)

### Loop through defined period ###
for(i in 1:length(yrs_ind)){
  
  ### Calculate total annual precipitation in each year
  tot_prec_df$MAP[i] <-
    sum(c(
      climdata$Jan.1[yrs_ind[i]],
      climdata$Feb.1[yrs_ind[i]],
      climdata$Mar.1[yrs_ind[i]],
      climdata$Apr.1[yrs_ind[i]],
      climdata$May.1[yrs_ind[i]],
      climdata$Jun.1[yrs_ind[i]],
      climdata$Jul.1[yrs_ind[i]],
      climdata$Aug.1[yrs_ind[i]],
      climdata$Sep.1[yrs_ind[i]],
      climdata$Oct.1[yrs_ind[i]],
      climdata$Nov.1[yrs_ind[i]],
      climdata$Dec.1[yrs_ind[i]]
    ))
  
}

## Plot time series of total precipitation ## 
png(filename = paste0(gsub('.{0,4}$', '', climfile_name),
                      "_MAP.png"),
    width = 1100,
    height = 600)
par(mar = c(5, 6, 2, 2))
# Plot #
plot(tot_prec_df,
     ylim = c(0,max(tot_prec_df$MAP)+100),
     cex=2, cex.axis= 1.6, cex.lab = 2.5,
     type="o", pch=21, col="black", bg="#1261A0",
     xlab="Simulation year",
     ylab="Total annual precipitation [mm]")

# Add vertical line to seperate scenario and spinup phases
abline(v = yrs_spinup, col = "black")
mtext("Spin-up phase", side=3, adj=0, cex=1.5)#x = start_year_plot_spinup+10, y = 0+max(anpp_transient$Total)*0.05,
#"Spin-up phase")
mtext("Scenario phase", side=3, adj=1, cex=1.5)

# close png
dev.off()


#------------------------------------------------------------------------------------------------------------#
############################################## Plot CO2 ######################################################
#------------------------------------------------------------------------------------------------------------#

## Plot time series of ambient CO2 concentration ## 
png(filename = paste0(gsub('.{0,4}$', '', climfile_name),
                      "_CO2.png"),
    width = 1100,
    height = 600)
par(mar = c(5, 6, 2, 2))
# Create plot #
plot(yrs_plot, climdata$CO2.ppm[yrs_ind],
     ylim = c(200,1000),
     cex=2, cex.axis= 1.6, cex.lab = 2.5,
     type="o", pch=21, col="black", bg="#3A3B3C",
     xlab="Simulation year",
     ylab=expression("Atmospheric CO"[2]*" [ppm]"))

# Add vertical line to seperate scenario and spinup phases
abline(v = yrs_spinup, col = "black")
mtext("Spin-up phase", side=3, adj=0, cex=1.5)#x = start_year_plot_spinup+10, y = 0+max(anpp_transient$Total)*0.05,
#"Spin-up phase")
mtext("Scenario phase", side=3, adj=1, cex=1.5)

# close png
dev.off()


#------------------------------------------------------------------------------------------------------------#
############################################### END ##########################################################
#------------------------------------------------------------------------------------------------------------#
