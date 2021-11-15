#------------------------------------------------------------------------------------------------------------#
############################## Create plots indicating vegetation changes ####################################
#------------------------------------------------------------------------------------------------------------#
##### Author: Lasse T. Keetz
##### Date:   2021-10-12 [YYYY-MM-DD]
##### Description: Reads in LPJ-GUESS Education output (.out) files and creates various output variable plots
#####              to analyze vegetation dynamics in the grid cell
#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

### Set working directory to where you stored the ".out" files ### 
# OBS! Do not remove the r"()" ('r' and parentheses)! --> R syntax for "raw strings" #
setwd(
  r"(C:\Users\lassetk\OneDrive - Universitetet i Oslo\Teaching\GEO9915_DM_DGVM\DGVM\Lab\Guesswork\Testrun_tropics\)"
  )

### Set number of spinup years you used in your simulation ###
yrs_spinup <- 1000

### Now you can run the entire script using CTRL+ALT+A, or click on "Code"->"Run Region"->"Run All"

#------------------------------------------------------------------------------------------------------------#
########################################### Data preparation #################################################
#------------------------------------------------------------------------------------------------------------#

### Read output files ###
# Annual net primary production #
anpp <- read.table("anpp.out", header=T, sep="", dec=".") 
# Leaf area index #
lai <- read.table("lai.out", header=T, sep="", dec=".") 
### Repeat for different output files if you want! ### 
#...#


### Extract data for transient period ###
# Years from transient period representing "real years" for plotting #
yrs_transient <- anpp$Year[anpp$Year > yrs_spinup]
# yrs_corresp <- (yrs_transient - yrs_spinup) + 2000 
# This is assuming the spin-up phase represents "modern-day"

# Subset data, omitting spin up phase #
anpp_transient <- subset(anpp,
                         anpp$Year>yrs_spinup) 

# Subset data, 200 years before transient phase #
start_year_plot_spinup <- yrs_spinup - 200
anpp_spinup_subset <- subset(anpp,
                             anpp$Year>=(start_year_plot_spinup) & anpp$Year<=yrs_spinup) 

#------------------------------------------------------------------------------------------------------------#
############################################ Plot annual NPP #################################################
#------------------------------------------------------------------------------------------------------------#

## Plot time series of Total ANPP ## 
# Plot margins #
png(filename = "Simulation_TotalNPP.png",
    width = 600,
    height = 600)
par(mar = c(5, 6, 2, 2))
# Plot #
col_scenario = "#50C87890"
plot(yrs_transient, anpp_transient$Total,
     cex=2,cex.axis= 1.6, cex.lab = 2.5,
     xlim=c(start_year_plot_spinup, max(yrs_transient)),
     ylim=c(0, max(anpp_transient$Total)),
     type="p",pch=21,
     col="#00000060", bg=col_scenario,
     xlab="Simulation year",
     ylab=expression(paste("Total annual NPP [Kg C m"^"-2","]")))

col_spinup = "#3A3B3C30"
points(start_year_plot_spinup:yrs_spinup, anpp_spinup_subset$Total,
       cex=2, type="p", pch=21,
       col="#00000030", bg=col_spinup)

### Add linear regression to investigate trends ### 
lm <- lm(anpp_transient$Total~yrs_transient)
col_trendline = "#DC143C70"
abline(lm, col=col_trendline, lwd=2.5)

# Explore the 'lm' object for additional statistics by removing the comment
# in the next line!
# summary(lm)

# Add vertical line to seperate scenario and spinup phases
abline(v = yrs_spinup, col = "black")
mtext("Spin-up phase", side=3, adj=0, cex=1.5)#x = start_year_plot_spinup+10, y = 0+max(anpp_transient$Total)*0.05,
     #"Spin-up phase")
mtext("Scenario phase", side=3, adj=1, cex=1.5)

legend("topleft", 
       legend=c("ANPP - scenario", "ANPP - spinup", "Lin. trend for scenario"),
       col=c("black", "black", col_trendline),
       pch=c(21,21,NA),lty=c(NA,NA,1), pt.bg = c(col_scenario,col_spinup),
       cex = 1.5)

# close png
dev.off()


#------------------------------------------------------------------------------------------------------------#
##################################### Calculate total LAI 30 year averages ###################################
#------------------------------------------------------------------------------------------------------------#

### Which PFTs are present? ###
# Calculate average PFT LAI at beginning and end of transient period #
n_years <- 10

beg <- (min(yrs_transient)-n_years):min(yrs_transient)
end <- yrs_transient[yrs_transient>=max(yrs_transient)-n_years]

# Determine corresponding indices in lai.out file #
beg.ind <- which(lai$Year %in% beg)
end.ind <- which(lai$Year %in% end)

pft_lai_beg <- list()
pft_lai_end <- list()

### Loop through all the individual PFTs in file ###
for(i in 1:length(lai)){
  if(names(lai[i])!="Lon" &&
     names(lai[i])!="Lat" && 
     names(lai[i])!="Year" && 
     names(lai[i])!="Total"){
    
    ### Test which PFTs are present in beginning and end period ###
    # LAI values for current pft... #
    cur = lai[i][,1]
    
    # Calculate mean values at previously determined indices #
    if(mean(cur[beg.ind]) > 0){
      pft_lai_beg[[length(pft_lai_beg)+1]] = mean(cur[beg.ind])
      # Set name of PFT #
      names(pft_lai_beg)[length(pft_lai_beg)] = names(lai[i])
    }
    
    # Same for end period... #
    if(mean(cur[end.ind]) > 0){
      pft_lai_end[[length(pft_lai_end)+1]] = mean(cur[end.ind])
      # Set name of PFT #
      names(pft_lai_end)[length(pft_lai_end)] = names(lai[i])
    }
    
  }
}

##############################################

# Store the names of all PFTs that are contained in either df
all_PFTS = union(names(pft_lai_end), names(pft_lai_beg))
# Name for the time periods to group by
periods = c(paste0("Last ",n_years," years of spinup phase"),
            paste0("Last ",n_years," years of scenario phase"))

grouped_matrix = matrix(data=NA,
                        nrow=length(periods),
                        ncol=length(all_PFTS),
                        dimnames=list(
                          periods, all_PFTS
                        ))

### Fill matrix ###

for(pft in all_PFTS){
  
  # Start with spinup phase
  if(pft %in% names(pft_lai_beg)){
    grouped_matrix[periods[1], pft] = as.numeric(pft_lai_beg[pft])
  } else {
    grouped_matrix[periods[1], pft] = 0
  }
  
  # Now the scenario phase
  if(pft %in% names(pft_lai_end)){
    grouped_matrix[periods[2], pft] = as.numeric(pft_lai_end[pft])
  } else {
    grouped_matrix[periods[2], pft] = 0
  }
}

#------------------------------------------------------------------------------------------------------------#
########### Plot mean LAI for end of the simulation phases, averaged over n_years ############################
#------------------------------------------------------------------------------------------------------------#

png(filename = "compare_phases_mean_lai.png",
    width = 600,
    height = 600)

par(mar = c(5, 5, 2, 2))

barplot(height = grouped_matrix,
        xlab = "Plant Functional Type",
        ylab = paste0("Mean LAI - ", n_years, " year average"),
        cex.axis= 1.5, cex.names= 1.5, cex.lab = 2, 
        args.legend = list(cex = 1.5),
        col = c("#3A3B3C20", "#1261A0"),
        legend.text = rownames(grouped_matrix),
        beside = TRUE,
        ylim = c(0, max(grouped_matrix)+1))
dev.off()

## Compare corresponding output files from your other simulations. ## 

## Explore the files further with plot(), hist() etc. ## 


#------------------------------------------------------------------------------------------------------------#
############################################### END ##########################################################
#------------------------------------------------------------------------------------------------------------#