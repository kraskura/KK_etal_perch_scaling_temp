if (!require("rspatial")) remotes::install_github('rspatial/rspatial')

# install and load the needed packages
# install.packages("pacman")
pacman::p_load(ggplot2, raster, usmap, tidyverse)
library(rspatial)

# Santa Barbara county : 06083 (3 islands)
SB.fip<-"06083" 
# Ventura county: 06111 (2 islands)
Ventura.fip<-"06111"
# LA county: 06037 (1 island)
LA.fip<-"06037"

# function to make adjustable alpha color level to regular plot function. 
add.alpha <- function(col=NULL, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


# required data to plot
cn <- sp_data('pt_countries') # for this to work, use remote github install of a package ("R version 4.2.0 (2022-04-22)")
v<-read.csv("./Range-map/Map.csv")
v2<-read.csv("./Range-map/Collection_sites.csv")

dir.create(path = "./Range-map/Figures", recursive = T)
# 
png(paste( "./Range-map/Figures/Fig1A_rangeMAP",Sys.Date(),".png",sep=""), width = 8, height = 8, units = "in", res = 300)
  plot(cn, xlim=c(-125, -120), ylim=c(25,45), axes=TRUE)
  points(y = v$Center.Lat[v$Overall.Probability>0.90], x = v$Center.Long[v$Overall.Probability>0.90], cex=1.5, col=add.alpha("#11973A",0.8), pch = 19)
  # points(y = v2$Lat[2:4], x = v2$Lon[2:4], cex=1.5, col="brown", pch = "*")
  points(y = v$Center.Lat[v$Overall.Probability<0.90 & v$Overall.Probability>0.50], x = v$Center.Long[v$Overall.Probability<0.90 & v$Overall.Probability>0.50], cex=1.5, col=add.alpha("#ADB00B",0.8), pch = 19)
  points(y = v$Center.Lat[ v$Overall.Probability<0.50], x = v$Center.Long[v$Overall.Probability<0.50 ], cex=1.5, col=add.alpha("#ECB51F",0.8), pch = 19)
  points(y = v2$Lat[1], x = v2$Lon[1], cex=2, lwd=2, col="black", bg = "red", pch = 23)
dev.off()
# ggsave(filename = paste( "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Perch /Respiromtery/Figs/Fig1_CAmap",Sys.Date(),".png",sep=""), width = 3, height = 3, units = "in")


# counties plot
df <- data.frame(
  fips = c(SB.fip, Ventura.fip, LA.fip), 
  values = c(1,1,1)
)

png(paste( "./Range-map/Figures/Figure1A_insetCA",Sys.Date(),".png",sep=""), width = 5, height = 5, units = "in", res = 300)
  plot_usmap(data = df, include =  "CA",  alpha = 0.25, color = "black")+
  theme(legend.position = "right")+
  scale_fill_continuous(
    low = "white", high = "white")
dev.off()  # theme(legend.position = "none")


