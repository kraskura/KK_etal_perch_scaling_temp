
library(tidyverse)
library(ggplot2)
library(corrplot)
library(ggridges)
library(car)
library(reshape)
library(gridExtra)
library(ggpubr)
library(lubridate)


# USED SBC LTER DATA CITATION:
# Santa Barbara Coastal LTER, L. Washburn, C. Gotschalk, and D. Salazar. 2022. SBC LTER: Ocean: Currents and Biogeochemistry: Moored CTD and ADCP data from Naples Reef Mooring (NAP), ongoing since 2001 ver 30. Environmental Data Initiative. https://doi.org/10.6073/pasta/11ce26bb3ab7afa41a2f4bfb7836cbca (Accessed dataset: naples_mooring_nap_20220330.csv).

# 1. Set up work space and read in data: -----
# for perch: 
setwd(dir = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/data-LTER/Temp_Naples_Perch/")

dir.create(path = "../SBC-LTER-NAPL-temp/Figures/RateChange/", showWarnings = TRUE, recursive = TRUE)
dir.create(path = "../SBC-LTER-NAPL-temp/ExportedData/RateChange/", showWarnings = TRUE, recursive = TRUE)

data.n<-read.csv("../Data/naples_mooring_nap_20220330.csv")

# 2. Dataset organisation --------------
data.n<-data.n[, c("year", "month", "day","decimal_time","Temp_adcp", "Temp_top", "Temp_mid", "Temp_bot", "Temperature")]

str(data.n)

data.n <- data.n[!(data.n$Temp_top == 9999 | 
                     data.n$Temp_mid == 9999 |
                     data.n$Temp_bot == 9999 |
                     data.n$Temperature == 9999), ]
data.n<-data.n[, c("year", "month", "day", "decimal_time", "Temp_top", "Temp_mid", "Temp_bot", "Temperature")]
# data.h <- data.h[!(data.h$temp_c == 9999), ]


length(which(data.n[]==9999)) # == 0 (good)

# change fornat date, add timezone correctly 
x.n<-as.POSIXct(as.Date(paste(data.n$year, "-", data.n$month, "-", data.n$day, sep = ""), tzone = "PDT")) +3600*5 + 3600*24*data.n$decimal_time
x.n<-force_tz(x.n, tzone = "America/Los_Angeles")
TIME.n<-strptime(x.n, "%Y-%m-%d %H:%M:%S")
data.n$TIME<-TIME.n
data.n$SITE<-as.factor("NAPL")


# for exploring the data; useful time-stamps as factorial variables .. 
data.n$YEAR<-as.numeric(data.n$year)
data.n$MONTH<-as.numeric(data.n$month)
data.n$DAY<-as.numeric(data.n$day)
data.n$HR<-as.numeric(as.character(substr(data.n$TIME, start = 12, stop = 13)))
data.n$MIN<-as.numeric(as.character(substr(data.n$TIME, start = 15, stop = 16)))

data.n$MO_DAY<-substr(data.n$TIME, start=6, stop=10)
data.n$MO_DAY_HR<-paste(data.n$MO_DAY,"_", data.n$HR, sep="")
data.n$DATE_LOCAL<-paste(data.n$YEAR,"-", data.n$MONTH, "-", data.n$DAY, sep="")

data.n$TEMP_C<-data.n$Temperature
data.n$TEMP_C_top<-data.n$Temp_top
data.n$TEMP_C_mid<-data.n$Temp_mid
data.n$TEMP_C_bot<-data.n$Temp_bot

# location of temp mooring 
data.n_top<-data.n[, c("DATE_LOCAL", "SITE", "YEAR", "MONTH", "DAY", "HR", "MIN", "MO_DAY", "MO_DAY_HR", "TEMP_C_top", "TIME")]
data.n_bot<-data.n[, c("DATE_LOCAL", "SITE", "YEAR", "MONTH", "DAY", "HR", "MIN", "MO_DAY", "MO_DAY_HR", "TEMP_C_bot", "TIME")]
data.n_mid<-data.n[, c("DATE_LOCAL", "SITE", "YEAR", "MONTH", "DAY", "HR", "MIN", "MO_DAY", "MO_DAY_HR", "TEMP_C_mid", "TIME")]
data.n<-data.n[, c("DATE_LOCAL", "SITE", "YEAR", "MONTH", "DAY", "HR", "MIN", "MO_DAY", "MO_DAY_HR", "TEMP_C", "TIME")]

# Explore, summarize data by day, by site, all years together
temp_sum.n <- data.n %>%
  dplyr:::group_by(DATE_LOCAL) %>%
  dplyr::summarize(mean_temp = mean(TEMP_C), min_temp = min(TEMP_C), max_temp = max(TEMP_C), var_temp = var(TEMP_C),
                   sd_temp = sd(TEMP_C),
                   mean_YEAR = mean (YEAR), mean_DAY = mean (DAY), mean_MO = mean(MONTH), mean_HR = mean (HR) )

# by year and date
temp_sum_Y.n <- data.n %>%
  dplyr:::group_by(DATE_LOCAL, YEAR) %>%
  dplyr:::summarize(mean_temp = mean(TEMP_C), min_temp = min(TEMP_C), max_temp = max(TEMP_C), var_temp = var(TEMP_C),
                   sd_temp = sd(TEMP_C),
                   mean_YEAR = mean (YEAR), mean_DAY = mean (DAY), mean_MO = mean(MONTH), mean_HR = mean (HR) , .groups = "keep")

temp_sum.n$MO_DAY<-as.factor(substr(temp_sum.n$DATE_LOCAL, start=6, stop=10))
temp_sum_Y.n$MO_DAY<-as.factor(substr(temp_sum_Y.n$DATE_LOCAL, start=6, stop=10))
temp_sum.n<-as.data.frame(temp_sum.n)
temp_sum_Y.n<-as.data.frame(temp_sum_Y.n)

# all data that are within 12 C +/- 0.5 of selected temperature treatments 
data.n20<-data.n[c(data.n$TEMP_C >= 19.5 & data.n$TEMP_C <= 20.5),] # 20ºC
data.n20<-data.n[c(data.n$TEMP_C >= 15.5 & data.n$TEMP_C <= 16.5),] # 16ºC
data.n12<-data.n[c(data.n$TEMP_C <= 12.5),] # 12ºC
data.n22<-data.n[c(data.n$TEMP_C >= 21.5),] # 22ºC

# note temp category
data.n12$TEMP_categ<-"12C"
data.n20$TEMP_categ<-"20C"
data.n22$TEMP_categ<-"22C"
data.n$TEMP_categ<-"all"
# data9999$TEMP_categ<-"15-17C" # control/ acclim

# dataset with selected treatment categories
data_T_cat.n<-rbind(data.n12, data.n20, data.n22, data.n)


# 3. Function to estimate temp rate changes ------

rate_T_change<-function(data, site , year, month, dataset_name, save.plot = TRUE, save.data=TRUE){
  
  test_data0 <- data[c(data$SITE == site & data$YEAR == year & data$MONTH == month), ]
  test_data0$SITE<-factor(test_data0$SITE)
  
  d1_temp<-as.data.frame(matrix(nrow=0, ncol=17))
  colnames(d1_temp)<-c("SITE", "SERIAL", "DATE_LOCAL", "TIME_LOCAL", "TEMP_C", "YEAR", "MONTH", "DAY", "MO_DAY", "HR", "MO_DAY_HR", "MIN", "TIME", "TIME_dec", "diff_time", "diff_temp_c", "diff_hourly_total")
  
  
  d1 <- test_data0
  d1$TIME_dec <- as.numeric(d1$TIME, origin="1970-01-01",tz="PDT")
  d1<-d1[order(d1$TIME_dec),]
  
  d1$diff_time<-NA
  d1$diff_temp_c<-NA
  d1$changeRate_CperH<-NA
  
  for(j in 2:nrow(d1)){
    d1$diff_time[j]<-(as.numeric(d1$TIME_dec[j]) - as.numeric(d1$TIME_dec[j-1])) / 60  /60 # h
    d1$diff_temp_c[j]<-as.numeric(d1$TEMP_C[j] - d1$TEMP_C[j-1]) # temp change
    d1$changeRate_CperH[j]<- d1$diff_temp_c[j]/ d1$diff_time[j]
    # d1$diff_hourly[j]<- d1$diff_temp_c[j] * (d1$diff_time[j]/60)
  }
  
  # **********************************************************
  # **********************************************************
  # if (i == n-h_delta){
  g2Peaks<-0 # peak counter = runnig count of censecutive peaks that are abeve the threshold
  peakThresh<-2
  peakStreak<-3 
  inPeak<-FALSE
  streaksVec<-c()
  thisStreak<-c()
  maxStreak<-0
  maxStreakInd<-NA
  
  for(i in 1:nrow(d1)){
    if(is.na(d1$changeRate_CperH[i])){next}
    # print(d1$changeRate_CperH[i])
    if(d1$changeRate_CperH[i]>=peakThresh){
      g2Peaks<-g2Peaks+1
    }else{
      g2Peaks<-0
    }
    if(g2Peaks>maxStreak){
      maxStreak<-g2Peaks
      maxStreakInd<-i
    }
    
    if(g2Peaks>=peakStreak){
      # print('peaked')
      if(inPeak==FALSE){
        thisStreak<-c(i-1)
        # print('initlist')
      }
      inPeak<-TRUE
      thisStreak<-append(thisStreak,i)
    }
    if(g2Peaks==0){
      # print('outofpeak')
      if(inPeak==TRUE){
        # print('save this peak')
        streaksVec<-append(streaksVec,thisStreak)
      }
      inPeak<-FALSE
    }
  }
  
  d1_maxStreak<-d1[(maxStreakInd-maxStreak-10) : (maxStreakInd+10),]
  # **********************************************************
  # **********************************************************
  
  # if (i == n-h_delta){
  g2Peaks<-0 # peak counter = runnig count of censecutive peaks that are abeve the threshold
  peakThresh<--2
  peakStreak<-3
  inPeak<-FALSE
  streaksVec<-c()
  thisStreak<-c()
  maxStreakDrop<-0
  maxStreakDropInd<-NA
  
  for(i in 1:nrow(d1)){
    if(is.na(d1$changeRate_CperH[i])){next}
    # print(d1$changeRate_CperH[i])
    if(d1$changeRate_CperH[i]<=peakThresh){
      g2Peaks<-g2Peaks+1
    }else{
      g2Peaks<-0
    }
    if(g2Peaks>maxStreakDrop){
      maxStreakDrop<-g2Peaks
      maxStreakDropInd<-i
    }
    
    if(g2Peaks>=peakStreak){
      # print('peaked')
      if(inPeak==FALSE){
        thisStreak<-c(i-1)
        # print('initlist')
      }
      inPeak<-TRUE
      thisStreak<-append(thisStreak,i)
    }
    if(g2Peaks==0){
      # print('outofpeak')
      if(inPeak==TRUE){
        # print('save this peak')
        streaksVec<-append(streaksVec,thisStreak)
      }
      inPeak<-FALSE
    }
  }
  
  d1_maxStreakDrop<-d1[(maxStreakDropInd-maxStreakDrop-10) : (maxStreakDropInd+10),]  
  # **********************************************************
  # **********************************************************
  
  d1_maxStreak$Streak<-"Increase"
  d1_maxStreak$Streak[d1_maxStreak$changeRate_CperH>=1.8]<-"targetUP"
  d1_maxStreakDrop$Streak<-"Decrease"
  d1_maxStreakDrop$Streak[d1_maxStreakDrop$changeRate_CperH<=-1.8]<-"targetDOWN"
  d1_max_streaks<-rbind(d1_maxStreak, d1_maxStreakDrop)
  # **********************************************************
  # **********************************************************
  
  # plot<-ggplot(d1, aes(changeRate_CperH))+
  #   geom_histogram()+
  #   scale_x_continuous(limits = c(-3, 3), labels = scales::comma, name = expression(Temperature~change~(degree*C~hour^-1)))+
  #   # xlim(0.05, 2.5)+
  #   theme_pubr()+
  #   # xlab(expression(Temperature~change~(degree*C~hour^-1)))+
  #   ylab(expression(N~datapoints~recorded))+
  #   theme(axis.title.x = element_text(size = 9), 
  #         axis.title.y = element_text(size = 9),
  #         axis.text.y = element_blank(),
  #         axis.text.x = element_text(size = 7, face = "bold"), 
  #         plot.title = element_text(size = 6))+
  #   ggtitle(paste("YEAR:", year, 
  #                 "MONTH:", month, 
  #                 "SITE: ", site,
  #                 "DATASET: ", dataset_name))
  # 
  plot2<-ggplot(d1_max_streaks, aes(x = as.POSIXct(TIME), y = TEMP_C, color = Streak))+
    geom_point(pch = 19, size=3)+
    geom_line( size=0.5, color= "grey")+
    scale_x_datetime(name = "Time h")+
    # xlim(0.05, 2.5)+
    theme_pubr()+
    scale_colour_viridis_d()+
    facet_wrap(.~factor(DAY), nrow = 6, ncol = 6, scales = "free")+
    # xlab(expression(Temperature~change~(degree*C~hour^-1)))+
    ylab("Temp ºC")+
    theme(axis.title.x = element_text(size = 9), 
          axis.title.y = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 7, face = "bold"), 
          plot.title = element_text(size = 6))+
    ggtitle(paste("YEAR:", year, 
                  "MONTH:", month, 
                  "SITE: ", site,
                  "DATASET: ", dataset_name))
  
  # plot3<-ggplot(d1, aes(x = as.POSIXct(TIME), y = TEMP_C))+
  #   geom_hline(yintercept = 12, color = "blue")+
  #   geom_hline(yintercept = 16, color = "green")+
  #   geom_hline(yintercept = 20, color = "red")+
  #   geom_line( size=0.5, color= "black")+
  #   scale_x_datetime()+
  #   # xlim(0.05, 2.5)+
  #   scale_x_time(name = "Time h", )+
  #   theme_pubr()+
  #   scale_colour_viridis_d()+
  #   facet_wrap(.~factor(DAY), nrow = 6, ncol = 6, scales = "free")+
  #   # xlab(expression(Temperature~change~(degree*C~hour^-1)))+
  #   ylab("Temp ºC")+
  #   theme(axis.title.x = element_text(size = 5), 
  #         axis.title.y = element_text(size = 5),
  #         axis.text.y = element_text(size = 5),
  #         axis.text.x = element_text(size = 5), 
  #         plot.title = element_text(size = 5))+
  #   ggtitle(paste("YEAR:", year, 
  #                 "MONTH:", month, 
  #                 "SITE: ", site,
  #                 "DATASET: ", dataset_name))
  # 
  # ggsave(plot = plot, filename = paste("./Figures/RateChange/plotRATES",dataset_name, year, month, site, ".png", sep="_"), width = 3, height = 3, units = "in")
  ggsave(plot = plot2, filename = paste("./Figures/RateChange/plotSTREAKS",dataset_name, year, month, site, ".png", sep="_"), width = 7, height = 3, units = "in")
  # ggsave(plot = plot3, filename = paste("./Figures/RateChange/plotTEMPS",dataset_name, year, month, site, ".png", sep="_"), width = 10, height = 10, units = "in")
  write.csv( d1, file = paste("./ExportedData/RateChange/data",dataset_name, year, month, site, ".csv", sep="_"), row.names = FALSE)
  write.csv( d1_max_streaks, file = paste("./ExportedData/RateChange/dataSTREAKS_",dataset_name, year, month, site,  ".csv", sep="_"), row.names = FALSE)
  
  
  # return(assign(x = filename, d1_temp))
  # }
  # }  
  
  # }
  # return(d1_temp) 
}

## Apply the function for specific months (numeric) and years (numeric)
rate_T_change(data = data.n, site = "NAPL", year = 2016,month = 6, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2016,month = 8, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2016,month = 5, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2016,month = 7, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2016,month = 9, dataset_name = "SBCLTER_NAPL")

rate_T_change(data = data.n, site = "NAPL", year = 2017,month = 6, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2017,month = 8, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2017,month = 5, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2017,month = 7, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2017,month = 9, dataset_name = "SBCLTER_NAPL")

rate_T_change(data = data.n, site = "NAPL", year = 2018,month = 6, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2018,month = 8, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2018,month = 5, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2018,month = 7, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n, site = "NAPL", year = 2018, month = 9, dataset_name = "SBCLTER_NAPL")

rate_T_change(data = data.n,  site = "NAPL", year = 2019, month = 5, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2019, month = 6, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2019, month = 7, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2019, month = 8, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2019, month = 9, dataset_name = "SBCLTER_NAPL")

rate_T_change(data = data.n,  site = "NAPL", year = 2020, month = 5, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2020, month = 6, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2020, month = 7, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2020, month = 8, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2020, month = 9, dataset_name = "SBCLTER_NAPL")

rate_T_change(data = data.n,  site = "NAPL", year = 2021, month = 5, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2021, month = 6, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2021, month = 7, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2021, month = 8, dataset_name = "SBCLTER_NAPL")
rate_T_change(data = data.n,  site = "NAPL", year = 2021, month = 9, dataset_name = "SBCLTER_NAPL")

# select data to be used:

# Temperature acute drop data
data.list<-c("./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2020_6_NAPL_.csv",
             "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2021_5_NAPL_.csv", 
             "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2021_6_NAPL_.csv", 
             "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2019_5_NAPL_.csv", 
             "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2020_5_NAPL_.csv", 
             "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2018_6_NAPL_.csv",
             "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2017_7_NAPL_.csv",
             "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2017_6_NAPL_.csv",
             "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2016_6_NAPL_.csv")

# Temperature acute rise data
data.listIN<-c("./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2016_6_NAPL_.csv", 
               "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2019_9_NAPL_.csv",
               "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2017_8_NAPL_.csv",
               "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2018_8_NAPL_.csv",
               "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2018_9_NAPL_.csv",
               "./ExportedData/RateChange/dataSTREAKS__SBCLTER_NAPL_2016_7_NAPL_.csv")


for (i in 1:length(data.list)){
  dd1<-read.csv(data.list[i])
  dd1<-dd1[c(dd1$Streak == "Decrease" | dd1$Streak == "targetDOWN"), ]
  dd1$numTime<-c(1:nrow(dd1))
  dd1$group<-i
  
  if(i == 1){
    dataDECREASE<-dd1
  }else{
    dataDECREASE<-rbind(dataDECREASE, dd1)
  }
}


for (i in 1:length(data.listIN)){
  dd1<-read.csv(data.listIN[i])
  dd1<-dd1[c(dd1$Streak == "Increase" | dd1$Streak == "targetUP"), ]
  dd1$numTime<-c(1:nrow(dd1))
  dd1$group<-i
  
  if(i == 1){
    dataINCREASE<-dd1
  }else{
    dataINCREASE<-rbind(dataINCREASE, dd1)
  }
}

# 4. Figures: Thesis and JEB -------

# warm: #FFB4A9, #E9A099, #DE4949, #900013, "grey"
# cold: #7EB0E1, #578BBA, #2E6894, #00284E
cols.4<-c("#3596B5","#49817B","#A2416E", "#CD6C95")
vlines <- data.frame(xint = c(12, 16, 20, 22), grp = c("#3596B5", "#49817B", "#A2416E", "#CD6C95"))


# The closest site to the shore
plot_site2<-ggplot(data = data.n[c(data.n$MONTH == 5 | data.n$MONTH == 6 | data.n$MONTH == 7 | data.n$MONTH == 8 | data.n$MONTH == 9), ],
                   aes(x = TEMP_C, group = MONTH, y = factor(MONTH)))+
  geom_density_ridges(fill = "grey", alpha = 1)+
  geom_vline(data = vlines, aes(xintercept = xint, colour = grp), linetype = "solid", alpha = 1, size =3, show.legend = FALSE) +
  scale_color_manual(values =cols.4)+
  theme_ridges(grid = FALSE)+
  scale_x_continuous(breaks = c(12,16, 20, 22), labels = c("12", "16", "20", "22 ºC"))+
  scale_y_discrete(breaks =c("5", "6", "7", "8", "9"), labels = c("May", "Jun", "Jul", "Aug", "Sept"))+
  theme(axis.title.x =element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15))




up1<-ggplot(data = dataINCREASE, aes(y=TEMP_C, x = numTime*20/60, group = group ))+
  annotate("rect", xmin = 3, xmax = 5, ymin = -Inf, ymax = Inf,
           alpha = .2)+
  geom_line()+
  geom_hline(yintercept = 16, color = cols.4[2], alpha = 1, size = 2)+
  geom_hline(yintercept = 20, color = cols.4[3], alpha = 1, size = 2)+
  ylab(expression(Temperature~degree*C))+
  xlab(Relative~Time~(h))+
  theme_pubr()+
  ylim(10, 22)+
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.y =  element_text(size = 15),
        axis.text.x = element_text(size = 15))

down1<-ggplot(data = dataDECREASE, aes(y=TEMP_C, x = numTime*20/60, group = group ))+
  annotate("rect", xmin = 3, xmax = 5, ymin = -Inf, ymax = Inf,
           alpha = .2)+
  geom_line()+
  geom_hline(yintercept = 16, color = cols.4[2], alpha = 1, size = 2)+
  geom_hline(yintercept = 12, color = cols.4[1], alpha = 1, size = 2)+
  ylab(expression(Temperature~degree*C))+
  xlab(Relative~Time~(h))+
  ylim(10, 22)+
  theme_pubr()+
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.y =  element_text(size = 15),
        axis.text.x = element_text(size = 15))

ggsave(plot =  plot_site2, filename = paste( "./Figures/Figure2B_ridges",Sys.Date(),".png",sep=""), width = 4, height = 3, units = "in")
ggsave(plot =  up1, filename = paste( "./Figures/Figure2C_increase",Sys.Date(),".png",sep=""), width = 3, height = 3, units = "in")
ggsave(plot =  down1, filename = paste( "./Figures/Figure2C_decrease",Sys.Date(),".png",sep=""), width = 3, height = 3, units = "in")


