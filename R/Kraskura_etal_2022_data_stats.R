
library(emmeans)
library(lme4)
library(ggsci)
library(cowplot)
library(ggpubr)
library(grid)
library(gridExtra)
library(gtable)

# Custom functions used in the script: -----

# function to order model selection based on the lowest BIC score
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return( BIC.t)
}

# Function to format ggplot. 
# Args:
#   plot: saved ggplot, variable name, must provide 
#   title: the title of the plot; default= "" (no title)
#	  y_title: the title of y axis, must provide 
#	  x_title: the title of x axis, must provide
#	  print: logical, TRUE = print (default), FALSE = do not print the plot
#
# Returns: formatted ggplot. Formatting themes include: axis, text size, label sizes, no gridlines, panel border around the plot


ggformat<-function(plot, title="", y_title, x_title, print=TRUE, size_text = 15){
  
  plot_name<-deparse(substitute(plot))
  
  plot<- plot +
    theme_classic()+
    ggtitle(title)+							
    ylab(y_title)+ 							
    xlab(x_title)+ 							
    theme(axis.text.y=element_text(size=size_text, colour= 'black'),
          axis.text.x=element_text(size=size_text, colour= 'black'),
          axis.line.y=element_line(colour = 'black',size=0.5),
          axis.line.x=element_line(colour = 'black',size=0.5),
          axis.ticks.y=element_line(size=0.5),
          # axis.ticks.x=element_line(size=0.5), 
          axis.ticks.x.bottom = element_line(size=0.5, colour = "black"))+
    theme(axis.title.y=element_text(size=size_text),
          axis.title.x=element_text(size=size_text),
          panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))
  
  if (print==TRUE){	
    print(plot)	
  }
  
  assign(plot_name, plot, envir=parent.frame())
  
}

# color pallete:
cols.4<-c("#3596B5", "#49817B", "#CD6C95", "#A2416E")
# dark:
col.d<-c("#002640")



# 1. Setup and read in data -----

data<-read.csv("./Data/Kraskura_etal_perch_temp_size_RESP.csv")
data.abt<-read.csv("./Data/Kraskura_etal_perch_temp_size_ABT.csv")
data.abtID<-read.csv("./Data/Kraskura_etal_perch_temp_size_CardTempTol.csv")

dir.create("./Figures", recursive = TRUE)
# convert needed columns to factor
factor_cols.data = c("FishID", "ExperimentID", "treatm", "Tank", "RespoID", "sex", "Ch", "treatm", "origin", "sizeClass", "pregnant")
data = data %>% 
  mutate(across(all_of(factor_cols.data), factor))

# convert needed columns to factor
data.abt = data.abt %>% 
  mutate(across(all_of(c("FishID", "treatm", "sizeClass")), factor))
data.abtID = data.abtID %>% 
  mutate(across(all_of(c("FishID", "sizeClass", "origin")), factor))

# 1.2. temperature specific datasets (all fish)
data12<-data[data$treatm=="12",]
data20<-data[data$treatm=="20",]
data16<-data[data$treatm=="16",]
data22<-data[data$treatm=="22",]
data24<-data[data$treatm=="24",]

# 1.3. temp specific ABT data (4 fish have two measurements at 16 C, start of the test delayed, others all good )
abt.16.0<-data.abt[(data.abt$TEMP == 16),] # n = 31
data.abt16<-abt.16.0[!duplicated(abt.16.0$FishID),] # n = 27; taling the first 16 C HR measurement
data.abt20<-data.abt[(data.abt$TEMP == 20),]
data.abt22<-data.abt[(data.abt$TEMP == 22),]
data.abt24<-data.abt[(data.abt$TEMP == 24),]

# 1.4. lmer datasets; exclude NAs 
data.fas.mod<-data[!is.na(data$FAS),]
data.as.mod<-data[!is.na(data$AS),]
data.mmr.mod<-data[!is.na(data$mmr),]
data.rmr.mod<-data[!is.na(data$rmr),]

# 2. Linear models  --------
## 2.1. lme4::lmer cont temp metabolism -------
regrmrO<-lmer(log(rmr) ~ log(BW) + t_mean + origin + (1|FishID), data.rmr.mod, REML = F)
regrmrSO<-lmer(log(rmr) ~ log(BW) + t_mean + sex + origin + (1|FishID), data.rmr.mod, REML = F)
regrmrS<-lmer(log(rmr) ~ log(BW) + t_mean + sex + (1|FishID), data.rmr.mod, REML = F)
regrmr<-lmer(log(rmr) ~ log(BW) + t_mean + (1|FishID), data.rmr.mod, REML = F)
# regrmr2<-lmer(log(rmr) ~ log(BW) * t_mean + sex + (1|FishID), data.rmr.mod, REML = F)

regmmrO<-lmer(log(mmr) ~ log(BW) + t_mean + origin + (1|FishID), data.mmr.mod, REML = F)
regmmrSO<-lmer(log(mmr) ~ log(BW) + t_mean + sex + origin + (1|FishID), data.mmr.mod, REML = F)
regmmrS<-lmer(log(mmr) ~ log(BW) + t_mean + sex + (1|FishID), data.mmr.mod, REML = F)
regmmr<-lmer(log(mmr) ~ log(BW) + t_mean + (1|FishID), data.mmr.mod, REML = F)
# regmmr2<-lmer(log(mmr) ~ log(BW) * t_mean + sex + (1|FishID), data.mmr.mod, REML = F)

regasO<-lmer(log(AS) ~ log(BW) + t_mean + origin + (1|FishID), data.as.mod, REML = F)
regasSO<-lmer(log(AS) ~ log(BW) + t_mean + sex + origin + (1|FishID), data.as.mod, REML = F)
regasS<-lmer(log(AS) ~ log(BW) + t_mean + sex + (1|FishID), data.as.mod, REML = F)
regas<-lmer(log(AS) ~ log(BW) + t_mean + (1|FishID), data.as.mod, REML = F)
# regas2<-lmer(log(AS) ~ log(BW) * t_mean + (1|FishID), data.as.mod, REML = F)

BICdelta(BIC(regmmrS, regmmr, regmmrSO, regmmrO)) # missing sex on few fish, see dataset
BICdelta(BIC(regrmrS, regrmr, regrmrSO, regrmrO)) # missing sex on few fish, see dataset
BICdelta(BIC(regas, regasO, regasSO, regasS)) # missing sex on few fish, see dataset

# residuals 
plot(resid(regrmrO))
plot(resid(regmmr))
plot(resid(regas))

summary(regrmrO)
summary(regmmr)
summary(regas)


## 2.2. lme4::lmer discrete temp metabolism -------
disc.regrmrO<-lmer(log(rmr) ~ log(BW) + treatm + origin + (1|FishID), data.rmr.mod, REML = F)
disc.regrmrSO<-lmer(log(rmr) ~ log(BW) + treatm + sex + origin + (1|FishID), data.rmr.mod, REML = F)
disc.regrmrS<-lmer(log(rmr) ~ log(BW) + treatm + sex + (1|FishID), data.rmr.mod, REML = F)
disc.regrmr<-lmer(log(rmr) ~ log(BW) + treatm + (1|FishID), data.rmr.mod, REML = F)
disc.regrmrX<-lmer(log(rmr) ~ log(BW) * treatm + (1|FishID), data.rmr.mod, REML = F)
# disc.regrmr2<-lmer(log(rmr) ~ log(BW) * treatm + sex + (1|FishID), data.rmr.mod, REML = F)

disc.regmmrO<-lmer(log(mmr) ~ log(BW) + treatm + origin + (1|FishID), data.mmr.mod, REML = F)
disc.regmmrSO<-lmer(log(mmr) ~ log(BW) + treatm + sex + origin + (1|FishID), data.mmr.mod, REML = F)
disc.regmmrS<-lmer(log(mmr) ~ log(BW) + treatm + sex + (1|FishID), data.mmr.mod, REML = F)
disc.regmmr<-lmer(log(mmr) ~ log(BW) + treatm + (1|FishID), data.mmr.mod, REML = F)
disc.regmmrX<-lmer(log(mmr) ~ log(BW) * treatm + (1|FishID), data.mmr.mod, REML = F)
# disc.regmmr2<-lmer(log(mmr) ~ log(BW) * treatm + sex + (1|FishID), data.mmr.mod, REML = F)

disc.regasO<-lmer(log(AS) ~ log(BW) + treatm + origin + (1|FishID), data.as.mod, REML = F)
disc.regasSO<-lmer(log(AS) ~ log(BW) + treatm + sex + origin + (1|FishID), data.as.mod, REML = F)
disc.regasS<-lmer(log(AS) ~ log(BW) + treatm + sex + (1|FishID), data.as.mod, REML = F)
disc.regas<-lmer(log(AS) ~ log(BW) + treatm + (1|FishID), data.as.mod, REML = F)
# disc.regas2<-lmer(log(AS) ~ log(BW) * treatm + (1|FishID), data.as.mod, REML = F)

# scaling with discrete temp
disc.regfasO<-lmer(log(FAS) ~ log(BW) + treatm + origin + (1|FishID), data.fas.mod, REML = F)
disc.regfasSO<-lmer(log(FAS) ~ log(BW) + treatm + sex + origin + (1|FishID), data.fas.mod, REML = F)
disc.regfasS<-lmer(log(FAS) ~ log(BW) + treatm + sex + (1|FishID), data.fas.mod, REML = F)
disc.regfas<-lmer(log(FAS) ~ log(BW) + treatm + (1|FishID), data.fas.mod, REML = F)
disc.regfas2<-lmer(log(FAS) ~ log(BW) * treatm + (1|FishID), data.fas.mod, REML = F)

BICdelta(BIC(disc.regmmrS, disc.regmmr, disc.regmmrSO, disc.regmmrO, disc.regmmrX))
BICdelta(BIC(disc.regrmrS, disc.regrmr, disc.regrmrSO, disc.regrmrO, disc.regrmrX))
BICdelta(BIC(disc.regas, disc.regasO, disc.regasSO, disc.regasS))
BICdelta(BIC(disc.regfas, disc.regfasO, disc.regfasSO, disc.regfasS, disc.regfas2))

plot(resid(disc.regmmr))
plot(resid(disc.regrmrO))
plot(resid(disc.regas))
plot(resid(disc.regfasO))

summary(disc.regmmr)
summary(disc.regrmrO)
summary(disc.regas)
summary(disc.regfasO)



## 2.3. glm temp-specific metabolism----------
reg.12.rmr<-glm(log((rmr)) ~ log(BW), data12, family = "gaussian")
reg.16.rmr<-glm(log((rmr)) ~ log(BW), data16, family = "gaussian")
reg.20.rmr<-glm(log((rmr)) ~ log(BW), data20, family = "gaussian")
reg.22.rmr<-glm(log((rmr)) ~ log(BW), data22, family = "gaussian")

reg.12.mmr<-glm(log((mmr)) ~ log(BW), data12, family = "gaussian")
reg.16.mmr<-glm(log((mmr)) ~ log(BW), data16, family = "gaussian")
reg.20.mmr<-glm(log((mmr)) ~ log(BW), data20, family = "gaussian")
reg.22.mmr<-glm(log((mmr)) ~ log(BW), data22, family = "gaussian")

reg.12.as<-glm(log(AS) ~ log(BW), data12, family = "gaussian")
reg.16.as<-glm(log(AS) ~ log(BW), data16, family = "gaussian")
reg.20.as<-glm(log(AS) ~ log(BW), data20, family = "gaussian")
reg.22.as<-glm(log(AS) ~ log(BW), data22, family = "gaussian")

reg.12.fas<-glm(log(FAS) ~ log(BW), data12, family = "gaussian")
reg.16.fas<-glm(log(FAS) ~ log(BW), data16, family = "gaussian")
reg.20.fas<-glm(log(FAS) ~ log(BW), data20, family = "gaussian")
reg.22.fas<-glm(log(FAS) ~ log(BW), data22, family = "gaussian")

# no biological difference
# reg.12.mmrChase<-glm(log((mmrChase)) ~ log(BW), data12, family = "gaussian")
# reg.16.mmrChase<-glm(log((mmrChase)) ~ log(BW), data16, family = "gaussian")
# reg.20.mmrChase<-glm(log((mmrChase)) ~ log(BW), data20, family = "gaussian")
# reg.22.mmrChase<-glm(log((mmrChase)) ~ log(BW), data22, family = "gaussian")




## 2.4. lme4::lmer beat per minute. scaling -----------

# size discrete. param
modbpm<-lmer(log(bpm) ~  log(BW) + treatm + (1|FishID), data.abt, REML = F)
modbpm1<-lmer(log(bpm) ~  log(BW) * treatm + (1|FishID), data.abt, REML = F)
modbpm0<-lmer(log(bpm) ~  log(BW) + (1|FishID), data.abt, REML = F)

BICdelta(BIC(modbpm, modbpm0, modbpm1))

plot(resid(modbpm))
summary(modbpm)

## 2.5. lm beat per minute. temp specific scaling -----------
m16<-lm(log(bpm)~log(BW), data.abt16)
m20<-lm(log(bpm)~log(BW), data.abt20)
m22<-lm(log(bpm)~log(BW), data.abt22)
m24<-lm(log(bpm)~log(BW), data.abt24)


## 2.6. lm cardiac thermal tolerance: scaling ------

mod.arr0<-lm(log(temp_ARRH) ~ 1, data = data.abtID)
mod.arr1<-lm(log(temp_ARRH) ~ log(BW), data = data.abtID)
mod.arrO<-lm(log(temp_ARRH) ~ log(BW) + origin, data = data.abtID)
mod.arrS<-lm(log(temp_ARRH) ~ log(BW) + sex, data = data.abtID)
mod.arrSO<-lm(log(temp_ARRH) ~ log(BW) + sex + origin, data = data.abtID)

mod.peak0<-lm(log(Tpeak) ~ 1, data = data.abtID)
mod.peak1<-lm(log(Tpeak) ~ log(BW), data = data.abtID)
mod.peakO<-lm(log(Tpeak) ~ log(BW) + origin, data = data.abtID)
mod.peakS<-lm(log(Tpeak) ~ log(BW) + sex, data = data.abtID)
mod.peakSO<-lm(log(Tpeak) ~ log(BW) + sex + origin, data = data.abtID)

mod.bp0<-lm(log(breakpoint_Cels) ~ 1, data = data.abtID)
mod.bp1<-lm(log(breakpoint_Cels) ~ log(BW), data = data.abtID)
mod.bpO<-lm(log(breakpoint_Cels) ~ log(BW) + origin, data = data.abtID)
mod.bpS<-lm(log(breakpoint_Cels) ~ log(BW) + sex, data = data.abtID)
mod.bpSO<-lm(log(breakpoint_Cels) ~ log(BW) + sex + origin, data = data.abtID)

mod.HRp0<-lm(log(HRpeak) ~ 1, data = data.abtID)
mod.HRp1<-lm(log(HRpeak) ~ log(BW), data = data.abtID)
mod.HRpO<-lm(log(HRpeak) ~ log(BW) + origin, data = data.abtID)
mod.HRpS<-lm(log(HRpeak) ~ log(BW) + sex, data = data.abtID)
mod.HRpSO<-lm(log(HRpeak) ~ log(BW) + sex + origin, data = data.abtID)


BICdelta(BIC(mod.arr0, mod.arr1, mod.arrO, mod.arrS, mod.arrSO)) 
BICdelta(BIC(mod.peak0, mod.peak1, mod.peakO, mod.peakS, mod.peakSO)) 
BICdelta(BIC(mod.HRp0, mod.HRp1, mod.HRpO, mod.HRpS, mod.HRpSO))
BICdelta(BIC(mod.bp0, mod.bp1, mod.bpO, mod.bpS, mod.bpSO))

# check resid
plot(mod.arr1)
plot(mod.peak1)
plot(mod.bpO)
plot(mod.HRp0)

summary(mod.arr1) 
summary(mod.peak1)
summary(mod.HRp0)
summary(mod.bpO)
# summary(mod.bp1)



## 2.7. metabolism correlations -------
modmrcor1<-lm(mmr ~ rmr + treatm, data = data)
modmrcor0<-lm(mmr ~ rmr , data = data)

BICdelta(BIC(modmrcor1, modmrcor0))
summary(modmrcor1)

## 2.8. heart size scaling  ---------
heart.scale <- lm(log(heart.mass/1000)~log(BW), data = data.abtID)
summary(heart.scale)




# 3.Mass specific measurements metabolism, STATS and estimated model means ----------
## 3.1. LMER reference grid metabolism data with discrete temperature -----------
pred.means.mmr1kg<-as.data.frame(summary(ref_grid(disc.regmmr, at = list(BW= 1),calc = c(n = ".wgt."))))
pred.means.rmr1kg<-as.data.frame(summary(ref_grid(disc.regrmrO, at = list(BW= 1),calc = c(n = ".wgt."))))
pred.means.as1kg<-as.data.frame(summary(ref_grid(disc.regas, at = list(BW= 1),calc = c(n = ".wgt."))))
pred.means.fas1kg<-as.data.frame(summary(ref_grid(disc.regfasO, at = list(BW= 1),calc = c(n = ".wgt."))))
pred.means.bpm1kg<-as.data.frame(summary(ref_grid(modbpm, at = list(BW = 1),calc = c(n = ".wgt."))))

# check estimated marginal means: BW = 1 is for body mass 1 kg 
# emmeans(disc.regmmr, c("treatm"), at = list(BW = 1), calc = c(n = ".wgt."),) # n =8
# emmeans(disc.regrmrO, c("treatm","origin"), at = list(BW = 1), calc = c(n = ".wgt.")) # n =8
# emmeans(disc.regas, c("treatm"), at = list(BW = 1), calc = c(n = ".wgt.")) # n =4
# emmeans(disc.regfasO, c("treatm","origin"), at = list(BW = 1), calc = c(n = ".wgt.")) # n=8
# emmeans(modbpm, c("treatm"), at = list(BW = 1), calc = c(n = ".wgt."))  # 10 g fish 


## 3.2. new data grid with each individual, predict values -------
data.rmr.pred <- data.rmr.mod[, c("FishID", "origin", "treatm", "BW", "rmr", "mmr", "AS", "FAS")]
data.rmr.pred$BW <- 1
data.mmr.pred <- data.mmr.mod[, c("FishID", "origin", "treatm", "BW", "rmr", "mmr", "AS", "FAS")]
data.mmr.pred$BW <- 1
data.as.pred <- data.as.mod[, c("FishID", "origin", "treatm", "BW", "rmr", "mmr", "AS", "FAS")]
data.as.pred$BW <- 1
data.fas.pred <- data.fas.mod[, c("FishID", "origin", "treatm", "BW", "rmr", "mmr", "AS", "FAS")]
data.fas.pred$BW <- 1
data.bpm.pred <- data.abt[, c("FishID",  "treatm", "BW", "bpm")]
data.bpm.pred$BW <- 1
# pred 2 used to get only temp specific data
# data.bpm.pred2 <-as.data.frame(expand.grid(treatm = c("16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28"),
#                                            BW = c(0.01, 0.35)))
data.bpm.pred2 <-as.data.frame(expand.grid(treatm = c("16", "20", "22",  "24"),
                                           BW = c(0.01, 0.35)))
data.bpm.pred2$pre <-predict(modbpm, re.form = NA, newdata = data.bpm.pred2)
data.bpm.pred2$TEMP<-as.numeric(as.character(data.bpm.pred2$treatm))

# using a loop below adding individual specific values. "full data"
data.rmr.pred$rmr_pred <- NA
data.mmr.pred$mmr_pred <- NA
data.as.pred$as_pred <- NA
data.rmr.pred$fas_pred <- NA
data.bpm.pred$bpm_pred <- NA

data.rmr.pred$mean_rmr<-NA
data.mmr.pred$mean_mmr<-NA
data.as.pred$mean_as<-NA
data.fas.pred$mean_fas<-NA
data.bpm.pred$mean_bpm <- NA

data.rmr.pred$mean_rmrSE<-NA
data.mmr.pred$mean_mmrSE<-NA
data.as.pred$mean_asSE<-NA
data.fas.pred$mean_fasSE<-NA
data.bpm.pred$mean_bpmSE <- NA

data.rmr.pred$mean_rmrn<-NA
data.mmr.pred$mean_mmrn<-NA
data.as.pred$mean_asn<-NA
data.fas.pred$mean_fasn<-NA
data.bpm.pred$mean_bpmn <- NA

data.rmr.pred$q10_rmr<-NA
data.mmr.pred$q10_mmr<-NA
data.as.pred$q10_as<-NA
data.fas.pred$q10_fas<-NA
data.bpm.pred$q10_bpm <- NA

## 3.3. Q10 calculations: use estimated means ---------
# https://www.physiologyweb.com/calculators/q10_calculator.html

pred.means.bpm1kg$q10<-NA
for(i in 2:nrow(pred.means.bpm1kg)){
  
  R2<-exp(pred.means.bpm1kg$prediction[i])
  R1<-exp(pred.means.bpm1kg$prediction[i-1])
  T2<-as.numeric(as.character(pred.means.bpm1kg$treatm[i]))
  T1<-as.numeric(as.character(pred.means.bpm1kg$treatm[i-1]))
  
  pred.means.bpm1kg$q10[i]<-(R2/R1)^(10/(T2-T1))
  
}

pred.means.rmr1kg$q10<-NA
for(i in 2:nrow(pred.means.rmr1kg)){
  
  R2<-exp(pred.means.rmr1kg$prediction[i])
  R1<-exp(pred.means.rmr1kg$prediction[i-1])
  T2<-as.numeric(as.character(pred.means.rmr1kg$treatm[i]))
  T1<-as.numeric(as.character(pred.means.rmr1kg$treatm[i-1]))
  pred.means.rmr1kg$q10[i]<-(R2/R1)^(10/(T2-T1))
}

pred.means.mmr1kg$q10<-NA
for(i in 2:nrow(pred.means.mmr1kg)){
  
  R2<-exp(pred.means.mmr1kg$prediction[i])
  R1<-exp(pred.means.mmr1kg$prediction[i-1])
  T2<-as.numeric(as.character(pred.means.mmr1kg$treatm[i]))
  T1<-as.numeric(as.character(pred.means.mmr1kg$treatm[i-1]))
  pred.means.mmr1kg$q10[i]<-(R2/R1)^(10/(T2-T1))
}

pred.means.as1kg$q10<-NA
for(i in 2:nrow(pred.means.as1kg)){
  
  R2<-exp(pred.means.as1kg$prediction[i])
  R1<-exp(pred.means.as1kg$prediction[i-1])
  T2<-as.numeric(as.character(pred.means.as1kg$treatm[i]))
  T1<-as.numeric(as.character(pred.means.as1kg$treatm[i-1]))
  pred.means.as1kg$q10[i]<-(R2/R1)^(10/(T2-T1))
}

pred.means.fas1kg$q10<-NA
for(i in 2:nrow(pred.means.fas1kg)){
  
  R2<-exp(pred.means.fas1kg$prediction[i])
  R1<-exp(pred.means.fas1kg$prediction[i-1])
  T2<-as.numeric(as.character(pred.means.fas1kg$treatm[i]))
  T1<-as.numeric(as.character(pred.means.fas1kg$treatm[i-1]))
  pred.means.fas1kg$q10[i]<-(R2/R1)^(10/(T2-T1))
}


for (i in 1:nrow(data.rmr.pred)){
  data.rmr.pred$rmr_pred[i] <- pred.means.rmr1kg[which(pred.means.rmr1kg$treatm == data.rmr.pred$treatm[i] & pred.means.rmr1kg$origin == data.rmr.pred$origin[i]), "prediction"] + residuals(disc.regrmrO)[i] 
  data.rmr.pred$mean_rmr[i]<-pred.means.rmr1kg[which(pred.means.rmr1kg$treatm == data.rmr.pred$treatm[i] & pred.means.rmr1kg$origin == data.rmr.pred$origin[i]), "prediction"]
  data.rmr.pred$mean_rmrSE[i]<-pred.means.rmr1kg[which(pred.means.rmr1kg$treatm == data.rmr.pred$treatm[i] & pred.means.rmr1kg$origin == data.rmr.pred$origin[i]), "SE"]
  data.rmr.pred$mean_rmrn[i]<-pred.means.rmr1kg[which(pred.means.rmr1kg$treatm == data.rmr.pred$treatm[i] & pred.means.rmr1kg$origin == data.rmr.pred$origin[i]), "n"]
  data.rmr.pred$q10_rmr[i]<-pred.means.rmr1kg[which(pred.means.rmr1kg$treatm == data.rmr.pred$treatm[i] & pred.means.rmr1kg$origin == data.rmr.pred$origin[i]), "q10"]
}
for (i in 1:nrow(data.mmr.pred)){
  data.mmr.pred$mmr_pred[i] <- pred.means.mmr1kg[which(pred.means.mmr1kg$treatm == data.mmr.pred$treatm[i]), "prediction"] + residuals(disc.regmmrO)[i] 
  data.mmr.pred$mean_mmr[i]<-pred.means.mmr1kg[which(pred.means.mmr1kg$treatm == data.mmr.pred$treatm[i]), "prediction"]
  data.mmr.pred$mean_mmrSE[i]<-pred.means.mmr1kg[which(pred.means.mmr1kg$treatm == data.mmr.pred$treatm[i]), "SE"]
  data.mmr.pred$mean_mmrn[i]<-pred.means.mmr1kg[which(pred.means.mmr1kg$treatm == data.mmr.pred$treatm[i]), "n"]
  data.mmr.pred$q10_mmr[i]<-pred.means.mmr1kg[which(pred.means.mmr1kg$treatm == data.mmr.pred$treatm[i]), "q10"]
}
for (i in 1:nrow(data.as.pred)){
  data.as.pred$as_pred[i] <- pred.means.as1kg[which(pred.means.as1kg$treatm == data.as.pred$treatm[i]), "prediction"] + residuals(disc.regasO)[i] 
  data.as.pred$mean_as[i]<-pred.means.as1kg[which(pred.means.as1kg$treatm == data.as.pred$treatm[i]), "prediction"]
  data.as.pred$mean_asSE[i]<-pred.means.as1kg[which(pred.means.as1kg$treatm == data.as.pred$treatm[i]), "SE"]
  data.as.pred$mean_asn[i]<-pred.means.as1kg[which(pred.means.as1kg$treatm == data.as.pred$treatm[i]), "n"]
  data.as.pred$q10_as[i]<-pred.means.as1kg[which(pred.means.as1kg$treatm == data.as.pred$treatm[i]), "q10"]
}
for (i in 1:nrow(data.fas.pred)){
  data.fas.pred$fas_pred[i] <- pred.means.fas1kg[which(pred.means.fas1kg$treatm == data.fas.pred$treatm[i] & pred.means.fas1kg$origin == data.fas.pred$origin[i]), "prediction"] + residuals(disc.regfasO)[i] 
  data.fas.pred$mean_fas[i]<-pred.means.fas1kg[which(pred.means.fas1kg$treatm == data.fas.pred$treatm[i] & pred.means.fas1kg$origin == data.fas.pred$origin[i]), "prediction"]
  data.fas.pred$mean_fasSE[i]<-pred.means.fas1kg[which(pred.means.fas1kg$treatm == data.fas.pred$treatm[i] & pred.means.fas1kg$origin == data.fas.pred$origin[i]), "SE"]
  data.fas.pred$mean_fasn[i]<-pred.means.fas1kg[which(pred.means.fas1kg$treatm == data.fas.pred$treatm[i] & pred.means.fas1kg$origin == data.fas.pred$origin[i]), "n"]
  data.fas.pred$q10_fas[i]<-pred.means.fas1kg[which(pred.means.fas1kg$treatm == data.fas.pred$treatm[i]), "q10"]
}

for (i in 1:nrow(data.bpm.pred)){
  data.bpm.pred$bpm_pred[i] <- pred.means.bpm1kg[which(pred.means.bpm1kg$treatm == data.bpm.pred$treatm[i]), "prediction"] + residuals(modbpm)[i] 
  data.bpm.pred$mean_bpm[i]<-pred.means.bpm1kg[which(pred.means.bpm1kg$treatm == data.bpm.pred$treatm[i]), "prediction"]
  data.bpm.pred$mean_bpmSE[i]<-pred.means.bpm1kg[which(pred.means.bpm1kg$treatm == data.bpm.pred$treatm[i]), "SE"]
  data.bpm.pred$mean_bpmn[i]<-pred.means.bpm1kg[which(pred.means.bpm1kg$treatm == data.bpm.pred$treatm[i]), "n"]
  data.bpm.pred$mean_bpmn[i]<-pred.means.bpm1kg[which(pred.means.bpm1kg$treatm == data.bpm.pred$treatm[i]), "n"]
  data.bpm.pred$q10_bpm[i]<-pred.means.bpm1kg[which(pred.means.bpm1kg$treatm == data.bpm.pred$treatm[i]), "q10"]
}

data.bpm.pred.plot<-data.bpm.pred[!data.bpm.pred$treatm=="28", ]
pred.means.bpm1kg<-pred.means.bpm1kg[!pred.means.bpm1kg$treatm=="28", ]




# 4. ANOVA and Post-Hoc: -----
## 4.1. Type II Anovas --------
car:::Anova(modbpm, type = "II")
car:::Anova(disc.regfasO, type = "II")
car:::Anova(disc.regas, type = "II")
car:::Anova(disc.regmmr, type = "II")
car:::Anova(disc.regrmrO, type = "II")

# LM
car:::Anova(mod.arr1, type = "II")
car:::Anova(mod.peak1, type = "II")
# car:::Anova(mod.HRp0, type = "II") # size not sign. 
car:::Anova(mod.bpO, type = "II")


## 4.2. Post hocs ---------
pairs(emmeans::emmeans(modbpm, "treatm"))
pairs(emmeans::emmeans(disc.regfasO, c("treatm")))
pairs(emmeans::emmeans(disc.regas, c("treatm")))
pairs(emmeans::emmeans(disc.regmmr, c("treatm")))
pairs(emmeans::emmeans(disc.regrmrO, c("treatm")))





# 5. Obtain estimates: b (slope), a (intercepts), and CIs -----
## 5.1. metabolic rates and bpms: -------
# intercepts
scaling.int.mmr<-fixef(disc.regmmr)[1] # at 12 C
scaling.int.rmr<-fixef(disc.regrmrO)[1] # at 12 C and in field
scaling.int.as<-fixef(disc.regas)[1] # at 12 C
scaling.int.fas<-fixef(disc.regfasO)[1] # at 12 C and in field

# Scaling slopes:
mmr.b<-fixef(disc.regmmr)[2]
rmr.b<-fixef(disc.regrmrO)[2]
as.b<-fixef(disc.regas)[2]
fas.b<-fixef(disc.regfasO)[2]
bpm.b<-fixef(modbpm)[2]

# 2.5 % CIs meatabolism 
mmr.ciL<-as.data.frame(confint.merMod(disc.regmmr))[4,1]
rmr.ciL<-as.data.frame(confint.merMod(disc.regrmrO))[4,1]
as.ciL<-as.data.frame(confint.merMod(disc.regas))[4,1]
fas.ciL<-as.data.frame(confint.merMod(disc.regfasO))[4,1]

# 97.5 % CIs meatabolism 
mmr.ciH<-as.data.frame(confint.merMod(disc.regmmr))[4,2]
rmr.ciH<-as.data.frame(confint.merMod(disc.regrmrO))[4,2]
as.ciH<-as.data.frame(confint.merMod(disc.regas))[4,2]
fas.ciH<-as.data.frame(confint.merMod(disc.regfasO))[4,2]

# 2.5 and 97.5% CIs heart rates, bpm 
bpm.ciL<-as.data.frame(confint.merMod(modbpm))[4,1]
bpm.ciH<-as.data.frame(confint.merMod(modbpm))[4,2]



## 5.2. temp specific scaling, metabolism --------
# slopes, metabolism
slope12.rmr<-coef(reg.12.rmr)[2]
slope16.rmr<-coef(reg.16.rmr)[2]
slope20.rmr<-coef(reg.20.rmr)[2]
slope22.rmr<-coef(reg.22.rmr)[2]
slope12.mmr<-coef(reg.12.mmr)[2]
slope16.mmr<-coef(reg.16.mmr)[2]
slope20.mmr<-coef(reg.20.mmr)[2]
slope22.mmr<-coef(reg.22.mmr)[2]
slope12.as<-coef(reg.12.as)[2]
slope16.as<-coef(reg.16.as)[2]
slope20.as<-coef(reg.20.as)[2]
slope22.as<-coef(reg.22.as)[2]
slope12.fas<-coef(reg.12.fas)[2]
slope16.fas<-coef(reg.16.fas)[2]
slope20.fas<-coef(reg.20.fas)[2]
slope22.fas<-coef(reg.22.fas)[2]

int12.rmr<-coef(reg.12.rmr)[1]
int16.rmr<-coef(reg.16.rmr)[1]
int20.rmr<-coef(reg.20.rmr)[1]
int22.rmr<-coef(reg.22.rmr)[1]
int12.mmr<-coef(reg.12.mmr)[1]
int16.mmr<-coef(reg.16.mmr)[1]
int20.mmr<-coef(reg.20.mmr)[1]
int22.mmr<-coef(reg.22.mmr)[1]

# 2.5 % CI limits, metabolism
ciL.12.rmr<-as.data.frame(confint(reg.16.rmr))[2,1]
ciL.16.rmr<-as.data.frame(confint(reg.16.rmr))[2,1]
ciL.20.rmr<-as.data.frame(confint(reg.20.rmr))[2,1]
ciL.22.rmr<-as.data.frame(confint(reg.22.rmr))[2,1]
ciL.12.mmr<-as.data.frame(confint(reg.12.mmr))[2,1]
ciL.16.mmr<-as.data.frame(confint(reg.16.mmr))[2,1]
ciL.20.mmr<-as.data.frame(confint(reg.20.mmr))[2,1]
ciL.22.mmr<-as.data.frame(confint(reg.22.mmr))[2,1]
ciL.12.as<-as.data.frame(confint(reg.12.as))[2,1]
ciL.16.as<-as.data.frame(confint(reg.16.as))[2,1]
ciL.20.as<-as.data.frame(confint(reg.20.as))[2,1]
ciL.22.as<-as.data.frame(confint(reg.22.as))[2,1]
ciL.12.fas<-as.data.frame(confint(reg.12.fas))[2,1]
ciL.16.fas<-as.data.frame(confint(reg.16.fas))[2,1]
ciL.20.fas<-as.data.frame(confint(reg.20.fas))[2,1]
ciL.22.fas<-as.data.frame(confint(reg.22.fas))[2,1]

# 97.5 % CI limits metabolism 
ciH.12.rmr<-as.data.frame(confint(reg.16.rmr))[2,2]
ciH.16.rmr<-as.data.frame(confint(reg.16.rmr))[2,2]
ciH.20.rmr<-as.data.frame(confint(reg.20.rmr))[2,2]
ciH.22.rmr<-as.data.frame(confint(reg.22.rmr))[2,2]
ciH.12.mmr<-as.data.frame(confint(reg.12.mmr))[2,2]
ciH.16.mmr<-as.data.frame(confint(reg.16.mmr))[2,2]
ciH.20.mmr<-as.data.frame(confint(reg.20.mmr))[2,2]
ciH.22.mmr<-as.data.frame(confint(reg.22.mmr))[2,2]
ciH.12.as<-as.data.frame(confint(reg.12.as))[2,2]
ciH.16.as<-as.data.frame(confint(reg.16.as))[2,2]
ciH.20.as<-as.data.frame(confint(reg.20.as))[2,2]
ciH.22.as<-as.data.frame(confint(reg.22.as))[2,2]
ciH.12.fas<-as.data.frame(confint(reg.12.fas))[2,2]
ciH.16.fas<-as.data.frame(confint(reg.16.fas))[2,2]
ciH.20.fas<-as.data.frame(confint(reg.20.fas))[2,2]
ciH.22.fas<-as.data.frame(confint(reg.22.fas))[2,2]

## 5.3. temp specific scaling, heart rates --------
# slopes
hr16.b<-coef(m16)[2]
hr20.b<-coef(m20)[2]
hr22.b<-coef(m22)[2]
hr24.b<-coef(m24)[2]

# 2.5 % CI limits 
hr16.ciL<-as.data.frame(confint(m16))[2,1]
hr20.ciL<-as.data.frame(confint(m20))[2,1]
hr22.ciL<-as.data.frame(confint(m22))[2,1]
hr24.ciL<-as.data.frame(confint(m24))[2,1]

# 97.5 % CI limits 
hr16.ciH<-as.data.frame(confint(m16))[2,2]
hr20.ciH<-as.data.frame(confint(m20))[2,2]
hr22.ciH<-as.data.frame(confint(m22))[2,2]
hr24.ciH<-as.data.frame(confint(m24))[2,2]



## 5.4. lms; cardiac thermal tolerance metrics, heart mass ------------
# model slopes 
arr.b<-coef(mod.arr1)[2]
Tpeak.b<-coef(mod.peak1)[2]
HR.b<-coef(mod.HRp0)[2]
bp.b<-coef(mod.bpO)[2]
heart.b<-coef(heart.scale)[2]

heart.a<-coef(heart.scale)[1]


# 97.5 % CI limits 
arr.ciH<-as.data.frame(confint(mod.arr1))[2,2]
Tpeak.ciH<-as.data.frame(confint(mod.peak1))[2,2]
HR.ciH<-as.data.frame(confint(mod.HRp0))[2,2]
bp.ciH<-as.data.frame(confint(mod.bpO))[2,2]
heart.ciH<-as.data.frame(confint(heart.scale))[2,2]

# 2.5 % CI limits 
arr.ciL<-as.data.frame(confint(mod.arr1))[2,1]
Tpeak.ciL<-as.data.frame(confint(mod.peak1))[2,1]
HR.ciL<-as.data.frame(confint(mod.HRp0))[2,1]
bp.ciL<-as.data.frame(confint(mod.bpO))[2,1]
heart.ciL<-as.data.frame(confint(heart.scale))[2,1]


# 6. Plotting metrics, setup data frames. -----------

# ******
names<-c("mmr", "rmr", "as", "bpm",
         "mmr", "mmr", "mmr", "mmr", 
         "rmr", "rmr", "rmr", "rmr", 
         "as", "as", "as", "as", 
         "fas", "fas", "fas", "fas",
         "fas", "arr", "Tpeak",
         "HR", "bp", "heart",
         "bpm", "bpm", "bpm", "bpm")

temps<-c("all", "all", "all", "all",
         '12', '16', "20", "22",
         '12', '16', "20", "22",
         '12', '16', "20", "22",
         '12', '16', "20", "22",
         "all", "all", "all",
         "all", "all", "all",
         "16", "20", "22","24")

slopes<-c(mmr.b, rmr.b, as.b, bpm.b,
          slope12.mmr, slope16.mmr, slope20.mmr, slope22.mmr,
          slope12.rmr, slope16.rmr, slope20.rmr, slope22.rmr,
          slope12.as, slope16.as, slope20.as, slope22.as,
          slope12.fas, slope16.fas, slope20.fas, slope22.fas,
          fas.b,  arr.b,  Tpeak.b,
          HR.b,  bp.b,  heart.b,
          hr16.b, hr20.b, hr22.b, hr24.b)

ciL<-c(mmr.ciL, rmr.ciL, as.ciL, bpm.ciL,
       NA, NA, NA, NA,
       NA, NA, NA, NA,
       NA, NA, NA, NA,
       NA, NA, NA, NA,
       fas.ciL,  arr.ciL,  Tpeak.ciL,
       HR.ciL,  bp.ciL,  heart.ciL,
       NA, NA, NA, NA)

ciH<-c(mmr.ciH, rmr.ciH, as.ciH, bpm.ciH,
       NA, NA, NA, NA,
       NA, NA, NA, NA,
       NA, NA, NA, NA,
       NA, NA, NA, NA,
       fas.ciH,  arr.ciH,  Tpeak.ciH,
       HR.ciH,  bp.ciH,  heart.ciH,
       NA, NA, NA, NA)

ciL.all<-c(mmr.ciL, rmr.ciL, as.ciL, bpm.ciL,
           ciL.12.mmr, ciL.16.mmr, ciL.20.mmr, ciL.22.mmr,
           ciL.12.rmr, ciL.16.rmr, ciL.20.rmr, ciL.22.rmr,
           ciL.12.as, ciL.16.as, ciL.20.as, ciL.22.as,
           ciL.12.fas, ciL.16.fas, ciL.20.fas, ciL.22.fas,
           fas.ciL,  arr.ciL,  Tpeak.ciL,
           HR.ciL,  bp.ciL,  heart.ciL,
           hr16.ciL, hr20.ciL, hr22.ciL, hr24.ciL)

ciH.all<-c(mmr.ciH, rmr.ciH, as.ciH, bpm.ciH,
           ciH.12.mmr, ciH.16.mmr, ciH.20.mmr, ciH.22.mmr,
           ciH.12.rmr, ciH.16.rmr, ciH.20.rmr, ciH.22.rmr,
           ciH.12.as, ciH.16.as, ciH.20.as, ciH.22.as,
           ciH.12.fas, ciH.16.fas, ciH.20.fas, ciH.22.fas,
           fas.ciH,  arr.ciH,  Tpeak.ciH,
           HR.ciH,  bp.ciH,  heart.ciH,
           hr16.ciH, hr20.ciH, hr22.ciH, hr24.ciH)

color<-c("grey", "grey", "grey", "grey", 
         "#3596B5", "#49817B", "#CD6C95", "#A2416E",
         "#3596B5", "#49817B", "#CD6C95", "#A2416E",
         "#3596B5", "#49817B", "#CD6C95", "#A2416E",
         "#3596B5", "#49817B", "#CD6C95", "#A2416E",
         "grey", "grey", "grey",
         "grey", "grey", "grey",
         "#49817B", "#CD6C95", "#A2416E", "red")

fill<-c("black", "black", "black", "black", 
        "#3596B5", "#49817B", "#CD6C95", "#A2416E",
        "#3596B5", "#49817B", "#CD6C95", "#A2416E",
        "#3596B5", "#49817B", "#CD6C95", "#A2416E",
        "#3596B5", "#49817B", "#CD6C95", "#A2416E",
        "black", "black", "black",
        "black", "black", "black",
        "#49817B", "#CD6C95", "#A2416E", "red")


size <- c("B", "B","B","B",
          "A", "A", "A", "A",
          "A", "A", "A", "A",
          "A", "A", "A", "A",
          "A", "A", "A", "A",
          "B", "B","B",
          "B", "B","B",
          "A", "A", "A", "A")


fig7data<-data.frame(matrix(nrow = 30, ncol = 8),
                     row.names = 1:30)
colnames(fig7data)<-c("Performance","slopes", "CI.L", "CI.H", "color", "fill", "size", "temps")
fig7data$Performance<-names
fig7data$slopes<-slopes
fig7data$CI.L<-ciL
fig7data$CI.H<-ciH
fig7data$color<-color
fig7data$size<-size
fig7data$fill<-fill
fig7data$temps<-temps

fig7data <- fig7data[order(fig7data$slopes), ]
fig7data$orderPerf<-fct_reorder(fig7data$Performance, fig7data$slopes)

fig7data.top<-fig7data[fig7data$slopes>0.5, ]
fig7data.bottom<-fig7data[fig7data$slopes<0.5, ]
fig7data.top<-fig7data.top[!is.na(fig7data.top$slopes),]
fig7data.bottom<-fig7data.bottom[!is.na(fig7data.bottom$slopes),]
fig7data.top[fig7data.top$size== "B","slopes"]<-round(fig7data.top[fig7data.top$size== "B","slopes"], 3)
fig7data.bottom[fig7data.bottom$size== "B","slopes"]<-round(fig7data.bottom[fig7data.bottom$size== "B","slopes"], 3)

# create lists for plotting all temp specific plots:
slope.mmr.list<-list(slope12.mmr, slope16.mmr, slope20.mmr, slope22.mmr)
slope.rmr.list<-list(slope12.rmr, slope16.rmr, slope20.rmr, slope22.rmr)
int.mmr.list<-list(int12.mmr, int16.mmr, int20.mmr, int22.mmr)
int.rmr.list<-list(int12.rmr, int16.rmr, int20.rmr, int22.rmr)
data.frame.list<-list(data12, data16, data20, data22)



# size summary of data, (fish were morphometrically measured multiple times, get means)
size.sum<-data %>% 
  dplyr::group_by(FishID, sex, sizeClass, origin) %>% 
  dplyr:: summarize(BWg_m= mean(BW*1000, na.rm = T),
            BWg_sd= sd(BW*1000, na.rm = T),
            SL_m= mean(SL.cm, na.rm = T),
            SL_sd= sd(SL.cm, na.rm = T))


slopesCIdata<-data.frame(matrix(nrow = 30, ncol = 5),
                         row.names = 1:30)
colnames(slopesCIdata)<-c("Performance","slopes", "CI.L", "CI.H",  "temps")
slopesCIdata$Performance<-names
slopesCIdata$slopes<-slopes
slopesCIdata$CI.L<-ciL.all
slopesCIdata$CI.H<-ciH.all
slopesCIdata$temps<-temps


# means for final figure 7 
f.rmr<-pred.means.rmr1kg[pred.means.rmr1kg$origin == "field", ]
f.mmr<-pred.means.mmr1kg
f.as<-pred.means.as1kg
f.fas<-pred.means.fas1kg[pred.means.fas1kg$origin == "field",]
f.rmr$p<-"rmr"
f.mmr$p<-"mmr"
f.as$p<-"as"
f.fas$p<-"fas"
f.bpm<-pred.means.bpm1kg
f.bpm$p<-"bpm"

plot.final.mr<-rbind(f.rmr[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")],
                     f.mmr[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")],
                     f.as[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")],
                     f.bpm[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")],
                     f.fas[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")])

plot.final.mr$q10<-as.numeric(plot.final.mr$q10)
plot.final.mr$treatm<-as.numeric(as.character(plot.final.mr$treatm))
plot.final.mr2<-rbind(f.rmr[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")],
                      f.mmr[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")],
                      f.as[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")],
                      f.fas[, c("BW", "treatm","prediction", "SE", "df", "n", "q10", "p")])
plot.final.mr2$treatm<-as.factor(as.character(plot.final.mr2$treatm))



# Field fish only, for plotting purposes:
dataF<-data[c(data$origin == "field"),]

dataF12<-dataF[dataF$treatm=="12",]
dataF20<-dataF[dataF$treatm=="20",]
dataF16<-dataF[dataF$treatm=="16",]
dataF22<-dataF[dataF$treatm=="22",]
dataF24<-dataF[dataF$treatm=="24",]




# 7. Figures ---------

## Figure 3 AB ------
thesis2A<-ggformat(plot = (
  ggplot(data=data, aes(y=log(mmr), x=log(BW), shape = pregnant))+
    geom_point( color="black", fill = cols.4[4], size=2, show.legend = F)+
    geom_point(data=data20, aes(y=log(mmr), x=log(BW)),  color="black", fill = cols.4[3], size=2, show.legend = F)+
    geom_point(data=data16, aes(y=log(mmr), x=log(BW)),  color="black", fill = cols.4[2], size=2, show.legend = F)+
    geom_point(data=data12, aes(y=log(mmr), x=log(BW)),  color="black", fill = cols.4[1], size=2, show.legend = F)+
    scale_shape_manual(values = c(21, 22))+
    theme_classic()+
    ylim(-4, 2)+
    xlim(-6, 0)+
    geom_abline(slope = mmr.b, intercept = scaling.int.mmr, color = cols.4[1], lty="solid")+
    geom_abline(slope = mmr.b, intercept = scaling.int.mmr + fixef(disc.regmmr)[3], color = cols.4[2], lty="solid")+
    geom_abline(slope = mmr.b, intercept = scaling.int.mmr + fixef(disc.regmmr)[4], color = cols.4[3], lty="solid")+
    geom_abline(slope = mmr.b, intercept = scaling.int.mmr + fixef(disc.regmmr)[5], color = cols.4[4], lty="solid")+
    annotate(geom = "text",x = -5.6, y = 1.9, hjust = 0,
             label = deparse(bquote(~MMR == italic(ln) ~"("* .(round(as.numeric(scaling.int.mmr),2))*")"~ BM^.(round(as.numeric(mmr.b),2)))),
             parse = T, color = cols.4[1])+
    annotate(geom = "text",x = -5.6, y = 1.5, hjust = 0,
             label = deparse(bquote(~MMR == italic(ln) ~"("* .(round(as.numeric(scaling.int.mmr + fixef(disc.regmmr)[3]),2))*")"~ BM^.(round(as.numeric(mmr.b),2)))),
             parse = T, color = cols.4[2])+
    annotate(geom = "text",x = -5.6, y = 1.1, hjust = 0,
             label = deparse(bquote(~MMR == italic(ln) ~"("* .(round(as.numeric(scaling.int.mmr + fixef(disc.regmmr)[4]),2))*")"~ BM^.(round(as.numeric(mmr.b),2)))),
             parse = T, color = cols.4[3])+
    annotate(geom = "text",x = -5.6, y = 0.7, hjust = 0,
             label = deparse(bquote(~MMR == italic(ln) ~"("* .(round(as.numeric(scaling.int.mmr + fixef(disc.regmmr)[5]),2))*")"~ BM^.(round(as.numeric(mmr.b),2)))),
             parse = T, color = cols.4[4])
),y_title = expression(italic(ln)~MMR), x_title = expression(italic(ln)~Body~mass~(kg)))

thesis2B<-ggformat(plot = (
  ggplot(data=data, aes(y=log(rmr), x=log(BW), shape = pregnant))+
    geom_point( color="black", fill = cols.4[4], size=2, show.legend = F)+
    geom_point(data=data20, aes(y=log(rmr), x=log(BW)),  color="black", fill = cols.4[3], size=2, show.legend = F)+
    geom_point(data=data16, aes(y=log(rmr), x=log(BW)),  color="black", fill = cols.4[2], size=2, show.legend = F)+
    geom_point(data=data12, aes(y=log(rmr), x=log(BW)),  color="black", fill = cols.4[1], size=2, show.legend = F)+
    scale_shape_manual(values = c(21, 22))+
    theme_classic()+
    ylim(-4, 2)+
    xlim(-6, 0)+
    geom_abline(slope = rmr.b, intercept = scaling.int.rmr, color = cols.4[1], lty="solid")+
    geom_abline(slope = rmr.b, intercept = scaling.int.rmr + fixef(disc.regrmr)[3], color = cols.4[2], lty="solid")+
    geom_abline(slope = rmr.b, intercept = scaling.int.rmr + fixef(disc.regrmr)[4], color = cols.4[3], lty="solid")+
    geom_abline(slope = rmr.b, intercept = scaling.int.rmr + fixef(disc.regrmr)[5], color = cols.4[4], lty="solid")+
    annotate(geom = "text",x = -5.6, y = 1.9, hjust = 0,
             label = deparse(bquote(~RMR == italic(ln) ~"("* .(round(as.numeric(scaling.int.rmr),2))*")"~ BM^.(round(as.numeric(rmr.b),2)))),
             parse = T, color = cols.4[1])+
    annotate(geom = "text",x = -5.6, y = 1.5, hjust = 0,
             label = deparse(bquote(~RMR == italic(ln) ~"("* .(round(as.numeric(scaling.int.rmr + fixef(disc.regrmr)[3]),2))*")"~ BM^.(round(as.numeric(rmr.b),2)))),
             parse = T, color = cols.4[2])+
    annotate(geom = "text",x = -5.6, y = 1.1, hjust = 0,
             label = deparse(bquote(~RMR == italic(ln) ~"("* .(round(as.numeric(scaling.int.rmr + fixef(disc.regrmr)[4]),2))*")"~ BM^.(round(as.numeric(rmr.b),2)))),
             parse = T, color = cols.4[3])+
    annotate(geom = "text",x = -5.6, y = 0.7, hjust = 0,
             label = deparse(bquote(~RMR == italic(ln) ~"("* .(round(as.numeric(scaling.int.rmr + fixef(disc.regrmr)[5]),2))*")"~ BM^.(round(as.numeric(rmr.b),2)))),
             parse = T, color = cols.4[4])
),y_title = expression(italic(ln)~RMR), x_title = expression(italic(ln)~Body~mass~(kg)))




## Figure 3 CDEF------
thesis2C<-ggformat(
  plot = (
    ggplot(data=data12, aes(y=log(rmr), x=log(BW)))+
      geom_point(pch=23, color=cols.4[1], fill = "white", size=2, alpha = 1)+
      geom_point(data = data12, mapping = (aes(y=log(mmr), x=log(BW))),  pch=21, color=cols.4[1], fill = "white", size=2, alpha = 1)+
      geom_point(data = dataF12, mapping = (aes(y=log(rmr), x=log(BW))),  pch=23, color=cols.4[1], fill = cols.4[1], size=2, alpha = 1)+
      geom_point(data = dataF12, mapping = (aes(y=log(mmr), x=log(BW))),  pch=21, color=cols.4[1], fill = cols.4[1], size=2, alpha = 1)+
      theme_classic()+
      ylim(-4, 2)+
      xlim(-6, 1)+
      geom_abline(slope = slope12.rmr, intercept = int12.rmr, color = cols.4[1])+
      geom_abline(slope = slope12.mmr, intercept = int12.mmr, color = cols.4[1])+
      annotate("text", x = -6, y = 1.6, hjust = 0, label = expression(bold(C)), size = 5, parse = T)+
      annotate("text", x = -0.75, y = -3.5, label = expression(12*degree*C), color = cols.4[1], size = 5)+
      annotate("text", x = -5.25, y = 1.6, hjust = 0, label = deparse(bquote(~italic(b)[MMR] ~"="~ .(round(as.numeric(slope12.mmr), 3)))), size = 4, parse = T)+
      annotate("text", x = -5.25, y = 1, hjust = 0, label = deparse(bquote(~italic(b)[RMR] ~"="~ .(round(as.numeric(slope12.rmr), 3)))), size = 4, parse = T)
  ), y_title = expression(italic(ln)~MR), x_title = expression(italic(ln)~Body~mass~(kg)), print = T)

thesis2D<-ggformat(
  plot = (
    ggplot(data=data16, aes(y=log(rmr), x=log(BW)))+
      geom_point(pch=23, color=cols.4[2], fill = "white", size=2, alpha = 1)+
      geom_point(data = data16, mapping = (aes(y=log(mmr), x=log(BW))),  pch=21, color=cols.4[2], fill = "white", size=2, alpha = 1)+
      geom_point(data = dataF16, mapping = (aes(y=log(rmr), x=log(BW))),  pch=23, color=cols.4[2], fill = cols.4[2], size=2, alpha = 1)+
      geom_point(data = dataF16, mapping = (aes(y=log(mmr), x=log(BW))),  pch=21, color=cols.4[2], fill = cols.4[2], size=2, alpha = 1)+
      theme_classic()+
      ylim(-4, 2)+
      xlim(-6, 0)+
      geom_abline(slope = slope16.rmr, intercept = int16.rmr, color = cols.4[2])+
      geom_abline(slope = slope16.mmr, intercept = int16.mmr, color = cols.4[2])+
      annotate("text", x = -6, y = 1.6, hjust = 0, label = expression(bold(D)), size = 5, parse = T)+
      annotate("text", x = -0.75, y = -3.5, label = expression(16*degree*C), color = cols.4[2], size = 5)+
      annotate("text", x = -5.25, y = 1.6, hjust = 0, label = deparse(bquote(~italic(b)[MMR] ~"="~ .(round(as.numeric(slope16.mmr), 3)))), size = 4, parse = T)+
      annotate("text", x = -5.25, y = 1, hjust = 0, label = deparse(bquote(~italic(b)[RMR] ~"="~ .(round(as.numeric(slope16.rmr), 3)))), size = 4, parse = T)
  ), y_title = expression(italic(ln)~MR), x_title = expression(italic(ln)~Body~mass~(kg)), print = F)


thesis2E<-ggformat(
  plot = (
    ggplot(data=data20, aes(y=log(rmr), x=log(BW)))+
      geom_point(pch=23, color=cols.4[3], fill = "white", size=2, alpha = 1)+
      geom_point(data = data20, mapping = (aes(y=log(mmr), x=log(BW))),  pch=21, color=cols.4[3], fill = "white", size=2, alpha = 1)+
      geom_point(data = dataF20, mapping = (aes(y=log(rmr), x=log(BW))),  pch=23, color=cols.4[3], fill = cols.4[3], size=2, alpha = 1)+
      geom_point(data = dataF20, mapping = (aes(y=log(mmr), x=log(BW))),  pch=21, color=cols.4[3], fill = cols.4[3], size=2, alpha = 1)+
      theme_classic()+
      ylim(-4, 2)+
      xlim(-6, 0)+
      geom_abline(slope = slope20.rmr, intercept = int20.rmr, color = cols.4[3])+
      geom_abline(slope = slope20.mmr, intercept = int20.mmr, color = cols.4[3])+
      annotate("text", x = -6, y = 1.6, hjust = 0, label = expression(bold(E)), size = 5, parse = T)+
      annotate("text", x = -0.75, y = -3.5, label = expression(20*degree*C), color = cols.4[3], size = 5)+
      annotate("text", x = -5.25, y = 1.6, hjust = 0, label = deparse(bquote(~italic(b)[MMR] ~"="~ .(round(as.numeric(slope20.mmr), 3)))), size = 4, parse = T)+
      annotate("text", x = -5.25, y = 1, hjust = 0, label = deparse(bquote(~italic(b)[RMR] ~"="~ .(round(as.numeric(slope20.rmr), 3)))), size = 4, parse = T)
  ), y_title = expression(italic(ln)~MR), x_title = expression(italic(ln)~Body~mass~(kg)), print = F)

thesis2F<-ggformat(
  plot = (
    ggplot(data=data22, aes(y=log(rmr), x=log(BW)))+
      geom_point(pch=23, color=cols.4[4], fill = "white", size=2, alpha = 1)+
      geom_point(data = data22, mapping = (aes(y=log(mmr), x=log(BW))),  pch=21, color=cols.4[4], fill = "white", size=2, alpha = 1)+
      geom_point(data = dataF22, mapping = (aes(y=log(rmr), x=log(BW))),  pch=23, color=cols.4[4], fill = cols.4[4], size=2, alpha = 1)+
      geom_point(data = dataF22, mapping = (aes(y=log(mmr), x=log(BW))),  pch=21, color=cols.4[4], fill = cols.4[4], size=2, alpha = 1)+
      theme_classic()+
      ylim(-4, 2)+
      xlim(-6, 0)+
      geom_abline(slope = slope22.rmr, intercept = int22.rmr, color = cols.4[4])+
      geom_abline(slope = slope22.mmr, intercept = int22.mmr, color = cols.4[4])+
      annotate("text", x = -6, y = 1.6, hjust = 0, label = expression(bold(F)), size = 5, parse = T)+
      annotate("text", x = -0.75, y = -3.5, label = expression(22*degree*C), color = cols.4[4], size = 5)+
      annotate("text", x = -5.25, y = 1.6, hjust = 0, label = deparse(bquote(~italic(b)[MMR] ~"="~ .(round(as.numeric(slope22.mmr), 3)))), size = 4, parse = T)+
      annotate("text", x = -5.25, y = 1, hjust = 0, label = deparse(bquote(~italic(b)[RMR] ~"="~ .(round(as.numeric(slope22.rmr), 3)))), size = 4, parse = T)
  ), y_title = expression(italic(ln)~MR), x_title = expression(italic(ln)~Body~mass~(kg)), print = F)

thesis2C <- thesis2C + theme( axis.title.y = element_blank(), 
                              axis.title.x = element_blank())
thesis2D <- thesis2D + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                             axis.title.x = element_blank())
thesis2E <- thesis2E + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                             axis.title.x = element_blank())
thesis2F <- thesis2F + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                             axis.title.x = element_blank())


## Figure 4 ABCD ------
plot.size <- ggplot(data.abt, aes(x = temp_mean, y = bpm, group=FishID, size = mass.g/1000, color = mass.g/1000, fill = mass.g/1000))+
  geom_errorbar(aes(ymin = bpm - bpm.SD, ymax = bpm + bpm.SD, group=FishID), size=0.5)+
  geom_point(pch=21, size=1)+
  geom_line(lwd=0.5)+
  annotate("segment",  x = 16, xend = 16, color=cols.4[2], y = 60, yend = 70, size=1.5, arrow = arrow( angle = 40, length = unit(0.25,"cm")))+
  annotate("segment",  x = 20, xend = 20, color=cols.4[3], y = 60, yend = 70, size=1.5, arrow = arrow( angle = 40, length = unit(0.25,"cm")))+
  annotate("segment",  x = 22, xend = 22, color=cols.4[4], y = 60, yend = 70, size=1.5, arrow = arrow( angle = 40, length = unit(0.25,"cm")))+
  annotate("segment",  x = 24, xend = 24, color="red", y = 60, yend = 70, size=1.5, arrow = arrow( angle = 40, length = unit(0.25,"cm")))+
  scale_color_gradient(low = "#97AFB9", high = "#002D47", name = "kg")+
  scale_fill_gradient(low = "#97AFB9", high = "#002D47", name = "kg")+
  scale_x_continuous(limits = c(14, 28), breaks = c(14, 16, 18, 20, 22, 24, 26, 28))+
  scale_y_continuous(limits = c(58, 180), breaks = c(60, 80, 100, 120, 140, 160))+
  annotate("text", x = -4.2, y = 4.35, hjust = 0, label = deparse(bquote(~italic(b) ~"="~ .(round(as.numeric(hr16.b), 3)))), size = 4, parse = T, color = cols.4[2])
ggformat(plot.size, print = F, y_title = expression(italic(f)[Hmax]~(beats~min^-1)), x_title = expression(Temperature~degree*C), size_text = 15)

plot.size<-plot.size+theme(legend.position = c(0.32, 0.81),
                           legend.direction = "horizontal", 
                           legend.key.width = unit(0.7,"cm"), 
                           legend.key.height = unit(0.3,"cm"), 
                           legend.margin = margin(0, 0, 0, 0, "cm"))


plot.bpm.lmer<-ggplot()+
  geom_label(pred.means.bpm1kg, mapping = aes(x = treatm, y = exp(prediction), label = round(q10,2), fill = NULL), label.size = NA, nudge_x = 0, nudge_y = -5.5, size = 3)+
  geom_line(pred.means.bpm1kg, mapping = aes(y=exp(prediction), x=treatm, group= factor(BW)), alpha=1, show.legend = F, color = "black", position = position_dodge(width = 0.4))+
  geom_errorbar(pred.means.bpm1kg, mapping = aes(ymax = exp(prediction+SE), ymin = exp(prediction-SE), x=treatm), width = 0.5, inherit.aes = F, color = "black", position = position_dodge(width = 0.4), show.legend = F)+
  geom_point(pred.means.bpm1kg, mapping = aes(y=exp(prediction), x=treatm),
             show.legend = F, alpha = 1, size = 3, fill = "grey70", color = "black", pch=21, position = position_dodge(width = 0.4))+
  geom_point(pred.means.bpm1kg[pred.means.bpm1kg$treatm == "16",], mapping = aes(y=exp(prediction), x=treatm),
             show.legend = F, alpha = 1, size = 3, fill = cols.4[2], color = "black", pch=21, position = position_dodge(width = 0.4))+
  geom_point(pred.means.bpm1kg[pred.means.bpm1kg$treatm == "20",], mapping = aes(y=exp(prediction), x=treatm),
             show.legend = F, alpha = 1, size = 3, fill = cols.4[3], color = "black", pch=21, position = position_dodge(width = 0.4))+
  geom_point(pred.means.bpm1kg[pred.means.bpm1kg$treatm == "22",], mapping = aes(y=exp(prediction), x=treatm),
             show.legend = F, alpha = 1, size = 3, fill = cols.4[4], color = "black", pch=21, position = position_dodge(width = 0.4))+
  geom_point(pred.means.bpm1kg[pred.means.bpm1kg$treatm == "24",], mapping = aes(y=exp(prediction), x=treatm),
             show.legend = F, alpha = 1, size = 3, fill = "red", color = "black", pch=21, position = position_dodge(width = 0.4))+
  annotate("text", x = 9, y = 75, fontface = "italic", label = "Mass-adjusted \n for 1kg fish")+
  scale_x_discrete( breaks = c("14", "16", "18", "20", "22", "24", "26", "28"))+
  theme_classic()
ggformat(plot.bpm.lmer, y_title = expression(italic(f)[Hmax]~(beats~min^-1)), x_title = "Temperature ºC", size_text = 15)

plot.abt.scale<-ggplot(data.abt, aes(x = log(mass.g/1000), y = log(bpm), fill = TEMP, color = TEMP, group = TEMP))+
  scale_y_continuous(limits = c(log(65), log(160)))+
  geom_point(data.abt[data.abt$TEMP == 16,], mapping = aes(x = log(mass.g/1000), y = log(bpm), group = TEMP), fill = cols.4[2], pch=21, size=2, color="black", show.legend = F)+
  geom_smooth(data.abt[data.abt$TEMP == 16,], mapping = aes(x = log(mass.g/1000), y = log(bpm), group = TEMP), method = "lm", se = FALSE, color = cols.4[2], show.legend = F)+
  geom_point(data.abt[data.abt$TEMP == 20,], mapping = aes(x = log(mass.g/1000), y = log(bpm), group = TEMP), fill = cols.4[3], pch=21,size=2, color="black", show.legend = F)+
  geom_smooth(data.abt[data.abt$TEMP == 20,], mapping = aes(x = log(mass.g/1000), y = log(bpm), group = TEMP), method = "lm", se = FALSE, color = cols.4[3], show.legend = F)+
  geom_point(data.abt[data.abt$TEMP == 24,], mapping = aes(x = log(mass.g/1000), y = log(bpm), group = TEMP), fill = "red", pch=21,size=2, color="black", show.legend = F)+
  geom_smooth(data.abt[data.abt$TEMP == 24,], mapping = aes(x = log(mass.g/1000), y = log(bpm), group = TEMP), method = "lm", se = FALSE, color = "red", show.legend = F)+
  geom_point(data.abt[data.abt$TEMP == 22,], mapping = aes(x = log(mass.g/1000), y = log(bpm), group = TEMP), fill = cols.4[4], pch=21,size=2, color="black", show.legend = F)+
  geom_smooth(data.abt[data.abt$TEMP == 22,], mapping = aes(x = log(mass.g/1000), y = log(bpm), group = TEMP), method = "lm", se = FALSE, color = cols.4[4], show.legend = F)+
  xlim(-5, -1)+
  annotate("text", x = -5, hjust=0, y = 4.35, label = "16ºC", color = cols.4[2], size = 4)+
  annotate("text", x = -5, hjust=0, y = 4.30, label = "20ºC", color = cols.4[3], size = 4)+
  annotate("text", x = -5, hjust=0, y = 4.25, label = "22ºC", color = cols.4[4], size = 4)+
  annotate("text", x = -5, hjust=0, y = 4.2, label = "24ºC", color = "red", size = 4)+
  annotate("text", x = -4.2, y = 4.35, hjust = 0, label = deparse(bquote(~italic(b) ~"="~ .(round(as.numeric(hr16.b), 3)))), size = 4, parse = T, color = cols.4[2])+
  annotate("text", x = -4.2, y = 4.30, hjust = 0, label = deparse(bquote(~italic(b) ~"="~ .(round(as.numeric(hr20.b), 3)))), size = 4, parse = T, color = cols.4[3])+
  annotate("text", x = -4.2, y = 4.25, hjust = 0, label = deparse(bquote(~italic(b) ~"="~ .(round(as.numeric(hr22.b), 3)))), size = 4, parse = T, color = cols.4[4])+
  annotate("text", x = -4.2, y = 4.2, hjust = 0, label = deparse(bquote(~italic(b) ~"="~ .(round(as.numeric(hr24.b), 3)))), size = 4, parse = T, color = "red")
ggformat(plot.abt.scale, y_title = expression(italic(f)[Hmax]~(beats~min^-1)),
         x_title = expression(italic(ln)~Body~mass~(kg)),
         size_text = 15, print = F)

plot2<-ggplot(data.abt, aes(x = log(BW), y = log(bpm), fill = TEMP, color = TEMP, group = TEMP))+
  geom_point(pch=21, color = "black", size=2, alpha = 0.7)+
  xlim(-5, -1)+
  geom_line(data = data.bpm.pred2, mapping = aes(x = log(BW), y = pre, color = treatm), size = 0.8, show.legend = F)+
  scale_color_manual(values =c(cols.4[2:4], "red"))+
  scale_fill_gradient(low = "grey90", high = "grey35", name = "ºC")+
  scale_y_continuous(limits = c(log(65), log(160)))+
  annotate("text", x = -4.5, y = 4.35, hjust = 0, label = deparse(bquote(~italic(b) ~"="~ .(round(as.numeric(bpm.b), 3)))), size = 4, parse = T, color = "black")
ggformat(plot2, y_title = expression(italic(f)[Hmax]~(beats~min^-1)),
         x_title = expression(italic(ln)~Body~mass~(kg)), size_text = 15, print = F)

plot2<-plot2+theme(legend.position = c(0.32, 0.08),
                   legend.direction = "horizontal", 
                   legend.key.width = unit(0.7,"cm"), 
                   legend.key.height = unit(0.3,"cm"))




## Figure 5 ABCD: ------

temparrh2<-ggformat(
  plot = (ggplot(data.abtID, aes(y = log(temp_ARRH), x = log(BW), fill= origin))+
            geom_abline(intercept = coef(mod.arr1)[1], slope = coef(mod.arr1)[2], color = "black", size = 1)+
            geom_point( color = "black", show.legend = F, pch = 21, size = 2)+
            annotate(geom = "text",x = -3.8, y = 3, hjust = 0,
                     label = deparse(bquote(~italic(ln)~T[ARR] == italic(ln) ~"("* .(round(as.numeric(coef(mod.arr1)[1]),2))*")"~ BM^.(round(as.numeric(coef(mod.arr1)[2]),2)))),
                     parse = T, color = "black")+
            annotate(geom = "text",x = -3.8, y = 2.96, hjust = 0,
                     label = " n = 29", color = "black")+
            # annotate(geom = "text", x = -3.8, y = 2.96, hjust = 0, label = "*origin: p<0.001", size = 2.5)+
            xlim(-6, -0.5)+
            ylim(2.89, 3.4)+
            scale_fill_manual(values = c("black", "white"))
  ),
  x_title = expression(italic(ln)~Body~mass*~(kg)),
  y_title = expression(italic(ln)~T[ARR]*~(degree*C)), size_text = 15, print = F)

tpeak2<-ggformat(
  plot = (ggplot(data.abtID, aes(y = log(Tpeak), x = log(BW), fill= origin))+
            geom_abline(intercept = coef(mod.peak1)[1], slope = coef(mod.peak1)[2], color = "black", size = 1)+
            geom_point( color = "black", show.legend = F, pch = 21, size = 2)+
            annotate(geom = "text",x = -3.8, y = 3, hjust = 0,
                     label = deparse(bquote(~italic(ln)~T[PEAK] == italic(ln) ~"("* .(round(as.numeric(coef(mod.peak1)[1]),2))*")"~ BM^.(round(as.numeric(coef(mod.peak1)[2]),2)))),
                     parse = T, color = "black")+
            annotate(geom = "text",x = -3.8, y = 2.96, hjust = 0,
                     label = " n = 30", color = "black")+
            xlim(-6, -0.5)+
            ylim(2.89, 3.4)+
            scale_fill_manual(values = c("black", "white"))
  ),
  x_title = expression(italic(ln)~Body~mass*~(kg)),
  y_title = expression(italic(ln)~T[PEAK]*~(degree*C)), size_text = 15, print = F)

tabt2<-ggformat(
  plot = (ggplot(data.abtID, aes(y = log(breakpoint_Cels), x = log(BW), fill= origin))+
            geom_abline(intercept = coef(mod.bpO)[1], slope = coef(mod.bpO)[2], color = "black", size = 1)+
            geom_point( color = "black", show.legend = F, pch = 21, size = 2)+
            annotate(geom = "text",x = -3.8, y = 3.35, hjust = 0,
                     label = deparse(bquote(~italic(ln)~T[AB] == italic(ln) ~"("* .(round(as.numeric(coef(mod.bpO)[1]),2))*")"~ BM^.(round(as.numeric(coef(mod.bpO)[2]),2)))),
                     parse = T, color = "black")+
            annotate(geom = "text", x = -3.8, y = 3.31, hjust = 0, label = "*origin: p<0.05", size = 3)+
            annotate(geom = "text",x = -3.8, y = 3.27, hjust = 0,
                     label = " n = 27", color = "black")+
            xlim(-6, -0.5)+
            ylim(2.89, 3.4)+
            scale_fill_manual(values = c("black", "white"))
  ),
  x_title = expression(italic(ln)~Body~mass*~(kg)),
  y_title = expression(italic(ln)~T[AB]*~(degree*C)), size_text = 15, print = F)

hrmax2<-ggformat(
  plot = (ggplot(data.abtID, aes(y = log(HRpeak), x = log(BW), fill= origin))+
            # geom_abline(intercept = coef(mod.HRp0)[1], slope = coef(mod.HRp0)[2], color = "black", size = 1)+
            geom_point( color = "black", show.legend = F, pch = 21, size = 2)+
            # annotate(geom = "text",x = -4.25, y = 4.74, hjust = 0,
            #          label = deparse(bquote(~italic(ln)~f[Hmax] == italic(ln) ~"("* .(round(as.numeric(coef(mod.HRp0)[1]),2))*")"~ BM^.(round(as.numeric(coef(mod.HRp0)[2]),2)))),
            #          parse = T, color = "black")+
            annotate(geom = "text", x = -1.9, y = 5.12, hjust = 0, label = "n = 30")+
            xlim(-6, -0.5)+
            ylim(4.71,5.13)+
            scale_fill_manual(values = c("black", "white"))),
  x_title = expression(italic(ln)~Body~mass~(kg)),
  y_title = expression(italic(ln)~PEAK[italic(f)*Hmax]~(beats~min^-1)))

## Figure 6 ABCD ------ 

plot.rmr.lmer<-ggplot(data=data.rmr.pred, aes(y=exp(rmr_pred), x=as.numeric(as.character(treatm)), shape = origin, fill = factor(treatm), group=factor(origin), color = factor(treatm)))+
  geom_line(pred.means.rmr1kg, mapping = aes(y=exp(prediction), x=as.numeric(as.character(treatm)), linetype = factor(origin), group=factor(origin)),
            alpha=1, show.legend = F, color = "black", position = position_dodge(width = 0.4))+
  # geom_boxplot(mapping= aes(group = factor(treatm)), alpha=0.4, show.legend = F)+
  geom_point(show.legend = F, alpha = 0.4, size = 2, position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2))+
  geom_errorbar(pred.means.rmr1kg, mapping = aes(ymax = exp(prediction+SE), ymin = exp(prediction-SE), group = factor(origin), x=as.numeric(as.character(treatm))), width = 0.2, inherit.aes = F, color = "black", position = position_dodge(width = 0.4), show.legend = F)+
  geom_point(pred.means.rmr1kg, mapping = aes(y=exp(prediction), x=as.numeric(as.character(treatm)), shape = origin, fill = factor(treatm), group=factor(origin), color = factor(treatm)),
             show.legend = F, alpha = 1, size = 4, color = "black", position = position_dodge(width = 0.4))+
  scale_shape_manual(values = c(21,23))+
  scale_x_continuous(breaks = c(12, 14, 16,18,20, 22)) +
  scale_fill_manual(values = c(cols.4))+
  scale_color_manual(values = c(cols.4))+
  ylim(0,4.1)+
  theme_classic()
ggformat(print =T, plot.rmr.lmer, y_title = expression(RMR~(mgO[2]~min^-1~kg^-1)), x_title = "Temperature ºC", size_text = 15)
# plot.rmr.lmer<-plot.rmr.lmer+theme(legend.position = c(0.8, 0.2))

plot.fas.lmer<-ggplot(data=data.fas.pred, aes(y=exp(fas_pred), x=as.numeric(as.character(treatm)), shape = origin, fill = factor(treatm), group=factor(origin), color = factor(treatm)))+
  # geom_line(size=0.2, alpha=0.6)+
  # geom_boxplot(mapping= aes(group = factor(treatm)), alpha=0.4, show.legend = F)+
  geom_line(pred.means.fas1kg, mapping = aes(y=exp(prediction), x=as.numeric(as.character(treatm)), linetype = factor(origin), group=factor(origin)),
            alpha=1, show.legend = F, color = "black", position = position_dodge(width = 0.4))+
  geom_point(show.legend = F, alpha = 0.4, size = 2, position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2))+
  geom_errorbar(pred.means.fas1kg, mapping = aes(ymax = exp(prediction+SE), ymin = exp(prediction-SE), group = factor(origin), x=as.numeric(as.character(treatm))), width = 0.2, inherit.aes = F, color = "black", position = position_dodge(width = 0.4), show.legend = F)+
  geom_point(pred.means.fas1kg, mapping = aes(y=exp(prediction), x=as.numeric(as.character(treatm)), shape = origin, fill = factor(treatm), group=factor(origin), color = factor(treatm)),
             show.legend = F, alpha = 1, size = 4, color = "black", position = position_dodge(width = 0.4))+
  scale_shape_manual(values = c(21,23))+
  scale_fill_manual(values = c(cols.4))+
  scale_x_continuous(breaks = c(12, 14, 16,18,20, 22)) +
  scale_color_manual(values = c(cols.4))+
  ylim(1,3.1)+
  theme_classic()
ggformat(plot.fas.lmer, y_title = "FAS (MMR / RMR)", x_title = "Temperature ºC", size_text = 15)


plot.AS.lmer<-ggplot(data=data.as.pred, mapping = aes(y=exp(as_pred), x=as.numeric(as.character(treatm)), fill = factor(treatm), color = factor(treatm), group =  factor(treatm)))+
  # geom_line(size=0.2, alpha=0.6)+
  # geom_boxplot(mapping= aes(group = factor(treatm)), alpha=0.4, show.legend = F)+
  geom_line(pred.means.as1kg, mapping = aes(y=exp(prediction), x=as.numeric(as.character(treatm)), group= factor(BW)), alpha=1, show.legend = F, color = "black", position = position_dodge(width = 0.4))+
  geom_jitter(show.legend = F, alpha = 0.4, size = 2,  pch=21, width = 0.05)+
  geom_errorbar(pred.means.as1kg, mapping = aes(ymax = exp(prediction+SE), ymin = exp(prediction-SE), x=as.numeric(as.character(treatm))), width = 0.2, inherit.aes = F, color = "black", position = position_dodge(width = 0.4), show.legend = F)+
  geom_point(pred.means.as1kg, mapping = aes(y=exp(prediction), x=as.numeric(as.character(treatm)),  fill = factor(treatm), color = factor(treatm)),
             show.legend = F, alpha = 1, size = 4, color = "black", pch=21, position = position_dodge(width = 0.4))+
  # scale_shape_manual(values = c(21, 21, 21, 21))+
  scale_fill_manual(values = c(cols.4))+
  scale_color_manual(values = c(cols.4))+
  scale_x_continuous(breaks = c(12, 14, 16,18,20, 22)) +
  theme_classic()+
  ylim(0,4.1)
ggformat(plot.AS.lmer, y_title = expression(AAS~(mgO[2]~min^-1~kg^-1)), x_title = "Temperature ºC", size_text = 15)

plot.mmr.lmer<-ggplot(data=data.mmr.pred, aes(y=exp(mmr_pred), x=as.numeric(as.character(treatm)), fill = factor(treatm),color = factor(treatm)))+
  # geom_line(size=0.2, alpha=0.6)+
  # geom_boxplot(mapping= aes(group = factor(treatm)), alpha=0.4, show.legend = F)+
  geom_line(pred.means.mmr1kg, mapping = aes(y=exp(prediction), x=as.numeric(as.character(treatm)), group= factor(BW)), alpha=1, show.legend = F, color = "black", position = position_dodge(width = 0.4))+
  geom_jitter(show.legend = F, alpha = 0.4, size = 2, pch=21, width = 0.05)+
  geom_errorbar(pred.means.mmr1kg, mapping = aes(ymax = exp(prediction+SE), ymin = exp(prediction-SE), x=as.numeric(as.character(treatm))), width = 0.2, inherit.aes = F, color = "black", position = position_dodge(width = 0.4), show.legend = F)+
  geom_point(pred.means.mmr1kg, mapping = aes(y=exp(prediction), x=as.numeric(as.character(treatm)),  fill = factor(treatm), color = factor(treatm)),
             show.legend = F, alpha = 1, size = 4, color = "black", pch=21, position = position_dodge(width = 0.4))+
  # scale_shape_manual(values = c(21, 21, 21, 21))+
  scale_fill_manual(values = c(cols.4))+
  scale_color_manual(values = c(cols.4))+
  scale_x_continuous(breaks = c(12, 14, 16,18,20, 22)) +
  ylim(0,6)+
  theme_classic()
ggformat(plot.mmr.lmer, y_title = expression(MMR~(mgO[2]~min^-1~kg^-1)), x_title = "Temperature ºC", size_text = 15)




## Figure 7 Final ABCD----- 

p1<-ggformat(plot=
               ggplot(fig7data.top, aes(x = as.numeric(slopes),  y = orderPerf, fill = color,  size = factor(size), shape = factor(size)))+
               geom_errorbarh(mapping = aes(xmin = CI.L, xmax = CI.H,),
                              color = "black", size = 0.6, height = 0.1)+
               xlim(0.66, 1)+
               geom_point( color = "grey30", show.legend = F)+
               scale_size_manual(values = c(3,4))+
               geom_text(size=4, color = "black", nudge_y = 0.3,  data =  fig7data.top[fig7data.top$size== "B",] ,
                         aes(x = slopes, label = slopes,  y = orderPerf))+
               scale_shape_manual(values = c(21, 2))+
               scale_fill_manual(values = c(cols.4[3], cols.4[1], cols.4[2], cols.4[4], "black"))+
               scale_color_manual(values = c(cols.4[3], cols.4[1], cols.4[2], cols.4[4], "black"))+
               scale_y_discrete(labels = c(expression(RMR),
                                           expression(MMR),
                                           expression(VM),
                                           expression(AAS)))+
               theme_classic()
             , y_title = "", x_title = "Scaling Slopes", size_text = 13)
p1<-p1+theme(plot.margin = margin(0,1,-0.5,-0.3, "cm"))


p2<-ggformat(plot=
               ggplot(fig7data.bottom, aes(x = as.numeric(slopes),  y = orderPerf, fill = color,  size = factor(size), shape = factor(size)))+
               geom_vline(xintercept = 0, color = "grey70", linetype = "dashed", size=0.5)+
               geom_errorbarh(mapping = aes(xmin = CI.L, xmax = CI.H),
                              color = "black", size = 0.6, height = 0.1)+
               xlim(-0.1, 0.1)+
               geom_point( color = "grey30", show.legend = F)+
               scale_size_manual(values = c(3, 4))+
               geom_text(size=4, color = "black", nudge_y = 0.3,
                         data =  fig7data.bottom[fig7data.bottom$size== "B",] , aes(x = slopes, label = slopes,  y = orderPerf))+
               scale_shape_manual(values = c(21, 2))+
               scale_fill_manual(values = c(cols.4[1], cols.4[2], cols.4[4], cols.4[3], "black", "red"))+
               # scale_color_manual(values = c(cols.4[3], cols.4[1], cols.4[2], cols.4[4], "black","red"))+
               scale_y_discrete(labels = c(expression(italic(f)[Hmax]),
                                           expression(T[AB]),
                                           expression(T[ARR]),
                                           expression(T[PEAK]),
                                           expression(FAS)))+
               theme_classic()
             , y_title = "", x_title = "Scaling Slopes", size_text = 13)
p2<-p2+theme(plot.margin = margin(0,1,1,-0.5, "cm"))

p3.q10<-
  ggplot(plot.final.mr, mapping = aes(x = treatm, y=q10, group=factor(p),shape = factor(p), 
                                      fill = treatm, label = p))+
  geom_line(alpha=1, show.legend = F,size=1,  color = "grey30")+
  geom_vline(show.legend = F, xintercept = 24, size=0.5, color = "red")+# geom_point(data = f.fas, size = 3,show.legend = F,
  geom_point(show.legend = F,  alpha = 1, size = 2)+
  scale_shape_manual(values = c(22, 21, 24, 23, 25))+
  scale_x_continuous(breaks = c(16, 20, 22, 24), labels = c("16", "20", "22", "24"))+
  scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5), limits = c(0, 2.7))+
  scale_color_gradient(low = "grey70", high = "grey35", name = "ºC")+
  scale_fill_gradient(low = "grey90", high = "grey35", name = "ºC")+
  annotate("segment",  x = 16, xend = 16, color=cols.4[2], y = 0, yend = 0.3, size=1.2, arrow = arrow(angle = 40, length = unit(0.25,"cm")))+
  annotate("segment",  x = 20, xend = 20, color=cols.4[3], y = 0, yend = 0.3, size=1.2, arrow = arrow(angle = 40, length = unit(0.25,"cm")))+
  annotate("segment",  x = 22, xend = 22, color=cols.4[4], y = 0, yend = 0.3, size=1.2, arrow = arrow(angle = 40, length = unit(0.25,"cm")))+
  # geom_errorbar(pred.means.rmr1kg, mapping = aes(ymax = exp(prediction+SE)*10, ymin = exp(prediction-SE)*10, group = factor(BW), x=treatm), width = 0.5, inherit.aes = F, color = "black", show.legend = F)+
  annotate("text", x = 15.5, y = 1.7, hjust = 1, label = expression(bold(MMR)), size = 4, parse = T)+
  annotate("text", x = 15.5, y = 1.1, hjust = 1, label = expression(bold(paste("AAS"))),  size = 4, parse = T)+
  annotate("text", x = 15.5, y = 0.75, hjust = 1, label = expression(bold(paste("FAS"))),  size = 4, parse = T)+
  annotate("text", x = 15.5, y = 2.28, hjust = 1, label = expression(bold(RMR)), size = 4, parse = T)+
  annotate("text", x = 16.6, y = 2, hjust = 1, label = expression(bold(italic(f)[Hmax])), size = 4, parse = T)+
  theme_classic()
ggformat(p3.q10, x_title = expression(Acute~treatments~(degree*C)),
         y_title = expression(Q[10]), size_text = 13, print = T)
# 


p3<-
  ggplot(plot.final.mr2, mapping = aes(x = as.numeric(as.character(treatm)), y=exp(prediction), group=factor(p), 
                                       fill = factor(treatm),
                                       color = factor(treatm),
                                       shape = factor(p)))+
  # geom_line(data = f.fas,
  #           mapping = aes(x = treatm, y=exp(prediction)-1, group=factor(p), fill = factor(treatm), shape = factor(p)),
  #           show.legend = F,size=1,  color = "grey70")+
  geom_line(alpha=1, show.legend = F,size=1,  color = "grey70")+
  geom_line(pred.means.bpm1kg,  inherit.aes = F,
            mapping = aes(y=exp(prediction)/25, x=as.numeric(as.character(treatm)), group= factor(BW)),
            alpha=1, show.legend = F, size=1,  color = "grey70")+
  geom_vline(xintercept = 24, size=0.5, color = "red")+# geom_point(data = f.fas, size = 2,show.legend = F,
  #            mapping = aes(x = treatm, y=exp(prediction)-1, group=factor(p), fill = factor(treatm), shape = factor(p)))+
  geom_point(show.legend = F, alpha = 1, size = 2)+
  scale_fill_manual(values = c(cols.4))+
  scale_color_manual(values = c(cols.4))+
  scale_shape_manual(values = c(22,  24, 23, 25))+
  scale_x_continuous(breaks = c( 12, 16, 20, 22, 24), limits= c(7, 28))+
  scale_y_continuous(sec.axis = sec_axis( trans=~.*25, name=expression(italic(f)[Hmax]~(beats~min^-1))))+
  # geom_errorbar(pred.means.bpm1kg, mapping = aes(ymax = exp(prediction+SE), ymin = exp(prediction-SE)/100, x=treatm),
  #               width = 0.5, inherit.aes = F, color = "black", position = position_dodge(width = 0.4), show.legend = F)+
  geom_point(pred.means.bpm1kg, mapping = aes(y=exp(prediction)/25, x=as.numeric(as.character(treatm))),
             show.legend = F, alpha = 1, size = 2, fill = "grey70", color = "grey70", pch=21, inherit.aes = F)+
  geom_point(pred.means.bpm1kg[pred.means.bpm1kg$treatm == "16",], mapping = aes(y=exp(prediction)/25, x=as.numeric(as.character(treatm))),
             show.legend = F, alpha = 1, size = 2, fill = cols.4[2], color = cols.4[2], pch=21, inherit.aes = F)+
  geom_point(pred.means.bpm1kg[pred.means.bpm1kg$treatm == "20",], mapping = aes(y=exp(prediction)/25, x=as.numeric(as.character(treatm))),
             show.legend = F, alpha = 1, size = 2, fill = cols.4[3], color = cols.4[3], pch=21, inherit.aes = F)+
  geom_point(pred.means.bpm1kg[pred.means.bpm1kg$treatm == "22",], mapping = aes(y=exp(prediction)/25, x=as.numeric(as.character(treatm))),
             show.legend = F, alpha = 1, size = 2, fill = cols.4[4], color = cols.4[4], pch=21, inherit.aes = F)+
  geom_point(pred.means.bpm1kg[pred.means.bpm1kg$treatm == "24",], mapping = aes(y=exp(prediction)/25, x=as.numeric(as.character(treatm))),
             show.legend = F, alpha = 1, size = 2, fill = "red", color = "red", pch=21, inherit.aes = F)+
  # geom_errorbar(pred.means.rmr1kg, mapping = aes(ymax = exp(prediction+SE)*10, ymin = exp(prediction-SE)*10, group = factor(BW), x=treatm), width = 0.5, inherit.aes = F, color = "black", show.legend = F)+
  annotate("text", x = 11, y = 2.8, hjust = 1, label = expression(bold(MMR)), size = 4, parse = T)+
  annotate("text", x = 11, y = 1.65, hjust = 1, label = expression(bold(paste("AAS"))),  size = 4, parse = T)+
  annotate("text", x = 11, y = 2.2, hjust = 1, label = expression(bold(paste("FAS"))),  size = 4, parse = T)+
  annotate("text", x = 11, y = 1.2, hjust = 1, label = expression(bold(RMR)), size = 4, parse = T)+
  annotate("text", x = 21, y = 3, hjust = 1, label = expression(bold(italic(f)[Hmax])), size = 4, parse = T)+
  theme_classic()
ggformat(p3, x_title = expression(Temperature~(degree*C)),
         y_title = expression(MR~(mgO[2]~min^-1~kg^-1)), size_text = 13, print = T)


## Figure S1 ----------

supl.mmr.plot<-gather(data[, c("FishID", "treatm", "mmr", "mmrChase")], 
                      performance, value, mmr:mmrChase, factor_key=TRUE)


MMRcomp1<-ggformat( plot =
            (ggplot(data=supl.mmr.plot, aes(x=treatm, y=value, fill = performance))+
               geom_boxplot(alpha = 0.3, outlier.colour = NULL, show.legend = F)+
               geom_point(pch=21, position = position_dodge(width = 0.7), show.legend = F)+
               scale_color_manual(values = c("dodgerblue", "orange"))+
               scale_fill_manual(values = c("dodgerblue", "orange"))+
               theme_classic()+
               geom_point(mapping = aes(y = 3.1, x = 0.8), pch=19, size=3, color = "dodgerblue", show.legend = F)+
               geom_point(mapping = aes(y = 2.9, x = 0.8), pch=19, size=3, color = "orange", show.legend = F)+
               annotate(geom = "text", x = 1, y = 3.1, label = "MMR", size = 3.5, hjust =0)+
               annotate(geom = "text", x = 1, y = 2.9, label = expression(MMR[CHASE]), size = 3.5, hjust =0)),
          y_title = expression(MMR~(mgO[2]~min^-1)), x_title = "Treatment" ,
          title = "", print = F )

temp_names <- c(
  `12` = "12ºC",
  `16` = "16ºC",
  `20` = "20ºC",
  `22` = "22ºC"
)

MMRcomp2<-ggformat( plot =
            (ggplot(data=data, aes(x=mmr, y=mmrChase, fill = treatm, group = treatm))+
               geom_point(pch=21, alpha = 1, position = position_dodge(width = 0.5), size=2, show.legend = F)+
               scale_color_manual(values = cols.4)+
               scale_fill_manual(values = cols.4)+
               theme_classic()+
               geom_abline(slope = 1, intercept = 0, color = "black", lty = 2)+
               facet_wrap(.~treatm, nrow = 2, labeller = as_labeller(temp_names))),
          y_title = expression(MMR~(mgO[2]~min^-1)), x_title = expression(MMR[CHASE]~(mgO[2]~min^-1)),
          title = "", print =F)


## Figure S2 ------

temparrh2S<-ggformat(
  plot = (ggplot(data.abtID, aes(y = (temp_ARRH), x = (BW), fill= origin))+
            geom_point( color = "black", show.legend = F, pch = 21, size = 2)+
            scale_fill_manual(values = c("black", "white"))+
            ylim(17, 30)+
            xlim(0, 0.25)
  ),
  x_title = expression((Body~mass*~kg)),
  y_title = expression((T[ARR]*~degree*C)), size_text = 15, print = F)

tpeak2S<-ggformat(
  plot = (ggplot(data.abtID, aes(y = (Tpeak), x = (BW), fill= origin))+
            geom_point( color = "black", show.legend = F, pch = 21, size = 2)+
            ylim(17, 30)+
            xlim(0, 0.25)+
            scale_fill_manual(values = c("black", "black"))
  ),
  x_title = expression((Body~mass*~kg)),
  y_title = expression((T[PEAK]*~degree*C)), size_text = 15, print = F)

tabt2S<-ggformat(
  plot = (ggplot(data.abtID, aes(y = (breakpoint_Cels), x = (BW), fill= origin))+
            geom_point( color = "black", show.legend = F, pch = 21, size = 2)+
            ylim(17, 30)+
            xlim(0, 0.25)+
            scale_fill_manual(values = c("black", "white"))
  ),
  x_title = expression((Body~mass*~kg)),
  y_title = expression((T[AB]*~degree*C)), size_text = 15, print = F)

hrmax2S<-ggformat(
  plot = (ggplot(data.abtID, aes(y = (HRpeak), x = (BW), fill= origin))+
            geom_point( color = "black", show.legend = F, pch = 21, size = 2)+
            # ylim(21, 30)+
            xlim(0, 0.25)+
            scale_fill_manual(values = c("black", "white"))
  ),
  x_title = expression((Body~mass~kg)),
  y_title = expression((PEAK[italic(f)*Hmax]~(beats~min^-1))))


## Figure S3 -------------
plotheart1<-ggformat(
  plot = (ggplot(data.abt, aes(x = log((mass.g/1000)), y = log(heart.mass/1000)))+
            geom_point(pch=21, color = "black", fill = "grey", size = 3)+
            ylim(-12.5, -8)+
            xlim(-6, -1)+
            geom_smooth(method="lm", se=FALSE, color= "black")+
            annotate(geom = "text",x = -5.8, y = -8.3, hjust = 0,
                     label = deparse(bquote(~italic(ln)~VM == italic(ln) ~"("* .(round(as.numeric(heart.a),2))*")"~ BM^.(round(as.numeric(heart.b),2)))),
                     parse = T, color = "black")),
  y_title = expression(italic(ln)~Ventricle~mass~(kg)),
  x_title = expression(italic(ln)~Body~mass~(kg)), size_text = 15)





# 9. SAVING FIGURES ------
## Save Fig 3 AB --------
cowplot::plot_grid(thesis2A, thesis2B, align = "hv",
                   nrow = 1, ncol = 2,
                   labels = c("A", "B"), 
                   label_x = c(0.16, 0.16), 
                   label_y = c(0.9, 0.9)) %>% 
ggsave(filename = paste( "./Figures/Fig2_AB",Sys.Date(),".png",sep=""), width = 8, height = 4)

## Save Fig 3 CDEF --------
grid.arrange(grobs = list(thesis2C, thesis2D,thesis2E, thesis2F),
             ncol = 4, nrow =1, widths=c(1.08,1,1,1),
             bottom = textGrob(expression(italic(ln)~Body~mass~(kg)),rot = 0, gp = gpar(fontsize = 18)),
             left = textGrob(expression(italic(ln)~MR),rot = 90, gp = gpar(fontsize = 18))) %>% 
ggsave(filename = paste( "./Figures/Fig2_CDEF",Sys.Date(),".png",sep=""), width = 9.5, height = 3)

## Save Fig 4 ABCD ---------
cowplot:::plot_grid(plot.size, plot.bpm.lmer, plot2, plot.abt.scale,
                    nrow = 2, 
                    ncol =2, labels = "AUTO", align = "vh",
                    label_x = c(0.21, 0.21),
                    label_y = c(0.9, 0.9), label_size = 14) %>%
ggsave(filename = paste( "./Figures/Fig3_ABCD",Sys.Date(),".png",sep=""),
         width = 8, height = 7.8)

## Save Fig 5 ABCD ---------
cowplot::plot_grid(tpeak2, temparrh2, tabt2, hrmax2, 
                   align = "hv", label_x = c(0.2, 0.2), label_y = c(0.9, 0.9),labels = "AUTO" ) %>% 
ggsave(filename = paste( "./Figures/Fig6_ABCD",Sys.Date(),".png",sep=""), width = 8.2, height = 8.2)

## Save Fig 6 ABCD ---------
cowplot:::plot_grid(plot.mmr.lmer, plot.rmr.lmer, plot.AS.lmer, plot.fas.lmer,
                    nrow = 2, 
                    ncol =2, labels = "AUTO", align = "vh",
                    label_x = c(0.25, 0.25),
                    label_y = c(0.9, 0.9), label_size = 14) %>%
ggsave(filename = paste( "./Figures/Fig4_ABCD",Sys.Date(),".png",sep=""),
         width = 8.5, height = 7.8)

## Save Fig 7 ABCD -------
cowplot:::plot_grid(p1, p2, p3, p3.q10, nrow =2, ncol=2, labels = "AUTO", align = "hv",label_x = c(0.2, 0.2), label_y = c(0.9, 0.9)) %>% 
  ggsave(filename = paste( "./Figures/Fig7_FINAL",Sys.Date(),".png",sep=""), width = 10, height = 8)


## Save Fig S3 -----

ggsave(plot =plotheart1, filename = paste( "./Figures/Fig5_",Sys.Date(),"_rvmscale", Sys.Date(),".png",sep=""),
       height = 3.5, width = 3.5)

## Save Fig S2 --------
cowplot::plot_grid(tpeak2S, temparrh2S, tabt2S, hrmax2S, 
                   align = "hv", label_x = c(0.2, 0.2), label_y = c(0.9, 0.9),labels = "AUTO" ) %>% 
ggsave(filename = paste( "./Figures/Fig6SUPPL_ABCD",Sys.Date(),".png",sep=""), width = 8.2, height = 8.2)

## Save Fig S1 ----------
cowplot::plot_grid(MMRcomp1, MMRcomp2, align = "hv",
                   nrow = 1, ncol = 2,
                   labels = c("A", "B"), 
                   label_y = c(0.95, 0.95)) %>% 
ggsave(filename = paste( "./Figures/FigS1_MMRs",Sys.Date(),".png",sep=""), width = 6, height = 3)

# 10. Summaries ----

# mean diff between TARR and TPEAK
summary(data.abtID$temp_ARRH - data.abtID$Tpeak, na.rm = T)


# cardiac thermal tolerance indices
means.heartID<-data.abtID %>% 
  dplyr:::summarise(min.PEAKfhmax = min(HRpeak, na.rm = T),
                    max.PEAKfhmax = max(HRpeak, na.rm = T),
                    mean.PEAKfhmax = mean(HRpeak, na.rm = T),
                    sd.PEAKfhmax = sd(HRpeak, na.rm = T),
                    n.PEAKfhmax = length(HRpeak)- sum(is.na(HRpeak)),
                    
                    min.breakpoint_Cels = min(breakpoint_Cels, na.rm = T),
                    max.breakpoint_Cels = max(breakpoint_Cels, na.rm = T),
                    mean.breakpoint_Cels = mean(breakpoint_Cels, na.rm = T),
                    sd.breakpoint_Cels = sd(breakpoint_Cels, na.rm = T),
                    n.breakpoint_Cels = length(breakpoint_Cels)- sum(is.na(breakpoint_Cels)),
                    
                    min.Tpeak = min(Tpeak, na.rm = T),
                    max.Tpeak = max(Tpeak, na.rm = T),
                    mean.Tpeak = mean(Tpeak, na.rm = T),
                    sd.Tpeak = sd(Tpeak, na.rm = T),
                    n.Tpeak = length(Tpeak)- sum(is.na(Tpeak)),
                    
                    min.temp_ARRH = min(temp_ARRH, na.rm = T),
                    max.temp_ARRH = max(temp_ARRH, na.rm = T),
                    mean.temp_ARRH = mean(temp_ARRH, na.rm = T),
                    sd.temp_ARRH = sd(temp_ARRH, na.rm = T),
                    n.temp_ARRH = length(temp_ARRH)- sum(is.na(temp_ARRH)),
                    
                    min.BW = min(BW, na.rm = T),
                    max.BW = max(BW, na.rm = T),
                    mean.BW = mean(BW, na.rm = T),
                    sd.BW = sd(BW, na.rm = T),
                    n.BW = length(BW)- sum(is.na(BW)))


means.heartID.SizeClass<-data.abtID %>% 
  group_by(sizeClass) %>% 
  dplyr:::summarise(min.PEAKfhmax = min(HRpeak, na.rm = T),
                    max.PEAKfhmax = max(HRpeak, na.rm = T),
                    mean.PEAKfhmax = mean(HRpeak, na.rm = T),
                    sd.PEAKfhmax = sd(HRpeak, na.rm = T),
                    n.PEAKfhmax = length(HRpeak)- sum(is.na(HRpeak)),
                    
                    min.breakpoint_Cels = min(breakpoint_Cels, na.rm = T),
                    max.breakpoint_Cels = max(breakpoint_Cels, na.rm = T),
                    mean.breakpoint_Cels = mean(breakpoint_Cels, na.rm = T),
                    sd.breakpoint_Cels = sd(breakpoint_Cels, na.rm = T),
                    n.breakpoint_Cels = length(breakpoint_Cels)- sum(is.na(breakpoint_Cels)),
                    
                    min.Tpeak = min(Tpeak, na.rm = T),
                    max.Tpeak = max(Tpeak, na.rm = T),
                    mean.Tpeak = mean(Tpeak, na.rm = T),
                    sd.Tpeak = sd(Tpeak, na.rm = T),
                    n.Tpeak = length(Tpeak)- sum(is.na(Tpeak)),
                    
                    min.temp_ARRH = min(temp_ARRH, na.rm = T),
                    max.temp_ARRH = max(temp_ARRH, na.rm = T),
                    mean.temp_ARRH = mean(temp_ARRH, na.rm = T),
                    sd.temp_ARRH = sd(temp_ARRH, na.rm = T),
                    n.temp_ARRH = length(temp_ARRH)- sum(is.na(temp_ARRH)),
                    
                    min.BW = min(BW, na.rm = T),
                    max.BW = max(BW, na.rm = T),
                    mean.BW = mean(BW, na.rm = T),
                    sd.BW = sd(BW, na.rm = T),
                    n.BW = length(BW)- sum(is.na(BW)))

mean(data.abtID[data.abtID$sizeClass == "YOY", "breakpoint_Cels"], na.rm = T) - mean(data.abtID[data.abtID$sizeClass == "> YOY", "breakpoint_Cels"], na.rm =T)
mean(data.abtID[data.abtID$sizeClass == "YOY", "Tpeak"]) - mean(data.abtID[data.abtID$sizeClass == "> YOY", "Tpeak"], na.rm =T)
mean(data.abtID[data.abtID$sizeClass == "YOY", "temp_ARRH"]) - mean(data.abtID[data.abtID$sizeClass == "> YOY", "temp_ARRH"], na.rm =T)

