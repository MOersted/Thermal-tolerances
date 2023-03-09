##########
#Title: A unifying model to estimate thermal tolerance limits in ectotherms across static, dynamic and fluctuating exposures to thermal stress 
#Authors: Lisa Bjerregaard Jørgensen, Hans Malte, Michael Ørsted, Nikolaj Andreasen Klahn & Johannes Overgaard
#Date: 20 April 2021
#Version: 1.1
##########

# This script derives TDT parameters from static experiments, where time to failure (tcoma) is measured at one or more constant temperatures. An input data template is provided (static_input.csv)

# Install packages not yet installed
packages <- c("rootSolve")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

########## INPUT STATIC DATA ##########
static = read.csv("static_input.csv")

########## SELECT TYPES OF OUTPUT DATA ##########

# 1.1	Output : Tolerable temperature at a given exposure time
# duration (t_coma, in minutes) for which sCTmax should be predicted. Note 1 day is 1440 minutes
# extra_t_coma = c(1*1440,3*1440,7*1440) #example e.g. extrapolate to 1, 3, and 7 days
extra_t_coma = c(5,60,120) #example

# 1.2	Output: Knockdown time at a given temperature
# temperature (sCTmax, in degC) for which t_coma should be predicted
# extra_sCTmax = c(35,42) #example

# 1.3	Output: Dynamic CTmax in ramping assays
#Input the ramprates (in degC/min) for which dCTmax should be predicted
ramprates = c(0.05,0.10,0.25)

# 1.4	Output: Knockdown time in dynamic fluctuating temperatures
#If t_coma under a fluctuating temperature profile is wanted, input the temperature profile
fluc_temp = read.csv("fluctuating_temperature_profile.csv",as.is = T)

#loop over group
groups = unique(as.character(static$group))

if(exists("extra_t_coma")){sCTmax_extra_t_coma = data.frame()}
if(exists("extra_sCTmax")){sCTmax_extra_sCTmax = data.frame()}
if(exists("ramprates")){dCTmax.predictions = data.frame()}
if(exists("fluc_temp")){fluctemp.predictions = data.frame(group = groups,t_coma=NA)}

#Parameters
# z: temperature sensitivity coefficient [dimesionless]
# sCTmax: static CTmax [deg Celsius]
# tLs: static CTmaxtime [min]
# To: Ramp start temp [deg Celsius]
# Tc: temperature where damage accumulation starts [deg Celsius]
Tc = 19 #example
To = 19 #example

pdf(file="TDT_from_static.pdf",width = 12,height = 8)
for(i in 1:length(groups))
  {
  tmp = static[static$group == groups[i],]
 
  #if extra sCTmax values are desired at other exposure durations (extra_t_coma), calculate these from the TDT
  if(exists("extra_t_coma"))
    {
    if(nrow(tmp)>1)
      {
      #In the case of two or more static temperatures the parameters of the TDT can be extracted from a linear regression
      fit = lm(log10(t_coma)~assay_temp,data=tmp)
      intercept = summary(fit)$coefficient[1,1]
      slope = summary(fit)$coefficient[2,1]
      z = -1/slope
      
      sCTmaxtmp = data.frame(group=rep(groups[i],times=length(extra_t_coma)),t_coma=extra_t_coma,assay_temp = NA,z_pred=z)
      for(k in 1:length(extra_t_coma))
        {
        sCTmaxtmp$assay_temp[k] = (log10(extra_t_coma[k])-intercept)/slope  
        }
      #plot TDT to assess goodness of fit
      if(!exists("extra_sCTmax")){
      plot(tmp$assay_temp,log10(tmp$t_coma),type='n',
           ylim=c(0,max(c(log10(extra_t_coma),log10(tmp$t_coma)))),
           xlim=c(min(c(sCTmaxtmp$assay_temp,static$assay_temp)),max(c(sCTmaxtmp$assay_temp,static$assay_temp))),
           xlab="Temperature (°C)",
           ylab=expression('log'[10]*' knockdown time (min)'),
           main=paste0(unique(tmp$group),", z = ",signif(z,3))
          )
      abline(fit,lwd=2,col="grey")
      points(tmp$assay_temp,log10(tmp$t_coma),pch=21,bg="black",cex=1.2)
      text(x=max(c(sCTmaxtmp$assay_temp,static$assay_temp))*0.99,y=max(c(log10(extra_t_coma),log10(tmp$t_coma)))*0.95,bquote(paste('R'^'2'*' = ',.(signif(summary(fit)$r.squared,digits=2)))),cex=1.5)
      points(sCTmaxtmp$assay_temp,log10(sCTmaxtmp$t_coma),pch=21,bg="white",cex=1.2)
      legend("bottomleft",legend = c("Observed","Predicted","TDT fit"),pch=c(21,21,NA),pt.bg = c("black","white",NA),lwd = c(NA,NA,2),col=c("black","black","grey"),bty='n')
      }
      tmp$z_pred = z
      }else{
      #in the case where t_coma has been measured at only one static temperature, one must supply an estimated z value (see main text)
      z = 2.5
      slope = -1/z
      intercept = -slope*tmp$assay_temp+log10(tmp$t_coma)
      
      sCTmaxtmp = data.frame(group=rep(groups[i],times=length(extra_t_coma)),t_coma=extra_t_coma,assay_temp = NA,z_pred=z)
      for(k in 1:length(extra_t_coma))
        {
        sCTmaxtmp$assay_temp[k] = (log10(extra_t_coma[k])-intercept)/slope  
        }
      #plot TDT (no fit possible with only one measurement)
      plot(tmp$assay_temp,log10(tmp$t_coma),type='n',
           ylim=c(0,max(c(log10(sCTmaxtmp$t_coma),log10(static$t_coma)))),
           xlim=c(min(c(sCTmaxtmp$assay_temp,static$assay_temp)),max(c(sCTmaxtmp$assay_temp,static$assay_temp))),
           xlab="Temperature (°C)",
           ylab=expression('log'[10]*' knockdown time (min)'),
           main=paste0(unique(tmp$group),", z = ",signif(z,3)," (estimated)")
          )
      abline(intercept,slope,lwd=2,col="grey")
      points(tmp$assay_temp,log10(tmp$t_coma),pch=21,bg="black",cex=1.2)
      text(x=max(c(sCTmaxtmp$assay_temp,static$assay_temp))*0.99,y=max(c(log10(extra_t_coma),log10(tmp$t_coma)))*0.95,bquote('R'^'2'*' = NA'),cex=1.5)
      points(sCTmaxtmp$assay_temp,log10(sCTmaxtmp$t_coma),pch=21,bg="white",cex=1.2)
      legend("bottomleft",legend = c("Observed","Predicted","TDT fit"),pch=c(21,21,NA),pt.bg = c("black","white",NA),lwd = c(NA,NA,2),col=c("black","black","grey"),bty='n')
      tmp$z_pred = z
      }
    sCTmax_extra_t_coma = rbind(sCTmax_extra_t_coma,tmp,sCTmaxtmp)
    }
  
  #if extra t_coma values are desired at other exposure temperatures (extra_sCTmax), calculate these from the TDT
  if(exists("extra_sCTmax"))
    {
    if(nrow(tmp)>1)
      {
      #In the case of two or more static temperatures the parameters of the TDT can be extracted from a linear regression
      fit = lm(log10(t_coma)~assay_temp,data=tmp)
      intercept = summary(fit)$coefficient[1,1]
      slope = summary(fit)$coefficient[2,1]
      z = -1/slope
      
      sCTmaxtmp = data.frame(group=rep(groups[i],times=length(extra_sCTmax)),t_coma=NA,assay_temp = extra_sCTmax,z_pred=z)
      for(k in 1:length(extra_sCTmax))
        {
        sCTmaxtmp$t_coma[k] = 10^(intercept+(slope*extra_sCTmax[k]))
        }
      #plot TDT to assess goodness of fit
      plot(tmp$assay_temp,log10(tmp$t_coma),type='n',
           ylim=c(min(c(log10(sCTmaxtmp$t_coma),log10(static$t_coma))),max(c(log10(sCTmaxtmp$t_coma),log10(static$t_coma)))),
           xlim=c(min(c(sCTmaxtmp$assay_temp,static$assay_temp)),max(c(sCTmaxtmp$assay_temp,static$assay_temp))),
           xlab="Temperature (°C)",
           ylab=expression('log'[10]*' knockdown time (min)'),
           main=paste0(unique(tmp$group),", z = ",signif(z,3))
          )
      abline(fit,lwd=2,col="grey")
      points(tmp$assay_temp,log10(tmp$t_coma),pch=21,bg="black",cex=1.2)
      text(x=max(c(sCTmaxtmp$assay_temp,static$assay_temp))*0.99,y=max(c(log10(sCTmaxtmp$t_coma),log10(static$t_coma)))*0.95,bquote(paste('R'^'2'*' = ',.(signif(summary(fit)$r.squared,digits=2)))),cex=1.5)
      points(sCTmaxtmp$assay_temp,log10(sCTmaxtmp$t_coma),pch=21,bg="white",cex=1.2)
      legend("bottomleft",legend = c("Observed","Predicted","TDT fit"),pch=c(21,21,NA),pt.bg = c("black","white",NA),lwd = c(NA,NA,2),col=c("black","black","grey"),bty='n')
      tmp$z_pred = z
      }else{
      #in the case where t_coma has been measured at only one static temperature, one must supply an estimated z value (see main text)
      z = 2.5
      slope = -1/z
      intercept = -slope*tmp$assay_temp+log10(tmp$t_coma)
      
      sCTmaxtmp = data.frame(group=rep(groups[i],times=length(extra_sCTmax)),t_coma=NA,assay_temp = extra_sCTmax,z_pred=z)
      for(k in 1:length(extra_sCTmax))
        {
        sCTmaxtmp$t_coma[k] = 10^(intercept+(slope*extra_sCTmax[k]))
        }
      #plot TDT (no fit possible with only one measurement)
      plot(tmp$assay_temp,log10(tmp$t_coma),type='n',
           ylim=c(min(c(log10(sCTmaxtmp$t_coma),log10(static$t_coma))),max(c(log10(sCTmaxtmp$t_coma),log10(static$t_coma)))),
           xlim=c(min(c(sCTmaxtmp$assay_temp,static$assay_temp)),max(c(sCTmaxtmp$assay_temp,static$assay_temp))),
           xlab="Temperature (°C)",
           ylab=expression('log'[10]*' knockdown time (min)'),
           main=paste0(unique(tmp$group),", z = ",signif(z,3))
      )
      abline(intercept,slope,lwd=2,col="grey")
      points(tmp$assay_temp,log10(tmp$t_coma),pch=21,bg="black",cex=1.2)
      text(x=max(c(sCTmaxtmp$assay_temp,static$assay_temp))*0.99,y=max(c(log10(sCTmaxtmp$t_coma),log10(static$t_coma)))*0.95,bquote('R'^'2'*' = NA'),cex=1.5)
      points(sCTmaxtmp$assay_temp,log10(sCTmaxtmp$t_coma),pch=21,bg="white",cex=1.2)
      legend("bottomleft",legend = c("Observed","Predicted","TDT fit"),pch=c(21,21,NA),pt.bg = c("black","white",NA),lwd = c(NA,NA,2),col=c("black","black","grey"),bty='n')
      tmp$z_pred = z
      }
    sCTmax_extra_sCTmax = rbind(sCTmax_extra_sCTmax,tmp,sCTmaxtmp)
    }
  
  if(exists("ramprates"))
    {
    #Predict dCTmax from sCTmax
    if(nrow(tmp)>1)
      {
      #In the case of two or more static temperatures dCTmax can be predicted from the parameters of the TDT
      fit = lm(log10(t_coma)~assay_temp,data=tmp)
      intercept = summary(fit)$coefficient[1,1]
      slope = summary(fit)$coefficient[2,1]
      z = -1/slope
      
      #plot TDT to assess goodness of fit
      if(!any(exists("extra_t_coma") | exists("extra_sCTmax"))){
        plot(tmp$assay_temp,log10(tmp$t_coma),type='n',
             ylim=c(0,max(log10(static$t_coma))),
             xlim=c(min(static$assay_temp),max(static$assay_temp)),
             xlab="Temperature (°C)",
             ylab=expression('log'[10]*' knockdown time (min)'),
             main=paste0(unique(tmp$group),", z = ",signif(z,3))
        )
        abline(fit,lwd=2,col="grey")
        points(tmp$assay_temp,log10(tmp$t_coma),pch=21,bg="black",cex=1.2)
        text(x=max(static$assay_temp)*0.99,y=max(log10(static$t_coma))*0.95,bquote(paste('R'^'2'*' = ',.(signif(summary(fit)$r.squared,digits=2)))),cex=1.5)
      }
      dCTmax.predictionstmp = data.frame("group" = rep(unique(tmp$group),times=length(ramprates)),ramprate=rep(ramprates,times=length(unique(tmp$group))),dCTmax_pred = NA,z_pred=z)
      tLs = max(tmp$t_coma)
      #for each desired ramprate, predict dCTmax
      for(j in 1:length(ramprates))
        {
        sCTmax = tmp$assay_temp[tmp$t_coma == tLs] #sCTmax = (log10(1*1440)-intercept)/slope   #tmp$t_coma == tLs
        dCTmax = To+(z/log(10))*log((log(10)*ramprates[j]*tLs/z)*exp((log(10)/z)*(sCTmax-To))+exp((log(10)/z)*(Tc-To)))
        dCTmax.predictionstmp$dCTmax_pred[j] = dCTmax
        }
      }else{
      #in the case where t_coma has been measured at only one static temperature, you must supply an estimated z value (see main text)
      z = 2.5
      slope = -1/z
      intercept = -slope*tmp$assay_temp+log10(tmp$t_coma)
      
      #plot TDT (no fit possible with only one measurement)
      if(!any(exists("extra_t_coma") | exists("extra_sCTmax"))){
        plot(tmp$assay_temp,log10(tmp$t_coma),type='n',
             ylim=c(0,max(log10(static$t_coma))),
             xlim=c(min(static$assay_temp),max(static$assay_temp)),
             xlab="Temperature (°C)",
             ylab=expression('log'[10]*' knockdown time (min)'),
             main=paste0(unique(tmp$group),", z = ",signif(z,3)," (estimated)")
        )
        abline(intercept,slope,lwd=2,col="grey")
        points(tmp$assay_temp,log10(tmp$t_coma),pch=21,bg="black",cex=1.2)
        text(x=max(static$assay_temp)*0.99,y=max(log10(static$t_coma))*0.95,bquote('R'^'2'*' = NA',),cex=1.5)
      }
      dCTmax.predictionstmp = data.frame("group" = rep(unique(tmp$group),times=length(ramprates)),ramprate=rep(ramprates,times=length(unique(tmp$group))),dCTmax_pred = NA,z_pred=z)
      tLs = max(tmp$t_coma)
      for(j in 1:length(ramprates))
      {
        sCTmax = tmp$assay_temp[tmp$t_coma == tLs] #sCTmax = (log10(1*1440)-intercept)/slope   #tmp$t_coma == tLs
        dCTmax = To+(z/log(10))*log((log(10)*ramprates[j]*tLs/z)*exp((log(10)/z)*(sCTmax-To))+exp((log(10)/z)*(Tc-To)))
        dCTmax.predictionstmp$dCTmax_pred[j] = dCTmax
      }
    }
    dCTmax.predictions = rbind(dCTmax.predictions,dCTmax.predictionstmp)
  }
  
  if(exists("fluc_temp"))
    {
    #The TDT is built from static temperatures and t_coma predicted from dynamic data. Additive injury is determined per time interval (in minutes) 
    fluc_temp$add_inj = NA
    for(time in 1:nrow(fluc_temp))
      {
      fluc_temp$add_inj[time] = ifelse(any(fluc_temp$temperature[time+1]>=Tc,fluc_temp$temperature[time]>=Tc),(100*(fluc_temp$time_min[time+1]-fluc_temp$time_min[time]))/(10^(slope*max(fluc_temp$temperature[time+1],fluc_temp$temperature[time])+intercept)),0)
      }
    fluc_temp$cum_inj = cumsum(fluc_temp$add_inj)
    
    #Select the percentage accumulated thermal injury for which tolerable duration should be predicted
    acc_injury = 100
    fluctemp.predictions$t_coma[i] = fluc_temp$time_min[which.min(abs(fluc_temp$cum_inj - acc_injury))]  
    }
  print(paste("Group",i,"out of",length(groups)))
}
dev.off()

if(exists("ramprates")){write.csv(dCTmax.predictions,file="dCTmax_predictions.csv")}
if(exists("extra_t_coma")){write.csv(sCTmax_extra_t_coma,file="extra_t_coma.csv",row.names = FALSE)}
if(exists("extra_sCTmax")){write.csv(sCTmax_extra_sCTmax,file="extra_sCTmax.csv",row.names = FALSE)}
if(exists("fluc_temp")){write.csv(fluctemp.predictions,file="fluctemp_predictions.csv",row.names = FALSE)}
