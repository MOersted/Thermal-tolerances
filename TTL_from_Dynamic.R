##########
#Title: A unifying model to estimate thermal tolerance limits in ectotherms across static, dynamic and fluctuating exposures to thermal stress 
#Authors: Lisa Bjerregaard Jørgensen, Hans Malte, Michael Ørsted, Nikolaj Andreasen Klahn & Johannes Overgaard
#Date: 5 March 2021
#Version: 1.0
##########

# This script derives TTL parameters from dynamic experiments, where the maximal temperatures tolerated (dynamic CTmax, dCTmax) are measured using one or more ramping rates. An input data template is provided (dynamic_input.csv). 

# Install packages not yet installed
packages <- c("rootSolve")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

########## INPUT DYNAMIC DATA ##########
dynamic = read.csv("dynamic_input.csv")

########## SELECT TYPES OF OUTPUT DATA ##########

# 2.1	Output : Tolerable temperature at a given exposure time
###IMPORTANT: this is required and not optional, as the other output data types are calculated from estimated sCTmax at these given t_coma exposure times
#Input the t_coma times (in minutes) for which sCTmax should be predicted
t_coma = c(1,10,60) #example

# 2.2	Output: Knockdown time at a given temperature 
#If additional/specific temperatures for which knockdown time (t_coma) should be predicted input them here (in °C)
extra_sCTmax = c(39,40,41) #example

# 2.3	Output: Dynamic CTmax in ramping assays
#Input other ramprates (in °C/min) for which dCTmax should be predicted
extra_ramprates = c(0.2,0.15,0.025) #example

# 2.4	Output: Knockdown time in dynamic fluctuating temperatures
#If t_coma under a fluctuating temperature profile is wanted, input the temperature profile
fluc_temp = read.csv("fluctuating_temperature_profile.csv",as.is = T)

#loop over group
groups = unique(as.character(dynamic$group))
sCTmax.predictions = NULL
if(exists("extra_sCTmax")){sCTmax.predictions_extra = NULL}
if(exists("extra_ramprates")){dCTmax.predictions = data.frame()}
if(exists("fluc_temp")){fluctemp.predictions = data.frame(group = groups,t_coma=NA)}

if(!exists("t_coma")){stop("t_coma times (in minutes) for which sCTmax should be predicted must be supplied")}

pdf(file="TTL_from_dynamic.pdf",width = 12,height = 8)
for(i in 1:length(groups))
  {
  dtmp = dynamic[dynamic$group == groups[i],]
  sCTmax.predictionstmp = data.frame(group = rep(groups[i],times=length(t_coma)),t_coma=t_coma,sCTmax_pred=NA,z_pred=NA)  
    
  #Parameters
  # z: temperature sensitivity coefficient [dimesionless]
  # sCTmax: static CTmax [deg Celsius]
  # tLs: static CTmaxtime [min]
  # To: Ramp start temp [deg Celsius]
  # Tc: temperature where damage accumulation starts [deg Celsius]
  
  Tc = 19 #example. This can be group specific
  To = 19 #example. This can be group specific
  
  ramprates = unique(dtmp$ramprate)
  if(!exists("extra_ramprates")){plot(dtmp$ramprate,dtmp$dCTmax,xlim=c(0,max(dynamic$ramprate)+0.1),ylim=c(min(dynamic$dCTmax),max(dynamic$dCTmax)),xlab="Ramp rate (°C/min)",ylab="dCTmax (°C)",main=groups[i],type="n")}
  # xpred = seq(from=min(dtmp$ramprate)*0.2,to=max(dtmp$ramprate)*1.2,length.out=500)
  for(j in 1:length(t_coma))
    {
    tLs =  t_coma[j]
    if(length(ramprates)>2)
      {
      ###In the case there are three or more ramp rates, we can fit a non-linear model, as we have >2 observations to solve 2 unknowns
      #for each desired t_coma, predict sCTmax (in °C)
      dCTmaxmod = nls(dCTmax ~ To+(z/log(10))*log((log(10)*ramprate*tLs/z)*exp((log(10)/z)*(sCTmax-To))+exp((log(10)/z)*(Tc-To))),start = list(sCTmax=35,z=2.5),data = dtmp)
      sCTmax.predictionstmp$z_pred[j] = as.numeric(dCTmaxmod$m$getPars()[2])
      sCTmax.predictionstmp$sCTmax_pred[j] = as.numeric(dCTmaxmod$m$getPars()[1])
      
      eq = function(x){To+(as.numeric(dCTmaxmod$m$getPars()[2])/log(10))*log((log(10)*x*tLs/as.numeric(dCTmaxmod$m$getPars()[2]))*exp((log(10)/as.numeric(dCTmaxmod$m$getPars()[2]))*(as.numeric(dCTmaxmod$m$getPars()[1])-To))+exp((log(10)/as.numeric(dCTmaxmod$m$getPars()[2]))*(Tc-To)))}
      if(j==length(t_coma) & !exists("extra_ramprates")){curve(eq, from=min(dtmp$ramprate)*0.2,to=max(dtmp$ramprate)*1.2, add=TRUE,lwd=2,col = "grey")}
      }
    if(length(ramprates)==2)
      {
      ####In the case there are only two ramp rates, we have to solve 2 equations with 2 unknowns; z and sCTmax (x[1] and x[2], respectively, in the below function)
      fn <- function(x)(c(F1=(x[1]/log(10))*log((log(10)*dtmp$ramprate[1]*tLs/x[1])*exp((log(10)/x[1])*(x[2]-To))+exp((log(10)/x[1])*(Tc-To)))+To-dtmp$dCTmax[1],
                          F2=(x[1]/log(10))*log((log(10)*dtmp$ramprate[2]*tLs/x[1])*exp((log(10)/x[1])*(x[2]-To))+exp((log(10)/x[1])*(Tc-To)))+To-dtmp$dCTmax[2]))
      roots = multiroot(fn,start = c(2.5,35),maxiter = 1000,positive=TRUE)
      sCTmax.predictionstmp$z_pred[j] = roots[[1]][1]
      sCTmax.predictionstmp$sCTmax_pred[j] = roots[[1]][2]
      eq = function(x){To+(roots[[1]][1]/log(10))*log((log(10)*x*tLs/roots[[1]][1])*exp((log(10)/roots[[1]][1])*(roots[[1]][2]-To))+exp((log(10)/roots[[1]][1])*(Tc-To)))}
      if(j==length(t_coma) & !exists("extra_ramprates")){curve(eq, from=min(dtmp$ramprate)*0.2,to=max(dtmp$ramprate)*1.2, add=TRUE,lwd=2,col = "grey")}
      }
    if(length(ramprates)==1)
      {
      ####In the case dCTmax has only been assessed for one ramp rates, one must supply an estimated z value (see main text)
      z = 2.5
      sCTmax = To+(z/log(10))*log((z/(log(10)*dtmp$ramprate[1]*tLs))*(exp((log(10)/z)*(dtmp$dCTmax[1]-To))-exp((log(10)/z)*(Tc-To))))
      sCTmax.predictionstmp$z_pred[j] = z
      sCTmax.predictionstmp$sCTmax_pred[j] = sCTmax
      eq = function(x){To+(z/log(10))*log((log(10)*x*tLs/z)*exp((log(10)/z)*(sCTmax-To))+exp((log(10)/z)*(Tc-To)))}
      if(j==length(t_coma) & !exists("extra_ramprates")){curve(eq, from=min(dtmp$ramprate)*0.2,to=max(dtmp$ramprate)*1.2, add=TRUE,lwd=2,col = "grey")}
      }
    }
  if(!exists("extra_ramprates")){points(dtmp$ramprate,dtmp$dCTmax,pch=21,bg="black",cex=1.2)}
  
  if(exists("extra_sCTmax"))
    {
    #x axis is assay temperature
    #y axis is log10(knockdown time) at constant temperature
    TTLfit = lm(log10(t_coma)~sCTmax_pred,data=sCTmax.predictionstmp)
    new = data.frame(sCTmax_pred = extra_sCTmax)
    
    sCTmax.predictions_extratmp = data.frame(group=rep(groups[i],times=length(extra_sCTmax)),extra_sCTmax=NA,extra_t_coma_pred=NA)
    sCTmax.predictions_extratmp$extra_sCTmax = extra_sCTmax
    sCTmax.predictions_extratmp$extra_t_coma_pred = 10^(predict(TTLfit, newdata = new))
    sCTmax.predictions_extra = rbind(sCTmax.predictions_extra,sCTmax.predictions_extratmp)
    rm(sCTmax.predictions_extratmp)
    }
  
  if(exists("extra_ramprates"))
    {
    dCTmax.predictionstmp = data.frame("group" = rep(unique(dtmp$group),times=length(extra_ramprates)),extra_ramprate=rep(extra_ramprates,times=length(unique(dtmp$group))),dCTmax_extra = NA,z_pred=sCTmax.predictionstmp$z_pred)
    tLs = max(sCTmax.predictionstmp$t_coma)
    #for each desired ramprate, predict dCTmax
    for(j in 1:length(extra_ramprates))
      {
      sCTmax = sCTmax.predictionstmp$sCTmax_pred[sCTmax.predictionstmp$t_coma == tLs] #sCTmax = (log10(1*1440)-intercept)/slope   #tmp$t_coma == tLs
      z = sCTmax.predictionstmp$z_pred[sCTmax.predictionstmp$t_coma == tLs]
      dCTmax = To+(z/log(10))*log((log(10)*extra_ramprates[j]*tLs/z)*exp((log(10)/z)*(sCTmax-To))+exp((log(10)/z)*(Tc-To)))
      dCTmax.predictionstmp$dCTmax_extra[j] = dCTmax
      }
    plot(dtmp$ramprate,dtmp$dCTmax,xlim=c(0,max(c(dynamic$ramprate,extra_ramprates))+0.1),ylim=c(min(c(dynamic$dCTmax,dCTmax.predictionstmp$dCTmax_extra)),max(c(dynamic$dCTmax,dCTmax.predictionstmp$dCTmax_extra))),xlab="Ramp rate (°C/min)",ylab="dCTmax (°C)",main=paste(groups[i]," "),type="n")
    curve(eq, from=min(c(dynamic$ramprate,extra_ramprates))*0.2,to=max(c(dynamic$ramprate,extra_ramprates))*1.2, add=TRUE,lwd=2,col = "grey")
    points(dtmp$ramprate,dtmp$dCTmax,pch=21,bg="black",cex=1.2)
    points(extra_ramprates,dCTmax.predictionstmp$dCTmax_extra,pch=21,bg="white",cex=1.2)
    legend("topleft",legend = c("Observed","Predicted","Dynamic TTL fit"),pch=c(21,21,NA),pt.bg = c("black","white",NA),lwd = c(NA,NA,2),col=c("black","black","grey"),bty='n')
    
    dCTmax.predictions = rbind(dCTmax.predictions,dCTmax.predictionstmp)
    rm(dCTmax.predictionstmp)
    }
  
  if(exists("fluc_temp"))
    {
    #The TTL is built from static temperatures and t_coma predicted from dynamic data. Additive injury is determined per time interval (in minutes) 
    TTLfit = lm(log10(t_coma)~sCTmax_pred,data=sCTmax.predictionstmp)
    intercept = summary(TTLfit)$coefficient[1,1]
    slope = summary(TTLfit)$coefficient[2,1]
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
  
  sCTmax.predictions = rbind(sCTmax.predictions,sCTmax.predictionstmp)
  rm(sCTmax.predictionstmp)
  print(paste("Group",i,"out of",length(groups)))
  }
dev.off()

write.csv(sCTmax.predictions,file="sCTmax_t_coma.csv")
if(exists("extra_sCTmax")){write.csv(sCTmax.predictions_extra,file="extra_sCTmax.csv",row.names = FALSE)}
if(exists("extra_ramprates")){write.csv(dCTmax.predictions,file="dCTmax_extra_ramprates.csv",row.names = FALSE)}
if(exists("fluc_temp")){write.csv(fluctemp.predictions,file="fluctemp_predictions.csv",row.names = FALSE)}
