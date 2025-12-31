# This script import data, report, config from jack mackerel stock assessment model, h1_1.07 version
# To create the operating model, historical abundance
# The is the first version on 9-may-2025 for Fernando Espíndola 
# RF == Robin Forrest add improvements to original code
# The last version on 9-nov-2025 for Fernando Espíndola with improvements made by RF
options(warn=-1);rm(list=ls())
library(ggplot2);library("dplyr")

# Read all JJM object --- 
setwd("model")
mod_current <- jjmR::readJJM("h1_1.07",path="config",input="input")
setwd("..")

#SET SAME PARAMETERS
nsim <- 20
fleet_name <- mod_current[[1]]$data$Fnames
nf <- length(fleet_name)
np <- 1 # Single sex model

# Active sd and autocorrelation for deviations of recruitment for projection
Perr <- c()
AC <- c()

#STOCK ID
Name = "Jack Mackerel"
proyears = 30
LowerTri = 0     # MSEtool definition: Integer. The number of recent years for which model estimates of recruitment are ignored (not reliably estimated by the assessment)

#Bacu_sd = 0.15   # standard deviation of acustic index observation error
#Bacu_ac = -0.14  # specify the autocorrelation of the acoustic index observation error

#Bmph_sd = 0.15   # standard deviation of mdph index observation error
#Bmph_ac = -0.14  # specify the autocorrelation of the mdph index observation error

#
# Historical parameters ----
rng_years <- mod_current[[1]]$data$years 
years <- seq(rng_years[1],rng_years[2],1)
nyears <- length(years)

# Operating model parameters ----
rng_age <- mod_current[[1]]$data$ages #1:12
age <- seq(rng_age[1],rng_age[2])
nage <- length(age)
nage_om <- max(age) + 1

spawn_time_frac <- 0.83 # Fraction of year when spawning occurs
SRrel <- 1 # (1 = Beverton-Holt, 2 = Ricker)

# Parameters scale, Warning!!! Check the length (years start 1:66) the rec_dev (Stcok-Recruitment relationship)
# Check R0 in tpl file => rec_dev==deviation of log mean history, N[1,2,,1]
# REvisar si es en escala logaritmica
#rec_dev <- mod_current[[1]]$output$Stock_1$Stock_Rec[,4]-mod_current[[1]]$output$Stock_1$Stock_Rec[,3] 
#dev_R0 <- mod_current[[1]]$output$Stock_1$rec_dev[12:66,2]

SRR <- mod_current[[1]]$output$Stock_1$Stock_Rec

# RF changes:
#  Changed years to 12:66 (1970:2024)
#  Changed EPR0 to EPR0=SRR[1,2]/SRR[1,3] (instead of EPR0=SRR[1,2]/SRR[1,4])
#  This makes R0 consistent with JJM
R0_start <- log(mean(SRR[, 4]))
h_start <- log(0.7 - 0.2)
opt <- nlminb(c(R0_start,h_start),SAMtool:::get_SR,E=SRR[12:66,2],R=SRR[12:66,4],EPR0=SRR[1,2]/SRR[1,3],type="BH")
SR_par <- SAMtool:::get_SR(opt$par,E=SRR[12:66,2],R=SRR[12:66,4],EPR0=SRR[1,2]/SRR[1,3],opt=F,figure=F,type="BH")

R0 <- SR_par$R0
SB0 <- SR_par$E0
phi0 <- SR_par$E0/SR_par$R0
h <- SR_par$h

# RF I don't think these get used but FYI they are backwards, should be
# log(SRR[,4])-log(SR_par$Rpred)

dev_SRR <- log(SR_par$Rpred)-log(SRR[,4])

# Look at values from SAMtool:::get_SR()
# JJM values: R0=16199; SB0=E0=25883; phi0=EPR0=25883/16199 = 1.598; steepness = 0.65

# Make a plot of the SR curve to see where R0 and SSB0 are coming from
jjmYPR <- mod_current[[1]]$output[[1]]$YPR

# RF add plot SR curve
png("figuras/SR_curve.png",height=3,width=4,res=400,units="in")
par(mar=c(2.9,2.8,1.7,1),mgp=c(1.8,0.5,0))
plot(jjmYPR[,2],jjmYPR[,4],pch=19,xlab="SSB",ylab="Recruits",main="JJM SRR from YPR object",xlim=c(0,31000),ylim=c(2000,38000), type="l",lty=1, lwd=2)
abline(v=jjmYPR[1,2], lty=2)
abline(h=jjmYPR[1,4], lty=2, col=2)
points(SRR[12:66,2],SRR[12:66,4],cex=0.7,col=4)
legend("topleft",legend=c("JJM SRR from output$Stock_Rec[,3]","JJM SSB0 from output$YPR[1,2] and output$Stock_Rec[1,2]","JJM R0 from output$YPR[1,4] and output$Stock_Rec[1,3]"),lty=c(1,2,2),col=c(1,1,2),
lwd=c(2,1,1),cex=0.6,bty="n")
dev.off()


# Natural mortality ----
# First age-class is zero, set mortality, maturity, weight = 0 to ensure spawners per
# recruit is equal between models
M <- local({
  maa <- array(NA, c(nsim,nage_om,nyears,np))
  for(x in 1:nsim){
    maa[x, -1, , ] <- mod_current[[1]]$parameters$N_Mort
  }
  maa[, 1, , ] <- 1e-15 + .Machine$double.eps
  maa
})
# Fishing mortality ----
# 4 fleets: 1 = N_Chile, 2 = SC_Chile, 3 = FarNorth, 4 = Offshore_Trawl
F_fleet <- sapply(1:4,function(i) mod_current[[1]]$output[[1]][[paste0("F_age_",i)]][,-1],simplify="array")
F_age <- apply(F_fleet, c(1,2), sum)
#
FM <- local({
  faa <- array(NA,c(nsim,nage_om,nyears,np,nf))
  for(x in 1:nsim) {
    for(f in 1:nf){faa[x,-1,,,f] <- t(F_fleet[,,f])}
  }
  faa[, 1, , , ] <- 0
  faa
})

# Sample random walk in selectivity in projection
V_fleet <- apply(F_fleet, c(1, 3), function(x) x/max(x))
V_all <- apply(F_age, 1, function(x) x/max(x))

# RF add this plot
#for(i in 1:nf){
#    matplot(1:nage,V_fleet[,,i], type="l", main = paste("Fleet",i))
#    lines(1:nage,V_fleet[,55,i], lwd=2)
#}
## matplot(t(mod_current$h1_1.07$output$Stock_1$sel_fsh_1[,3:14]), type="l", xlab="Age", ylab="mod_current sel_fsh_1", main="JJM output sel_fsh_1")
#matplot(V_fleet[,,1], type="l", xlab="Age", ylab="V_fleet fleet 1", main="Sim Vfleet Fleet 1")
## matplot(t(mod_current$h1_1.07$output$Stock_1$sel_fsh_2[,3:14]), type="l", xlab="Age", ylab="mod_current sel_fsh_2", main="JJM output sel_fsh_2")
#matplot(V_fleet[,,2], type="l", xlab="Age", ylab="V_fleet fleet 2", main="Sim Vfleet Fleet 2")
## matplot(t(mod_current$h1_1.07$output$Stock_1$sel_fsh_3[,3:14]), type="l", xlab="Age", ylab="mod_current sel_fsh_3", main="JJM output sel_fsh_3")
#matplot(V_fleet[,,3], type="l", xlab="Age", ylab="V_fleet fleet 3", main="Sim Vfleet Fleet 3")
## matplot(t(mod_current$h1_1.07$output$Stock_1$sel_fsh_4[,3:14]), type="l", xlab="Age", ylab="mod_current sel_fsh_4", main="JJM output sel_fsh_4")
#matplot(V_fleet[,,4], type="l", xlab="Age", ylab="V_fleet fleet 4", main="Sim Vfleet Fleet 4")
#

# RF cahnged nyears to proyears
# VPROJ <- array(NA,c(nage,nyears,nsim,nf))
VPROJ <- array(NA,c(nage,proyears,nsim,nf))


set.seed(324)
for(f in 1:nf){
  Vproj <- sapply(1:nsim, function(...) {
    delta <- rnorm(proyears, 0, 0.025)
    delta2 <- rnorm(proyears * length(age), 0, 0.005) %>% matrix(proyears, length(age))
    #Vlast <- F_age[nyears, ]/max(F_age[nyears, ])
    Vlast <- F_fleet[nyears,,f]/max(F_fleet[nyears,,f]) 
    Vout <- array(NA, c(length(age), proyears))
    Vout[, 1] <- log(Vlast) + delta[1] + delta2[1, ]
    for(y in 2:proyears) Vout[, y] <- Vout[, y-1] + delta[y] + delta2[y, ]
      V <- exp(Vout)
      V[V > 1] <- 1
      #V[11:12, ] <- 1 # Team commented this out on day 3
      return(V)
    },
    simplify = "array") #%>% aperm(c(3,1,2))
    VPROJ[,,,f] <- Vproj
}

# Selectividad
year_proy <- seq(years[length(years)]+1,length=proyears)

edad <- c(2,4,6,8,10,12);nedad<-length(edad);cols<-1:nedad
png("figuras/selectividad_age.png",height=4,width=5,res=400,units="in")
layout(matrix(c(1:4,5,5),ncol=2,byrow=T),widths=c(1,1),heights=c(5,5,1))
par(mar=c(5,4,2,2),cex.axis=0.8,cex.lab=0.9)
for(j in 1:4){
  matplot(year_proy,t(VPROJ[edad, ,1 ,j]),typ="l",xlab="Años",yaxt="n",
  ylab="Selectividad",col=0)
  for(i in 1:nedad){lines(year_proy,t(VPROJ[edad[i], ,1 , j]),lty=1,lwd=2,col=cols[i])}
  axis(2,las=1)
  mtext(fleet_name[j],side=3,line=0.7,cex=0.8)
}
par(mar=c(0.1,0.1,0.1,0.1))
plot(1:10,1:10,xlab="",ylab="",xaxt="n",yaxt="n",axes=FALSE,col=0)
legend("center",legend=edad,title="Edades",horiz=TRUE,col=cols,lty=1,lwd=2,cex=0.8,box.col=0)
#plot(age, Maturity[1,,1],typ='o',xlab="Age",ylab="Maturity")
dev.off()



# Abundance at age ----
N <- local({
  FM_sum <- apply(FM, 1:4, sum)
  naa <- array(NA,c(nsim, nage_om, nyears,np))
  for(x in 1:nsim) {
    naa[x, age + 1, , ] <- t(mod_current[[1]]$output$Stock_1$N[,2:nage_om])
  }
  # Backwards calculation for age 0 abundance
  for(a in 1) {
    naa[, a, 2:nyears - 1, ] <- naa[, a+1, 2:nyears, ] /
      exp(-1*(FM_sum[, a, 2:nyears - 1, ] + M[, a, 2:nyears - 1, ])) # F+M should be numerically zero!
  }
  naa
})


#plot_age<-2:6
#matplot(years, t(FM[1,plot_age,,1]),typ="l",ylab="F",col=plot_age,lty=plot_age-1)
#legend("topleft",legend=plot_age,col=plot_age,lty=plot_age-1)

#persp(x=0:max(age),y=years,z=FM[1,,1],xlab="Age",ylab="Year",zlab="F",theta=45 + 180)


# Maturity ogive ----
Maturity <- local({
  mat <- array(NA,c(nsim,nage_om,nyears,np))
  for(x in 1:nsim) {
    mat[x, -1, , ] <- t(mod_current[[1]]$output$Stock_1$mature_a)
  }
  mat[, 1, , ] <- 0 + .Machine$double.eps
  mat
})


png("figuras/madurez_age.png",height=4,width=5,res=400,units="in")
par(mar=c(5,4,2,2),cex.axis=0.8,cex.lab=0.8)
matplot(c(age,13),Maturity[1, ,1 ,1],typ="l",xlab="Edades",ylab="Madurez sexual",col=1:13)
for(i in 1:13){lines(c(age,13),Maturity[1, , i, 1],lty=1,lwd=2,col=i)}
legend("topright",legend=1:(max(age)+1),col=1:13,lty=1:13,cex=0.6)
#plot(age, Maturity[1,,1],typ='o',xlab="Age",ylab="Maturity")
dev.off()


# Stock weight at age ----
Weight_age <- local({
  wt <- array(NA,c(nsim,nage_om,nyears,np))
  for(x in 1:nsim) {
    wt[x, -1, , ] <- t(mod_current[[1]]$output$Stock_1$wt_a_pop)
  }
  wt[, 1, , ] <- 0
  wt
})
#


png("figuras/pesos_age.png",height=4,width=5,res=400,units="in")
par(mar=c(5,4,2,2),cex.axis=0.8,cex.lab=0.8)
matplot(years,t(Weight_age[1, , ,1]),typ="l",yaxt="n",xlab="Años",ylab="Pesos a la edad",col=1:13)
#for(i in 1:13){lines(c(age,13),Maturity[1, , i, 1],lty=1,lwd=2,col=i)}
legend("topright",legend=1:(max(age)+1),title="Edades",col=1:13,lty=1:13,cex=0.6)
axis(2,las=1)
#plot(age, Maturity[1,,1],typ='o',xlab="Age",ylab="Maturity")
dev.off()

Fec_age <- local({
  fec <- array(NA,c(nsim,nage_om,nyears,np))
  # Sum over length bins to find fecundity at age
  for(x in 1:nsim) {
    for(y in 1:nyears) {
      fec_y <- t(mod_current[[1]]$output$Stock_1$mature_a*mod_current[[1]]$output$Stock_1$wt_a_pop) # Matrix of fecundity by length and age
      fec[x,-1 , y, 1] <- fec_y
    }
  fec[, 1, , ] <- 0
  }
  fec
})


# Length at age 
# Growth parameters
Linf <- mod_current[[1]]$control$Linf[1,1]
K <- mod_current[[1]]$control$K[1,1]
L0 <- mod_current[[1]]$control$Lo_len[1,1]
t0 <- log(1 - L0/Linf)/K
#
# The current ADMB model does not use length data (Peru have length data)
# We need a placeholder for the operating model but we will have to ignore the length modeling
Length_age <- local({
  laa <- array(NA, c(nsim, nage_om, nyears,np))
  for(x in 1:nsim) {
    #laa[x,-1, , ] <- Linf * (1 - exp(-K * (age - t0)))
    #laa[x, 1, , ] <- L0
    laa[x, , , 1] <- 1:nage_om # ignore the length modelling 
  }
  laa
})


#
# sigma_edad - variability in length at age #linea 2104 esta ALK == Sigma_len:0.09-0.1-(-4)
#LatASD <- admb_rep[["alfa_edad"]] + admb_rep[["beta_edad"]] * Length_age
# Prob_age2len => mod:current[[1]]$output$Stock_1$P_age2len_1 

# Recruitment sampling
#log_Rdev <- admb_rep$Reclutas[3, ]
if (!missing(Perr)) Perr <- sd(dev_SRR,na.rm=TRUE) # 0.632
if (!missing(AC)) AC <- acf(dev_SRR[!is.na(dev_SRR)],plot=FALSE)$acf[2] # 0.827

#plot(dev_R0, panel.first = {grid(); abline(h = 0, lty = 2)}, typ = 'o')


# Make operating model ----
MOM <- MSEtool::Assess2MOM(
  Name = Name,
  proyears = proyears,
  CurrentYr = max(years),
  h = h,
  naa = N,
  faa = FM,
  waa = Weight_age,
  Mataa = Maturity,
  Maa = M,
  laa = Length_age,
  fecaa = Fec_age,
  nyr_par_mu = 1,
  LowerTri = LowerTri,
  R0 = R0,
  phi0 = phi0,
  SRrel = SRrel,
  spawn_time_frac = spawn_time_frac,
  Perr = Perr,
  AC = AC,
  Obs = MSEtool::Precise_Unbiased, #RF added these. Default was Imprecise_Unbiased
  Imp = MSEtool::Perfect_Imp
)


# Condition operating model

MOM@Stocks[[1]]@Name <- Name
MOM@Fleets[[1]][[1]]@Name <- fleet_name[1]
MOM@Fleets[[1]][[2]]@Name <- fleet_name[2]
MOM@Fleets[[1]][[3]]@Name <- fleet_name[3]
MOM@Fleets[[1]][[4]]@Name <- fleet_name[4]

# This is a placeholder for growth parameters that we will not use in the operating model ----
MOM@cpars[[1]][[1]]$Linf <- MOM@cpars[[1]][[1]]$K <- MOM@cpars[[1]][[1]]$t0 <- rep(1, nsim)

index_name <- mod_current[[1]]$data$Inames

# Revisar la estructura de los indices de abundancia, ACUSPeru, ACusChile
# Ver el peso de los indices

for (f in 1:nf) {#nf

  # Add acoustic index and fishery CPUE to operating model ----
  # See class?Data
  RealData <- new("Data")

  # Add catch to operating model ----

  RealData@Cat <- matrix(mod_current[[1]]$data$Fcaton[,f], 1, nyears)
  RealData@CV_Cat <- matrix(mod_current[[1]]$data$Fcatonerr[,f], 1, nyears)

  # Add catch at age to operating model ----  
  # Check number or proportion at age ---
  RealData@CAA <- array(NA, c(1, nyears, nage_om))
  if(f!=3){
    RealData@CAA[1, , -1] <- mod_current[[1]]$data$Fagecomp[,,f]
  } else {
    year_peru <- mod_current[[1]]$output$Stock_1$phat_len_fsh_3[,1] 
    tallasperu <- mod_current[[1]]$output$Stock_1$phat_len_fsh_3[,2:42]
    claveperu<-mod_current[[1]]$output$Stock_1$P_age2len_1
    # Redimencionar a nyears y strart in 11 postition
    cedadperu<-matrix(NA,nyears,12)
    for(i in 11:nyears){
      for (j in 1:12){
      cedadperu[i,j]=sum(tallasperu[i-10,]*claveperu[j,])
      }
    }
    #ceadadperuanios=cbind(mod_current[[1]]$output$Stock_1$phat_len_fsh_3[,1],cedadperu)
    RealData@CAA[1, , -1] <- cedadperu
  }
  # tallasperu <- mod_current[[1]]$output$Stock_1$phat_len_fsh_3[2:42] # 2:42 bin length
  # length => 20:50 (1 step)
  # mod:current[[1]]$output$Stock_1$P_age2len_1
  
  # block, which will continue in the projection
  #y_f <- c(2011, 2011, 2010) # Starting year of most recent catchability block

  # Add fishery CPUE to operating model
  # openMSE assumes indices of abundance have constant catchability
  # If you want the most recent catchability, make NA the years old, like exemple...
  #RealData@VInd[1, years < y_f[f]] <- RealData@CV_VInd[1, years < y_f[f]] <- NA
  #RealData@VInd[RealData@VInd == 0] <- RealData@CV_VInd[RealData@VInd == 0] <- NA
  #RealData@CV_AddInd <- ifelse(admb_dat$indices[, 4] > 0, 0.15, NA) %>% array(c(1, 1, OM@nyears))

  if(f==2){#CPUE CS Chile
    RealData@VInd <- matrix(mod_current[[1]]$data$Index[,3], 1, nyears)
    RealData@CV_VInd <- ifelse(mod_current[[1]]$data$Index[,3] > 0, 0.15, NA) %>% matrix(1,nyears) 
    #matrix(mod_current[[1]]$data$Indexerr[,3], 1, nyears)
  }
  if(f==3){#CPUE FarNorth
    RealData@VInd <- matrix(mod_current[[1]]$data$Index[,6], 1, nyears)
    RealData@CV_VInd <- ifelse(mod_current[[1]]$data$Index[,6] > 0, 0.20, NA) %>% matrix(1,nyears) 
    #matrix(mod_current[[1]]$data$Indexerr[,6], 1, nyears)
  }
  if(f==4){#CPUE Offshore Trawl
    RealData@VInd <- matrix(mod_current[[1]]$data$Index[,7], 1, nyears)
    RealData@CV_VInd <- ifelse(mod_current[[1]]$data$Index[,7] > 0, 0.20, NA) %>% matrix(1,nyears)  
    #matrix(mod_current[[1]]$data$Indexerr[,7], 1, nyears)
  }
 
  # Add ACOUSTIC index to operating model - only values since 2019 for the most recent
  # catchability and selectivity block => y_i<-2019 => make NA value

  if(f==1){#_CHILE_NORTE==N_Chile
    RealData@AddInd <- array(mod_current[[1]]$data$Index[,2], c(1, 1, nyears))
    RealData@CV_AddInd <- ifelse(mod_current[[1]]$data$Index[,2] > 0, 0.50, NA) %>% array(c(1,1,nyears))  
    #array(mod_current[[1]]$data$Indexerr[,2], c(1, 1, nyears))
    #RealData@AddInd[1, 1, years < y_i] <- RealData@CV_AddInd[1, 1, years < y_i] <- NA
    #RealData@AddInd[RealData@AddInd == 0] <- RealData@CV_AddInd[RealData@AddInd == 0] <- NA
  
    # Add selectivity fron acustic index, ojo ver si son dos cero y no uno???
    RealData@AddIndV <- array(c(0,mod_current[[1]]$output$Stock_1$sel_ind_2[nyears,3:14]),c(1,1,nage_om))
    RealData@AddIunits <- 1  # Biomass or Abundance
    RealData@AddIndType <- 1
  
    # Sample observation error of the index for the projection
    # By default, this code calculates the standard error and autocorrelation from
    # the residuals in 2007-2022 => y_i<2019 => make NA value
    Yrs <- mod_current[[1]]$output$Stock_1$Obs_Survey_2[,1]
    Bacu_obs <- mod_current[[1]]$output$Stock_1$Obs_Survey_2[,2] #_2 column with observate data
    #Bacu_obs[Bacu_obs == 0] <- NA
    Bacu_pred <- mod_current[[1]]$output$Stock_1$Obs_Survey_2[,3] #_3 column with estimate data
    logi_year_f1 <- years %in% Yrs
    BAP_full_f1 <- vector("numeric",nyears); BAP_full_f1[] <- NA; k <- 1
    for(i in 1:nyears){if(logi_year_f1[i]=="TRUE"){BAP_full_f1[i]<-Bacu_pred[k];k<-k+1}}
    #plot(admb_rep$Yrs, Bacu_obs, typ = 'o')
    #lines(admb_rep$Yrs, Bacu_pred, col = 2, typ = 'o')
    #log_resid <- log(Bacu_obs[admb_rep$Yrs %in% c(2007:2022)]/Bacu_pred[admb_rep$Yrs %in% c(2007:2022)])
    log_resid_f1 <- log(Bacu_obs/Bacu_pred)
    Bacu_sd1 <- sd(log_resid_f1, na.rm = TRUE) # 0.81
    Bacu_ac1 <- acf(log_resid_f1,plot=FALSE)$acf[2]       # 0.02
  
    MOM@cpars[[1]][[f]]$AddIbeta <- matrix(1, nsim, 1)

    Ierr_f1 <- array(NA,c(nsim,1,nyears + proyears))
    Ierr_f1[, 1, 1:nyears] <- matrix(RealData@AddInd[1, 1, ]/BAP_full_f1,nsim,nyears,byrow = TRUE)
    res_ind_f1<-RealData@AddInd[1,1,]/BAP_full_f1
    res_ind_f1<-res_ind_f1[which(is.na(res_ind_f1)=="FALSE")]
    # Comment: correlation between source index, e.g. CPUE and index acustic 
    # We sample for the last historical year as well! (NA in the conditioning model)
    #log_Ierr_proj <- MSEtool:::sample_recruitment(
    #  log(Ierr[, 1, seq(1, nyears - 1)]),
    #  proyears = proyears + 1, procsd = Bacu_sd, AC = Bacu_ac, seed = 29
    #)
    log_Ierr_proj <- MSEtool:::sample_recruitment(matrix(log(res_ind_f1),nsim,length(res_ind_f1),byrow=TRUE),
      proyears = proyears, procsd = Bacu_sd1, AC = Bacu_ac1
    )
    Ierr_f1[, 1, seq(nyears+1, nyears + proyears)] <- exp(log_Ierr_proj)
    MOM@cpars[[1]][[f]]$AddIerr <- Ierr_f1
  }
  if(f==2){#_CHILE_CENTRO_SUR==SC_Chile_PS
    RealData@AddInd <- array(mod_current[[1]]$data$Index[,1], c(1, 1, nyears))
    RealData@CV_AddInd <- ifelse(mod_current[[1]]$data$Index[,1] > 0, 0.20, NA) %>% array(c(1,1,nyears))  
    #array(mod_current[[1]]$data$Indexerr[,1], c(1, 1, nyears))
    #RealData@AddInd[1, 1, years < y_i] <- RealData@CV_AddInd[1, 1, years < y_i] <- NA
    #RealData@AddInd[RealData@AddInd == 0] <- RealData@CV_AddInd[RealData@AddInd == 0] <- NA
  
    # Add selectivity fron acustic index
    RealData@AddIndV <- array(c(0,mod_current[[1]]$output$Stock_1$sel_ind_1[nyears,3:14]),c(1,1,nage_om))
    RealData@AddIunits <- 1  # Biomass
    RealData@AddIndType <- 1
  
    # Sample observation error of the index for the projection
    # By default, this code calculates the standard error and autocorrelation from
    # the residuals in 2007-2022 => y_i<2019 => make NA value
    Yrs <- mod_current[[1]]$output$Stock_1$Obs_Survey_1[,1]
    Bacu_obs <- mod_current[[1]]$output$Stock_1$Obs_Survey_1[,2] #_2 column with observate data
    #Bacu_obs[Bacu_obs == 0] <- NA
    Bacu_pred <- mod_current[[1]]$output$Stock_1$Obs_Survey_1[,3] #_3 column with estimate data
    logi_year_f2 <- years %in% Yrs
    BAP_full_f2 <- vector("numeric",nyears); BAP_full_f2[] <- NA; k <- 1
    for(i in 1:nyears){if(logi_year_f2[i]=="TRUE"){BAP_full_f2[i] <- Bacu_pred[k];k<-k+1}}
    #plot(admb_rep$Yrs, Bacu_obs, typ = 'o')
    #lines(admb_rep$Yrs, Bacu_pred, col = 2, typ = 'o')
    #log_resid <- log(Bacu_obs[admb_rep$Yrs %in% c(2007:2022)]/Bacu_pred[admb_rep$Yrs %in% c(2007:2022)])
    log_resid_f2 <- log(Bacu_obs/Bacu_pred)
    Bacu_sd2 <- sd(log_resid_f2,na.rm = TRUE)       # 0.81
    Bacu_ac2 <- acf(log_resid_f2,plot=FALSE)$acf[2]  # 0.02

    MOM@cpars[[1]][[f]]$AddIbeta <- matrix(1, nsim, 1)
    
    Ierr_f2 <- array(NA, c(nsim, 1, nyears + proyears))
    Ierr_f2[, 1, 1:nyears] <- matrix(RealData@AddInd[1, 1, ]/BAP_full_f2,nsim,nyears,byrow = TRUE)
    res_ind_f2<-RealData@AddInd[1,1,]/BAP_full_f2
    res_ind_f2<-res_ind_f2[which(is.na(res_ind_f2)=="FALSE")]
    # We sample for the last historical year as well! (NA in the conditioning model)
    #log_Ierr_proj <- MSEtool:::sample_recruitment(
    #  log(Ierr[, 1, seq(1, nyears - 1)]),
    #  proyears = proyears + 1, procsd = Bacu_sd, AC = Bacu_ac, seed = 29
    #)
    log_Ierr_proj <- MSEtool:::sample_recruitment(matrix(log(res_ind_f2),nsim,length(res_ind_f2),byrow=TRUE),
      proyears = proyears, procsd = Bacu_sd2, AC = Bacu_ac2
    )
    Ierr_f2[, 1, seq(nyears+1, nyears + proyears)] <- exp(log_Ierr_proj)
    MOM@cpars[[1]][[f]]$AddIerr <- Ierr_f2
    #
    # Projected error structure of the MPH index ---
    #
    RealData@SpInd <- matrix(mod_current[[1]]$data$Index[,4], c(1, nyears))
    RealData@CV_SpInd <- ifelse(mod_current[[1]]$data$Index[,4] > 0, 0.50, NA) %>% matrix(1,nyears)
    #matrix(mod_current[[1]]$data$Indexerr[,4], c(1, nyears)) # Is decimal o real value??

    Yrs_mph <- mod_current[[1]]$output$Stock_1$Obs_Survey_4[,1]
    Bmph_obs <- mod_current[[1]]$output$Stock_1$Obs_Survey_4[,2] #_2 column with observate data
    #Bacu_obs[Bacu_obs == 0] <- NA
    Bmph_pred <- mod_current[[1]]$output$Stock_1$Obs_Survey_4[,3] #_3 column with estimate data
    logi_year_mph <- years %in% Yrs_mph
    BMP_full=BMO_full <- vector("numeric",nyears); BMP_full[] <- NA; k <- 1; r <- 1
    for(i in 1:nyears){if(logi_year_mph[i]=="TRUE"){BMP_full[i] <- Bmph_pred[k];k<-k+1}}
    for(i in 1:nyears){if(logi_year_mph[i]=="TRUE"){BMO_full[i] <- Bmph_obs[r];r<-r+1}}
    #plot(admb_rep$Yrs, Bacu_obs, typ = 'o')
    #lines(admb_rep$Yrs, Bacu_pred, col = 2, typ = 'o')
    #log_resid <- log(Bacu_obs[admb_rep$Yrs %in% c(2007:2022)]/Bacu_pred[admb_rep$Yrs %in% c(2007:2022)])
    log_resid_mph <- log(Bmph_obs/Bmph_pred)
    Bmph_sd <- sd(log_resid_mph, na.rm = TRUE)       # 0.81
    Bmph_ac <- acf(log_resid_mph,plot=FALSE)$acf[2]  # 0.02

    MOM@cpars[[1]][[f]]$SpI_beta <- rep(1, nsim) 

    Ierr_mph <- array(NA, c(nsim, 1, nyears + proyears))
    Ierr_mph[, 1, 1:nyears] <- matrix(BMO_full/BMP_full, nsim, nyears, byrow = TRUE)
    #res_ind<-/BMP_full
    #res_ind<-res_ind[which(is.na(res_ind)=="FALSE")]
    # We sample for the last historical year as well! (NA in the conditioning model)
    # Warning: check dimension of Ierr_mph
    log_Ierr_proj_mph <- MSEtool:::sample_recruitment(
      matrix(log_resid_mph,nsim,length(log_resid_mph),byrow=TRUE),
      proyears = proyears, procsd = Bmph_sd, AC = Bmph_ac
    )
    Ierr_mph[, 1, seq(nyears+1, nyears + proyears)] <- exp(log_Ierr_proj_mph)
    MOM@cpars[[1]][[f]]$SpIerr_y <- Ierr_mph[, 1,]
  }
  if(f==3){#_PERU==FarNorth
    RealData@AddInd <- array(mod_current[[1]]$data$Index[,5], c(1, 1, nyears))
    RealData@CV_AddInd <-ifelse(mod_current[[1]]$data$Index[,5] > 0, 0.20, NA) %>% array(c(1,1,nyears))  
    #array(mod_current[[1]]$data$Indexerr[,5], c(1, 1, nyears))
    #RealData@AddInd[1, 1, years < y_i] <- RealData@CV_AddInd[1, 1, years < y_i] <- NA
    #RealData@AddInd[RealData@AddInd == 0] <- RealData@CV_AddInd[RealData@AddInd == 0] <- NA
  
    # Add selectivity fron acustic index
    RealData@AddIndV <- array(c(0,mod_current[[1]]$output$Stock_1$sel_ind_5[nyears,3:14]),c(1,1,nage_om))
    RealData@AddIunits <- 1  # Biomass
    RealData@AddIndType <- 1
  
    # Sample observation error of the index for the projection
    # By default, this code calculates the standard error and autocorrelation from
    # the residuals in 2007-2022 => y_i<2019 => make NA value
    Yrs <- mod_current[[1]]$output$Stock_1$Obs_Survey_5[,1]
    Bacu_obs <- mod_current[[1]]$output$Stock_1$Obs_Survey_5[,2] #_2 column with observate data
    #Bacu_obs[Bacu_obs == 0] <- NA
    Bacu_pred <- mod_current[[1]]$output$Stock_1$Obs_Survey_5[,3] #_3 column with estimate data
    logi_year_f3 <- years %in% Yrs
    BAP_full_f3 <- vector("numeric",nyears); BAP_full_f3[] <- NA; k <- 1
    for(i in 1:nyears){if(logi_year_f3[i]=="TRUE"){BAP_full_f3[i] <- Bacu_pred[k];k<-k+1}}
    #plot(admb_rep$Yrs, Bacu_obs, typ = 'o')
    #lines(admb_rep$Yrs, Bacu_pred, col = 2, typ = 'o')
    #log_resid <- log(Bacu_obs[admb_rep$Yrs %in% c(2007:2022)]/Bacu_pred[admb_rep$Yrs %in% c(2007:2022)])
    log_resid_f3 <- log(Bacu_obs/Bacu_pred)
    Bacu_sd3 <- sd(log_resid_f3, na.rm = TRUE)       # 0.81
    Bacu_ac3 <- acf(log_resid_f3,plot=FALSE)$acf[2]  # 0.02
    
    MOM@cpars[[1]][[f]]$AddIbeta <- matrix(1, nsim, 1)
    
    Ierr_f3 <- array(NA, c(nsim, 1, nyears + proyears))
    Ierr_f3[, 1, 1:nyears] <- matrix(RealData@AddInd[1, 1, ]/BAP_full_f3,nsim,nyears,byrow = TRUE)
    res_ind_f3<-RealData@AddInd[1,1,]/BAP_full_f3
    res_ind_f3<-res_ind_f3[which(is.na(res_ind_f3)=="FALSE")]
    # We sample for the last historical year as well! (NA in the conditioning model)
    #log_Ierr_proj <- MSEtool:::sample_recruitment(
    #  log(Ierr[, 1, seq(1, nyears - 1)]),
    #  proyears = proyears + 1, procsd = Bacu_sd, AC = Bacu_ac, seed = 29
    #)
    log_Ierr_proj <- MSEtool:::sample_recruitment(
      matrix(log(res_ind_f3),nsim,length(res_ind_f3),byrow=TRUE),
      proyears = proyears, procsd = Bacu_sd3, AC = Bacu_ac3
    )
    Ierr_f3[, 1, seq(nyears+1, nyears + proyears)] <- exp(log_Ierr_proj)
    MOM@cpars[[1]][[f]]$AddIerr <- Ierr_f3

  }


  # Add fishery weight at age to the operating model
  # wt_f <- paste("Peso", c("arrastre", "palangre", "espinel"))
  wt_f <- paste0("wt_fsh_",1:4)
  Cobs <- sapply(1:4, function(i) mod_current[[1]]$output[[1]][[paste0("Obs_catch_", i)]])
  Z <- F_age + 0.28
  CN <- sapply(1:4, function(f) F_fleet[, , f]/Z * (1 - exp(-Z)) * t(N[1, -1, ,1]), simplify = "array")
  CB <- CN * mod_current[[1]]$data$Fwtatage
  Wt_age_C <- apply(CB, 1:3, sum)/apply(CN, 1:3, sum)

  out <- array(0, c(nsim, nage_om, nyears + proyears))
  out[, -1, 1:nyears] <- array(Wt_age_C[,,f],c(nyears,nage,nsim)) %>% aperm(3:1)
  out[, -1, nyears + 1:proyears] <- array(Wt_age_C[nyears,,f],c(nage,proyears,nsim)) %>% aperm(c(3,1,2))
  #out
  wt_age_c <- local({
    wt <- array(0, c(nsim, nage_om, nyears + proyears))
    for (x in 1:nsim) {
      wt[x, -1, 1:nyears] <- t(mod_current[[1]]$output$Stock_1[[wt_f[f]]][1:nyears,2:13])
      wt[x, , nyears + 1:proyears] <- matrix(wt[x, , nyears], nage_om, proyears)
    }
    wt
  })

  MOM@cpars[[1]][[f]]$Wt_age_C <- wt_age_c 

  # RF Update selectivity in projection
  for(k in 1:nsim){
    MOM@cpars[[1]][[f]]$V[k, -1, nyears + 1:proyears] <- VPROJ[,,k,f] 
  }
  # Pass data to OM ----
  MOM@cpars[[1]][[f]]$Data <- RealData

  # Specify that the index is proportional to the vulnerable biomass ----
  MOM@cpars[[1]][[f]]$VI_beta <- rep(1, nsim)

  # Specify the multinomial sample size of the fishery catch at age
  MOM@Obs[[1]][[f]]@CAA_ESS <- c(100, 100)

  # Specify the nominal sample size (expands values returned by dmultinom function)
  MOM@Obs[[1]][[f]]@CAA_nsamp <- c(100, 100)

}

# End MOM configuration

# Resample recruitment for projection ----
Rec_his <- mod_current[[1]]$output$Stock_1$Stock_Rec[,4]

drec   <- dim(mod_current[[1]]$output$Stock_1$Stock_Rec)[1] # 66
irec   <- drec-nyears+1  # 12

recl   <- mod_current[[1]]$output$Stock_1$Stock_Rec[irec:drec,4]
mrec   <- mean(recl)
sd_rec <- sd(log(recl/mrec))
ac_rec <- acf(log(recl/mrec),plot=FALSE)


set.seed(28)
dev_rec <- MSEtool:::sample_recruitment(
  Perr_hist = log(recl/mrec) %>% matrix(nsim, 1),
  proyears = proyears,
  procsd = sd_rec,
  AC = ac_rec$acf[2]
)
Perr_y <- exp(dev_rec)
R0_OM <- MOM@cpars[[1]][[1]]$R0[1]
Perr_sy <- Perr_y*(mrec/R0_OM)
Perr_y_orig <- Perr_sy
Rec_ori <- Perr_y_orig*R0_OM

MOM@cpars[[1]][[1]]$Perr_y[, length(age) + nyears + 1:proyears] <- Perr_y_orig

MOM1 <- MOM2 <- MOM3 <- MOM4 <- MOM

# Save model ----
saveRDS(MOM1,file="MOM/jurel_MOM1.rds")

#============================================================
# RF update Perr_y for alternative projections
# Could refine these in many ways, including having a regime shift during proj period
# e.g., https://cdnsciencepub.com/doi/10.1139/F09-142, https://academic.oup.com/icesjms/article/71/8/2208/744434
#============================================================

# 1. Resampled Perr_sy with low recruitment

recl   <- mod_current[[1]]$output$Stock_1$Stock_Rec[45:58,4] # 2003 al 2016 
mrec   <- mean(recl)
sd_rec <- sd(log(recl/mrec))
ac_rec <- acf(log(recl/mrec),plot=FALSE)


set.seed(28)
dev_rec <- MSEtool:::sample_recruitment(
  Perr_hist = log(recl/mrec) %>% matrix(nsim, 1),
  proyears = proyears,
  procsd = sd_rec,
  AC = ac_rec$acf[2]
)
Perr_y <- exp(dev_rec)
R0_OM <- MOM@cpars[[1]][[1]]$R0[1]
Perr_sy <- Perr_y*(mrec/R0_OM)
Perr_y_low <- Perr_sy # Scenario 1: Object to add to cpars
Rec_low <- Perr_y_low*R0_OM

MOM2@cpars[[1]][[1]]$Perr_y[, length(age) + nyears + 1:proyears] <- Perr_y_low

saveRDS(MOM2,file="MOM/jurel_MOM2.rds")

# 2. Only sample from historical period with good recruitment - include AC
# Pick 1985:1995
# RF recommend pick periodo, 1985:1995 [27:37]
recl   <- mod_current[[1]]$output$Stock_1$Stock_Rec[27:43,4] # 1985 al 2001 
mrec   <- mean(recl)
sd_rec <- sd(log(recl/mrec))
ac_rec <- acf(log(recl/mrec),plot=FALSE)

set.seed(28)
dev_rec <- MSEtool:::sample_recruitment(
    Perr_hist = log(recl/mrec) %>% matrix(nsim, 1),
    proyears = proyears,
    procsd = sd_rec,    # 0.4752
    AC = ac_rec$acf[2]  # 0.5474
)
Perr_y <- exp(dev_rec)
R0_OM <- MOM@cpars[[1]][[1]]$R0[1]
Perr_sy <- Perr_y*(mrec/R0_OM)
Perr_y_high <- Perr_sy  # Scenario 3: Object to add to cpars
Rec_high <- Perr_y_high*R0_OM

MOM4@cpars[[1]][[1]]$Perr_y[, length(age) + nyears + 1:proyears] <- Perr_y_high

# Save model ----
saveRDS(MOM4,file="MOM/jurel_MOM3.rds")

# Add to plot and add legend
ini_years <- years[1] - irec + 1
year_ini <- seq(ini_years,years[1]-1,by=1)
year_his <- mod_current[[1]]$output$Stock_1$Stock_Rec[,1]
year_sta <- years[length(years)]+1
year_pro <- seq(year_sta,length=proyears)
year_tot <- c(year_his,year_pro)
year_plu <- year_tot[length(year_tot)]+1


ay1<-seq(0,71000,by=10000)
png("figuras/reclutas_tipos_moms.png",height=4,width=6,res=400,units="in")
par(mar=c(3.3,3.3,2,1.9),mgp=c(1.8,0.7,0))
#par(mar=c(5,4.5,2,2))
plot(year_tot,c(Rec_his,Rec_ori[1,]),typ="l",ylab="Reclutamientos x 10000",yaxt="n",ylim=c(0,max(ay1)),
xlab="Años",main="Simulación 1",col="white",cex=0.8)
#abline(h=mean(Rec_his),lty=3)
abline(h=R0,lty=3)
abline(v=2024,lty=2)
lines(year_his,Rec_his,lwd=1,col=1)
# Add Rec Scenarios
points(year_pro,Rec_ori[1,],pch=21,col=1,bg=1,cex=0.5)
points(year_pro,Rec_low[1,],pch=21,col=2,bg=2,cex=0.5)
points(year_pro,Rec_high[1,],pch=21,col=3,bg=3,cex=0.5)
#
text(year_tot[length(year_tot)]+4.2,R0,expression(R[paste("0")]),cex=0.9,pos=4,xpd=T,col=1,font=2)
axis(2,at=ay1,labels=sprintf("%1.0f",ay1/10000),las=1)
legend("topleft",legend=c("Reclutamientos históricos","Re-muestreado Perr_y con serie histórica",
"Re-muestreado Perr_y con AC (bajo)","Re-muestreado Perr_y con AC (alto)"),
pch=c(NA,21,21,21),pt.bg=c(NA,1,2,3),lty=c(1,NA,NA,NA),
col=c(1,1:3),bty="n",cex=0.8)
dev.off()



# Compare operating model historical dynamics and ADMB assessment ----
MHist <- MSEtool:::SimulateMOM(MOM,parallel=FALSE,silent=FALSE)

MHist[[1]][[1]]@SampPars$Obs$CAA_nsamp <-
MHist[[1]][[2]]@SampPars$Obs$CAA_nsamp <-
MHist[[1]][[3]]@SampPars$Obs$CAA_nsamp <-
MHist[[1]][[4]]@SampPars$Obs$CAA_nsamp <- rep(400, 100)

for(i in 1:3){
  MHist <- MSEtool:::SimulateMOM(get(paste0("MOM",i)),parallel=TRUE,silent=FALSE)
  MHist[[1]][[1]]@SampPars$Obs$CAA_nsamp <-
  MHist[[1]][[2]]@SampPars$Obs$CAA_nsamp <-
  MHist[[1]][[3]]@SampPars$Obs$CAA_nsamp <-
  MHist[[1]][[4]]@SampPars$Obs$CAA_nsamp <- rep(400, 100)

  saveRDS(MHist,file=paste0("MOM/jurel_Hist",i,".rds"))
}

#saveRDS(MHist,file="MOM/jurel_Hist1.rds")
# Recruitment cases
cat("\n")
cat("\n","Reclutamientos altos : ",sprintf("%1.2f",median(Rec_high)/R0_OM),"\n")
cat("\n","Reclutamientos medios: ",sprintf("%1.2f",median(Rec_ori)/R0_OM),"\n")
cat("\n","Reclutamientos bajos : ",sprintf("%1.2f",median(Rec_low)/R0_OM),"\n")
cat("\n")


# Compare historical SSB ----

y_ind <- mod_current[[1]]$output[[1]]$SSB[,1] %in% mod_current[[1]]$output[[1]]$msy_mt[,1]
sb_jj <- mod_current[[1]]$output[[1]]$SSB[y_ind,2]
 
SB_OM <- MHist[[1]][[1]]@TSdata$SBiomass[1, , ] %>% rowSums()   # operating model
SB_AD <- mod_current[[1]]$output$Stock_1$Stock_Rec[12:66,2]    # JJM-12:66

png("figuras/comparison_jjm_sb_mom.png",height=3,width=4,res=400,units="in")
par(mar=c(2.9,2.8,1,1),mgp=c(1.8,0.5,0))
plot(years,sb_jj,xlab="Años",ylab="Biomasa desovante",typ="o",ylim=c(0, 28000)) # ADMB output
graphics::grid()
lines(years,SB_OM,col="red") # OM output
legend("topleft",c("JJM","MOM"),col=1:2,pch=c(1,NA),lty=1,cex=0.8)
dev.off()


# Compare historical abundance at age ----

NOM1 <- MHist[[1]][[1]]@AtAge$Number[1,2:13 , , ] %>%
  apply(1:2, sum) %>%
  structure(dimnames = list(Age = 1:12, Year = years)) %>%
  reshape2::melt() %>%
  mutate(type = "MOM")
NAD1 <- mod_current[[1]]$output$Stock_1$N[,2:13] %>%
  t() %>%
  structure(dimnames = list(Age = 1:12, Year = years)) %>%
  reshape2::melt() %>%
  mutate(type = "JJM")

g <- ggplot(rbind(NOM1,NAD1), aes(Year, value, colour = type)) +
  facet_wrap(vars(Age), scales = "free_y") +
  geom_point() +
  geom_line() +
  scale_x_continuous(limits=c(1970,2024),breaks=seq(1970,2024,by=13),
  labels=sprintf("%1.0f",seq(1970,2024,by=13))) +
  scale_color_manual(name="Tipo",values=c("#FF9F00","#00B7EB"),labels=c("MOM","JJM")) +
  labs(x = "Años", y = "Abundancia")+
  theme(plot.title=element_text(family="Calibri",color="black",size=21,face="bold"),
  legend.position="top",
  legend.text=element_text(size=9))
ggsave(paste0("figuras/comparison_jjm_age_mom.png"),g,height=5,width=8)

NOM0age <- MHist[[1]][[1]]@AtAge$Number[1,1:13 , , ] %>% apply(1:2, sum)

NAD1age <- mod_current[[1]]$output$Stock_1$N[,2:13] %>% t() 
mat <- t(rbind(NOM0age[1,],NAD1age[1,],years))

png("figuras/comparison_age_0.png",height=3,width=4,res=400,units="in")
par(mar=c(2.9,2.8,1,1),mgp=c(1.8,0.5,0))
plot(mat[1:54,3],mat[1:54,1],xlab="Años",ylab="Abundancia edad 0",
typ="o",ylim=c(0,58000),xlim=c(1970,2024)) # 
graphics::grid()
lines(mat[1:54,3],mat[2:55,2],col="red") # OM output
legend("topright",c("JJM","MOM"),col=1:2,pch=c(1,NA),lty=1,cex=0.8)
dev.off()

TBOM <- 2*MHist$`Stock 1`$`Fleet 1`@TSdata$Biomass[1,,1]
TBAD <- mod_current$h1_1.07$output$Stock_1$TotBiom[,2]

png("figuras/comparison_total_biomass.png",height=3,width=4,res=400,units="in")
par(mar=c(2.9,2.8,1,1),mgp=c(1.8,0.5,0))
plot(years,TBAD,xlab="Años",ylab="Biomasa total",typ="o",ylim=c(0, 39000)) # 
graphics::grid()
lines(years,TBOM,col="red") # OM output
legend("topleft",c("JJM","MOM"),col=1:2,pch=c(1,NA),lty=1,cex=0.8)
dev.off()


# Plot HCR
source("99-functions-mmp.R")

CBA_max <- 1200
Index <- seq(0,2100,50)
RC2_20mt <- sapply(Index,ramp_control_rule,c(500,1300),c(270,2000))
RC2_15mt <- sapply(Index,ramp_control_rule,c(500,1300),c(270,1500))
RC5_20mt <- empirical_HCR3(Index,Cmin=270,Cref=2000,Ilim=500,Iref=1300)
RC5_15mt <- empirical_HCR3(Index,Cmin=270,Cref=1500,Ilim=500,Iref=1300)
RC4_20mt <- empirical_HCR4(Index,Cmin=270,Cref=2000,Ilim=500,Iref=1300)
RC4_15mt <- empirical_HCR4(Index,Cmin=270,Cref=1500,Ilim=500,Iref=1300)

#
#library(svglite)
#name148<-paste(getwd(),"/figuras/RCC_all.svg",sep="")
#svglite(file=name148,width=13,height=10)
png("figuras/RCC_all.png",width=10,height=8,res=400,units="in")
layout(matrix(c(1:3,3),ncol=2,byrow=T),widths=c(1,1),heights=c(1,1))
#par(mar=c(4.2,4.5,2.9,2.5),mgp=c(2.5,0.7,0),cex.axis=1.7,cex.lab=2.0)
par(mar=c(4.2,4.2,2.5,1.5),mgp=c(2.5,0.7,0),cex.axis=1.5,cex.lab=1.5)
plot(Index,RC2_20mt,typ='l',panel.first = grid(),col=0,
xlab="Índice acústico (I)",ylab="CBA",ylim=c(0,2600),xlim=c(0,2100),lwd=2)
lines(Index,RC2_20mt,col="red",lwd=2)
lines(Index,RC2_15mt,col="blue",lwd=2)
abline(v=1300,lty=3,col=8)
abline(v=500,lty=3,col=8)
text(1300,3030,expression(I[paste("ref")]),cex=1.4,pos=1,xpd=T)
text(500,3030,expression(I[paste("lim")]),cex=1.4,pos=1,xpd=T)
legend("bottomright",c("Empírica RCC2 20mt","Empírica RCC2 15mt"),
col=c("red","blue"),lwd=c(2,2),lty=c(1,1),bty="n",cex=1.3)
legend("topleft",legend="a)",bty="n",cex=1.4)
#
plot(Index,RC5_20mt,typ='l',panel.first = grid(),col=0,
xlab="Índice acústico (I)",ylab="CBA",ylim=c(0,2600),xlim=c(0,2100),lwd=2)
lines(Index,RC5_20mt,col="red",lwd=2)
lines(Index,RC5_15mt,col="blue",lwd=2)
abline(v=1300,lty=3,col=8)
abline(v=500,lty=3,col=8)
text(1300,3030,expression(I[paste("ref")]),cex=1.4,pos=1,xpd=T)
text(500,3030,expression(I[paste("lim")]),cex=1.4,pos=1,xpd=T)
legend("bottomright",c("Empírica RCC3 20mt","Empírica RCC3 15mt"),
col=c("red","blue"),lwd=c(2,2),lty=c(1,1),bty="n",cex=1.3)
legend("topleft",legend="b)",bty="n",cex=1.4)
#
par(mar=c(4.2,20,2.5,15),mgp=c(2.5,0.7,0),cex.axis=1.5,cex.lab=1.5)
plot(Index,RC4_20mt,typ='l',panel.first = grid(),col=0,
xlab="Índice acústico (I)",ylab="CBA",ylim=c(0,2600),xlim=c(0,2100),lwd=2)
lines(Index,RC4_20mt,col="red",lwd=2)
lines(Index,RC4_15mt,col="blue",lwd=2)
abline(v=1300,lty=3,col=8)
abline(v=500,lty=3,col=8)
text(1300,3030,expression(I[paste("ref")]),cex=1.4,pos=1,xpd=T)
text(500,3030,expression(I[paste("lim")]),cex=1.4,pos=1,xpd=T)
legend("bottomright",c("Empírica RCC4 20mt","Empírica RCC4 15mt"),
col=c("red","blue"),lwd=c(2,2),lty=c(1,1),bty="n",cex=1.3)
legend("topleft",legend="c)",bty="n",cex=1.4)
#
dev.off()

# Plot survey index north Chile and SSB JJM
BT_AC<-as.numeric(mod_current[[1]]$data$Index[,2])
dat.ind<-as.data.frame(cbind(years,BT_AC,SB_AD))
names(dat.ind)<-c("years","survey_north","ssb_jjm")
a0<-which(dat.ind$years>1999)
dat.end<-dat.ind[a0,]
a1<-which(is.na(dat.end$survey_north)=="FALSE")
ind.dat<-dat.end[a1,]

# Regression model
mod.ind<-lm(survey_north~ssb_jjm,data=ind.dat)
eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(mod.ind)[1]), digits = 2),
              b = format(unname(coef(mod.ind)[2]), digits = 2),
             r2 = format(summary(mod.ind)$r.squared, digits = 2)))
dd<-dim(ind.dat)
ff<-seq(1,dd[1],by=2)
png("figuras/survey_ssb.png",width=4,height=3,res=400,units = "in")
par(mar=c(2.9,2.8,1.5,1.5),mgp=c(1.8,0.5,0),cex.axis=0.8,cex.lab=0.9)
plot(ind.dat$ssb_jjm,ind.dat$survey_north,typ='p',panel.first = grid(),
     ylab="Índice acústico norte",xlab="Biomasa desovante JJM",col=1,pch=21,
     ylim=c(0,3000),xlim=c(0,20000),bg=4,cex=0.8)
for(i in ff){
  text(ind.dat$ssb_jjm[i],ind.dat$survey_north[i],label=ind.dat$years[i],pos=4,cex=0.5,col=4,font=2)
}
abline(a=mod.ind$coefficients[1],b=mod.ind$coefficients[2])
text(6300,2800,label=as.expression(eq),cex=0.8,parse=TRUE)
dev.off()



cat("\n")
cat("\n","<>< HASTA AQUI TODO BIEN <><","\n")
cat("\n")




