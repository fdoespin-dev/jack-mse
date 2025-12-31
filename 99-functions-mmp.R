# PM
Sin_Captura <- function(x, DataList, ...) {
  np <- length(DataList)
  nf <- length(DataList[[1]])

  RecList <- lapply(1:np, function(p) {
    lapply(1:nf, function(f) {
      rec <- new("Rec")
      rec@TAC <- 1e-15
      rec
    })
  })
  return(RecList)
}
class(Sin_Captura) <- "MMP"

# Perfect management is an effort-based MP
# MP specifies the fishing effort that corresponds to F40%
# In openMSE, the effort is relative to the last historical year of the operating model.
# Effort can also be more efficient in the projection
# This code ensures that the correct fishing mortality is specified regardless of the
# improvement in the fishing effort, compare with MSEtool::FMSYref
Manejo_Perfecto <- function(x, DataList, reps = 1, FMSY, ...) {

  np <- length(DataList)
  nf <- length(DataList[[1]])

  # F at age by fleet
  F_hist_age <- sapply(1:nf, function(f) {
    p <- 1
    Data <- DataList[[p]][[f]]
    nyears <- length(Data@Misc$FleetPars$Find[x, ])
    FM <- Data@Misc$FleetPars$qs[x] * Data@OM$FinF[x]
    V <- Data@Misc$FleetPars$V[x, , nyears]
    V * FM
  })

  # Apical F by fleet
  F_hist <- apply(F_hist_age, 2, max)

  # Apical F (population)
  F_hist_tot <- max(rowSums(F_hist_age))

  # Ratio of fleet F to population F
  F_hist_ratio <- F_hist/F_hist_tot

  # Population FMSY
  #FMSY <- DataList[[1]][[1]]@OM$FMSY[1]

  # Fleet F corresponding to population FMSY, preserving the ratio of F among fleets
  F_proj <- F_hist_ratio * FMSY

  RecList <- lapply(1:np, function(p) {
    lapply(1:nf, function(f) {
      Data <- DataList[[p]][[f]]

      y <- max(Data@Year) - Data@LHYear+1
      nyears <- length(Data@Misc$FleetPars$Find[x, ])

      q <- Data@Misc$FleetPars$qs[x]
      qvar <- Data@Misc$FleetPars$qvar[x,y] # future only, always equal 1
      if (length(qvar)<1) qvar <- 1
      qinc <- Data@Misc$FleetPars$qinc[x] # future only, equal zero
      qcur <- qvar * q*(1+qinc/100)^y # catchability this year always equal 1

      HistE <- Data@OM$FinF[x] # Last historical fishing effort
      MSYE <- F_proj[f]/qcur # effort for this year's FMSY

      Rec <- new('Rec')
      #Rec@Effort <- MSYE/F_hist[f]   # Is / or * OR only MSYE withot F_hist[f]
      Rec@Effort <- MSYE/HistE
      Rec
    })
  })
  return(RecList)
}
class(Manejo_Perfecto) <- "MMP"


# Acustic Biomass reference 1.27 mm ton and Cref = 1.0 mm ton
# Rampa Control Rule
# JJM bio <- seq(0,2000,by=50)
# JJM CBA <- ramp_control_rule(bio,control_point=c(500,1300),relF=c(200,1100))

ramp_control_rule <- function(B_B0,control_point=c(0.1, 0.2, 0.2),relF=c(0.5, 0.8, 1)) {
  ctl <- c(0, control_point, Inf)
  relF <- c(relF[1], relF, relF[length(relF)])
  ind <- findInterval(B_B0, ctl)
  slope <- (relF[ind + 1] - relF[ind])/(ctl[ind + 1] - ctl[ind]) 
  x_diff <- B_B0 - ctl[ind]
  alpha <- slope * x_diff + relF[ind]
  return(alpha)
}

# Empirical rules (jack)  c(0.0835,0.676,0.186,0.054)
# data_lag = 1 means one year ago, data_lag = 2 two years ago
# y_mean = 1 means last year assessment, y_mean = 3 means 3 years ago assessment  
Emp_RCC2_15mt <- function(x, DataList, reps = 1,frac = c(0.0835,0.676,0.186,0.054),
                          y_mean = 1, CBA_max = 1500, banda = c(0.20,0.10), data_lag = 1) {

  np <- length(DataList)
  nf <- length(DataList[[1]])

  Data <- DataList[[1]][[1]]
  ny <- length(Data@Year)

  # Create blank object for catch recommendation
  RecList <- lapply(1:np, function(p) {
    lapply(1:nf, function(f) {
      Rec <- new("Rec")
      Rec@TAC <- NA_real_
      return(Rec)
    })
  })

  # Retrieve the previous CBA
  if (max(Data@Year) == Data@LHYear) {
  # If we are at the beginning of the projection, we specify the real 2023 CBA
    CBA_previous <- 1200
  } else {
    CBA_previous <- as.numeric(Data@MPrec[x])
    for(f in 2:nf){
      CBA_previous <- CBA_previous + as.numeric(DataList[[1]][[f]]@MPrec[x])
    }
  }

  yind <- seq(ny - y_mean + 1, ny) - data_lag

  #Ind <- mean(Data@AddInd[x, 1, yind], na.rm = TRUE)
  Ind <- Data@AddInd[x, 1, yind]

  #RCC with clasic rampa rule
  CBA_now <- empirical_HCR2(Ind,Cmin=270,Cref=CBA_max,Ilim=500,Iref=1300)

  # Apply stabilization range
  if (CBA_now/CBA_previous > 1+banda[1]) {CBA_now <- CBA_previous*(1+banda[1])}  
  if (CBA_now/CBA_previous < 1-banda[2]) {CBA_now <- CBA_previous*(1-banda[2])}  

  for (f in 1:nf) {
    RecList[[1]][[f]]@TAC <- frac[f] * CBA_now
  }
  return(RecList)
}
class(Emp_RCC2_15mt) <- "MMP"

Emp_RCC2_20mt <- Emp_RCC2_15mt

formals(Emp_RCC2_20mt)$CBA_max <- 2000
class(Emp_RCC2_20mt) <- "MMP"


Emp_RCC3_15mt <- function(x, DataList, reps = 1,frac = c(0.0835,0.676,0.186,0.054),
                          y_mean = 1, CBA_max = 1500, banda = c(0.20,0.10), data_lag = 1) {

  np <- length(DataList)
  nf <- length(DataList[[1]])

  Data <- DataList[[1]][[1]]
  ny <- length(Data@Year)

  # Create blank object for catch recommendation
  RecList <- lapply(1:np, function(p) {
    lapply(1:nf, function(f) {
      Rec <- new("Rec")
      Rec@TAC <- NA_real_
      return(Rec)
    })
  })

  # Retrieve the previous CBA
  if (max(Data@Year) == Data@LHYear) {
  # If we are at the beginning of the projection, we specify the real 2023 CBA
    CBA_previous <- 1200
  } else {
    CBA_previous <- as.numeric(Data@MPrec[x])
    for(f in 2:nf){
      CBA_previous <- CBA_previous + as.numeric(DataList[[1]][[f]]@MPrec[x])
    }
  }

  yind <- seq(ny - y_mean + 1, ny) - data_lag

  #Ind <- mean(Data@AddInd[x, 1, yind], na.rm = TRUE)
  Ind <- Data@AddInd[x, 1, yind]

  # RCC with rampa ilimited and follow Cref up limit
  CBA_now <- empirical_HCR3(Ind,Cmin=270,Cref=CBA_max,Ilim=500,Iref=1300)

  # Apply stabilization range
  if (CBA_now/CBA_previous > 1+banda[1]) {CBA_now <- CBA_previous*(1+banda[1])}  
  if (CBA_now/CBA_previous < 1-banda[2]) {CBA_now <- CBA_previous*(1-banda[2])}  

  for (f in 1:nf) {
    RecList[[1]][[f]]@TAC <- frac[f] * CBA_now
  }
  return(RecList)
}
class(Emp_RCC3_15mt) <- "MMP"

Emp_RCC3_20mt <- Emp_RCC3_15mt

formals(Emp_RCC3_20mt)$CBA_max <- 2000
class(Emp_RCC3_20mt) <- "MMP"


Emp_RCC4_15mt <- function(x, DataList, reps = 1,frac = c(0.0835,0.676,0.186,0.054),
                          y_mean = 1, CBA_max = 1500, banda = c(0.20,0.10), data_lag = 1) {

  np <- length(DataList)
  nf <- length(DataList[[1]])

  Data <- DataList[[1]][[1]]
  ny <- length(Data@Year)

  # Create blank object for catch recommendation
  RecList <- lapply(1:np, function(p) {
    lapply(1:nf, function(f) {
      Rec <- new("Rec")
      Rec@TAC <- NA_real_
      return(Rec)
    })
  })

  # Retrieve the previous CBA
  if (max(Data@Year) == Data@LHYear) {
  # If we are at the beginning of the projection, we specify the real 2023 CBA
    CBA_previous <- 1200
  } else {
    CBA_previous <- as.numeric(Data@MPrec[x])
    for(f in 2:nf){
      CBA_previous <- CBA_previous + as.numeric(DataList[[1]][[f]]@MPrec[x])
    }
  }

  yind <- seq(ny - y_mean + 1, ny) - data_lag

  #Ind <- mean(Data@AddInd[x, 1, yind], na.rm = TRUE)
  Ind <- Data@AddInd[x, 1, yind]

  # RCC with rampa exponential
  CBA_now <- empirical_HCR4(Ind,Iref=1300,Cref=CBA_max,Cmin=270,Ilim=500)

  # Apply stabilization range
  if (CBA_now/CBA_previous > 1+banda[1]) {CBA_now <- CBA_previous*(1+banda[1])}  
  if (CBA_now/CBA_previous < 1-banda[2]) {CBA_now <- CBA_previous*(1-banda[2])}  

  for (f in 1:nf) {
    RecList[[1]][[f]]@TAC <- frac[f] * CBA_now
  }
  return(RecList)
}
class(Emp_RCC4_15mt) <- "MMP"

Emp_RCC4_20mt <- Emp_RCC4_15mt

formals(Emp_RCC4_20mt)$CBA_max <- 2000
class(Emp_RCC4_20mt) <- "MMP"


# RCC with rampa in Iref and Ilim and Cref-Cmin
empirical_HCR2 <- function(I_y, Iref = Iref, Cref = Cref, Cmin = Cmin, Ilim = Ilim) {
  n <- length(I_y)
  cba<-vector("numeric",length=n)
  slope <- (Cref - Cmin)/(Iref - Ilim)
  for(i in 1:n){
    val <- I_y[i]
    if(val < Iref & val > Ilim){cba[i] <- slope * (val - Ilim) + Cmin}
    if(val <= Ilim){cba[i]<-Cmin}
    if(val >= Iref){cba[i]<-Cref}

  }
  return(cba)
}


# RCC with rampa exponential
empirical_HCR4 <- function(I_y, Iref = Iref, Cref = Cref, Cmin = Cmin, Ilim = Ilim){
  n <- length(I_y)
  cba <- vector("numeric",length=n)
  val_lim <- Ilim/Iref
  if(Cref==1500){ax1<-0.114;ax2<-0.043}
  if(Cref==2000){ax1<-0.055;ax2<-0.023}
  for (t in 1:n){
    val <- I_y[t]/Iref
    if(val <= val_lim){cba[t] <- 270}
    if(val > val_lim & val < 1){cba[t] <- Cref*(ax1+(val-ax2)^2.5)}
    if(val >= 1){cba[t] <- Cref*(val)^0.6}
  }
  return(cba)
}

# RCC with rampa ilimited and follow Cref up limit
empirical_HCR3 <- function(I_y, Cmin = Cmin, Cref = Cref, Ilim = Ilim, Iref = Iref) {
  n <- length(I_y)
  cba<-vector("numeric",length=n)
  slope <- (Cref - Cmin)/(Iref - Ilim)
  for(i in 1:n){
    val <- I_y[i]
    if(val > Ilim){cba[i] <- slope * (val - Ilim) + Cmin}
    if(val <= Ilim){cba[i] <- Cmin}
  }
  return(cba)
}




