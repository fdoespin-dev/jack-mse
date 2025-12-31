options(warn=-1)
library(MSEtool)
library(dplyr)

start_time<-Sys.time()

# Run simulation ----
do_parallel <- TRUE

if (do_parallel) {
  ncores <- 10
  if (.Platform$OS.type == "unix") {
    source("99-sfInitFork.R")
    sfInitFork(TRUE, ncores)
  } else {
    MSEtool::setup(ncores)
  }
}

# Define management procedures ----
source("99-functions-mmp.R")

# Declare MP interval application ----
MP_test <- c("Emp_RCC2_15mt", "Emp_RCC2_20mt", "Emp_RCC3_15mt","Emp_RCC3_20mt",
             "Emp_RCC4_15mt", "Emp_RCC4_20mt")

MP_interval <- ifelse(grepl("Sin|Manejo", MP_test), 200, 1) %>%
  structure(names = MP_test)

# Run projection ----
if (do_parallel) {

  sfLibrary(dplyr)
  sfExport(list = c("empirical_HCR2","empirical_HCR3","empirical_HCR4"))

  parallel <- structure(as.list(grepl("actual|Emp", MP_test)), names = MP_test)
} else {
  parallel <- FALSE
}


# Run projections for i MOM ----
for (i in 1:3) {

  MOM <- readRDS(file.path("MOM", paste0("jurel_MOM", i, ".rds")))
  #MOM@nsim <- 20 # Set in 01 make mom 4 fleet point R source code

  Hist <- MSEtool::SimulateMOM(MOM)
  nyears <- length(Hist[[1]][[1]]@Data@Year)
  FMSY <- Hist[[1]][[1]]@Ref$ByYear$FMSY[1, nyears]
  formals(Manejo_Perfecto)$FMSY <- FMSY
  class(Manejo_Perfecto) <- "MMP"

  Hist[[1]][[1]]@Misc$MOM@interval <- MP_interval

  tictoc::tic()
  MMSE <- ProjectMOM(Hist, MPs = MP_test, checkMPs = FALSE, parallel = parallel)
  tictoc::toc()

  saveRDS(MMSE, file = file.path("MOM", paste0("jurel_MMSE", i, ".rds")))
  rm(Hist, MMSE)
}


if (do_parallel) sfStop()

#
end_time<-Sys.time()
diff_time<-unclass(end_time-start_time)
cat("\n")
print(paste0("Tiempo transcurrido = ",round(diff_time,4)," ",attributes(diff_time)$units))





