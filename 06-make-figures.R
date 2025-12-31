
source("99-functions-proj-figures.R")

library(dplyr)
library(openMSE)


flag<-1

OM_name <- c("(1) Reclutamientos históricos","(2) Reclutamientos bajos",
              "(3) Reclutamientos altos")

ini_mom <- 1
fin_mom <- 3 

MMSE_list <- lapply(c(ini_mom:fin_mom),function(i) {
  mmse <- readRDS(file = paste0("MOM/jurel_MMSE", i, ".rds"))

  hist <- readRDS(file = paste0("MOM/jurel_Hist", i, ".rds"))

  # Append FMSY to mmse object
  nyears <- length(hist[[1]][[1]]@Data@Year)
  mmse@RefPoint$ByYear$FMSY <- hist[[1]][[1]]@Ref$ByYear$FMSY[1,nyears]

  # Append historical object
  mmse@multiHist <- hist

  return(mmse)
})


proyears <- MMSE_list[[1]]@proyears

# OM names
#OM_name <- c("4 flota")

## Generate time series figures ----

# SSB ----
g <- plot_array(MMSE_list,OM_name,type="SB") +
  coord_cartesian(xlim=c(2010,2024+proyears),ylim=c(0,40000),expand=FALSE)
ggsave("figuras/projection_SB.png",g,height=4.1,width=6.3)


# B/BMSY ----
g <- plot_array(MMSE_list,OM_name,type="B_BMSY") +
geom_hline(yintercept=1.0,linetype=3) +
coord_cartesian(xlim=c(2010,2024+proyears),ylim=c(0,4.2),expand = FALSE)
ggsave("figuras/projection_BMSY.png",g,height=4.1,width=6.3)


# F/FMSY ----
g <- plot_array(MMSE_list,OM_name,type="F_FMSY") +
  geom_hline(yintercept=1.0,linetype= 3) +
  coord_cartesian(xlim=c(2010,2024+proyears),ylim=c(0,2.5),expand=FALSE)
ggsave("figuras/projection_FMSY.png",g,height=4.1,width=6.3)


# Catch ----
g <- plot_array(MMSE_list,OM_name,type="Catch") +
  coord_cartesian(xlim=c(2010,2024+proyears),ylim=c(0,4000),expand=FALSE)
ggsave("figuras/projection_Catch.png",g,height=4.1,width=6.3)


# Recruitments ----
g <- plot_array(MMSE_list,OM_name,type="R") +
  coord_cartesian(xlim=c(2010,2024+proyears),ylim=c(0,35000),expand=FALSE)
ggsave("figuras/projection_recruit.png",g,height=4.1,width=6.3)



# Acoustic Index

g <- plot_index(MMSE_list[[1]],OM_name,type="Acustico Chile") + # type=c("MPH","Acustico Chile","Acustico Peru")
  coord_cartesian(xlim=c(2006, 2024+proyears),ylim=c(0,4200),expand=FALSE)
ggsave("figuras/projection_index.png",g,height=4.1,width=7.5)



# Kobe time plot biomass averaged across the 4 reference OMs ----
Kavg <- lapply(ini_mom:fin_mom,function(x) {
  g <- plot_Kobe_time(MMSE_list[[x]],type="B")
  g$data %>% mutate(OM=paste("OM",x))
}) %>% bind_rows() %>%
  group_by(MP,Year,name) %>%
  summarise(value=mean(value))
g <- plot_Kobe_time(output=Kavg) +
  geom_hline(yintercept=0.5,linetype=2,col="grey40") +
  geom_hline(yintercept=0.25,linetype=3,col="grey40") +
  #geom_vline(xintercept=2022+c(2,7)+0.75,linetype=3,col="grey40") +
  ggtitle("MO promedio")
ggsave("figuras/projection_kobe_avg.png",g,height=4.5,width=6.5)



# Kobe time plot biomass averaged across the 4 reference OMs ----
Favg <- lapply(ini_mom:fin_mom, function(x) {
  g <- plot_Kobe_time(MMSE_list[[x]],type="F")
  g$data %>% mutate(OM = paste("OM",x))
}) %>% bind_rows() %>%
  group_by(MP, Year, name) %>%
  summarise(value = mean(value))
g <- plot_Kobe_time(output=Favg,type="F") +
  geom_hline(yintercept=0.5,linetype=2,col="grey40") +
  #geom_vline(xintercept = 2021 + c(5, 12, 24), linetype = 3, col = "grey40") +
  ggtitle("MO promedio")
ggsave("figuras/projection_kobe_F_avg.png",g,height=4.5,width=6.5)




cat("\n")
cat("\n","<>< HASTA AQUI TODO BIEN <><","\n")
cat("\n")


