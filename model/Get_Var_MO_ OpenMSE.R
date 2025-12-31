#  Ignacio Payá
#  Abril 30 2025

#install.packages("devtools")
#devtools::install_github("SPRFMO/jjmR")
#library(jjmR)
#wd=setwd('C:/Users/ignacio.paya/OneDrive - Instituto de Fomento Pesquero/JUREL/FIP_OPEN_MSE/MdelSUR_workshop-southern-hake/ADMB/jurel')
#getwd()
# need subfolders: config, input, results

# Chequeo
#Run without estimation and with pdf file
#jjmR::runit(mod="h1_1.07",est=FALSE,exec=NULL,path="config",input="input",output="results",
#version="2015MS",pdf=TRUE,portrait=TRUE)
#mod_current <- jjmR::readJJM("h1_1.07.ls", path = "config", input = "input")
mod_current <- jjmR::readJJM("h1_1.07",path="config",input="input")

#jjmR::kobe(mod_current)

# Find variables
#names(mod_current)
#names(mod_current[[1]])
#mod_current[[1]]$info
#mod_current[[1]]$data
#mod_current[[1]]$control$modelName
#mod_current[[1]]$output
# outputs
#names(mod_current[[1]]$output$Stock_1)

break
## Table function
library(tidyverse)
tabla.resumen <- function(modelo) {
  y_ind <- modelo[[1]]$output[[1]]$SSB[,1] %in% modelo[[1]]$output[[1]]$msy_mt[,1]
  #table_qtt <- data.frame(yr  = sc08.01[[1]]$output[[1]]$msy_mt[,1]) %>%
  table_qtt <- data.frame(yr  = modelo[[1]]$output[[1]]$msy_mt[,1]) %>%
    mutate(
      #yflo      = round(as.vector(rowSums(modelo[[1]]$data$Fcaton)), 0),
      ssb       = round(modelo[[1]]$output[[1]]$SSB[y_ind,2], 0),
      r         = round(modelo[[1]]$output[[1]]$R[,2], 0),
      ssb_msy   = round(modelo[[1]]$output[[1]]$msy_mt[,10], 0),
      fflo      = round(modelo[[1]]$output[[1]]$msy_mt[,6], 2),
      fmsy      = round(modelo[[1]]$output[[1]]$msy_mt[,5], 2),
      ssb0      = round(modelo[[1]]$output[[1]]$msy_mt[,11], 2),
      r_fmsy    = round(modelo[[1]]$output[[1]]$msy_mt[,4], 2),
      r_ssbmsy  = round(modelo[[1]]$output[[1]]$msy_mt[,13], 2),
      msy       = round(modelo[[1]]$output[[1]]$msy_mt[,8], 2),
      msyL      = round(modelo[[1]]$output[[1]]$msy_mt[,9], 2),
      spr_ratio = round(modelo[[1]]$output[[1]]$msy_mt[,7], 2),
      spr_ft    = round(modelo[[1]]$output[[1]]$msy_mt[,2], 2)
    )
  return(table_qtt)
}

# Table
x=data.frame(tabla.resumen(mod_current))
#
library("writexl")
write_xlsx(x, path=paste0(getwd(),"/Results/",paste0(mod_current[[1]]$control$modelName,"_Results.xlsx")))
## desde jjmR
get_msy_mt(mod_current)
#### Biologicas: , , ,
# mortalidad natural
M=mod_current[[1]]$output$Stock_1$M
# ojiva de madurez a la edad
mature_a=mod_current[[1]]$output$Stock_1$mature_a
plot(mature_a)
# pesos a la edad en la población
wt_a_pop=mod_current[[1]]$output$Stock_1$wt_a_pop
plot(wt_a_pop)
# longitud a la edad ??????

# pesqueras:
# año + mortalidad por pesca de las flotas
F_age_1=mod_current[[1]]$output$Stock_1$F_age_1
matplot(F_age_1[20:30,2:13])

F_age_2=mod_current[[1]]$output$Stock_1$F_age_2
matplot(F_age_2[1,2:13])

F_age_3=mod_current[[1]]$output$Stock_1$F_age_3
matplot(F_age_3[1,2:13])

F_age_4=mod_current[[1]]$output$Stock_1$F_age_4
matplot(F_age_4[1,2:13])

# stock + año + selectividades
sel_fsh_1=mod_current[[1]]$output$Stock_1$sel_fsh_1
matplot(sel_fsh_1[30:35,3:14])

sel_fsh_2=mod_current[[1]]$output$Stock_1$sel_fsh_2
matplot(sel_fsh_2[30,3:14])

sel_fsh_3=mod_current[[1]]$output$Stock_1$sel_fsh_3
matplot(sel_fsh_3[30,3:14])

sel_fsh_4=mod_current[[1]]$output$Stock_1$sel_fsh_4
matplot(sel_fsh_1[30,3:14])

#abundancia población
N=mod_current[[1]]$output$Stock_1$N
matplot(N[30:35,2:13])

# Año + índices
mod_current[[1]]$output$Stock_1$Index_names[1]
Obs_Survey_1=mod_current[[1]]$output$Stock_1$Obs_Survey_1
plot(Obs_Survey_1[,1:2])
matplot(Obs_Survey_1[,1],Obs_Survey_1[,2:4])

# Acustica Norte

Obs_Survey_2=mod_current[[1]]$output$Stock_1$Obs_Survey_2
matplot(Obs_Survey_2[,1],Obs_Survey_2[,2:4])

Obs_Survey_3=mod_current[[1]]$output$Stock_1$Obs_Survey_3
matplot(Obs_Survey_3[,1],Obs_Survey_3[,2:4])

Obs_Survey_4=mod_current[[1]]$output$Stock_1$Obs_Survey_4
matplot(Obs_Survey_4[,1],Obs_Survey_4[,2:4])

Obs_Survey_5=mod_current[[1]]$output$Stock_1$Obs_Survey_5
matplot(Obs_Survey_5[,1],Obs_Survey_5[,2:4])

Obs_Survey_6=mod_current[[1]]$output$Stock_1$Obs_Survey_6
matplot(Obs_Survey_6[,1],Obs_Survey_6[,2:4])

Obs_Survey_7=mod_current[[1]]$output$Stock_1$Obs_Survey_7
matplot(Obs_Survey_7[,1],Obs_Survey_7[,2:4])

#  estructuras
mod_current[[1]]$output$Stock_1$Index_names

# las series de captura
Obs_Survey_7=mod_current[[1]]$output$Stock_1$
  plot(Obs_Survey_7[,1:2])



#y estructuras de longitud por flota

wt_fsh_1=mod_current[[1]]$output$Stock_1$wt_fsh_1
matplot(wt_fsh_1[1,2:13])
wt_fsh_2=mod_current[[1]]$output$Stock_1$wt_fsh_2
matplot(wt_fsh_2[1,2:13])



