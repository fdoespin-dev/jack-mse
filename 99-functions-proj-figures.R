
library(tidyverse)

.get_array <- function(.mse, .names, type = c("SB", "B_BMSY", "F_FMSY", "Catch","R"),annual=FALSE) {
  type <- match.arg(type)

  nyears <- .mse@nyears
  proyears <- .mse@proyears
  fin_year <- 2024
  nf <- .mse@nfleets
  Year <- fin_year + seq(1,.mse@proyears)
  Year_hist <- fin_year - seq(.mse@nyears,1) + 1 
  x <- switch(
    type,
    "SB" =  .mse@SSB[, 1, ,],
    "B_BMSY" = .mse@SB_SBMSY[ ,1 , ,],
    "F_FMSY" = local({
      # Need to calculate F-at-age since selectivity varies by fleet
      V_xaf <- sapply(1:nf, function(f) {
        .mse@multiHist[[1]][[f]]@SampPars$Fleet$V[, , .mse@nyears + 1]
      }, simplify = "array")
      FM_xfmy <- .mse@FM[, 1, , , ]

      FM_xafmy <- array(NA,c(.mse@nsim,dim(V_xaf)[2],nf,length(.mse@MPs[[1]]),.mse@proyears))
      ind <- as.matrix(
        expand.grid(x=1:.mse@nsim,a=1:dim(V_xaf)[2],f=1:nf,m=1:length(.mse@MPs[[1]]),y=1:.mse@proyears)
      )
      FM_xafmy[ind] <- FM_xfmy[ind[,-2]] * V_xaf[ind[,c(1:3)]]
      FM_xamy <- apply(FM_xafmy,c(1, 2, 4, 5),sum)  # Sum across fleets
      FM_xmy <- apply(FM_xamy,c(1, 3, 4),max)
      F_FMSY <- FM_xmy/.mse@RefPoint$ByYear$FMSY
      return(F_FMSY)
    }),
    "CBA" = apply(.mse@TAC[, 1, , ,],c(1,3,4),sum),
    "Catch" = apply(.mse@Catch[, 1, , ,],c(1,3,4),sum),
    "R" = apply(.mse@N[, 1, 2, , ,],c(1,2,3),sum),
    NULL
  )

  if (!is.null(x)) {
    x <- x %>%
      structure(dimnames = list(Simulation = 1:.mse@nsim, MP = unlist(.mse@MPs), Year = Year)) %>%
      reshape2::melt()
    #}
  }

  xhist <- switch(
    type,
    "SB" = apply(.mse@multiHist[[1]][[1]]@TSdata$SBiomass,1:2,sum),
    "B_BMSY" = apply(.mse@multiHist[[1]][[1]]@TSdata$SBiomass,1:2,sum)/.mse@RefPoint$ByYear$SSBMSY[, 1,1 ,1:.mse@nyears], 
    "F_FMSY" = sapply(1:nf, function(f) .mse@multiHist[[1]][[f]]@AtAge$F.Mortality, simplify = "array") %>%
      apply(c(1, 2, 3, 4), sum) %>% # Sum across fleet
      apply(c(1, 3), max) %>% # Maximum across age
      `/`(.mse@RefPoint$ByYear$FMSY),
    "Catch" = lapply(1:4, function(i) .mse@multiHist[[1]][[i]]@Data@Cat) %>% simplify2array() %>% apply(1:2,sum),
    "R" = apply(.mse@multiHist[[1]][[1]]@AtAge$Number[, 1, ,],1:2,sum),
    NULL
  )
  if (!is.null(xhist)) {
    xhist <- xhist %>%
      array(c(.mse@nsim, .mse@nyears, .mse@nMPs)) %>%
      aperm(c(1, 3, 2)) %>%
      structure(dimnames = list(Simulation = 1:.mse@nsim, MP = unlist(.mse@MPs), Year = Year_hist )) %>%
      reshape2::melt()
  }

  if (type == "Index") {
    index_list <- lapply(1:.mse@nMPs, function(i) {
      Data <- .mse@PPD[[i]]

      Data@AddInd[, 1, ] %>%
        structure(dimnames = list(Simulation = 1:.mse@nsim, Year = Data@Year)) %>%
        reshape2::melt() %>%
        mutate(MP = .mse@MPs[i])
    })

    x <- do.call(rbind, index_list) %>%
      dplyr::filter(!is.na(value)) %>%
      mutate(MP = factor(MP, levels = .mse@MPs))
  }

  output <- rbind(x, xhist) %>% dplyr::mutate(OM = .names)

  return(output)
}

sea_green1 <- rgb(46,139,87,alpha=70,maxColorValue=255)
sea_green2 <- rgb(162,255,203,alpha=80,maxColorValue=255)



plot_array <- function(MMSE, names, sims, type = c("SB","B_BMSY","F_FMSY","Catch","CBA","R"),annual=FALSE) {
  type <- match.arg(type)

  if (is.list(MMSE)) {
     array_df <- Map(.get_array, .mse = MMSE, .name = names, type = type, annual = annual) %>%
      dplyr::bind_rows()

  } else if (is(MMSE, "MMSE")) {
    array_df <- .get_array(MMSE, names, type, annual = annual)
  } else {
    stop("MMSE list or object not found.")
  }

  #
  array_hist <- filter(array_df,Year %in% seq(1970,2024,1)) %>% mutate(Simulation=1)

  #
  ylab <- switch(type,
                 "SB" = "Biomasa desovante",
                 "B_BMSY" = expression(BD/BD[RMS]),
                 "F_FMSY" = expression(F/F[RMS]),
                 "Catch" = "Captura",
                 "R" = "Reclutamientos",
                 "CBA" = "CBA")

  if (!missing(sims)) {

    array_sims <- filter(array_df, Simulation %in% sims) %>%
      mutate(Simulation = factor(Simulation))
    g <- ggplot(array_sims, aes(Year, value, group = Simulation, linetype = Simulation)) +
      facet_grid(vars(MP), vars(OM)) +
      geom_line() +
      #geom_hline(yintercept = 1, linetype = 3) +
      expand_limits(y = 0) +
      theme(panel.spacing = unit(0.2, "in"),
            legend.position = "bottom",
            legend.text=element_text(size=8.5),
            axis.text.x = element_text(angle = 0, vjust = 0.5),
            strip.background = element_blank()) +
      coord_cartesian(expand = FALSE) +
      labs(x = "Años", y = ylab, linetype = "Simulación")

  } else {

    array_band <- array_df %>%
      group_by(Year, MP, OM) %>%
      summarise(med = median(value),
                lwr = quantile(value, 0.05),
                upr = quantile(value, 0.95),
                lwr2 = quantile(value, 0.25),
                upr2 = quantile(value, 0.75))


    g <- ggplot(array_band, aes(Year, med)) +
      facet_wrap(vars(MP)) +
      geom_ribbon(alpha = 0.1, aes(ymin = lwr,  ymax = upr,  fill = OM)) +
      geom_ribbon(alpha = 0.2, aes(ymin = lwr2, ymax = upr2, fill = OM)) + # Plot only the interquartile range
      geom_line(aes(colour = OM)) + 
         expand_limits(y = 0) +
      geom_line(data=array_hist,aes(Year,value),linewidth=0.5,color=1) +
      theme(panel.spacing = unit(0.2, "in"),
            axis.text.x = element_text(angle = 0, vjust = 0.5),
            strip.background = element_blank()) +
      coord_cartesian(expand = FALSE) +
      labs(x = "Años", y = ylab, fill = "MO", colour = "MO") +
      theme(legend.position = "bottom",
            legend.text=element_text(size=8.5)) +
      guides(fill = guide_legend(ncol = 2),
             colour = guide_legend(ncol = 2))
  }

  if (any(array_df$Year < 2024)) { # 2022 is the historical period
    g <- g + geom_vline(xintercept = 2024, linetype = 2)
  }

  g
}


plot_index <- function(MMSE, names, type = c("MPH", "Acustico Chile", "Acustico Peru"), sims) {

  index_text <- switch(
    type,
    "MPH" = "MMSE@PPD[[1]][[2]][[i]]@SpInd",
    "Acustico Chile" = "MMSE@PPD[[1]][[1]][[i]]@AddInd[, 1, ]",
    "Acustico Peru" = "MMSE@PPD[[1]][[3]][[i]]@AddInd[, 1, ]"
  )

  index_list <- lapply(1:MMSE@nMPs, function(i) {
    
    Data <- MMSE@PPD[[1]][[1]][[i]]
    ind <- parse(text = index_text) %>% eval()

    ind %>%
      structure(dimnames = list(Simulation = 1:MMSE@nsim, Year = Data@Year)) %>%
      reshape2::melt() %>%
      mutate(MP = MMSE@MPs[[1]][i])
  })

  index_df <- do.call(rbind, index_list)

  if (!missing(sims)) {

    #if (length(sims) > 2) stop("This function only supports two simulations")

    array_sim <- index_df %>%
      mutate(OM = "OM 1") %>%  # One OM
      filter(Simulation %in% sims) %>%
      filter(!is.na(value)) %>%
      mutate(Simulation = factor(Simulation),
             MP = factor(MP, levels = unlist(MMSE@MPs)))

    g <- ggplot(array_sim, aes(Year, value, group = Simulation, linetype = Simulation)) +
      #facet_grid(vars(OM), vars(MP)) + #the first one
      facet_wrap(vars(MP)) +
      geom_line() +
      expand_limits(y = 0) +
      theme(panel.spacing = unit(0.2, "in"),
            legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust = 0.5),
            strip.background = element_blank()) +
      coord_cartesian(expand = FALSE) +
      labs(y = type, x = "Años") +
      geom_vline(xintercept = MMSE@OM[[1]][[1]]$CurrentYr[1], linetype = 2) # Many OM

  } else {
    array_band <- index_df %>%
      mutate(OM = "OM 1") %>% # One OM
      group_by(Year, MP, OM) %>%
      summarise(med = median(value, na.rm = TRUE),
                lwr = quantile(value, 0.05, na.rm = TRUE),
                upr = quantile(value, 0.95, na.rm = TRUE),
                lwr2 = quantile(value, 0.25, na.rm = TRUE),
                upr2 = quantile(value, 0.75, na.rm = TRUE)) %>%
      mutate(MP = factor(MP, levels = unlist(MMSE@MPs))) %>%
      filter(!is.na(med))

    g <- ggplot(array_band, aes(Year, med)) +
      #facet_grid(vars(OM), vars(MP)) + #the first one
      facet_wrap(vars(MP)) +
      geom_ribbon(fill = sea_green2, aes(ymin = lwr, ymax = upr)) +  #grey90
      geom_ribbon(fill = sea_green1, aes(ymin = lwr2, ymax = upr2)) + #grey70
      geom_line(color="darkgreen") +
      expand_limits(y = 0) +
      theme(panel.spacing = unit(0.2, "in"),
            axis.text.x = element_text(angle = 0, vjust = 0.5),
            strip.background = element_blank()) +
      coord_cartesian(expand = FALSE) +
      labs(y = type, x = "Años") +
      #geom_vline(xintercept = MSE@OM$CurrentYr[1], linetype = 2) # One OM
      geom_vline(xintercept = MMSE@OM[[1]][[1]]$CurrentYr[1], linetype = 2) # Many OM

  }

  g
}



plot_HCR <- function(B_B0, F_FMSY) {
  Brel <- seq(0, 1, length.out = 200)
  Frel <- HCRlin(Brel, 0.2, 0.4, 0.75, 1)
  plot(
    Brel, Frel,
    xlab = expression(B/B[0]),
    ylab = expression(F/F[MSY]),
    ylim = c(0, 1.25), type = "l",
    lwd = 2,
    panel.first = {
      abline(v = c(0.2, 0.4), col = "red", lty = 3)
      abline(h = c(0.75, 1), col = "red", lty = 3)
    }
  )
  if (!missing(B_B0) && !missing(F_FMSY)) {
    matpoints(B_B0, F_FMSY, pch = 1, typ = "p")
    legend(
      "bottomright",
      c("HCR", "Operating model"),
      lwd = c(2, 0),
      pch = c(NA, 1),
      col = c(1, 4)
    )
  }

  invisible()
}




plot_Kobe_time <- function(MMSE, output, type = c("B", "F")) {
  type <- match.arg(type)

  if (!missing(MMSE) && missing(output)) {
    nsim <- MMSE@nsim
    nyears <- MMSE@nyears
    nf <- MMSE@nfleets
    proyears <- MMSE@proyears
    nMPs <- length(unlist(MMSE@MPs)) 
    if (type == "B") {
     #B_BMSY <- .get_array(MMSE, .name = "OM", type = "B_BMSY") %>% bind_rows()
     B_BMSY = 2 * (MMSE@SSB[, 1, ,]/MMSE@RefPoint$ByYear$SSBMSY[, 1, ,(nyears+1):(nyears+proyears)])
     Kgreen <- structure(
       B_BMSY >= 1.05,
       dimnames = list(Simulation=1:MMSE@nsim,MP=unlist(MMSE@MPs),Year=2024+seq(1,MMSE@proyears))) %>%
       apply(2:3, mean) %>%
       reshape2::melt() %>%
       mutate(val = "green", name = "Subexplotado")
      Kdarkgreen <- structure(
        B_BMSY < 1.05 & B_BMSY >= 0.95,
        dimnames = list(Simulation=1:MMSE@nsim,MP=unlist(MMSE@MPs),Year=2024+seq(1,MMSE@proyears))) %>%
        apply(2:3, mean) %>%
        reshape2::melt() %>%
        mutate(val = "orange", name = "Plena explotado")
      Kyellow <- structure(
        B_BMSY < 0.95 & B_BMSY >= 0.5,
       dimnames = list(Simulation=1:MMSE@nsim,MP=unlist(MMSE@MPs),Year=2024+seq(1,MMSE@proyears))) %>%
        apply(2:3, mean) %>%
        reshape2::melt() %>%
        mutate(val = "yellow", name = "Sobreexplotado")
      Kred <- structure(
        B_BMSY < 0.5,
        dimnames = list(Simulation=1:MMSE@nsim,MP=unlist(MMSE@MPs),Year=2024+seq(1,MMSE@proyears))) %>%
        apply(2:3, mean) %>%
        reshape2::melt() %>% 
        mutate(val = "red", name = "Agotado")

      output <- rbind(Kgreen, Kdarkgreen, Kyellow, Kred)

    } else {
    F_FMSY = local({
      # Need to calculate F-at-age since selectivity varies by fleet
      V_xaf <- sapply(1:nf, function(f) {
        MMSE@multiHist[[1]][[f]]@SampPars$Fleet$V[, , MMSE@nyears + 1]
      }, simplify = "array")
      FM_xfmy <- MMSE@FM[, 1, , , ]

      FM_xafmy <- array(NA,c(MMSE@nsim,dim(V_xaf)[2],nf,length(MMSE@MPs[[1]]),MMSE@proyears))
      ind <- as.matrix(
        expand.grid(x=1:MMSE@nsim,a=1:dim(V_xaf)[2],f=1:nf,m=1:length(MMSE@MPs[[1]]),y=1:MMSE@proyears)
      )
      FM_xafmy[ind] <- FM_xfmy[ind[,-2]] * V_xaf[ind[,c(1:3)]]
      FM_xamy <- apply(FM_xafmy,c(1, 2, 4, 5),sum)  # Sum across fleets
      FM_xmy <- apply(FM_xamy,c(1, 3, 4),max)
      F_FMSY <- FM_xmy/MMSE@RefPoint$ByYear$FMSY
      return(F_FMSY)
    })
   #F_FMSY = sapply(1:nf, function(f) MMSE@multiHist[[1]][[f]]@AtAge$F.Mortality, simplify = "array") %>%
   #  apply(c(1, 2, 3, 4), sum) %>% apply(c(1, 3), max) %>% `/`(MMSE@RefPoint$ByYear$FMSY)
#    F_FMSY = mapply("/",apply(MMSE@FM[, 1, , ,],c(1,3:4),sum),MMSE@RefPoint$ByYear$FMSY[, 1, ,(nyears+1):(nyears+proyears)]) %>% array(c(MMSE@nsim,length(unlist(MMSE@MPs)),MMSE@proyears))
    output <- structure(
      F_FMSY <= 1,
      dimnames = list(Simulation=1:MMSE@nsim,MP=unlist(MMSE@MPs),Year=2024+seq(1,MMSE@proyears))) %>%
      apply(2:3, mean) %>%
      reshape2::melt() %>%
      mutate(val = "grey60", name = "No sobrepesca")

    }
  }

  if (type == "B") {
    fill_values = c("Subexplotado" = "green",
                    "Plena explotado" = "orange",
                    "Sobreexplotado" = "yellow",
                    "Agotado" = "red")

  } else {
    fill_values = c("No sobrepesca" = "grey60")
  }

  g <- output %>%
    mutate(name = factor(name, levels = names(fill_values))) %>%
    ggplot(aes(Year, value)) +
    geom_col(aes(fill = name), width = 1, colour = NA) +
    facet_wrap(vars(MP)) +
    scale_fill_manual(values = fill_values) +
    labs(x = "Años", y = "Probabilidad", fill = NULL) +
    theme(legend.position = "bottom",
          strip.background = element_blank()) +
    guides(fill = guide_legend(ncol = 2)) +
    coord_cartesian(expand = FALSE)

  g
}


