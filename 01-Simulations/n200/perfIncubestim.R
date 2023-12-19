


perfIncubestim <- function(nsim = 10, n = 100, coarseness = 1, K = 10,
                           incubdist = c("LogNormal","Weibull","MixWeibull", "Gamma"),
                           niter = 1000, tmax = 20, tgridlen = 300){
  tic <- proc.time()
  tg <- seq(0, tmax, length = tgridlen)
  dtg <- tg[2] - tg[1]
  fI <- matrix(0, nrow = nsim, ncol = tgridlen)
  Mselection <- c()
  MCMCaccept <- c()
  pengpos <- c()
  penoptim <- c()
  PEMean <- c()
  PESD  <- c()
  PEq <- matrix(0, nrow = nsim, ncol = 19)
  meanexposure <- c()
  progbar <- utils::txtProgressBar(min = 1, max = nsim, initial = 1,
                                   style = 3, char =">")
  incubdist <- match.arg(incubdist)

  for(s in 1:nsim){

    # Data generation
    maxDobs <- 1e+10
    while(maxDobs > tmax){
      datgen <- incubsim(incubdist = incubdist,
                         coarsenes = coarseness,
                         n = n,
                         tmax = tmax, tgridlen = tgridlen)
      x <- datgen$Dobsincub
      maxDobs <- max(x)
    }
    meanexposure[s] <- datgen$Expomean

    # Fit the incubation density
    Ifit <- estimIncub(x = x, tmax = tmax, K = K, niter = niter,
                       tgridlen = tgridlen)

    # Record density estimate on grid tg
    fI[s, ] <- Ifit$ftg
    pengpos[s] <- Ifit$pengpos
    penoptim[s] <- Ifit$penval
    MCMCaccept[s] <- Ifit$mcmcrate
    Mselection[s] <- Ifit$modselect

    # Record point estimate of statistics
    stats <- as.data.frame(Ifit$stats)
    PEMean[s] <- stats$PE[1]
    PESD[s]   <- stats$PE[2]
    PEq[s, ] <- stats$PE[3:nrow(stats)]

    utils::setTxtProgressBar(progbar, s)
  }
  close(progbar)

  # Target features of incubation density
  if(incubdist == "LogNormal"){
    Incubmean <- datgen$Iparams[1]
    Incubsd <- datgen$Iparams[2]
    Ipdf <- dlnorm(tg, meanlog = Incubmean, sdlog = Incubsd)
    meanI <- exp(Incubmean + 0.5 * (Incubsd^2))
    sdI <- sqrt(exp(2*Incubmean+Incubsd^2)*(exp(Incubsd^2)-1))
    qgrid <- seq(0.05, 0.95, by = 0.05)
    qtrue <- qlnorm(p = qgrid, meanlog = Incubmean, sdlog = Incubsd)
    tarid <- "LN"
    taridfull <- "Log-Normal"
    destimcol <- grDevices::rgb(145,218,255, maxColorValue = 255)
    destimed <- grDevices::rgb(50, 55, 230, maxColorValue = 255)
    if (coarseness == 1) {
      cid <- "E1"
      if (n == 200) {
        scenarionum <- 1
      } else if (n == 100) {
        scenarionum <- NA
      }
    } else if (coarseness == 2) {
      cid <- "E2"
      if (n == 200){
        scenarionum <- 2
      } else if (n == 100){
        scenarionum <- NA
      }
    }
  } else if (incubdist == "Weibull"){
    Ishape <- datgen$Iparams[1]
    Iscale <- datgen$Iparams[2]
    Ipdf <- stats::dweibull(tg, shape = Ishape, scale = Iscale)
    meanI <- Iscale * gamma(1 + 1 / Ishape)
    sdI <- sqrt((Iscale ^ 2) * ((gamma(1 + 2 / Ishape) - (
      gamma(1 + 1 / Ishape) ^ 2))))
    qgrid <- seq(0.05, 0.95, by = 0.05)
    qtrue <- stats::qweibull(p = qgrid, shape = Ishape, scale = Iscale)
    tarid <- "WB"
    taridfull <- "Weibull"
    destimcol <- grDevices::rgb(145,255,187, maxColorValue = 255)
    destimed <- grDevices::rgb(52,158,92, maxColorValue = 255)
    if (coarseness == 1) {
      cid <- "E1"
      if (n == 200) {
        scenarionum <- 3
      } else if (n == 100) {
        scenarionum <- NA
      }
    } else if (coarseness == 2) {
      cid <- "E2"
      if (n == 200) {
        scenarionum <- 4
      } else if (n == 100) {
        scenarionum <- NA
      }
    }
  } else if (incubdist == "MixWeibull"){
    Ipdf <- datgen$Ipdf
    meanI <- datgen$meanmix
    sdI <- datgen$sdmix
    qtrue <- datgen$qmix
    tarid <- "MixWB"
    taridfull <- "Weibull-mixture"
    destimcol <- grDevices::rgb(221,192,124, maxColorValue = 255)
    destimed <- grDevices::rgb(87, 61, 0, maxColorValue = 255)
    if (coarseness == 1) {
      cid <- "E1"
      if (n == 200) {
        scenarionum <- 5
      } else if (n == 100) {
        scenarionum <- NA
      }
    } else if (coarseness == 2) {
      cid <- "E2"
      if (n == 200) {
        scenarionum <- 6
      } else if (n == 100) {
        scenarionum <- NA
      }
    }
  } else if (incubdist == "Gamma"){
    Ishape <- datgen$Iparams[1]
    Irate <- datgen$Iparams[2]
    Ipdf <- stats::dgamma(tg, shape = Ishape, rate = Irate)
    meanI <- Ishape / Irate
    sdI <- sqrt(Ishape / (Irate ^ 2))
    qgrid <- seq(0.05, 0.95, by = 0.05)
    qtrue <- stats::qgamma(p = qgrid, shape = Ishape, rate = Irate)
    tarid <- "GA"
    taridfull <- "Gamma"
    destimcol <- grDevices::rgb(255,198,224, maxColorValue = 255)
    destimed <- grDevices::rgb(242,0,103, maxColorValue = 255)
    if (coarseness == 1) {
      cid <- "E1"
      if (n == 200) {
        scenarionum <- 7
      } else if (n == 100) {
        scenarionum <- NA
      }
    } else if (coarseness == 2) {
      cid <- "E2"
      if (n == 200) {
        scenarionum <- 8
      } else if (n == 100) {
        scenarionum <- NA
      }
    }
  }

  # Summary table
  simsum <- matrix(0, nrow = length(qtrue) + 2, ncol = 4)
  rownames(simsum) <- c("Mean","SD",  paste0("q", seq(0.05, 0.95, by = 0.05)))
  colnames(simsum) <- c("True","Average","Bias","RMSE")

  simsum[, 1] <- c(meanI, sdI, qtrue)
  simsum[, 2] <- c(mean(PEMean), mean(PESD), colMeans(PEq))
  simsum[, 3] <- simsum[, 2] - simsum[, 1]
  simsum[, 4] <- c(sqrt(mean((PEMean - meanI)^2)),
                   sqrt(mean((PESD - sdI)^2)),
                   sqrt(colMeans((PEq - matrix(rep(qtrue, nsim), nrow = nsim,
                                               byrow = T))^2)))

  summarytab <- simsum[c(1,2,3,7,12,17,21),]
  
  twoDG <- function(x) sprintf("%.2f", x)

  df <- data.frame(t = tg, ftrue = Ipdf, fmed = apply(fI,2,"median"),
                   fIfestim = t(fI))
  df <- df[1:head(which(tg > 16), 1), ]
  datcolnames <- colnames(df)
  legpos <- c(0.85,0.92)
  if(incubdist == "MixWeibull"){
    legpos <- c(0.5,0.91)
  }

  colors <- c("Target"= "black", "Estim" = destimcol, "Median" = destimed)
  linetypes <- c("Target" = 1,  "Estim" = 1, "Median" = 4)

  DEplot <- ggplot(data = df, aes(x = t)) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    labs(x = "Time from infection to onset (days)",
         y = "Probability density function", color = "",linetype = "") +
    theme_classic() +
    scale_y_continuous(labels = twoDG) +
    theme(legend.position = legpos)
  for(s in 4:(ncol(df))){
    DEplot <- DEplot + eval(parse(
      text = paste("ggplot2::geom_line(ggplot2::aes(y=",datcolnames[s],
                   ", color = 'Estim', linetype = 'Estim'))",sep = "")))
  }

  DEplot <- DEplot + geom_line(aes(y=ftrue, color = "Target",
                                   linetype = "Target"),
                               linewidth = 0.8) +
    geom_line(aes(y=fmed, color = "Median", linetype = "Median"),
              linewidth = 0.9) + ggtitle(paste0("Scenario ",scenarionum))

  #---- Squared Hellinger distances
  H2 <- c()
  for(s in 1:nsim) {
    H2[s] <- 0.5 * sum((sqrt(fI[s, ]) - sqrt(Ipdf)) ^ 2 * dtg)
  }

  H2df <- data.frame(H2 = H2)
  H2histo <- ggplot(data = H2df, aes(x=H2)) +
    geom_histogram(color = "gray40", fill = "gray75") +
    theme_classic() +
    labs(x = "Squared Hellinger distance", y = "Frequency") +
    ggtitle(paste0("Scenario ",scenarionum))

  ModelSelect <- matrix(0, nrow = 4, ncol = 1)
  colnames(ModelSelect) <- c("Percentage")
  rownames(ModelSelect) <- c("LPS","Log-Normal","Gamma","Weibull")
  ModelSelect[1,1] <- sum(Mselection == 1)/nsim * 100
  ModelSelect[2,1] <- sum(Mselection == 2)/nsim * 100
  ModelSelect[3,1] <- sum(Mselection == 3)/nsim * 100
  ModelSelect[4,1] <- sum(Mselection == 4)/nsim * 100

  #---- Plot histogram of maximum a posteriori penalty values
  hist(penoptim, breaks = sqrt(nsim ) + 5,
       xlab = "Penalty MAP", ylab = "", main = "Optimal penalty",
       col = "orange")

  hist(MCMCaccept, breaks = 10,
       xlab = "MCMC-Langevin acceptance rate", ylab = "",
       main = "MCMC(Langevin) acceptance",
       xlim = c(0,100), col = "cornsilk3")
  abline(v=57, col = "red", lty = 2, lwd = 2)

  toc <- proc.time() - tic

  cat("-----------------------------------------------------------\n")
  cat("Scenario", scenarionum, "\n")
  cat("-----------------------------------------------------------\n")
  cat(paste0("Percentage of lambda optim. interior: ",
             round((sum(pengpos=="interior")/nsim) * 100,2), " %.\n"))
  cat(paste0("Average exposure duration: ", round(mean(meanexposure),3),
             " day(s).\n"))
  cat(paste0("Number of simulations: ", nsim,".","\n"))
  cat(paste0("Sample size n: ", n, ".\n"))
  cat(paste0(K, " B-splines specified in time domain: [",0,",",tmax,"].\n"))
  cat("Average MCMC acceptance rate: ", mean(MCMCaccept),"%.\n")
  cat(paste0("MALA chain length: ", niter,".\n"))
  cat(paste0("Target incubation density is ", taridfull, "(",
             datgen$Iparams[1],",",datgen$Iparams[2],")", ".\n"))
  cat("-----------------------------------------------------------\n")
  print(round(summarytab,3))
  cat("-----------------------------------------------------------\n")
  cat("Time elapsed: ", round(toc[3]/60,3), "minutes. \n")
  cat("-----------------------------------------------------------\n")

  outlist <- list(DEplot = DEplot, H2histo = H2histo,
                  ModelSelect = ModelSelect, summarytab = summarytab,
                  tarid = tarid, n = n, cid = cid,
                  penoptim = penoptim)

  return(outlist)

}
