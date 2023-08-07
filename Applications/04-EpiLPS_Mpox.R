#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation Mpox                              #
#           Copyright Oswaldo Gressani 2023. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original xls data downloaded from Miura et al. (2022)
# https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2022.27.24.2200448

library("xlsx")
library("EpiLPS")

## Load Data
data <- read.csv("Datasets/Mpox.csv")
n <- nrow(data)

tEL <- data$Start.date.of.exposure
tER <- data$End.date.of.exposure
tSO <- data$Symptom.onset

# Making data continuous
tL <- c()
tR <- c()
expowindow<- c()
set.seed(123)
for(i in 1:n){
  tSO_temp <- tSO[i] + runif(1)
  tEL_temp <- tEL[i] + runif(1)
  tER_temp <- tER[i] + runif(1)
  while(tSO_temp <= tER_temp | tER_temp <= tEL_temp){
    tSO_temp <- tSO[i] + runif(1)
    tEL_temp <- tEL[i] + runif(1)
    tER_temp <- tER[i] + runif(1)
  }
  tL[i] <- tSO_temp - tER_temp
  tR[i] <- tSO_temp - tEL_temp
  expowindow[i] <- tER_temp - tEL_temp
}

data <- data.frame(tL = tL, tR = tR)

# Fit with EpiLPS

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                       tmax = 30)
incubfit$mcmcrate
plot(incubfit$tg, incubfit$ftg, type = "l", col = "blue")
rug(c(data$tL,data$tR))

# Comparison with Miura et al. 2022 who fitted a Lognormal with mean
# 9.0 (6.6-10.9) and 95th percentile 17.3 (13.0-29.0)

df2  <- data.frame("Reference","Virus","Distibution","Mean incubation","95th percentile")
df2[1,] <- c("Miura et al. 2022","Monkeypox","Lognormal",
            "9.0 (6.6-10.9)","17.3 (13.0-29.0)")
df2[2,] <- c("EpiLPS", "Monkeypox", "Lognormal",
            paste0(round(incubfit$stats[1,1],1), " (",
                   round(incubfit$stats[1,4],1), "-",
                   round(incubfit$stats[1,5],1), ")"),
            paste0(round(incubfit$stats[21,1],1), " (",
                   round(incubfit$stats[21,4],1), "-",
                   round(incubfit$stats[21,5],1), ")")
)
colnames(df2) <- c("Reference","Virus","Distribution","Mean incubation","95th percentile")
df2

write.table(df2, file = "EpiLPS_Miura.txt")


# Plot  incubation bounds

dataseg <- cbind(seq_len(18), data)
colnames(dataseg) <- c("index","lower","upper")

incub <- ggplot2::ggplot(data = dataseg, ggplot2::aes(x=index)) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), 
                          color = "darkslateblue",
                          linewidth = 1) +
  ggplot2::xlab("Individual index number") +
  ggplot2::ylab("Incubation bound") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 9)
  )


# Plot pdf + 95%CI---------------------------------------------
tdom <- incubfit$tg
fhat <- incubfit$ftg
flow <- apply(incubfit$ftgMCMC, 2, quantile, probs = 0.025)
fup <- apply(incubfit$ftgMCMC, 2, quantile, probs = 0.975)
densdat <- data.frame(tdom = tdom, fhat = fhat, flow = flow, fup = fup)


densplot <-
  ggplot2::ggplot(data = densdat, ggplot2::aes(x = tdom)) +
  ggplot2::xlab("Incubation period (days)") +
  ggplot2::ylab("Probability density function") +
  ggplot2::ggtitle("Monkeypox") +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = flow, ymax = fup),
    alpha = 0.5, fill = "gray60") + 
  ggplot2::geom_line(ggplot2::aes(y = fhat), 
                     color = "dodgerblue2", linewidth = 0.8) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 9)
  )

densplot

# Plot cdf + 95%CI---------------------------------------------

tdom <- incubfit$tg
dt <- tdom[2] - tdom[1]
Fhat <- cumsum(fhat * dt)
FMCMC <- matrix(0, nrow = nrow(incubfit$ftgMCMC), ncol = ncol(incubfit$ftgMCMC))
for(j in 1:nrow(incubfit$ftgMCMC)){
  FMCMC[j,] <- cumsum(incubfit$ftgMCMC[j,] * dt)
}
Flow <- apply(FMCMC, 2, quantile, probs = 0.025)
Fup <- apply(FMCMC, 2, quantile, probs = 0.975)
cdfdat <- data.frame(tdom = tdom, Fhat = Fhat, Flow = Flow, Fup = Fup)

cdfplot <-
  ggplot2::ggplot(data = cdfdat, ggplot2::aes(x = tdom)) +
  ggplot2::xlab("Incubation period (days)") +
  ggplot2::ylab("Cumulative distribution function") +
  ggplot2::ggtitle("Monkeypox") +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = Flow, ymax = Fup),
    alpha = 0.5, fill = "gray60") + 
  ggplot2::geom_line(ggplot2::aes(y = Fhat),
                     color = "brown3", linewidth = 0.8) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 9)
  )

cdfplot


gridExtra::grid.arrange(incub, densplot, cdfplot, nrow = 1, ncol = 3)
df2



























