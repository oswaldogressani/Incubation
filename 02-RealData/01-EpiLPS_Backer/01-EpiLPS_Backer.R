#------------------------------------------------------------------------------#
#                     EpiLPS incubation estimation Wuhan                       #
#           Copyright Oswaldo Gressani 2023. All rights reserved.              #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("xlsx")
library(tidyverse)

# Code snippet from Backer et al. 2020
dataraw <- read_tsv(file = "Backer_dataset.tsv")

dataraw <- dataraw %>% 
  mutate(tReport = as.integer((`reporting date` %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tSymptomOnset = as.integer((symptom_onset %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tStartExposure = as.integer((exposure_start %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tEndExposure = as.integer((exposure_end %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))) %>%
  # if no start of exposure (i.e. Wuhan residents) use arbitrarily chosen start exposure in far away in past (here half December 2019)
  mutate(tStartExposure = ifelse(is.na(tStartExposure), min(tSymptomOnset)-8, tStartExposure)) %>%
  # if symptom onset in Wuhan, exposure ends at symptom onset
  mutate(tEndExposure = ifelse(tEndExposure >= tSymptomOnset, tSymptomOnset, tEndExposure))
# End of code snippet from Backer et al. 2020

tEL <- dataraw$tStartExposure[-which(is.na(dataraw$exposure_start))]
n <- length(tEL)
mintEL <- abs(min(tEL))

tSO <- dataraw$tSymptomOnset[-which(is.na(dataraw$exposure_start))] + mintEL 
tEL <- tEL + mintEL 
tER <- dataraw$tEndExposure[-which(is.na(dataraw$exposure_start))] + mintEL 

all(tSO >= tER) # Check if symptom onset is larger or equal to exposure end
all(tER >= tEL) # Check if exposure end is larger or equal to exposure start

# Making data continuous
datacts <- matrix(0, nrow = n, ncol = 6)
colnames(datacts) <- c("tEL", "tER", "Expowin", "tSO", "tIL", "tIR")
tL <- c()
tR <- c()
expowindow<- c()
set.seed(1234)
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
 datacts[i, ] <- c(tEL_temp, tER_temp, expowindow[i], tSO_temp, tL[i], tR[i])
}

data <- data.frame(tL = tL, tR = tR)
data <- data[-c(which(expowindow > 19)),]
datacts <- datacts[-c(which(expowindow > 19)),]
datacts <- as.data.frame(datacts)

# Check data constraints and export data
nsub <- nrow(datacts)
C1 <- sum(datacts$tEL >= 0) == nsub
C2 <- sum(datacts$tER > datacts$tEL) == nsub
C3 <- sum(datacts$tSO > datacts$tER) == nsub
if(C1 & C2 & C3){
  print("Data ok. All constraints satisfied.")
}
write.xlsx(round(datacts,3), file = "Backer_continuous.xls")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                       tmax = 20)
incubfit$mcmcrate

df  <- data.frame("Reference","Virus","Distibution","Mean incubation","95th percentile")
df[1,] <- c("Backer et al. 2020","SARS-CoV-2","Lognormal",
            "4.5 (3.7-5.6)","8.0 (6.3-11.8)")
df[2,] <- c("EpiLPS", "SARS-CoV-2", "LogNormal",
            paste0(round(incubfit$stats[1,1],1), " (",
                   round(incubfit$stats[1,4],1), "-",
                   round(incubfit$stats[1,5],1), ")"),
            paste0(round(incubfit$stats[21,1],1), " (",
                   round(incubfit$stats[21,4],1), "-",
                   round(incubfit$stats[21,5],1), ")")
            )
colnames(df) <- c("Reference","Virus","Distribution","Mean incubation","95th percentile")
write.table(df, file = "EpiLPS_Backer.txt")

twoDG <- function(x) sprintf("%.2f", x)

# Plot  incubation bounds
dataseg <- cbind(seq_len(24), data)
colnames(dataseg) <- c("index","lower","upper")

incub <- ggplot2::ggplot(data = dataseg, aes(x=index)) +
  ggplot2::geom_linerange(aes(ymin = lower, ymax = upper), 
                          color = "darkslateblue",
                          linewidth = 1) +
  ggplot2::scale_x_continuous(name = "Individual index number",
                              limits = c(0,25)) +
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
  ggplot2::ggtitle("SARS-CoV-2") +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = flow, ymax = fup),
    alpha = 0.5, fill = "gray60") + 
  ggplot2::geom_line(ggplot2::aes(y = fhat), 
                     color = "dodgerblue2", linewidth = 0.8) +
  ggplot2::theme_classic() +
  ggplot2::xlim(0,15) +
  ggplot2::scale_y_continuous(labels = twoDG) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 9)
  )


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
  ggplot2::ggtitle("SARS-CoV-2") +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = Flow, ymax = Fup),
    alpha = 0.5, fill = "gray60") + 
  ggplot2::geom_line(ggplot2::aes(y = Fhat), 
                     color = "brown3", linewidth = 0.8) +
  ggplot2::theme_classic() +
  ggplot2::xlim(0,15) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 9)
  )

# Extract plots
pdf(file = "Backer_Ibounds.pdf", width = 8, height = 4)
suppressWarnings(incub)
dev.off()

pdf(file = "Backer_pdf.pdf", width = 8, height = 4)
suppressWarnings(densplot)
dev.off()

pdf(file = "Backer_cdf.pdf", width = 8, height = 4)
suppressWarnings(cdfplot)
dev.off()

svg(file = "EpiLPS_Backer.svg",width = 13, height = 4.5)
gridExtra::grid.arrange(incub, densplot, cdfplot, nrow = 1, ncol = 3)
dev.off()

df










