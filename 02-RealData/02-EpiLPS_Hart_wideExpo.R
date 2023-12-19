#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation transmission pairs                #
#           Copyright Oswaldo Gressani 2023. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original xls data downloaded from the article of Hart et al. (2020)
# https://elifesciences.org/articles/65534/figures#content

library("xlsx")
library("EpiLPS")

## Load Data
data <- read.xlsx("Datasets/Hart_dataset.xlsx", sheetIndex = 1)
n <- nrow(data)

E1L <- data$t_i1.L # Left bound for the day on which infector was infected
# min1L <- min(E1L, na.rm = TRUE)
E2L <- data$t_i2.L # Left bound for the day on which infectee was infected

E1R <- data$t_i1.R # Right bound for the day on which infector was infected
E2R <- data$t_i2.R # Right bound for the day on which infectee was infected

SO1 <- data$t_s1   # Day of symptom onset for infector
SO2 <- data$t_s2   # Day of symptom onset for infectee

E1R[which(E1R > SO1)] <- SO1[which(E1R > SO1)]
E2R[which(E2R > SO2)] <- SO2[which(E2R > SO2)]


# Making data continuous
datacts1 <- matrix(0, nrow = n, ncol = 6)
colnames(datacts1) <- c("tEL", "tER", "Expowin", "tSO", "tIL", "tIR")
datacts2 <- matrix(0, nrow = n, ncol = 6)
colnames(datacts2) <- c("tEL", "tER", "Expowin", "tSO", "tIL", "tIR")
t1L <- c()
t1R <- c()
t2L <- c()
t2R <- c()
expowindow1 <- c()
expowindow2 <- c()
set.seed(1234)
for(i in 1:n){
  # Infector
  tSO1_temp <- SO1[i] + runif(1)
  tE1L_temp <- E1L[i] + runif(1)
  tE1R_temp <- E1R[i] + runif(1)
  # Infectee
  tSO2_temp <- SO2[i] + runif(1)
  tE2L_temp <- E2L[i] + runif(1)
  tE2R_temp <- E2R[i] + runif(1)
  if(!anyNA(c(tE1L_temp, tE1R_temp))) {
    while (tSO1_temp <= tE1R_temp | tE1R_temp <= tE1L_temp) {
      tSO1_temp <- SO1[i] + runif(1)
      tE1L_temp <- E1L[i] + runif(1)
      tE1R_temp <- E1R[i] + runif(1)
    }
  }
  if(!anyNA(c(tE2L_temp, tE2R_temp))) {
    while (tSO2_temp <= tE2R_temp | tE2R_temp <= tE2L_temp) {
      tSO2_temp <- SO2[i] + runif(1)
      tE2L_temp <- E2L[i] + runif(1)
      tE2R_temp <- E2R[i] + runif(1)
    }
  }
  t1L[i] <- tSO1_temp - tE1R_temp
  t1R[i] <- tSO1_temp - tE1L_temp
  expowindow1[i] <- tE1R_temp - tE1L_temp
  datacts1[i, ] <- c(tE1L_temp, tE1R_temp, expowindow1[i], tSO1_temp, t1L[i], t1R[i])
  
  t2L[i] <- tSO2_temp - tE2R_temp
  t2R[i] <- tSO2_temp - tE2L_temp
  expowindow2[i] <- tE2R_temp - tE2L_temp
  datacts2[i, ] <- c(tE2L_temp, tE2R_temp, expowindow2[i], tSO2_temp, t2L[i], t2R[i])
}


df <- cbind(c(t1L, t2L), c(t1R, t2R))

datacts <- rbind(datacts1, datacts2)

# df <- cbind(c(t1L, t2L), c(t1R, t2R))
colnames(df) <- c("tL","tR")

df <- df[-which(is.na(df[,1])),]
df <- df[-which(is.na(df[,2])),]
df <- as.data.frame(df)
datacts <- datacts[-which(is.na(datacts[,5])),]
datacts <- datacts[-which(is.na(datacts[,6])),]
datacts <- as.data.frame(datacts)

# Check data constraints and export data
nsub <- nrow(datacts)
C1 <- sum(datacts$tEL >= 0) == nsub
C2 <- sum(datacts$tER > datacts$tEL) == nsub
C3 <- sum(datacts$tSO > datacts$tER) == nsub
if(C1 & C2 & C3){
  print("Data ok. All constraints satisfied.")
}

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = df, K = 20, niter = 20000, verbose = TRUE,
                       tmax = 23)
incubfit$mcmcrate
plot(incubfit$tg, incubfit$ftg, type = "l", col = "blue")
rug(c(df$tL,df$tR))

# Comparison with Xia et al. 2020 who fitted a Weibull with mean
# 4.9 (4.4-5.4) and 95th percentile 9.9 (8.9-11.2)

df2  <- data.frame("Reference","Virus","Distibution","Mean incubation","95th percentile")
df2[1,] <- c("Xia et al. 2020","SARS-CoV-2","Weibull",
            "4.9 (4.4-5.4)","9.9 (8.9-11.2)")
df2[2,] <- c("EpiLPS", "SARS-CoV-2", "Weibull",
            paste0(round(incubfit$stats[1,1],1), " (",
                   round(incubfit$stats[1,4],1), "-",
                   round(incubfit$stats[1,5],1), ")"),
            paste0(round(incubfit$stats[21,1],1), " (",
                   round(incubfit$stats[21,4],1), "-",
                   round(incubfit$stats[21,5],1), ")")
)
colnames(df2) <- c("Reference","Virus","Distribution","Mean incubation","95th percentile")
df2

write.table(df2, file = "EpiLPS_Hart_wideExpo.txt")

# Plot  incubation bounds

dataseg <- cbind(seq_len(nrow(df)), df)
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
  ggplot2::ggtitle("SARS-CoV-2") +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = flow, ymax = fup),
    alpha = 0.5, fill = "gray60") + 
  ggplot2::geom_line(ggplot2::aes(y = fhat), 
                     color = "dodgerblue2", linewidth = 0.8) +
  ggplot2::theme_classic() +
  ggplot2::xlim(0,20) +
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
  ggplot2::ggtitle("SARS-CoV-2") +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = Flow, ymax = Fup),
    alpha = 0.5, fill = "gray60") + 
  ggplot2::geom_line(ggplot2::aes(y = Fhat), 
                     color = "brown3", linewidth = 0.8) +
  ggplot2::theme_classic() +
  ggplot2::xlim(0,20) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 9)
  )

cdfplot

gridExtra::grid.arrange(incub, densplot, cdfplot, nrow = 1, ncol = 3)

































