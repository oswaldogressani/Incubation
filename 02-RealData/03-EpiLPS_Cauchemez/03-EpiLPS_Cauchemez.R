#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation MERS                              #
#           Copyright Oswaldo Gressani 2023. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original data can be found in Cauchemez et al. (2014)
# https://www.sciencedirect.com/science/article/pii/S1473309913703049

library("xlsx")
library("EpiLPS")

## Data
tL <- c(6,1,9,3,4,3,4)
tR <- c(9,4,12,4,4,3,4)
n <- length(tL)

set.seed(1234)

for(i in 1:n){
  tL_temp <- tL[i] + runif(1)
  tR_temp <- tR[i] + runif(1)
  while(tL_temp >= tR_temp){
    tL_temp <- tL[i] + runif(1)
    tR_temp <- tR[i] + runif(1)
  }
  tL[i] <- tL_temp
  tR[i] <- tR_temp
}

data <- data.frame(tL = tL, tR = tR)
write.xlsx(round(data,3), file = "Cauchemez_continuous.xls")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                               tmax = 20)

incubfit$mcmcrate
plot(incubfit$tg, incubfit$ftg, type = "l", col = "blue")
rug(c(data$tL,data$tR))

# From Cauchemez et al. 2014
LNmean <- 5.5
LNSD <- 2.5

# sdlogLN <- sqrt(log(((2.5)^2)/(exp(2*log(5.5))) + 1))
# meanlogLN <- log(5.5) - 0.5 * sdlogLN^2

funsig2 <- function(sig2){
  mu <- 0.5 * (2 * log(5.5) - sig2)
  val <- exp(2 * mu + sig2) * (exp(sig2) - 1) - (2.5 ^ 2)
  return(val)
}
sig2dom <- seq(0,3, length = 100)
plot(sig2dom, sapply(sig2dom, funsig2), type = "l")
abline(h=0)
sig2hat <- uniroot(funsig2, interval = c(0,3))$root
sdlog <- sqrt(sig2hat)
meanlog <- 0.5 * (2 * log(5.5) - sig2hat)
round(qlnorm(p = 0.95, meanlog = meanlog, sdlog = sdlog),1)

df2  <- data.frame("Reference","Virus","Distibution","Mean incubation","95th percentile")
df2[1,] <- c("Cauchemez et al. 2014","MERS","LogNormal",
             "5.5 (3.6-10.2)","10.2 (NA)")
df2[2,] <- c("EpiLPS", "MERS", "LogNormal",
             paste0(round(incubfit$stats[1,1],1), " (",
                    round(incubfit$stats[1,4],1), "-",
                    round(incubfit$stats[1,5],1), ")"),
             paste0(round(incubfit$stats[21,1],1), " (",
                    round(incubfit$stats[21,4],1), "-",
                    round(incubfit$stats[21,5],1), ")")
)
colnames(df2) <- c("Reference","Virus","Distribution","Mean incubation","95th percentile")

write.table(df2, file = "EpiLPS_Cauchemez.txt")

# Plot  incubation bounds
dataseg <- cbind(seq_len(7), data)
colnames(dataseg) <- c("index","lower","upper")

incub <- ggplot2::ggplot(data = dataseg, ggplot2::aes(x=index)) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), 
                          color = "darkslateblue",
                          linewidth = 1) +
  ggplot2::scale_x_continuous(name = "Individual index number",
                              limits = c(0,8)) +
  ggplot2::scale_y_continuous(name="Incubation bound",
                              limits = c(0,15)) +
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
  ggplot2::ylab("Probability density function") +
  ggplot2::ggtitle("MERS-CoV") +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = flow, ymax = fup),
    alpha = 0.5, fill = "gray60") + 
  ggplot2::geom_line(ggplot2::aes(y = fhat), 
                     color = "dodgerblue2", linewidth = 0.8) +
  ggplot2::theme_classic() +
  ggplot2::scale_x_continuous(name = "Incubation period (days)",
                              breaks = c(0,5,10,15,20)) +
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
  ggplot2::ylab("Cumulative distribution function") +
  ggplot2::ggtitle("MERS-CoV") +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = Flow, ymax = Fup),
    alpha = 0.5, fill = "gray60") + 
  ggplot2::geom_line(ggplot2::aes(y = Fhat), 
                     color = "brown3", linewidth = 0.8) +
  ggplot2::theme_classic() +
  ggplot2::scale_x_continuous(name = "Incubation period (days)",
                              breaks = c(0,5,10,15,20)) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 9)
  )

# Extract plots
pdf(file = "Cauchemez_Ibounds.pdf", width = 8, height = 4)
suppressWarnings(incub)
dev.off()

pdf(file = "Cauchemez_pdf.pdf", width = 8, height = 4)
suppressWarnings(densplot)
dev.off()

pdf(file = "Cauchemez_cdf.pdf", width = 8, height = 4)
suppressWarnings(cdfplot)
dev.off()

svg(file = "EpiLPS_Cauchemez.svg",width = 13, height = 4.5)
gridExtra::grid.arrange(incub, densplot, cdfplot, nrow = 1, ncol = 3)
dev.off()


df2























