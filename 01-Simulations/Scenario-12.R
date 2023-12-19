#------------------------------------------------------------------------------#
# Scenario 12 (Weibmix incubation density, n = 100, coarseness = 2 day).        #
# Measuring the performance of incubation estimation with LPS.                 #
# Copyright 2023 Oswaldo Gressani. All rights reserved.                        #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
source("perfIncubestim.R")

set.seed(1234) # A simple sequence of numbers

# Run simulation and extract summary
simobj <- perfIncubestim(nsim = 1000, n = 100, incubdist = "MixWeibull",
                         coarseness = 2, K = 20)

write.table(round(simobj$summarytab, 3),
            file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, ".txt"))

write.table(simobj$ModelSelect,
            file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, "MS.txt"))

pdf(file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, ".pdf"),
    width = 6.5, height = 5)
simobj$DEplot
dev.off()

pdf(file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, "Hellinger.pdf"),
    width = 6.5, height = 5)
simobj$H2histo
dev.off()

save(simobj, file = "S12.RData")







