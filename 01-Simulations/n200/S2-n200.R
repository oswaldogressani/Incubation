#------------------------------------------------------------------------------#
# Scenario 2 (Log-Normal incubation density, n = 200, coarseness = 2 day).     #
# Measuring the performance of incubation estimation with LPS.                 #
# Copyright 2023 Oswaldo Gressani. All rights reserved.                        #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
source("perfIncubestim.R")

set.seed(1982) # Italy wins FIFA World Cup

# Run simulation and extract summary
simobj <- perfIncubestim(nsim = 1000, n = 200, incubdist = "LogNormal",
                         coarseness = 2, K = 10)

write.table(round(simobj$summarytab, 3),
            file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, ".txt"))

write.table(simobj$ModelSelect,
            file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, "MS.txt"))

svg(file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, ".svg"),
    width = 6.5, height = 5)
simobj$DEplot
dev.off()

pdf(file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, "Hellinger.pdf"),
    width = 6.5, height = 5)
simobj$H2histo
dev.off()

save(simobj, file = "S2.RData")







