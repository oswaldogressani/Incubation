#------------------------------------------------------------------------------#
# Sensitivity-K20 (Log-Normal incubation density, n = 100, coarseness = 1 day).#
# Measuring the performance of incubation estimation with LPS.                 #
# Copyright 2023 Oswaldo Gressani. All rights reserved.                        #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("ggplot2")
source("perfIncubestim.R")

set.seed(1990) # Germany wins FIFA world cup

# Run simulation and extract summary
simobj <- perfIncubestim(nsim = 1000, n = 100, incubdist = "LogNormal",
                         coarseness = 1, K = 20)

write.table(round(simobj$summarytab, 3),
            file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, ".txt"))

write.table(simobj$ModelSelect,
            file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, "MS.txt"))

svg(file = paste0(simobj$tarid, "n", simobj$n, simobj$cid, ".svg"),
    width = 6.5, height = 5)
simobj$DEplot
dev.off()

save(simobj, file = "SensitivityK20.RData")







