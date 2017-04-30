## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(fig.width = 5, fig.height = 4, echo = TRUE)
library(spant)

## ----useage, message = FALSE---------------------------------------------
library(spant)
fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")

# import raw data
mrs_data <- read_mrs(fname, format = "spar_sdat")

# output basic data structure
print(mrs_data)

# plot data in the frequency domain
plot(mrs_data, xlim = c(5, 0.5))

## ----processing, message = FALSE-----------------------------------------
# apply water filter and align to tNAA resonance
mrs_proc <- hsvd_filt(mrs_data)
mrs_proc <- align(mrs_proc, 2.01)
plot(mrs_proc, xlim = c(5, 0.5))

## ----basis_sim, message = FALSE------------------------------------------
# simulate a typical basis set for short TE brain analysis
basis <- sim_basis_1h_brain_press(mrs_proc)

# output basis info
print(basis)

# plot basis signals
stackplot(basis, xlim = c(4, 0.5))

## ----fitting, message = FALSE--------------------------------------------
# perform VARPRO fitting to processed data
fit_res <- fit_mrs(mrs_proc, basis)

# plot the fit estimate, residual and baseline
plot(fit_res)

