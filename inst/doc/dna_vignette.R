## ----setup, include = FALSE---------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "", 
  fig.width = 12, 
  fig.asp = 0.62
)
options(knitr.kable.NA = '')

## ---- echo = FALSE------------------------------------------------------------
STR_profile <- rbind(
  `**Locus:**` = c("vWA","D18","TH01","D2","D8","D3","FGA","D16","D21","D19"),
  `**Alleles:**` = c("15, 18", "14, 17",  "6, 9.3", "17, 23", "12, 15", "15, 15", "19, 23", "11, 12", "28, 28", "13, 14")
)
kable(STR_profile, caption = "A DNA profile from the SGM plus multiplex")

## ---- message=FALSE-----------------------------------------------------------
library(DNAtools)
browseVignettes(package = "DNAtools")

