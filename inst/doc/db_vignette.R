## ----setup, include = FALSE---------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "", 
  fig.width = 8, 
  fig.asp = 0.72
)
options(knitr.kable.NA = '')

## ---- message=FALSE-----------------------------------------------------------
library(DNAtools)

## -----------------------------------------------------------------------------
data(dbExample, package = "DNAtools")
knitr::kable(head(dbExample)[,1:9])

## -----------------------------------------------------------------------------
allele_freqs <- lapply(1:10, function(x){
  al_freq <- table(c(dbExample[[x*2]], dbExample[[1+x*2]]))/(2*nrow(dbExample))
  al_freq[sort.list(as.numeric(names(al_freq)))]
})
names(allele_freqs) <- sub("\\.1", "", names(dbExample)[(1:10)*2])

## -----------------------------------------------------------------------------
db_summary <- dbCompare(dbExample, hit = 6, trace = FALSE)

## ---- echo = FALSE, results='asis'--------------------------------------------
db_summary$m[!DNAtools:::up.tri(db_summary$m)] <- NA
rownames(db_summary$m) <- paste0("**",rownames(db_summary$m),"**")
kable(db_summary$m)

## -----------------------------------------------------------------------------
plot(db_summary)

## -----------------------------------------------------------------------------
db_expect <- dbExpect(allele_freqs, n = nrow(dbExample), theta = 0, vector = TRUE)
db_sd <- sqrt(diag(dbVariance(allele_freqs, n = nrow(dbExample), theta = 0)))
plot(db_summary)
points(db_expect, pch = 4)
for(x in seq_along(db_sd)) segments(x0 = x, y0 = db_expect[x]-2*db_sd[x], y1 = db_expect[x]+2*db_sd[x])
legend("topright", bty = "n", pch = c(1,4,NA), c("Observed", "Expected", "95%-CI"),
       lty = c(NA,NA,1))

## -----------------------------------------------------------------------------
a_expect <- dbCollapse(dbExpect(probs = allele_freqs,  n = nrow(dbExample)))
a_observed <- dbCollapse(dbCompare(dbExample, trace = FALSE)$m)
plot(seq_along(a_observed)-1, a_observed, log = "y",
     xlab = "Number of matching alleles", ylab = "Count")
points(seq_along(a_expect)-1, a_expect, pch = 4)

## -----------------------------------------------------------------------------
relatives <- list(
  UN = dbExpect(probs = allele_freqs,  k = c(0,0,1), collapse = TRUE),
  FS = dbExpect(probs = allele_freqs,  k = c(1,2,1)/4, collapse = TRUE),
  FC = dbExpect(probs = allele_freqs,  k = c(0,1,3)/4, collapse = TRUE),
  PC = dbExpect(probs = allele_freqs,  k = c(0,1,0), collapse = TRUE),
  AV = dbExpect(probs = allele_freqs,  k = c(0,1,1)/2, collapse = TRUE)
)

## ---- echo=FALSE--------------------------------------------------------------
pos_range <- unlist(relatives)
pos_range <- range(pos_range[pos_range>0])

plot(NA, xlim = c(0,20), ylim = pos_range, log = "y", type = "n", 
     ylab = "Probability", xlab = "Number of matching alleles")
for(i in seq_along(relatives)) points(0:20, relatives[[i]], col = i, type = "b")
legend("top", horiz = TRUE, col = seq_along(relatives), inset = c(0, -0.1),
       legend = names(relatives), lty = 1, bty = "n", xpd = TRUE)

