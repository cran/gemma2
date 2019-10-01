## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
readr::read_csv(system.file("extdata", "mouse100.geno.txt", package = "gemma2"), col_names = FALSE) -> geno
readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE) -> pheno
readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100] -> kinship


## ------------------------------------------------------------------------
# isolate first SNP
g1 <- geno[1, - c(1:3)] # first 3 columns are SNP annotations!
t(as.matrix(g1)) -> g1m
as.matrix(pheno[, c(1, 6)]) -> phe16

## ------------------------------------------------------------------------
library(gemma2)
e_out <- eigen2(kinship)

## ------------------------------------------------------------------------
t(g1m) %*% e_out$vectors -> X1U

## ------------------------------------------------------------------------
MphEM(eval = e_out$values, 
      X = X1U, 
      Y = t(phe16) %*% e_out$vectors, 
      V_g = matrix(c(1.91352, 0, 0, 0.530827), nrow = 2), 
      V_e = matrix(c(0.320028, 0, 0, 0.561589), nrow = 2)
      ) -> foo

## ------------------------------------------------------------------------
length(foo)
class(foo)

## ------------------------------------------------------------------------
sapply(FUN = length, X = foo)

## ------------------------------------------------------------------------
MphEM

## ------------------------------------------------------------------------
sapply(FUN = function(x)x[[1]], X = foo) -> loglik
plot(loglik)

## ------------------------------------------------------------------------
e_out <- eigen2(as.matrix(kinship))
center_kinship(as.matrix(kinship)) -> kinship_centered
ec_out <- eigen2(kinship_centered)

## ------------------------------------------------------------------------
## first two
MphEM(eval = ec_out$values, 
      X = t(g1m) %*% ec_out$vectors, 
      Y = t(phe16) %*% ec_out$vectors, 
      V_g = matrix(c(1.91352, 0, 0, 0.530827), nrow = 2), 
      V_e = matrix(c(0.320028, 0, 0, 0.561589), nrow = 2)
      ) -> bar00
MphEM(eval = e_out$values, 
      X = t(g1m) %*% e_out$vectors, 
      Y = t(phe16) %*% e_out$vectors, 
      V_g = matrix(c(1.91352, 0, 0, 0.530827), nrow = 2), 
      V_e = matrix(c(0.320028, 0, 0, 0.561589), nrow = 2)
      ) -> bar10
## second two
MphEM(eval = ec_out$values, 
      X = t(cbind(1, g1m)) %*% ec_out$vectors, 
      Y = t(phe16) %*% ec_out$vectors, 
      V_g = matrix(c(1.91352, 0, 0, 0.530827), nrow = 2), 
      V_e = matrix(c(0.320028, 0, 0, 0.561589), nrow = 2)
      ) -> bar01

MphEM(eval = e_out$values, 
      X = t(cbind(1, g1m)) %*% e_out$vectors, 
      Y = t(phe16) %*% e_out$vectors, 
      V_g = matrix(c(1.91352, 0, 0, 0.530827), nrow = 2), 
      V_e = matrix(c(0.320028, 0, 0, 0.561589), nrow = 2)
      ) -> bar11


## ------------------------------------------------------------------------
bar00[[length(bar00)]][[2]]
bar01[[length(bar01)]][[2]]
bar10[[length(bar10)]][[2]]
bar11[[length(bar11)]][[2]]

## ------------------------------------------------------------------------
t(cbind(1, g1m)) %*% e_out$vectors

## ------------------------------------------------------------------------
MphEM(eval = e_out$values, 
      X = t(cbind(1, g1m)) %*% e_out$vectors, 
      Y = t(phe16) %*% e_out$vectors, 
      V_g = matrix(c(2.73, - 0.39, - 0.39, 0.74), nrow = 2), 
      V_e = matrix(c(0.15, 0.27, 0.27, 0.50), nrow = 2)
      ) -> out1
out1[[1]][c(2:3)]

## ------------------------------------------------------------------------
MphEM(eval = ec_out$values, 
      X = t(cbind(1, g1m)) %*% ec_out$vectors, 
      Y = t(phe16) %*% ec_out$vectors, 
      V_g = matrix(c(2.731324, - 0.394865, - 0.394865, 0.737116), nrow = 2), 
      V_e = matrix(c(0.147467, 0.272090, 0.272090, 0.502031), nrow = 2)
      ) -> out1
out1[[1]][c(2:3)]


