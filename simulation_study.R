### Package installation
install.packages('readxl')
install.packages('zoo')
install.packages('glmnet')
install.packages('tseries')
install.packages('caret')
install.packages('Matrix')
install.packages('corrplot')
install.packages('snow')
install.packages('foreach')
install.packages('tictoc')
install.packages('VIF')
install.packages('tikzDevice')
install.packages('gridExtra')
install.packages('dplyr')

### Call packages
library('readxl')
library('zoo')
library('glmnet')
library('Matrix')
library('tseries')
library('caret')
library('Matrix')
library('corrplot')
library('parallel')
library('snow')
library('MASS')
library('tictoc')
library('foreach')
library('VIF')
library('ggplot2')
library('tikzDevice')
library('gridExtra')
library('dplyr')


source('myfun.R')
source('dgp.R')
#### Setting for simulation: ----
## Number of different settings: n = 40, 80, 120, 200, 500, 1000
## Number of model for evaluation: OLS, Oracle OLS, AR, Lasso, and Adaptive Lasso
## Evaluation is carried out by computing Mean Squared Prediction Error
## For each setting, the study is replicated 1000 times. Final MSPE is the average of all MSPE in all replications
## For Lasso and Alasso, lambda will be chosen by 10 fold cross validation using a separate exploratory sample size n = 200 with 100 replication.
## Then the final c_lambda would be the median of all MSPE-optimized c_lambda.
## Takes about 2-3 hours to run

######################################################################################################
######################## DATA GENERATING PROCESS 2 - DEGENERACY IN THE LIMIT #########################
######################################################################################################
#### Simulating data ----
tic("total")
set.seed(26)
# Simulating data for finding c_lambda. c_lambda is the median of all CV-MSPE-optimized lambdas in each replication.
dgp2.sim.optim.200 <- replicate(100, fun.dgp2(n.burn = 1000, n.obs = 200))
# Simulating data for 6 settings

dgp2.sim.40 <- replicate(1000, fun.dgp2(n.burn = 500, n.obs = 40))
dgp2.sim.80 <- replicate(1000, fun.dgp2(n.burn = 500, n.obs = 80))
dgp2.sim.160 <- replicate(1000, fun.dgp2(n.burn = 500, n.obs = 160))
dgp2.sim.250 <- replicate(1000, fun.dgp2(n.burn = 500, n.obs = 250))
dgp2.sim.500 <- replicate(1000, fun.dgp2(n.burn = 500, n.obs = 500))
dgp2.sim.1000 <- replicate(1000, fun.dgp2(n.burn = 500, n.obs = 1000))
# dgp2.sim.10000 <- replicate(1000, fun.dgp2(n.burn = 500, n.obs = 10000))

# n.settings <- c(40, 80, 160, 250, 500, 1000, 10000)
n.settings <- c(40, 80, 160, 250, 500, 1000)
sim.scale = FALSE

toc()
#### Alasso performance assessment ----
## Finding optimized c_lambda by 10-folds cross validation for each replication of simulated data.

tic()
dgp2.alasso.lambda.optim <- fun.optim.lambda(dgp2.sim.optim.200, start.val = 0.5, 
                                             pen.factor = TRUE, pen.type = 'lm',
                                             scale = sim.scale)
dgp2.alasso.c_lambda <- median(dgp2.alasso.lambda.optim[, 1])
toc()

tic()
dgp2.lasso.lambda.optim <- fun.optim.lambda(dgp2.sim.optim.200, start.val = 0.5, 
                                            pen.factor = FALSE,
                                            scale = sim.scale)
dgp2.lasso.c_lambda <- median(dgp2.lasso.lambda.optim[, 1])
toc()

## Running Alasso on simulated data using obtained c_lambda * sqrt(n) / log(log(n))
## Two main metrics for performance assessment is MSPE and SR (Screening Error)
##    1) For MSPE: we will compute the OOS one-step-ahead squared prediction error for each replication
## in each simulation setting. MSPE will be the mean across all replications in each simulation setting.
##    2) For SR, there are three sub-metrics: 
##      - SR: overall success rate of classification into zero coefficients and non-zero coefficients
##      - SR 1: percentage of the correct selection in the active set
##      - SR 2: percentage of correct elimination of the zero coefficients

## Fitting for each replication in a setting and analyse performance metrics

#### ALASSO ----

dgp2.alasso.mspe.loglogn <- rep(0, length(n.settings))
dgp2.alasso.beta.loglogn <- list()
dgp2.alasso.sr.loglogn <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.alasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.alasso.current.fit <- fun.fit.lasso.all(get(dgp2.alasso.current.setting),
                                               c_lambda = dgp2.alasso.c_lambda,
                                               b_n = 'log(log(n))',
                                               pen.factor = TRUE,
                                               scale = sim.scale)
  dgp2.alasso.mspe.loglogn[i] <- mean(dgp2.alasso.current.fit$mspe)
  dgp2.alasso.beta.loglogn[[i]] <- dgp2.alasso.current.fit$beta
  dgp2.alasso.sr.loglogn[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                                beta.est = dgp2.alasso.current.fit$beta, n.sim = 1000)
}

dgp2.alasso.mspe.logn <- rep(0, length(n.settings))
dgp2.alasso.beta.logn <- list()
dgp2.alasso.sr.logn <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.alasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.alasso.current.fit <- fun.fit.lasso.all(get(dgp2.alasso.current.setting),
                                               c_lambda = dgp2.alasso.c_lambda,
                                               b_n = 'log(n)',
                                               pen.factor = TRUE,
                                               scale = sim.scale)
  dgp2.alasso.mspe.logn[i] <- mean(dgp2.alasso.current.fit$mspe)
  dgp2.alasso.beta.logn[[i]] <- dgp2.alasso.current.fit$beta
  dgp2.alasso.sr.logn[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                             beta.est = dgp2.alasso.current.fit$beta, n.sim = 1000)
}

dgp2.alasso.mspe.logn2 <- rep(0, length(n.settings))
dgp2.alasso.beta.logn2 <- list()
dgp2.alasso.sr.logn2 <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.alasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.alasso.current.fit <- fun.fit.lasso.all(get(dgp2.alasso.current.setting),
                                               c_lambda = dgp2.alasso.c_lambda,
                                               b_n = 'log(n) ^ 2',
                                               pen.factor = TRUE,
                                               scale = sim.scale)
  dgp2.alasso.mspe.logn2[i] <- mean(dgp2.alasso.current.fit$mspe)
  dgp2.alasso.beta.logn2[[i]] <- dgp2.alasso.current.fit$beta
  dgp2.alasso.sr.logn2[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                              beta.est = dgp2.alasso.current.fit$beta, n.sim = 1000)
}

dgp2.alasso.mspe.logn3 <- rep(0, length(n.settings))
dgp2.alasso.beta.logn3 <- list()
dgp2.alasso.sr.logn3 <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.alasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.alasso.current.fit <- fun.fit.lasso.all(get(dgp2.alasso.current.setting),
                                               c_lambda = dgp2.alasso.c_lambda,
                                               b_n = 'log(n) ^ 3',
                                               pen.factor = TRUE,
                                               scale = sim.scale)
  dgp2.alasso.mspe.logn3[i] <- mean(dgp2.alasso.current.fit$mspe)
  dgp2.alasso.beta.logn3[[i]] <- dgp2.alasso.current.fit$beta
  dgp2.alasso.sr.logn3[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                              beta.est = dgp2.alasso.current.fit$beta, n.sim = 1000)
}

#### LASSO ----

dgp2.lasso.mspe.n <- rep(0, length(n.settings))
dgp2.lasso.beta.n <- list()
dgp2.lasso.sr.n <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.lasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.lasso.current.fit <- fun.fit.lasso.all(get(dgp2.lasso.current.setting),
                                              c_lambda = dgp2.lasso.c_lambda,
                                              b_n = 'n',
                                              pen.factor = FALSE,
                                              scale = sim.scale)
  dgp2.lasso.mspe.n[i] <- mean(dgp2.lasso.current.fit$mspe)
  dgp2.lasso.beta.n[[i]] <- dgp2.lasso.current.fit$beta
  dgp2.lasso.sr.n[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                         beta.est = dgp2.lasso.current.fit$beta, n.sim = 1000)
}

dgp2.lasso.mspe.2n <- rep(0, length(n.settings))
dgp2.lasso.beta.2n <- list()
dgp2.lasso.sr.2n <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.lasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.lasso.current.fit <- fun.fit.lasso.all(get(dgp2.lasso.current.setting),
                                              c_lambda = dgp2.lasso.c_lambda,
                                              b_n = 'n ^ (1 / 2)',
                                              pen.factor = FALSE,
                                              scale = sim.scale)
  dgp2.lasso.mspe.2n[i] <- mean(dgp2.lasso.current.fit$mspe)
  dgp2.lasso.beta.2n[[i]] <- dgp2.lasso.current.fit$beta
  dgp2.lasso.sr.2n[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                          beta.est = dgp2.lasso.current.fit$beta, n.sim = 1000)
}

dgp2.lasso.mspe.3n <- rep(0, length(n.settings))
dgp2.lasso.beta.3n <- list()
dgp2.lasso.sr.3n <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.lasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.lasso.current.fit <- fun.fit.lasso.all(get(dgp2.lasso.current.setting),
                                              c_lambda = dgp2.lasso.c_lambda,
                                              b_n = 'n ^ (1 / 3)',
                                              pen.factor = FALSE,
                                              scale = sim.scale)
  dgp2.lasso.mspe.3n[i] <- mean(dgp2.lasso.current.fit$mspe)
  dgp2.lasso.beta.3n[[i]] <- dgp2.lasso.current.fit$beta
  dgp2.lasso.sr.3n[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                          beta.est = dgp2.lasso.current.fit$beta, n.sim = 1000)
}

#### OLS ----
## OLS
dgp2.ols.mspe <- rep(0, length(n.settings))
dgp2.ols.beta <- list()
for (i in 1:length(n.settings)){
  dgp2.ols.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.ols.current.fit <- fun.fit.ols.all(get(dgp2.ols.current.setting), scale = sim.scale)
  dgp2.ols.mspe[i] <- dgp2.ols.current.fit$mspe
  dgp2.ols.beta[[i]] <- dgp2.ols.current.fit$beta
}
## Oracle OLS
dgp2.oracle.mspe <- rep(0, length(n.settings))
dgp2.oracle.beta <- list()
for (i in 1:length(n.settings)){
  dgp2.oracle.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.oracle.current.fit <- fun.fit.oracle.all(get(dgp2.oracle.current.setting),
                                                beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                                scale = sim.scale)
  dgp2.oracle.mspe[i] <- dgp2.oracle.current.fit$mspe
  dgp2.oracle.beta[[i]] <- dgp2.oracle.current.fit$beta
}


# Print results
dgp2.mspe.result.df <- t(data.frame(alasso.mspe = dgp2.alasso.mspe, lasso.mspe = dgp2.lasso.mspe, 
                                    ols.mspe = dgp2.ols.mspe, oracle.mspe = dgp2.oracle.mspe))
colnames(dgp2.mspe.result.df) <- paste(n.settings)

#### Preparing tables ----

table1.mspe <- data.frame(n = n.settings,
                          Oracle = dgp2.oracle.mspe,
                          OLS = dgp2.ols.mspe,
                          alasso.loglogn = dgp2.alasso.mspe.loglogn,
                          alasso.logn = dgp2.alasso.mspe.logn,
                          alasso.logn2 = dgp2.alasso.mspe.logn2,
                          alasso.logn3 = dgp2.alasso.mspe.logn3,
                          lasso.n = dgp2.lasso.mspe.n,
                          lasso.2n = dgp2.lasso.mspe.2n,
                          lasso.3n = dgp2.lasso.mspe.3n)

table2.alasso.sr <- data.frame(n = n.settings,
                               SR.loglogn = dgp2.alasso.sr.loglogn[, 1],
                               SR.logn = dgp2.alasso.sr.logn[, 1],
                               SR.logn2 = dgp2.alasso.sr.logn2[, 1],
                               SR.logn3 = dgp2.alasso.sr.logn3[, 1])

table2.alasso.sr1 <- data.frame(n = n.settings,
                                SR1.loglogn = dgp2.alasso.sr.loglogn[, 2],
                                SR1.logn = dgp2.alasso.sr.logn[, 2],
                                SR1.logn2 = dgp2.alasso.sr.logn2[, 2],
                                SR1.logn3 = dgp2.alasso.sr.logn3[, 2])

table2.alasso.sr2 <- data.frame(n = n.settings,
                                SR2.loglogn = dgp2.alasso.sr.loglogn[, 3],
                                SR2.logn = dgp2.alasso.sr.logn[, 3],
                                SR2.logn2 = dgp2.alasso.sr.logn2[, 3],
                                SR2.logn3 = dgp2.alasso.sr.logn3[, 1])

table3.lasso.sr <- data.frame(n = n.settings,
                              SR.n = dgp2.lasso.sr.n[, 1],
                              SR.2n = dgp2.lasso.sr.2n[, 1],
                              SR.3n = dgp2.lasso.sr.3n[, 1])

table3.lasso.sr1 <- data.frame(n = n.settings,
                               SR1.n = dgp2.lasso.sr.n[, 2],
                               SR1.2n = dgp2.lasso.sr.2n[, 2],
                               SR1.3n = dgp2.lasso.sr.3n[, 2])

table3.lasso.sr2 <- data.frame(n = n.settings,
                               SR2.n = dgp2.lasso.sr.n[, 3],
                               SR2.2n = dgp2.lasso.sr.2n[, 3],
                               SR2.3n = dgp2.lasso.sr.3n[, 3])

#### Displaying results ----
print(table1.mspe)
print(cbind(table2.alasso.sr, table2.alasso.sr1[, 2:5], table2.alasso.sr2[, 2:5]))
print(cbind(table3.lasso.sr, table3.lasso.sr1[, 2:4], table3.lasso.sr2[, 2:4]))


