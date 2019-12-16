source('myfun.R')
source('dgp.R')
#### Setting for simulation: ----
## Number of different settings: n = 40, 80, 120, 200, 500, 1000
## Number of model for evaluation: OLS, Oracle OLS, AR, Lasso, and Adaptive Lasso
## Evaluation is carried out by computing Mean Squared Prediction Error
## For each setting, the study is replicated 100 times. Final MSPE is the average of all MSPE in all replications
## For Lasso and Alasso, lambda will be chosen by 10 fold cross validation in each replication at n = 200.
## Then the final lambda would be the median of all MSPE-optimized lambda.
sim.scale = FALSE
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

## Running Alasso on simulated data using obtained c_lambda*sqrt(n/log(log(n)))
## Two main metrics for performance assessment is MSPE and SR (Screening Error)
##    1) For MSPE: we will compute the OOS one-step-ahead squared prediction error for each replication
## in each simulation setting. MSPE will be the mean across all replications in each simulation setting.
##    2) For SR, there are three sub-metrics: 
##      - SR: overall success rate of classification into zero coefficients and non-zero coefficients
##      - SR 1: percentage of the correct selection in the active set
##      - SR 2: percentage of correct elimination of the zero coefficients

## Fitting for each replication in a setting and analyse perfomance metrics

# Adaptive Lasso
tic()
dgp2.alasso.mspe <- rep(0, length(n.settings))
dgp2.alasso.beta <- list()
dgp2.alasso.success.rate <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.alasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.alasso.current.fit <- fun.fit.lasso.all(get(dgp2.alasso.current.setting),
                                               c_lambda = dgp2.alasso.c_lambda, 
                                               pen.factor = TRUE,
                                               scale = sim.scale)
  dgp2.alasso.mspe[i] <- mean(dgp2.alasso.current.fit$mspe)
  dgp2.alasso.beta[[i]] <- dgp2.alasso.current.fit$beta
  dgp2.alasso.success.rate[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                                  beta.est = dgp2.alasso.current.fit$beta, n.sim = 1000)
}
toc()

# Plain Lasso
tic()
dgp2.lasso.mspe <- rep(0, length(n.settings))
dgp2.lasso.beta <- list()
dgp2.lasso.success.rate <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp2.lasso.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.lasso.current.fit <- fun.fit.lasso.all(get(dgp2.lasso.current.setting),
                                              c_lambda = dgp2.lasso.c_lambda, 
                                              pen.factor = FALSE, 
                                              scale = sim.scale)
  dgp2.lasso.mspe[i] <- mean(dgp2.lasso.current.fit$mspe)
  dgp2.lasso.beta[[i]] <- dgp2.lasso.current.fit$beta
  dgp2.lasso.success.rate[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                                 beta.est = dgp2.lasso.current.fit$beta, n.sim = 1000)
}
toc()



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

#### AR ----
dgp2.ar.mspe <- rep(0, length(n.settings))
dgp2.ar.order <- rep(0, length(n.settings))
for (i in 1:length(n.settings)){
  dgp2.ar.current.setting <- paste('dgp2.sim.', n.settings[i], sep = '')
  dgp2.ar.current.fit <- fun.fit.ar.all(get(dgp2.ar.current.setting), horizon = 1)
  dgp2.ar.mspe[i] <- dgp2.ar.current.fit$mspe
  dgp2.ar.order[i] <- dgp2.ar.current.fit$order
}

# Print results
dgp2.mspe.result.df <- t(data.frame(alasso.mspe = dgp2.alasso.mspe, lasso.mspe = dgp2.lasso.mspe, 
                                    ols.mspe = dgp2.ols.mspe, oracle.mspe = dgp2.oracle.mspe))
colnames(dgp2.mspe.result.df) <- paste(n.settings)


######################################################################################################
######################## DATA GENERATING PROCESS 3 - HIGH COLLINEARITY ###############################
######################################################################################################
#### Simulating data ----
tic("total")
set.seed(26)
# Simulating data for finding c_lambda. c_lambda is the median of all CV-MSPE-optimized lambdas in each replication.
dgp3.sim.optim.200 <- replicate(100, fun.dgp3(n.burn = 1000, n.obs = 200))
# Simulating data for 6 settings
dgp3.sim.40 <- replicate(1000, fun.dgp3(n.burn = 500, n.obs = 40))
dgp3.sim.80 <- replicate(1000, fun.dgp3(n.burn = 500, n.obs = 80))
dgp3.sim.160 <- replicate(1000, fun.dgp3(n.burn = 500, n.obs = 160))
dgp3.sim.250 <- replicate(1000, fun.dgp3(n.burn = 500, n.obs = 250))
dgp3.sim.500 <- replicate(1000, fun.dgp3(n.burn = 500, n.obs = 500))
dgp3.sim.1000 <- replicate(1000, fun.dgp3(n.burn = 500, n.obs = 1000))
# dgp3.sim.10000 <- replicate(1000, fun.dgp3(n.burn = 500, n.obs = 10000))

# n.settings <- c(40, 80, 160, 250, 500, 1000, 10000)
n.settings <- c(40, 80, 160, 250, 500, 1000)
toc()
#### Alasso performance assessment ----
## Finding optimized c_lambda by 10-folds cross validation for each replication of simulated data.

tic()
dgp3.alasso.lambda.optim <- fun.optim.lambda(dgp3.sim.optim.200, start.val = 0.5, 
                                             pen.factor = TRUE, pen.type = 'lm',
                                             scale = sim.scale)
dgp3.alasso.c_lambda <- median(dgp3.alasso.lambda.optim[, 1])
toc()

tic()
dgp3.lasso.lambda.optim <- fun.optim.lambda(dgp3.sim.optim.200, start.val = 0.5, 
                                            pen.factor = FALSE,
                                            scale = sim.scale)
dgp3.lasso.c_lambda <- median(dgp3.lasso.lambda.optim[, 1])
toc()

## Fitting for each replication in a setting and analyse perfomance metrics

# Adaptive Lasso
tic()
dgp3.alasso.mspe <- rep(0, length(n.settings))
dgp3.alasso.beta <- list()
dgp3.alasso.success.rate <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp3.alasso.current.setting <- paste('dgp3.sim.', n.settings[i], sep = '')
  dgp3.alasso.current.fit <- fun.fit.lasso.all(get(dgp3.alasso.current.setting),
                                               c_lambda = dgp3.alasso.c_lambda, 
                                               pen.factor = TRUE, pen.type = 'lm',
                                               scale = sim.scale)
  dgp3.alasso.mspe[i] <- mean(dgp3.alasso.current.fit$mspe)
  dgp3.alasso.beta[[i]] <- dgp3.alasso.current.fit$beta
  dgp3.alasso.success.rate[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                                  beta.est = dgp3.alasso.current.fit$beta, n.sim = 1000)
}
toc()

# Plain Lasso
tic()
dgp3.lasso.mspe <- rep(0, length(n.settings))
dgp3.lasso.beta <- list()
dgp3.lasso.success.rate <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp3.lasso.current.setting <- paste('dgp3.sim.', n.settings[i], sep = '')
  dgp3.lasso.current.fit <- fun.fit.lasso.all(get(dgp3.lasso.current.setting),
                                              c_lambda = dgp3.lasso.c_lambda, 
                                              pen.factor = FALSE,
                                              scale = sim.scale)
  dgp3.lasso.mspe[i] <- mean(dgp3.lasso.current.fit$mspe)
  dgp3.lasso.beta[[i]] <- dgp3.lasso.current.fit$beta
  dgp3.lasso.success.rate[i, ] <- fun.perf.lasso(beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0),
                                                 beta.est = dgp3.lasso.current.fit$beta, n.sim = 1000)
}
toc()



#### OLS ----
## OLS
dgp3.ols.mspe <- rep(0, length(n.settings))
dgp3.ols.beta <- list()
for (i in 1:length(n.settings)){
  dgp3.ols.current.setting <- paste('dgp3.sim.', n.settings[i], sep = '')
  dgp3.ols.current.fit <- fun.fit.ols.all(get(dgp3.ols.current.setting), 
                                          scale = sim.scale)
  dgp3.ols.mspe[i] <- dgp3.ols.current.fit$mspe
  dgp3.ols.beta[[i]] <- dgp3.ols.current.fit$beta
}

dgp3.oracle.mspe <- rep(0, length(n.settings))
dgp3.oracle.beta <- list()
for (i in 1:length(n.settings)){
  dgp3.oracle.current.setting <- paste('dgp3.sim.', n.settings[i], sep = '')
  dgp3.oracle.current.fit <- fun.fit.oracle.all(get(dgp3.oracle.current.setting),
                                                beta.true = c(0, - 0.1, 1 / sqrt(n.settings[i]), 0, 0.4, 0.2, - 0.2, 0, 0), 
                                                scale = sim.scale)
  dgp3.oracle.mspe[i] <- dgp3.oracle.current.fit$mspe
  dgp3.oracle.beta[[i]] <- dgp3.oracle.current.fit$beta
}

#### AR ----
dgp3.ar.mspe <- rep(0, length(n.settings))
dgp3.ar.order <- rep(0, length(n.settings))
for (i in 1:length(n.settings)){
  dgp3.ar.current.setting <- paste('dgp3.sim.', n.settings[i], sep = '')
  dgp3.ar.current.fit <- fun.fit.ar.all(get(dgp3.ar.current.setting), horizon = 1)
  dgp3.ar.mspe[i] <- dgp3.ar.current.fit$mspe
  dgp3.ar.order[i] <- dgp3.ar.current.fit$order
}

# Print results
dgp3.mspe.result.df <- t(data.frame(alasso.mspe = dgp3.alasso.mspe, lasso.mspe = dgp3.lasso.mspe, 
                                    ols.mspe = dgp3.ols.mspe, oracle.mspe = dgp3.oracle.mspe))
colnames(dgp3.mspe.result.df) <- paste(n.settings)



######################################################################################################
########################## DATA GENERATING PROCESS 4 - AS IN LEE(2018) ###############################
######################################################################################################
#### Simulating data ----
tic("total")
set.seed(26)
# Simulating data for finding c_lambda. c_lambda is the median of all CV-MSPE-optimized lambdas in each replication.
dgp4.sim.optim.200 <- replicate(100, fun.dgp4(n.burn = 500, n.obs = 100))
# Simulating data for 6 settings
dgp4.sim.40 <- replicate(1000, fun.dgp4(n.burn = 500, n.obs = 40))
dgp4.sim.80 <- replicate(1000, fun.dgp4(n.burn = 500, n.obs = 80))
dgp4.sim.160 <- replicate(1000, fun.dgp4(n.burn = 500, n.obs = 160))
dgp4.sim.250 <- replicate(1000, fun.dgp4(n.burn = 500, n.obs = 250))
dgp4.sim.500 <- replicate(1000, fun.dgp4(n.burn = 500, n.obs = 500))
dgp4.sim.1000 <- replicate(1000, fun.dgp4(n.burn = 500, n.obs = 1000))
# dgp4.sim.10000 <- replicate(1000, fun.dgp4(n.burn = 500, n.obs = 10000))

# n.settings <- c(40, 80, 160, 250, 500, 1000, 10000)
n.settings <- c(40, 80, 160, 250, 500, 1000)
toc()


# Scaling variables



#### Alasso performance assessment ----
## Finding optimized c_lambda by 10-folds cross validation for each replication of simulated data.

tic()
dgp4.alasso.lambda.optim <- fun.optim.lambda(dgp4.sim.optim.200, start.val = 0.5, 
                                             pen.factor = TRUE, pen.type = 'lm',
                                             scale = sim.scale)
dgp4.alasso.c_lambda <- median(dgp4.alasso.lambda.optim[, 1])
toc()

tic()
dgp4.lasso.lambda.optim <- fun.optim.lambda(dgp4.sim.optim.200, start.val = 0.5, 
                                            pen.factor = FALSE,
                                            scale = sim.scale)
dgp4.lasso.c_lambda <- median(dgp4.lasso.lambda.optim[, 1])
toc()

## Fitting for each replication in a setting and analyse perfomance metrics

# Adaptive Lasso
tic()
dgp4.alasso.mspe <- rep(0, length(n.settings))
dgp4.alasso.beta <- list()
dgp4.alasso.success.rate <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp4.alasso.current.setting <- paste('dgp4.sim.', n.settings[i], sep = '')
  dgp4.alasso.current.fit <- fun.fit.lasso.all(get(dgp4.alasso.current.setting),
                                               c_lambda = dgp4.alasso.c_lambda, 
                                               pen.factor = TRUE,
                                               scale = sim.scale)
  dgp4.alasso.mspe[i] <- mean(dgp4.alasso.current.fit$mspe)
  dgp4.alasso.beta[[i]] <- dgp4.alasso.current.fit$beta
  dgp4.alasso.success.rate[i, ] <- fun.perf.lasso(beta.true = c(0.4, 0, 1 / sqrt(n.settings[i]), 0, 0.3, - 0.3, 0, 0),
                                                  beta.est = dgp4.alasso.current.fit$beta, n.sim = 1000)
}
toc()

# Plain Lasso
tic()
dgp4.lasso.mspe <- rep(0, length(n.settings))
dgp4.lasso.beta <- list()
dgp4.lasso.success.rate <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)){
  dgp4.lasso.current.setting <- paste('dgp4.sim.', n.settings[i], sep = '')
  dgp4.lasso.current.fit <- fun.fit.lasso.all(get(dgp4.lasso.current.setting),
                                              c_lambda = dgp4.lasso.c_lambda, 
                                              pen.factor = FALSE, 
                                              scale = sim.scale)
  dgp4.lasso.mspe[i] <- mean(dgp4.lasso.current.fit$mspe)
  dgp4.lasso.beta[[i]] <- dgp4.lasso.current.fit$beta
  dgp4.lasso.success.rate[i, ] <- fun.perf.lasso(beta.true = c(0.4, 0, 1 / sqrt(n.settings[i]), 0, 0.3, - 0.3, 0, 0),
                                                 beta.est = dgp4.lasso.current.fit$beta, n.sim = 1000)
}
toc()



#### OLS ----
## OLS
dgp4.ols.mspe <- rep(0, length(n.settings))
dgp4.ols.beta <- list()
for (i in 1:length(n.settings)){
  dgp4.ols.current.setting <- paste('dgp4.sim.', n.settings[i], sep = '')
  dgp4.ols.current.fit <- fun.fit.ols.all(get(dgp4.ols.current.setting), 
                                          scale = sim.scale)
  dgp4.ols.mspe[i] <- dgp4.ols.current.fit$mspe
  dgp4.ols.beta[[i]] <- dgp4.ols.current.fit$beta
}

dgp4.oracle.mspe <- rep(0, length(n.settings))
dgp4.oracle.beta <- list()
for (i in 1:length(n.settings)){
  dgp4.oracle.current.setting <- paste('dgp4.sim.', n.settings[i], sep = '')
  dgp4.oracle.current.fit <- fun.fit.oracle.all(get(dgp4.oracle.current.setting),
                                                beta.true = c(0.4, 0, 1 / sqrt(n.settings[i]), 0, 0.3, - 0.3, 0, 0), 
                                                scale = sim.scale)
  dgp4.oracle.mspe[i] <- dgp4.oracle.current.fit$mspe
  dgp4.oracle.beta[[i]] <- dgp4.oracle.current.fit$beta
}

#### AR ----
dgp4.ar.mspe <- rep(0, length(n.settings))
dgp4.ar.order <- rep(0, length(n.settings))
for (i in 1:length(n.settings)){
  dgp4.ar.current.setting <- paste('dgp4.sim.', n.settings[i], sep = '')
  dgp4.ar.current.fit <- fun.fit.ar.all(get(dgp4.ar.current.setting), horizon = 1)
  dgp4.ar.mspe[i] <- dgp4.ar.current.fit$mspe
  dgp4.ar.order[i] <- dgp4.ar.current.fit$order
}

# Store results
dgp4.mspe.result.df <- t(data.frame(alasso.mspe = dgp4.alasso.mspe, lasso.mspe = dgp4.lasso.mspe, 
                                    ols.mspe = dgp4.ols.mspe, oracle.mspe = dgp4.oracle.mspe))
colnames(dgp4.mspe.result.df) <- paste(n.settings)

#### Display results
print(t(dgp2.mspe.result.df))
print(dgp2.alasso.success.rate)
print(dgp2.lasso.success.rate)

print(t(dgp3.mspe.result.df))
print(dgp3.alasso.success.rate)
print(dgp3.lasso.success.rate)

print(t(dgp4.mspe.result.df))
print(dgp4.alasso.success.rate)
print(dgp4.lasso.success.rate)


