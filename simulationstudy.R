source('myfun.R')

#### Setting for simulation: ----
## Number of different settings: n = 40, 80, 120, 200, 500, 1000
## Number of model for evaluation: OLS, Oracle OLS, AR, Lasso, and Adaptive Lasso
## Evaluation is carried out by computing Mean Squared Prediction Error
## For each setting, the study is replicated 100 times. Final MSPE is the average of all MSPE in all replications
## For Lasso and Alasso, lambda will be chosen by 10 fold cross validation in each replication at n = 200.
## Then the final lambda would be the median of all MSPE-optimized lambda.

#### Simulating data ----
tic("total")
set.seed(3000)
# Simulating data for finding c_lambda. c_lambda is the median of all CV-MSPE-optimized lambdas in each replication.
sim.optim.200 <- replicate(100, fun.dgp2(n.burn = 1000, n.obs = 200))
# Simulating data for 6 settings. The lambda for each setting is c_lambda*sqrt(n/log(log(n)))
sim.40 <- replicate(500, fun.dgp2(n.burn = 1000, n.obs = 40))
sim.80 <- replicate(500, fun.dgp2(n.burn = 1000, n.obs = 80))
sim.120 <- replicate(500, fun.dgp2(n.burn = 1000, n.obs = 120))
sim.200 <- replicate(500, fun.dgp2(n.burn = 1000, n.obs = 200))
sim.400 <- replicate(500, fun.dgp2(n.burn = 1000, n.obs = 400))
sim.1000 <- replicate(500, fun.dgp2(n.burn = 1000, n.obs = 1000))

n.settings <- c(40, 80, 120, 200, 400, 1000)
toc()
#### Alasso performance assessment ----

## Finding optimized c_lambda by 10-folds cross validation for each replication of simulated data.
tic()
lambda.optim.test <- fun.optim.lambda(sim.optim.200, start.val = 0.03)
c_lambda <- median(lambda.optim.test[, 1])
toc()
## Running Alasso on simulated data using obtained c_lambda*sqrt(n/log(log(n)))
## Two main metrics for performance assessment is MSPE and SR (Screening Error)
##  1) For MSPE: we will compute the OOS one-step-ahead squared prediction error for each replication
## in each simulation setting. MSPE will be the mean across all replications in each simulation setting.
##  2) For SR, there are three sub-metrics: 
##    - SR: overall success rate of classification into zero coefficients and non-zero coefficients
##    - SR 1: percentage of the correct selection in the active set
##    - SR 2: percentage of correct elimination of the zero coefficients

# Fitting Alasso for each replication in a setting


n.rep <- 500
current.setting <- 40
current.lambda <- c_lambda * sqrt(current.setting / log(log(current.setting)))
current.sim <- get(paste('sim.', current.setting, sep = ''))
current.list.y <- current.sim[seq(from = 1, to = length(current.sim), by = 2)]
current.list.x <- current.sim[seq(from = 0, to = length(current.sim), by = 2)]
current.index <- fun.train.test.split(1:current.setting, test.portion = .2)
current.list.fit <- list()
current.beta <- Matrix(0, nrow = n.rep, ncol = 9, sparse = T)
current.mspe <- rep(0, 500)
for (i in 1:n.rep){
  current.train.x <- current.list.x[[i]][current.index$train, ]
  current.train.y <- current.list.y[[i]][current.index$train]
  current.test.x <- t(as.matrix(current.list.x[[i]][current.index$test[1], ]))
  current.test.y <- current.list.y[[i]][current.index$test[1]]
  current.penalty <- 1/abs(lm(current.train.y ~ current.train.x)$coefficients[-1])
  current.fit <- fun.fit.lasso(current.train.y, current.train.x, lambda = current.lambda, pen.factor = current.penalty)
  current.beta[i, ] <- current.fit$beta
  current.mspe[i] <- fun.pred.lasso(current.test.y, current.test.x, current.fit)
  current.list.fit[[i]] <- current.fit
}
toc()
mean(current.mspe)
median(current.mspe)
min(current.mspe)


#### OLS ----
# A linear model will be fitted to all replications in each data setting.
# Data will be divided into 

# fun.cv.ols(var.y, var.x, kf = 10)
# split.index <- fun.ts.split(var.y, k = 10)
# fun.pred.ols(var.y[split.index$train[[2]]], var.x[split.index$train[[2]], ],
#              var.y[split.index$val[[2]]], var.x[split.index$val[[2]], ])
# fun.cv.lasso(var.y, var.x, lambda = 0.5, kf = 10)

# Now we will 

#### Testsbereich

ridge.coef <- glmnet(sim.200[[2]], sim.200[[1]], alpha = 0, lambda = 100)$beta@x

y.all <- sim.optim.200[seq(from = 1, to = 200, by = 2)]
x.all <- sim.optim.200[seq(from = 2, to = 200, by = 2)]
lambda.grid <- seq(from = 0, to = 0.5, by = 0.001)

test.apply <- sapply(lambda.grid, function(a) fun.cv.lasso(y.all[[1]], x.all[[1]], kf = 10, lambda = a))
test.optim <- optim(0.02, function(a) fun.cv.lasso(y.all[[1]], x.all[[1]], kf = 10, lambda = a),
                    method = "L-BFGS-B", lower = 0)
plot(lambda.grid, test.apply, type = 'l')

fun.cv.lasso(y.all[[1]], x.all[[1]], kf = 10, lambda = 0.3794631)
glmnet(x.all[[1]], y.all[[1]], lambda = 0.3794631)$beta
n.sim <- length(sim.40)/2
y.all <- sim.40[seq(from = 1, to = n.sim*2, by = 2)]
x.all <- sim.40[seq(from = 2, to = n.sim*2, by = 2)]

ridge.coef <- glmnet(x.all[[1]], y.all[[1]], alpha = 0, lambda = 0)$beta@x
ols.coef <- lm(y.all[[1]] ~ x.all[[1]])$coefficients[-1]
pen.factor <- 1/abs(ridge.coef)
test.fit <- fun.fit.lasso(y.all[[1]], x.all[[1]], lambda = seq(from = 0.01, to = 1, length.out = 200))
beta.mat <- as.matrix(test.fit$beta)

plot(function(a) 0.5 * sqrt(a/log(log(a))), from = 10, to = 1e10 )
plot(function(a) a^(1 / 2), from = 10, to = 1e10, add = T, col = 'red')
plot(function(a) a^(1 / 3), from = 10, to = 1e10, add = T, col = 'blue')
plot(function(a) a^(1 / 3) / a^(1 / 2), from = 10, to = 1e10, add = T, col = 'green')