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

### Functions:
## Please set working directory to the location of this file before running
source('myfun.R')
source('dgp.R')
### Import and clean Data
data.original <- read_excel("PredictorData2018.xlsx")

# Convert dates
data.original[[1]] <- as.Date(as.character(paste(data.original[[1]], '01', sep = '')), format = "%Y%m%d")
data <- subset(data.original, data.original$yyyymm >= as.Date('1927-01-01'))
date.start.year <- as.numeric(format(data[[1]][1], '%Y'))
date.start.month <- as.numeric(format(data[[1]][1], '%m'))
date.end.year <- as.numeric(format(data[[1]][length(data[[1]])], '%Y'))
date.end.month <- as.numeric(format(data[[1]][length(data[[1]])], '%m'))
# Convert data to time series
for (i in 2:length(data)) {
  data[[i]] <- ts(as.numeric(data[[i]]), start = c(date.start.year, date.start.month), frequency = 12)
}

### Preparing variables
lag = 1
# Dependent Var
var.premium <- data$CRSP_SPvw - data$Rfree
var.premium <- var.premium[- length(var.premium)]
n.obs <- length(var.premium)

# Regressors
var.dp <- log(data$D12) - log(data$Index) # d/p ratio
var.dy <- log(data$D12) - log(c(data$Index[-1], data$Index[n.obs])) # d/y ratio
var.ep <- log(data$E12) - log(data$Index) # e/p ratio
var.de <- log(data$D12) - log(data$E12) # d/e ratio
var.svar <- data$svar # stock variance
# var.csp <- data$csp # cross-sectional premium, not included, insufficient sample size
var.bm <- data$`b/m` # book to market ratio
var.ntis <- data$ntis # net equity expansion
# var.eqis # percent equity issuing # not included, insufficient sample size
var.tbl <- data$tbl # treasury bill/T-bill rate
var.lty <- data$lty # long term government bond yield
var.ltr <- data$ltr # long term government bond return
var.tms <- var.lty - var.tbl # term spread
var.dfy <- data$BAA - data$AAA # default yield spread
var.dfr <- data$corpr - data$ltr # default return spread
var.infl <- data$infl # inflation from Consumer Price Index

var.all <- cbind(var.dp, var.dy, var.ep, var.de, var.svar, var.bm, var.ntis, var.tbl,
                        var.lty, var.ltr, var.tms, var.dfy, var.dfr, var.infl) # design matrix of all regressors
var.all.lag <- var.all[- 1, ]

## Testing for stationarity/cointegration using ADF test
adf.pval.premium <- adf.test(var.premium, 'stationary', k = 1)
adf.pval.all <- apply(var.all, 2, function(x) adf.test(x, 'stationary', k = 1)$p.value)
ar1.coef.all <- apply(var.all, 2, function(x) ar(x, order.max = 1)$ar)
ar1.coef.premium <- ar(var.premium, order.max = 1)$ar
# Print results
print(adf.pval.all)
print(ar1.coef.all)
print(ar1.coef.premium)
# --> Test result is quite promising, as the null of unit root is rejected in 12/14 variables at 5%, and premium is stationary
adf.test(diff(var.tbl))
adf.test(diff(var.lty))
# Testing 1st difference of these 2 non-stationary time series confirms that both are of I(1)
# We will now also test for cointegration between the two: var.tbl and var.lty
adf.test(lm(var.tbl ~ var.lty)$residuals, k = 1)
# Residual series generated from simple regression between 2 var is stationary, implying cointegration property

## Detecting collinearity
fit.ols <- lm(var.premium ~ var.all.lag, singular.ok = TRUE)
summary(fit.ols)
# -> as expected from variable construction, var.de and var.tms is perfectly collinear with other variables
# We will try to compute VIF
var.all.vif.score <- foreach(i = 1:ncol(var.all), .combine = 'c') %do% {
  1 / (1 - summary(lm(var.all[, i] ~ var.all[, - i]))$r.squared)
}
print(var.all.vif.score)

# Compute eigenvalue and trace of design matrix C_n
C_n <- 1 / length(var.premium) * t(var.all) %*% var.all
var.all.eig <- eigen(C_n, only.values = TRUE)$values
print(min(abs(var.all.eig)))
print(sum(diag(C_n)))
# --> near-singularity detected

### Performing LASSO/ALASSO
## Prepare and scale variables for each settings
## DEBUGGING PART. PLEASE REFER TO TEST FILE. 
## RUN THE TEST FILE BEFORE RUNNING THIS TO SEE THE PROBLEMS WITH DATA
var.all.lag[186:245, 8]
tbl.alt<- rep(NA, 60)
for (i in 60:1){
  tbl.alt[i] <- mean(var.all.lag[(186 - i - 1):(245 - i - 1), 8])
}
var.all.lag[186:245, 8] <- rev(tbl.alt)
## END

# Prepare data (transform horizons) for estimation
y <- var.premium
n <- length(y)
y.1_12 <- var.premium
y.1_4 <- fun.horizon.transform(y.1_12, 1/4)
y.1_2 <- fun.horizon.transform(y.1_12, 1/2)
y.1 <- fun.horizon.transform(y.1_12, 1)
y.2 <- fun.horizon.transform(y.1_12, 2)
y.3 <- fun.horizon.transform(y.1_12, 3)

x <- var.all.lag
x.1_12 <- var.all.lag

x.1_4 <- x.1_12[1:length(y.1_4), ]
x.1_2 <- x.1_12[1:length(y.1_2), ]
x.1 <- x.1_12[1:length(y.1), ]
x.2 <- x.1_12[1:length(y.2), ]
x.3 <- x.1_12[1:length(y.3), ]

horizon.setting <- c('1_12', '1_4', '1_2', '1', '2', '3')

# Finding optimal lambda by 10-fold cross-validation
# Another consequence of singularity is that OLS cannot be used for estimation of weights for the adaptive LASSO.
# Therefore, ridge regression estimates will be used as a replacement due to them being good approximation
#   of non-regularized least square estimates in singular design with lambda near zero (Knight and Fu, 2000).
# Test ridge fit

coefficients(lm(y ~ x))
fit.ridge <- glmnet(x, y, alpha = 0)
coef.ridge <- coef(fit.ridge, s = 1e-10, exact = TRUE, x = x, y = y)[- 1] # omitting intercept
coef.ridge

# The estimated coefficients for ridge regression is stable and can be used to compute penalty factor for alasso
## Next we will use 10-fold rolling windows cross-validation to find optimal initial value of lambda for our problem at hand.
# As a first step, we will try to examine the path of cv.mspe given lambda path.

emp.scale <- FALSE # set scaling to FALSE

tic()
alasso.optim.mse <- sapply(seq(from = 0.0001, to = 0.5, length.out = 200), 
                           function(s) fun.cv.lasso(y.1, x.1, kf = 5, lambda = s, pen.factor = TRUE, pen.type = 'ridge', scale = emp.scale))
alasso.lambda.min <- seq(from = 0.0001, to = 0.5, length.out = 200)[which.min(alasso.optim.mse)]
plot(seq(from = 0.0001, to = 0.5, length.out = 200), alasso.optim.mse, type = 'l')
which.min(alasso.optim.mse)
print(alasso.lambda.min)
toc()

# As we can see, the mspe is decreasing very fast with increasing lambda and stay flat after certain values (of lambda)
# Therefore, optimization algorithms is very volatile and dependent on given initial value.
# After the examination above, we will set initial value for optimization to a value close to min (0.01).
# Reader can try and change the initial value below to see the volatility


alasso.optim <- optim(par = alasso.lambda.min, 
                      function(a) fun.cv.lasso(y.1, x.1, 
                                               kf = 5, lambda = a, 
                                               pen.factor = TRUE, pen.type = 'ridge',
                                               scale = emp.scale),
                      method = "L-BFGS-B", lower = 0.0001)
# toc()
print(alasso.optim$par)
print(alasso.optim$message)

## Getting to the real beef. Here, parallel computing is used to speen up the execution of codes

tic()
cl <- parallel::makeCluster(6) # default is 6 because there are 6 horizon setting that we want to run in parallel 
parallel::setDefaultCluster(cl = cl)

# exporting variables and packages to all nodes
envir <- environment(fun.lasso.predict)
parallel::clusterExport(cl, varlist = ls(envir), envir = envir)

parallel::clusterEvalQ(cl, { 
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
})

toc()

#### Start execution.
#### WARNING: This takes very long to complete: around 8 hours for 6 clusters running in parallel
tic()
par.alasso.result.10 <- parLapply(cl, horizon.setting,
                               function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                             pen.factor = TRUE, pen.type = 'ridge', 
                                                             scale = emp.scale))
toc()

tic()
par.lasso.result.10 <- parLapply(cl, horizon.setting,
                               function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                             pen.factor = FALSE, 
                                                             scale = emp.scale))
toc()

tic()
par.alasso.result.15 <- parLapply(cl, horizon.setting,
                                  function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                                pen.factor = TRUE, pen.type = 'ridge', window = 15,
                                                                scale = emp.scale))
toc()

tic()
par.lasso.result.15 <- parLapply(cl, horizon.setting,
                                 function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                               pen.factor = FALSE, window = 15,
                                                               scale = emp.scale))
toc()

tic()
par.alasso.result.20 <- parLapply(cl, horizon.setting,
                                  function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                                pen.factor = TRUE, pen.type = 'ridge', window = 20,
                                                                scale = emp.scale))
toc()

tic()
par.lasso.result.20 <- parLapply(cl, horizon.setting,
                                 function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                               pen.factor = FALSE, window = 20,
                                                               scale = emp.scale))
toc()
stopCluster(cl) # closing nodes

## Fetching neccesary resultsL MPSE
alasso.mspe.10 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.alasso.result.10[[i]]$mspe)
lasso.mspe.10 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.lasso.result.10[[i]]$mspe)
alasso.mspe.15 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.alasso.result.15[[i]]$mspe)
lasso.mspe.15 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.lasso.result.15[[i]]$mspe)
alasso.mspe.20 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.alasso.result.20[[i]]$mspe)
lasso.mspe.20 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.lasso.result.20[[i]]$mspe)



alasso.sel.rate.10 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.alasso.result.10[[i]]$mspe)


### Performing OLS


par.ols.result.10 <- lapply(horizon.setting, 
                            function(a) fun.ols.emp(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                    window = 10, horizon = 1,
                                                    scale = emp.scale))
par.ols.result.15 <- lapply(horizon.setting, 
                            function(a) fun.ols.emp(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                    window = 15, horizon = 1,
                                                    scale = emp.scale))
par.ols.result.20 <- lapply(horizon.setting, 
                            function(a) fun.ols.emp(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                    window = 20, horizon = 1,
                                                    scale = emp.scale))

ols.mspe.10 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.ols.result.10[[i]]$mspe)
ols.mspe.15 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.ols.result.15[[i]]$mspe)
ols.mspe.20 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.ols.result.20[[i]]$mspe)

### Performing RWwD

par.rw.result.10 <- lapply(horizon.setting, 
                           function(a) fun.rw.emp(get(paste('y.', a, sep = '')), window = 10, horizon = 1))
par.rw.result.15 <- lapply(horizon.setting, 
                           function(a) fun.rw.emp(get(paste('y.', a, sep = '')), window = 15, horizon = 1))
par.rw.result.20 <- lapply(horizon.setting, 
                           function(a) fun.rw.emp(get(paste('y.', a, sep = '')), window = 20, horizon = 1))

rw.mspe.10 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.rw.result.10[[i]]$mspe)
rw.mspe.15 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.rw.result.15[[i]]$mspe)
rw.mspe.20 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.rw.result.15[[i]]$mspe)


#### Preparing results
# window = 10
emp.result.10 <- data.frame(ols = ols.mspe.10, rwwd = rw.mspe.10,
                            alasso = alasso.mspe.10, lasso = lasso.mspe.10)
rownames(emp.result.10) <- paste(horizon.setting)

emp.result.percent.10 <- data.frame(rel.ols = ols.mspe.10 / ols.mspe.10, rel.rwwd = rw.mspe.10 / ols.mspe.10,
                                    rel.alasso = alasso.mspe.10 / ols.mspe.10, rel.lasso = lasso.mspe.10 / ols.mspe.10)
rownames(emp.result.percent.10) <- paste(horizon.setting)
# window = 15
emp.result.15 <- data.frame(ols = ols.mspe.15, rwwd = rw.mspe.15,
                            alasso = alasso.mspe.15, lasso = lasso.mspe.15)
rownames(emp.result.15) <- paste(horizon.setting)

emp.result.percent.15 <- data.frame(rel.ols = ols.mspe.15 / ols.mspe.15, rel.rwwd = rw.mspe.15 / ols.mspe.15,
                                    rel.alasso = alasso.mspe.15 / ols.mspe.15, rel.lasso = lasso.mspe.15 / ols.mspe.15)
rownames(emp.result.percent.15) <- paste(horizon.setting)
# window = 20
emp.result.20 <- data.frame(ols = ols.mspe.20, rwwd = rw.mspe.20,
                            alasso = alasso.mspe.20, lasso = lasso.mspe.20)
rownames(emp.result.20) <- paste(horizon.setting)

emp.result.percent.20 <- data.frame(rel.ols = ols.mspe.20 / ols.mspe.20, rel.rwwd = rw.mspe.20 / ols.mspe.20,
                                    rel.alasso = alasso.mspe.20 / ols.mspe.20, rel.lasso = lasso.mspe.20 / ols.mspe.20)
rownames(emp.result.percent.20) <- paste(horizon.setting)


#### Printing results
# Table 4: Estimated AR(1) coefficients
print(ar1.coef.all)

# Table 5: MPSE
print(data.frame(emp.result.10, emp.result.percent.10))
print(data.frame(emp.result.15, emp.result.percent.15))
print(data.frame(emp.result.20, emp.result.percent.20))


