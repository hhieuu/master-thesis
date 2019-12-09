#### Main working file
# ## Package installation
# install.packages('readxl')
# install.packages('zoo')
# install.packages('glmnet')
# install.packages('tseries')
# install.packages('caret') 
# install.packages('Matrix')
# install.packages('corrplot')
# install.packages('snow')

# ## Call packages
# library('readxl')
# library('zoo')
# library('glmnet')
# library('Matrix')
# library('tseries')
# library('caret')
# library('Matrix')
# library('corrplot')
# library(parallel)
# library(snow)
# library('MASS')

### Functions:
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
var.csp <- data$csp # cross-sectional premium
var.bm <- data$`b/m` # book to market ratio
var.ntis <- data$ntis # net equity expansion
# var.eqis # percent equity issuing
var.tbl <- data$tbl # treasury bill/T-bill rate from 1920 to 1933
var.lty <- data$lty # long term government bond yield, 1919 to 1925
var.ltr <- data$ltr # long term government bond return, 1926 to 2005
var.tms <- var.lty - var.tbl # term spread
var.dfy <- data$BAA - data$AAA # default yield spread
var.dfr <- data$corpr - data$ltr # default return spread
var.infl <- data$infl # inflation from Consumer Price Index, 1919 to 2005
# var.ik # 
var.all <- cbind(var.dp, var.dy, var.ep, var.de, var.svar, var.bm, var.ntis, var.tbl,
                        var.lty, var.ltr, var.tms, var.dfy, var.dfr, var.infl) # design matrix of all regressors
var.all.lag <- var.all[- 1, ]

## Testing for stationarity/cointegration using ADF test
adf.pval.premium <- adf.test(var.premium, 'stationary', k = 1)
adf.pval.all <- apply(var.all, 2, function(x) adf.test(x, 'stationary', k = 1)$p.value)
# --> Test result is quite promising, as the null of unit root is rejected in 12/14 variables at 5%
adf.test(diff(var.tbl))
adf.test(diff(var.lty))

# Testing 1st difference of these 2 non-stationary time series confirms that both are of I(1)
# We will now also test for cointegration between the two: var.tbl and var.lty
adf.test(lm(var.tbl ~ var.lty)$residuals)
# Residual series generated from simple regression between 2 var is stationary, implying cointegration property

## Detecting collinearity
fit.ols <- lm(var.premium ~ var.all.lag)
summary(fit.ols)
# -> as expected from variable construction, var.de and var.tms is perfectly collinear with other variables
# -> some other variables exhibit high correlation as well, demonstrated in correlogram:
corrplot(cor(var.all), 'number')
# This will, together with LASSO's problems, further invalidate 'standard' inference

## Fitting every series with AR(1) model
fit.ar.coef <- rep(NA, times = ncol(var.all))
for (i in 1:ncol(var.all)){
  fit.ar.coef[i] <- ar(var.all[ , i], order.max = 1)$ar
}


### Performing LASSO
## Prepare and scale variables for each settings


y <- var.premium
n <- length(y)
# y.1_12 <- var.premium[480:n]
y.1_12 <- var.premium[1:160]
y.1_4 <- fun.horizon.transform(y.1_12, 1/4)
y.1_2 <- fun.horizon.transform(y.1_12, 1/2)
y.1 <- fun.horizon.transform(y.1_12, 1)
y.2 <- fun.horizon.transform(y.1_12, 2)
y.3 <- fun.horizon.transform(y.1_12, 3)

x <- scale(var.all.lag, TRUE, TRUE)
# x.1_12 <- scale(var.all.lag, TRUE, TRUE)[480:n, ]
x.1_12 <- scale(var.all.lag, TRUE, TRUE)[1:160, ]
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

fit.ridge <- glmnet(x, y, alpha = 0)
coef.ridge <- coef(fit.ridge, s = 1e-10, exact = TRUE, x = x, y = y)[- 1] # omitting intercept
coef.ridge

# The estimated coefficients for ridge regression is stable and can be used to compute penalty factor for alasso
## Next we will use 10-fold rolling windows cross-validation to find optimal initial value of lambda for our problem at hand.
# As a first step, we will try to examine the path of cv.mspe given lambda path.

tic()
alasso.optim.mse <- sapply(seq(from = 0.0001, to = 0.3, length.out = 200), 
                           function(s) fun.cv.lasso(y.1, x.1, kf = 5, lambda = s, pen.factor = TRUE, pen.type = 'ridge'))
plot(1:length(alasso.optim.mse), alasso.optim.mse, type = 'l')
alasso.lambda.min <- seq(from = 0.0001, to = 0.3, length.out = 200)[which.min(alasso.optim.mse)]
which.min(alasso.optim.mse)
print(alasso.lambda.min)

# As we can see, the mspe is decreasing very fast with increasing lambda and stay flat after certain values (of lambda)
# Therefore, optimization algorithms is very volatile and dependent on given initial value.
# After the examination above, we will set initial value for optimization to a value close to min (0.01).
# Reader can try and change the initial value below to see the volatility

alasso.optim <- optim(par = alasso.lambda.min, 
                      function(a) fun.cv.lasso(y.1, x.1, kf = 5, lambda = a, pen.factor = TRUE, pen.type = 'ridge'),
                      method = "L-BFGS-B", lower = 0.0001)
toc()
print(alasso.optim$par)
print(alasso.optim$message)



# tic()
# test.pred <- fun.lasso.predict(y.1, x.1)
# toc() # 10395 secs to run for full data
# mean(test.pred$mspe) 
# length(test.pred$mspe)

# Using parallel apply to speed up. Still very long to run

cl <- parallel::makeCluster(detectCores() - 1)
# on.exit(parallel::stopCluster(cl), add = TRUE)
envir <- environment(fun.lasso.predict)
parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
parallel::setDefaultCluster(cl = cl)

tic()
par.alasso.result.10 <- parLapply(cl, horizon.setting,
                               function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                             init.val = 0.03))
toc()

tic()
par.lasso.result.10 <- parLapply(cl, horizon.setting,
                               function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                             pen.factor = FALSE, init.val = 0.04))
toc()

tic()
par.alasso.result.12 <- parLapply(cl, horizon.setting,
                                  function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                                window = 12, init.val = 0.03))
toc()

tic()
par.lasso.result.12 <- parLapply(cl, horizon.setting,
                                 function(a) fun.lasso.predict(get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                               pen.factor = FALSE, window = 12, init.val = 0.03))
toc()
stopCluster(cl)

alasso.mspe.10 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.alasso.result.10[[i]]$mspe)
lasso.mspe.10 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.lasso.result.10[[i]]$mspe)
alasso.mspe.12 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.alasso.result.12[[i]]$mspe)
lasso.mspe.12 <- foreach(i = 1:length(horizon.setting), .combine = c) %do% mean(par.lasso.result.12[[i]]$mspe)






