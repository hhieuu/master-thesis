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
var.premium.lag <- var.premium[-1]
var.premium.stdrd <- scale(var.premium)
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
var.all.lag <- var.all[-1, ]

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
fit.ols <- lm(var.premium.lag ~ var.all.lag)
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
## Preparing lambda path
# scale variables (normalize)
nml <- function(y) sqrt(sum((y-mean(y))^2)/length(y)) #normalization function
var.all.scale <- scale(var.all, scale = apply(var.all, 2, nml))
var.premium.scale <- scale(var.premium, scale = nml(var.premium))

# Calculate lambda path
lambda.max <- max(abs(colSums(t(var.all.scale)%*%var.premium.scale)))/n.obs
lambda.epsilon <- .01 # fraction of lambda max to be lambda min
lambda.n <- 100 # length of lambda path
lambda.path <- round(exp(seq(log(lambda.max), log(lambda.max*lambda.epsilon), 
                             length.out = lambda.n)), digits = 10)

## Cross-validation
Rprof()
cv.expr <- expression(source('myfun.R'))
n.cluster <- detectCores()
cl <- makeCluster(n.cluster)
clusterExport(cl, c('var.premium.lag', 'var.all.lag', 'lambda.path', 'cv.expr'))
clusterEvalQ(cl, eval(cv.expr))
cv.mspe <- parSapply(cl, lambda.path, function(x) fun.cv.lasso(y = var.premium.lag, x = var.all.lag, lambda = x, kf = 10))
stopCluster(cl)
summaryRprof()

# plot cv error
plot(lambda.path, cv.mspe, pch = 20)
lambda.min <- lambda.path[which.min(cv.mspe)]

# replicate weights used in CV
pen.factor <- (abs(coef(lm(var.premium.lag ~ var.all.lag))[-1]))^-1
pen.factor[which(is.na(pen.factor))] <- 1

# refit using optimal lambda for Alasso
lambda.optim <- lambda.min
fit.lasso <- glmnet(var.all.lag, var.premium.lag, lambda = lambda.optim, penalty.factor = pen.factor)
coef(fit.lasso)

coef(glmnet(var.all.lag, var.premium.lag, lambda = 0.0001))
