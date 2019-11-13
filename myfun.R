
fun.ts.split <- function(y, k){ # create appropriate index for cross validation
  n <- length(y)
  k1 <- k + 1
  partition <- split(1:n, sort(factor(rep(1:k1, length = n)))) # split ts into folds
  train <- list()
  val <- list()
  for (i in 1:k){
    train[[i]] <- unlist(partition[1:i], F, F)
    val[[i]] <- partition[[i + 1]]
  }
  return(list(train = train, val = val))
}

fun.mse <- function(y, y.hat){ # Calculating MSE, or MSPE
  n <- length(y)
  mse <- 1/n*sum((y - y.hat)^2)
  return(mse)
}

# fun.pred.lasso <- function(y.train, x.train, y.val, x.val, lambda, pen.factor = NULL){ # calculate mspe for given train set, validation set, lambda
#   x.train <- Matrix(x.train)
#   x.val <- Matrix(x.val)
#   if (is.null(pen.factor)){
#     pen.factor <- rep(1, ncol(x.train))
#   }
#   fit <- glmnet(x.train, y.train, lambda = lambda, penalty.factor = pen.factor)
#   y.val.hat <- predict.glmnet(fit, x.val)
#   mspe <- fun.mse(y.val, y.val.hat)
#   return(mspe)
# }

fun.fit.lasso <- function(y.train, x.train, lambda, pen.factor = NULL){
  x.train <- Matrix(x.train)
  if (is.null(pen.factor)){
    pen.factor <- rep(1, ncol(x.train))
  }
  fit <- glmnet(x.train, y.train, lambda = lambda, penalty.factor = pen.factor)
  return(fit)
}

fun.pred.lasso <- function(y.val, x.val, fitted){ # calculate mspe for given fit, validation set, lambda
  y.val.hat <- predict.glmnet(fitted, x.val)
  mspe <- fun.mse(y.val, y.val.hat)
  return(mspe)
}


fun.cv.lasso <- function(y, x, lambda, kf){ # bring it all to calculate k-fold CV mspe
  n.obs <- length(y)
  cv.index <- fun.ts.split(y, k = kf)
  n.trainset <- length(cv.index$train)
  n.valset <- length(cv.index$val)
  mspe.mat <- rep(1, n.valset)
  for (i in 1:n.trainset) {
    y.train <- y[cv.index$train[[i]]]
    y.val <- y[cv.index$val[[i]]]
    x.train <- x[cv.index$train[[i]], ]
    x.val <- x[cv.index$val[[i]], ]
    pen.factor <- (abs(coef(lm(y ~ x))[-1]))^-1
    pen.factor[which(is.na(pen.factor))] <- 1
    fitted <- fun.fit.lasso(y.train, x.train, lambda, pen.factor = pen.factor)
    mspe.mat[i] <- fun.pred.lasso(y.val, x.val, fitted)
  }
  cv.mspe <- mean(mspe.mat)
  return(cv.mspe)
}

fun.pred.ols <- function(y.train, x.train, y.val, x.val){
  x.val <- cbind(rep(1, length(y.val)), x.val)
  fit.coef <- matrix(coef(lm(y.train ~ x.train)))
  y.val.hat <- x.val%*%fit.coef
  mspe <- fun.mse(y.val, y.val.hat)
  return(mspe)
}

fun.cv.ols <- function(y, x, kf){
  n.obs <- length(y)
  cv.index <- fun.ts.split(y, k = kf)
  n.trainset <- length(cv.index$train)
  n.valset <- length(cv.index$val)
  mspe.mat <- rep(1, n.valset)
  for (i in 1:n.trainset) {
    y.train <- y[cv.index$train[[i]]]
    y.val <- y[cv.index$val[[i]]]
    x.train <- x[cv.index$train[[i]], ]
    x.val <- x[cv.index$val[[i]], ]
    mspe.mat[i] <- fun.pred.ols(y.train, x.train, y.val, x.val)
  }
  cv.mspe <- mean(mspe.mat)
  return(cv.mspe)
}

#### Functions for DGPs ----
fun.dgp1 <- function(n.pred = 2, rho = 1, n.obs = 500, n.burn = 500, intercept, beta){
  beta <- matrix(beta)
  n.tot <- n.burn + n.obs + 1
  x.mat <- matrix(1, ncol = n.pred, nrow = n.tot) # create matrix for storage
  x.mat[1, ] <- rnorm(n.pred, mean = 0, sd = 10) # initial value
  err.mat <- matrix(rnorm(n.tot*(n.pred + 1), mean = 0, sd = 1), ncol = n.pred + 1, nrow = n.tot) # create matrix for error term
  for (i in 2:n.tot){
    pred.mat[i, ] <- rho*x.mat[i - 1, ] + err.mat[i, -1]
  }
  x.mat.lag <- pred.mat[- n.tot, ]
  y <- intercept + x.mat.lag%*%beta + err.mat[-1, 1]
  return(list(y = y[n.burn:(n.tot - 1)], x = x.mat.lag[n.burn:(n.tot - 1), ]))
}

fun.dgp2 <- function(n.burn = 500, n.obs = 500){
  n.tot <- n.burn + n.obs
  
  ## dgp for independent AR processes (no collinearity)
  n.ar <- 3
  rho.ar <- matrix(c(0.99, 0.90, -0.5)) # vector of ar coefficent
  x.ar <- matrix(NA, ncol = n.ar, nrow = n.tot)
  x.ar[1, ] <- rnorm(n.ar) # initial value
  x.ar.err <- matrix(rnorm(n.ar*n.tot), ncol = n.ar, nrow = n.tot) # error matrix
  for (i in 2:n.tot){
    x.ar[i, ] <- rho.ar*x.ar[i - 1, ] + x.ar.err[i, ] # simulating AR(1) variables
  }
  
  ## dgp for 2 highly persistent & collinear AR(1) process
  n.coll <- 2
  rho.coll <- matrix(c(0.98, 0.85))
  x.coll <- matrix(NA, ncol = n.coll, nrow = n.tot)
  x.coll[1, ] <- rnorm(n.coll) # init values
  x.coll.err <- matrix(mvrnorm(n.coll*n.tot, mu = c(0, 0), Sigma = matrix(c(1, 0.95, 0.95, 1), 2, 2)),
                       ncol = n.coll, nrow = n.tot) # error matrix. Error term is highly correlated (0.95)
  for (i in 2:n.tot){
    x.coll[i, ] <- rho.coll*x.coll[i - 1, ] + x.coll.err[i, ] # simulating AR(1) variables
  }
  
  ## dgp for a VECM with 4 I(1) variables and cointegration rank 2
  n.vecm <- 4
  r.vecm <- 2
  rho.nu <- c(0.2, 0.2)
  mat.coin <- matrix(c(1, -1, 0, 0, 0, 0, 1, -1), nrow = 2, ncol = 4, byrow = T) # cointegrating matrix
  mat.load <- matrix(c(0, 1, 0, 0, 0, 0, 0, 1), nrow = 2, ncol = 4, byrow = T) # loading matrix
  x.vecm <- matrix(NA, ncol = n.vecm, nrow = n.tot) # init matrix
  x.vecm[1, ] <- rnorm(n.vecm) # init values
  x.vecm.err <- matrix(NA, ncol = n.vecm, nrow = n.tot) # init matrix for error terms
  x.vecm.err.eps <- matrix(rnorm((n.vecm - r.vecm)*n.tot), ncol = n.vecm - r.vecm, nrow = n.tot)
  x.vecm.err.nu <- matrix(NA, ncol = n.vecm - r.vecm, nrow = n.tot)
  x.vecm.err.nu[1, ] <- rnorm(r.vecm)
  x.vecm.err.nu.err <- matrix(rnorm(r.vecm*n.tot), ncol = r.vecm, nrow = n.tot)
  for (i in 2:n.tot){ # calculating matrix eps2 = eps1 - nu1 and eps4 = eps3 - nu2
    x.vecm.err.nu[i, ] <- rho.nu*x.vecm.err.nu[i - 1, ] + x.vecm.err.nu.err[i, ]
  }
  for (i in 0:(r.vecm - 1)){ # arrange matrix of all error terms
    x.vecm.temp1 <- cbind(x.vecm.err.eps[ , i + 1], x.vecm.err.eps[ , i + 1] - x.vecm.err.nu[ , i + 1])
    x.vecm.err[ , (2*i + 1):(2*i + 2)] <- x.vecm.temp1
  }
  for (i in 2:n.tot){ # simulating I(1) process with 4 elements and cointegration rank 2
    x.vecm[i, ] <- x.vecm[i - 1, ] + t(mat.load)%*%mat.coin%*%x.vecm[i - 1, ] + x.vecm.err[i, ]
  }
  
  ## bring it all together for a x matrix and create y
  x <- cbind(x.ar, x.coll, x.vecm)
  y.intercept <- 0.3
  y.coef <- matrix(c(0, -0.2, 1/sqrt(n.obs), 0.6, 0, 0.3, -0.3, 0, 0.5), nrow = 1)
  y.err <- rnorm(n.tot)
  y <- y.intercept + x%*%t(y.coef) + y.err
  
  return(list(y = y[(n.burn + 1):n.tot], x = x[(n.burn + 1):n.tot, ]))
}



