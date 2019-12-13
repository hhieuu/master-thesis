#### General functions: ----
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

fun.train.test.split <- function(y, test.portion){ # Create trainset and test (hold-out) set
  n <- length(y)
  n.train <- floor(n - test.portion * n)
  train <- seq(from = 1, to = n.train, by = 1)
  test <- seq(from = n.train + 1, to = n, by = 1)
  return(list(train = train, test = test))
}

fun.mse <- function(y, y.hat){ # Calculating MSE, or MSPE
  n <- length(y)
  mse <- mean((y - y.hat) ^ 2)
  return(mse)
}

fun.horizon.transform <- function(y, horizon = 1/12){
  h <- 12 * horizon - 1
  n.new <- length(y) - h
  y.new <- rep(NA, n.new)
  for (i in 1:n.new){
    y.new[i] <- sum(y[i:(i + h)])
  }
  return(y.new)
}

#### Lasso/alasso functions: ----

fun.cv.lasso <- function(y, x, lambda, kf, 
                         pen.factor = TRUE, pen.type,
                         scale = TRUE){ # calculate k-fold CV mspe
  n.obs <- length(y)
  p <- ncol(x)
  cv.index <- fun.ts.split(y, k = kf)
  n.trainset <- length(cv.index$train)
  n.valset <- length(cv.index$val)
  mspe.mat <- rep(1, n.valset)
  
  for (i in 1:n.trainset) {
    y.train <- y[cv.index$train[[i]]]
    y.val <- y[cv.index$val[[i]]]
    x.train <- x[cv.index$train[[i]], ]
    x.val <- x[cv.index$val[[i]], ]
    # Scale if applicable
    if (scale) {
      train.mean <- colMeans(x.train)
      train.sd <- foreach(j = 1:p, .combine = c) %do% sd(x.train[, j])
      x.train <- scale(x.train, TRUE, TRUE)
      x.val <- scale(x.val, train.mean, train.sd)
    }
    # Check penalty factor
    if (pen.factor){
      if ('lm' %in% pen.type) {
        penalty <- (1 / abs(lm(y.train ~ x.train)$coefficients[-1])) ^ 1
      } else if ('ridge' %in% pen.type) {
        ridge <- glmnet(x.train, y.train, alpha = 0, nlambda = 50, standardize = FALSE)
        ridge.coef <- coef(ridge, s = 1e-8)[- 1]
        penalty <- (1 / abs(ridge.coef)) ^ 1
      } else {
        stop('Error: imput either lm or ridge for penalty factor calculation')
      }
      
      penalty <- p * penalty / sum(penalty)
      x.train <- scale(x.train, FALSE, penalty)
      x.val <- scale(x.val, FALSE, penalty)
    } else {
      penalty <- NULL
    }
    fit <- glmnet(x.train, y.train, nlambda = 50, standardize = FALSE)
    y.hat <- predict.glmnet(fit, x.val, s = lambda, exact = TRUE, 
                            x = x.train, y = y.train)
    mspe.mat[i] <- fun.mse(y.val, y.hat)
  }
  return(mean(mspe.mat))
}



fun.optim.lambda <- function(sample, start.val, 
                             pen.factor = TRUE, pen.type = 'lm',
                             scale = TRUE){
  n.sim <- length(sample)/2
  y.all <- sample[seq(from = 1, to = length(sample), by = 2)]
  x.all <- sample[seq(from = 2, to = length(sample), by = 2)]
  optim.lambda <- matrix(0, ncol = 2, nrow = n.sim)
  for (i in 1:n.sim){
    optim.tmp <- optim(par = start.val, 
                       function(a) fun.cv.lasso(y.all[[i]], x.all[[i]], kf = 10, lambda = a, 
                                                pen.factor = pen.factor, pen.type = 'lm',
                                                scale = scale),
                       method = "L-BFGS-B", lower = 0.001)
    optim.lambda[i, 1] <- optim.tmp$par
    optim.lambda[i, 2] <- optim.tmp$value
  }
  return(optim.lambda)
}


fun.fit.lasso.all <- function(data, c_lambda, 
                              pen.factor = TRUE, pen.type = 'lm', 
                              scale = TRUE) {
  n.rep <- length(data) / 2
  setting <- length(data[[1]])
  p <- nrow(data[[2]])
  test.portion <- .1
  train.portion <- 1 - test.portion
  y.list <- data[seq(from = 1, to = length(data), by = 2)]
  x.list <- data[seq(from = 2, to = length(data), by = 2)]
  index <- fun.train.test.split(1:setting, test.portion = test.portion)
  fit.list <- list()
  beta <- Matrix(0, nrow = n.rep, ncol = ncol(x.list[[1]]), sparse = T)
  mspe <- rep(0, n.rep)
  
  # Lambda
  if (pen.factor) {
    lambda <- c_lambda * (setting * train.portion) ^ (1 / 2) / log(log(setting * train.portion))
    # lambda <- c_lambda * (setting * train.portion) ^ (1 / 2) / log(setting * train.portion)
  } else {
    lambda <- c_lambda * (setting * train.portion) ^ (1 / 2)
  }
  
  for (i in 1:n.rep) {
    y.train <- y.list[[i]][index$train]
    y.test <- y.list[[i]][index$test]
    x.train <- x.list[[i]][index$train, ]
    x.test <- as.matrix(x.list[[i]][index$test, ])
    # Scale if applicable
    if (scale) {
      train.mean <- colMeans(x.train)
      train.sd <- foreach(j = 1:p, .combine = c) %do% sd(x.train[, j])
      x.train <- scale(x.train, TRUE, TRUE)
      x.test <- scale(x.test, train.mean, train.sd)
    }
    # Check penalty factor
    
    if (pen.factor) {
      if ('lm' %in% pen.type) {
        penalty <- (1 / abs(lm(y.train ~ x.train)$coefficients[-1])) ^ 1
      } else if ('ridge' %in% pen.type) {
        ridge <- glmnet(x.train, y.train, alpha = 0, nlambda = 50, standardize = FALSE)
        ridge.coef <- coef(ridge, s = 1e-8)[- 1]
        penalty <- (1 / abs(ridge.coef)) ^ 1
      } else {
        stop('Error: imput either lm or ridge for penalty factor calculation')
      }
      penalty <- p * penalty / sum(penalty)
      x.train <- scale(x.train, FALSE, penalty)
      x.test <- scale(x.test, FALSE, penalty)
    } else {
      penalty <- NULL
    }
    
    fit <- glmnet(x.train, y.train, standardize = FALSE, nlambda = 50)
    y.hat <- predict.glmnet(fit, x.test, s = lambda, # / (setting * train.portion),
                            x = x.train, y = y.train, exact = TRUE)
    mspe[i] <- fun.mse(y.hat, y.test)
    beta[i, ] <- coef(fit, x.test, s = lambda, # / (setting * train.portion), 
                      x = x.train, y = y.train, exact = TRUE)[- 1]
  }
  return(list(mspe = mspe, beta = beta))
}

fun.perf.lasso <- function(beta.true, beta.est, n.sim){
  # Getting index of true zero and nonzero variables
  index.zero.true <- which(beta.true == 0)
  index.nonzero.true <- which(beta.true != 0)
  
  index.zero.est <- lapply(1:n.sim, function(i) which(beta.est[i, ] == 0))
  index.nonzero.est <- lapply(1:n.sim, function(i) which(beta.est[i, ] != 0))
  
  # SR1: Percentage of corrent selection in the active set
  sr1 <- sapply(1:n.sim, function(i) sum(index.nonzero.est[[i]] %in% 
                                           index.nonzero.true))
  
  # SR2: percentage of correct elimination of the zero coefficients
  sr2 <- sapply(1:n.sim, function(i) sum(index.zero.est[[i]] %in%
                                           index.zero.true))
  
  # SR: overall success rate of classification into zero coefficients AND non-zero coefficients
  sr <- sapply(1:n.sim, function(i) sum(index.zero.est[[i]] %in% index.zero.true,
                                        index.nonzero.est[[i]] %in% index.nonzero.true))
  return(matrix(c(mean(sr) / length(beta.true),
                  mean(sr1) / length(index.nonzero.true),
                  mean(sr2)) / length(index.zero.true),
                nrow = 1, ncol = 3))
}

#### OLS functions ----

fun.fit.ols.all <- function(data, scale = TRUE) {
  n.rep <- length(data) / 2
  setting <- length(data[[1]])
  y.list <- data[seq(from = 1, to = length(data), by = 2)]
  x.list <- data[seq(from = 0, to = length(data), by = 2)]
  index <- fun.train.test.split(1:setting, test.portion = .1)
  fit.list <- list()
  beta <- Matrix(0, nrow = n.rep, ncol = ncol(x.list[[1]]) + 1, sparse = T)
  mspe <- rep(0, n.rep)
  
  for (i in 1:n.rep) {
    x.train <- x.list[[i]][index$train, ]
    y.train <- y.list[[i]][index$train]
    y.test <- y.list[[i]][index$test]
    x.test <- cbind(1, x.list[[i]][index$test, ])
    # Scale if applicable
    if (scale) {
      train.mean <- colMeans(x.train)
      train.sd <- foreach(j = 1:p, .combine = c) %do% sd(x.train[, j])
      x.train <- scale(x.train, TRUE, TRUE)
      x.test <- scale(x.test, train.mean, train.sd)
    }
    
    fit <- lm(y.train ~ x.train)
    y.hat <- x.test %*% matrix(coef(fit))
    beta[i, ] <- matrix(fit$coefficients, nrow = 1)
    mspe[i] <- fun.mse(y.test, y.hat)
  }
  return(list(mspe = mean(mspe), beta = beta))
}

fun.fit.oracle.all <- function(data, beta.true) {
  
  n.rep <- length(data) / 2
  setting <- length(data[[1]])
  y.list <- data[seq(from = 1, to = length(data), by = 2)]
  x.list <- data[seq(from = 0, to = length(data), by = 2)]
  index <- fun.train.test.split(1:setting, test.portion = .1)
  index.beta.true <- which(beta.true != 0)
  fit.list <- list()
  beta <- Matrix(0, nrow = n.rep, ncol = length(index.beta.true), sparse = T)
  mspe <- rep(0, n.rep)
  
  for (i in 1:n.rep) {
    x.train <- x.list[[i]][index$train, index.beta.true]
    y.train <- y.list[[i]][index$train]
    y.test <- y.list[[i]][index$test]
    x.test <- cbind(1, x.list[[i]][index$test, index.beta.true])
    # Scale if applicable
    if (scale) {
      train.mean <- colMeans(x.train)
      train.sd <- foreach(j = 1:p, .combine = c) %do% sd(x.train[, j])
      x.train <- scale(x.train, TRUE, TRUE)
      x.test <- scale(x.test, train.mean, train.sd)
    }
    
    fit <- lm(y.train ~ x.train)
    y.hat <- x.test %*% matrix(coef(fit))
    beta[i, ] <- matrix(fit$coefficients[- 1], nrow = 1)
    mspe[i] <- fun.mse(y.test, y.hat)
  }
  return(list(mspe = mean(mspe), beta = beta))
}

#### AR functions ----
fun.fit.ar.all <- function(data, horizon = 1) {
  n.rep <- length(data) / 2
  setting <- length(data[[1]])
  y.list <- data[seq(from = 1, to = length(data), by = 2)]
  index <- fun.train.test.split(1:setting, test.portion = .2)
  fit.list <- list()
  mspe <- rep(0, n.rep)
  order <- rep(0, n.rep)
  
  for (i in 1:n.rep) {
    y.train <- y.list[[i]][index$train]
    y.test <- y.list[[i]][index$test][1:horizon]
    fit <- ar(y.train)
    y.hat <- predict(fit, n.ahead = horizon)$pred
    order[i] <- fit$order
    mspe[i] <- fun.mse(y.test, y.hat)
  }
  return(list(mspe = mean(mspe), order = mean(order)))
}


#### Empirical functions ----

fun.lasso.predict <- function(y, x, pen.factor = TRUE, pen.type = 'ridge', 
                              window = 10, horizon = 1,
                              scale = TRUE){
  ## Splitting data into rolling windows
  index.roll <- createTimeSlices(y, initialWindow = 12 * window, horizon = horizon)
  n.test <- length(index.roll$test)
  p <- ncol(x)
  mspe <- rep(NA, n.test)
  y.hat.mat <- rep(NA, n.test)
  coef <- matrix(NA, nrow = n.test, ncol = p)
  
  for (i in 1:n.test){
    x.train <- x[index.roll$train[[i]], ]
    y.train <- y[index.roll$train[[i]]]
    x.test <- matrix(x[index.roll$test[[i]], ], nrow = horizon)
    y.test <- y[index.roll$test[[i]]]
    # First find lambda by optimizing w.r.t  cv
    lambda.seq <- sapply(seq(from = 0.001, to = 0.3, length.out = 30),
                         function(s) fun.cv.lasso(y.train, x.train, kf = 5, lambda = s,
                                                  pen.factor = pen.factor, pen.type = pen.type,
                                                  scale = scale))
    lambda.start <- seq(from = 0.0001, to = 0.3, length.out = 30)[which.min(lambda.seq)]
    lambda.optim <- optim(par = lambda.start, function(s) fun.cv.lasso(y.train, x.train, lambda = s, kf = 5,
                                                                       pen.factor = pen.factor, pen.type = pen.type),
                          method = "L-BFGS-B", lower = 1e-4)$par
    
    
    
    # Scale if applicable
    if (scale) {
      train.sd <- foreach(i = 1:p, .combine = c) %do% sd(x.train[, i])
      train.mean <- colMeans(x.train)
      x.train <- scale(x.train, TRUE, TRUE)
      x.test <- scale(x.test, center = train.mean, scale = train.sd)
    }
    
    # Find the penalty factor if applicable
    if (pen.factor){
      if ('lm' %in% pen.type) {
        penalty <- (1 / abs(lm(y.train ~ x.train)$coefficients[-1])) ^ 1
      } else if ('ridge' %in% pen.type) {
        ridge <- glmnet(x.train, y.train, alpha = 0, nlambda = 50, standardize = FALSE)
        penalty <- (1 / abs(coef(ridge, s = 1e-10)[- 1])) ^ 1
      } else {
        stop('Error: input either lm or ridge for penalty factor calculation')
      }
      penalty <- p * penalty / sum(penalty)
      x.train <- scale(x.train, FALSE, penalty)
      x.test <- scale(x.test, FALSE, penalty)
    } else {
      penalty <- NULL
    }
    
    # Fit lasso/alasso
    fit <- glmnet(x.train, y.train, nlambda = 50, standardize = FALSE)
    y.hat <- predict.glmnet(fit, x.test, s = lambda.optim, exact = TRUE, x = x.train, y = y.train)
    mspe[i] <- fun.mse(y.hat, y.test)
    coef[i, ] <- coef(fit, s = lambda.optim)[- 1]
    y.hat.mat[i] <- y.hat
  }
  return(list(mspe = mspe, y.hat = y.hat.mat , coef = coef))
}

fun.ols.emp <- function(y, x, window = 10, horizon = 1,
                        scale = TRUE){
  ## Splitting data into rolling windows
  index.roll <- createTimeSlices(y, initialWindow = 12 * window, horizon = horizon)
  n.test <- length(index.roll$test)
  p <- ncol(x)
  mspe <- rep(NA, n.test)
  coef <- matrix(NA, nrow = n.test, ncol = p + 1)
  y.hat.mat <- rep(NA, n.test)
  
  ## Performing OLS
  for (i in 1:n.test) {
    x.train <- x[index.roll$train[[i]], ]
    x.test <- cbind(1, matrix(x[index.roll$test[[i]], ], nrow = horizon))
    y.train <- y[index.roll$train[[i]]]
    y.test <- y[index.roll$test[[i]]]
    # Scale if applicable
    if (scale) {
      train.sd <- foreach(i = 1:p, .combine = c) %do% sd(x.train[, i])
      train.mean <- colMeans(x.train)
      x.train <- scale(x.train, TRUE, TRUE)
      x.test <- scale(x.test, center = train.mean, scale = train.sd)
    }
    # Fit OLS
    fit <- lm(y.train ~ x.train)
    index.na <- which(is.na(coef(fit)))
    y.hat <- x.test[, - index.na] %*% matrix(coef(fit))[- index.na, ]
    coef[i, ] <- matrix(fit$coefficients, nrow = 1)
    mspe[i] <- fun.mse(y.test, y.hat)
    y.hat.mat[i] <- y.hat
  }
  return(list(mspe = mspe, y.hat = y.hat.mat, coef = coef))
}

fun.rw.emp <- function(y, window = 10, horizon = 1){
  ## Splitting data into rolling windows
  index.roll <- createTimeSlices(y, initialWindow = 12 * window, horizon = horizon)
  n.test <- length(index.roll$test)
  mspe <- rep(NA, n.test)
  y.hat.mat <- rep(NA, n.test)
  
  ## Performing RWwD
  for (i in 1:n.test){
    y.train <- y[index.roll$train[[i]]]
    y.test <- y[index.roll$test[[i]]]
    
    rw.mean <- mean(y.train)
    mspe[i] <- fun.mse(y.test, rw.mean)
    y.hat.mat[i] <- rw.mean
  }
  return(list(mspe = mspe, y.hat = y.hat.mat))
}


