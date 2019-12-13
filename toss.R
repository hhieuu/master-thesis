fun.cv.lasso.all <- function(data, c_lambda, kf, pen.factor = TRUE){ 
  
  n.rep <- length(data) / 2
  setting <- length(data[[1]])
  y.list <- data[seq(from = 1, to = length(data), by = 2)]
  x.list <- data[seq(from = 0, to = length(data), by = 2)]
  index <- fun.train.test.split(1:setting, test.portion = .2)
  n.col <- ncol(x.list[[1]])
  
  # Lambda
  if (pen.factor) {
    lambda <- c_lambda * sqrt(setting / log(log(setting)))
  } else {
    lambda <- c_lambda * sqrt(setting)
  }
  
  mspe <- rep(0, n.rep)
  success.rate <- matrix(0, nrow = n.rep, ncol = 3)
  for (i in 1:n.rep){
    mspe.cv <- rep(0, kf)
    beta.cv <- Matrix(0, nrow = kf, ncol = n.col, sparse = T)
    for (j in 1:kf){
      cv.index <- fun.ts.split(1:setting, k = kf)
      x.train <- x.list[[i]][cv.index$train[[j]], ]
      x.val <- matrix(x.list[[i]][cv.index$val[[j]][1], ], nrow = 1)
      y.train <- y.list[[i]][cv.index$train[[j]]]
      y.val <- y.list[[i]][cv.index$val[[j]][1]]
      # Penalty factor (weights)
      if (pen.factor){
        # pen.factor.cv <- 1/abs(lm(y.train ~ x.train)$coefficients[-1])
        pen.factor.cv <- 1/abs(glmnet(x.train, y.train, alpha = 0, lambda = 1e8)$beta@x)
      } else {
        pen.factor.cv <- NULL
      }
      fit <- fun.fit.lasso(y.train, x.train, lambda, pen.factor = pen.factor.cv)
      mspe.cv[j] <- fun.pred.lasso(y.val, x.val, fit)
      beta.cv[j, ] <- fit$beta
    }
    mspe[i] <- mean(mspe.cv)
    success.rate[i, ] <- fun.perf.lasso(beta.true = c(0, -0.2, 1/sqrt(n.settings[i]), 0.6, 0, 0.3, -0.7, 0, 0.5),
                                        beta.est = beta.cv,
                                        n.obs = 1, n.sim = kf)
  }
  return(list(mspe = mean(mspe), success.rate = colMeans(success.rate)))
}

## Test using cross-validation: NOT WORTH IT, ABORT. Well, I spent like two hours on this already.
# tic()
# alasso.mspe.cv <- rep(0, length(n.settings))
# alasso.success.rate.cv <- matrix(0, nrow = length(n.settings), ncol = 3)
# for (i in 1:length(n.settings)){
#   alasso.current.setting.cv <- paste('sim.', n.settings[i], sep = '')
#   alasso.current.fit.cv <- fun.cv.lasso.all(get(alasso.current.setting.cv), c_lambda = alasso.c_lambda, 
#                                             kf = 10, pen.factor = TRUE)
#   alasso.mspe.cv[i] <- alasso.current.fit.cv$mspe
#   alasso.success.rate.cv[i, ] <-  alasso.current.fit.cv$success.rate
# }
# toc()


fun.dgp2.alt <- function(n.burn = 500, n.obs = 500){
  n.tot <- n.burn + n.obs
  
  ## dgp for 3 independent AR processes (no collinearity)
  n.ar <- 3
  rho.ar <- matrix(c(0.99, 0.90, -0.5)) # vector of ar coefficent
  x.ar <- matrix(NA, ncol = n.ar, nrow = n.tot)
  x.ar[1, ] <- rnorm(n.ar) # initial value
  x.ar.err <- matrix(rnorm(n.ar*n.tot), ncol = n.ar, nrow = n.tot) # error matrix
  for (i in 2:n.tot){
    x.ar[i, ] <- rho.ar*x.ar[i - 1, ] + x.ar.err[i, ] # simulating AR(1) variables
  }
  
  ## dgp for 2 highly-collinear AR(1) process
  # n.coll <- 2
  # rho.coll <- matrix(c(0.8, 0.75))
  # x.coll <- matrix(NA, ncol = n.coll, nrow = n.tot)
  # x.coll[1, ] <- rnorm(n.coll) # init values
  # x.cov.err <- (n.obs - 2) / n.obs * (1 - rho.coll[1] * rho.coll[2]) / sqrt((1 - rho.coll[1] ^ 2) * (1 - rho.coll[2] ^ 2))
  # x.coll.err <- mvrnorm(n.coll*n.tot, mu = c(0, 0), Sigma = matrix(c(2, x.cov.err, x.cov.err, 1), 2, 2)) # error matrix. Error term is highly correlated (0.95)
  # for (i in 2:n.tot){
  #   x.coll[i, ] <- rho.coll*x.coll[i - 1, ] + x.coll.err[i, ] # simulating AR(1) variables
  # }
  
  n.coll <- 2
  rho.coll.ar <- matrix(c(0.2, 0.2))
  x.coll.ar <- matrix(NA, ncol = n.coll, nrow = n.tot)
  x.coll.ar[1, ] <- rnorm(n.coll) # init values
  x.cov.ar <- (n.obs - 2) / n.obs * (1 - rho.coll.ar[1] * rho.coll.ar[2]) / sqrt((1 - rho.coll.ar[1] ^ 2) * (1 - rho.coll.ar[2] ^ 2))
  x.coll.ar.err <- mvrnorm(n.coll*n.tot, mu = c(0, 0), Sigma = matrix(c(1, x.cov.ar, x.cov.ar, 1), 2, 2)) # error matrix. Error term is highly correlated (0.95)
  for (i in 2:n.tot){
    x.coll.ar[i, ] <- rho.coll.ar * x.coll.ar[i - 1, ] + x.coll.ar.err[i, ] # simulating AR(1) variables
  }
  x.coll <- matrix(NA, ncol = n.coll, nrow = n.tot)
  x.coll[1, ] <- rnorm(n.coll) # init values
  x.coll.const <- matrix(runif(2), ncol = n.coll)
  for (i in 2:n.tot){
    x.coll[i, ] <- i * x.coll.const + x.coll[i - 1, ] + x.coll.ar[i, ]
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
  # y.coef <- matrix(c(0, -0.2, 1/sqrt(n.obs), 0.6, 0, 0.3, -0.7, 0, 0.5), nrow = 1)
  # y.coef <- matrix(c(0, -1.5, 1/sqrt(n.obs), 0.8, 0, 0, -0.9, 1.2, 0), nrow = 1)
  y.coef <- matrix(c(0, - 0.6, 1/sqrt(n.obs), 0, 0, 0.3, - 0.2, 0.3, 0), nrow = 1)
  y.err <- rnorm(n.tot)
  y <- y.intercept + x %*% t(y.coef) + y.err
  
  return(list(y = y[(n.burn + 1):n.tot], x = x[(n.burn + 1):n.tot, ]))
}


# Parallel optim test

# tic()
# alasso.lambda.optim.test.prl <- fun.optim.lambda.parallel(dgp2.sim.optim.200, start.val = 0.1, pen.factor = TRUE)
# alasso.c_lambda.prl <- median(alasso.lambda.optim.test.prl[, 1])
# toc()
# 
# tic() # ~ 49 mins to run, careful
# lambda.optim.test.alt.parallel <- fun.optim.lambda.alt.parallel(dgp2.sim.optim.200, start.val = 0.03)
# c_lambda.alt.parallel <- median(lambda.optim.test.alt.parallel[, 1])
# toc()


## Parallel Test

fun.optim.lambda.parallel <- function(sample, start.val, pen.factor = TRUE){
  n.sim <- length(sample)/2
  y.all <- sample[seq(from = 1, to = length(sample), by = 2)]
  x.all <- sample[seq(from = 2, to = length(sample), by = 2)]
  optim.lambda <- matrix(0, ncol = 2, nrow = n.sim)
  
  cl <- parallel::makeCluster(4)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  envir <- environment(fun.cv.lasso)
  parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
  parallel::setDefaultCluster(cl = cl)
  
  for (i in 1:n.sim){
    optim.tmp <- optimParallel(par = start.val, 
                               function(a) fun.cv.lasso(y.all[[i]], x.all[[i]], kf = 10, lambda = a, pen.factor = pen.factor),
                               method = "L-BFGS-B", lower = 0.001)
    optim.lambda[i, 1] <- optim.tmp$par
    optim.lambda[i, 2] <- optim.tmp$value
  }
  return(optim.lambda)
  parallel::stopCluster(cl)
}

fun.optim.lambda.alt.parallel <- function(sample, start.val, pen.factor = TRUE){
  n.sim <- length(sample)/2
  y.all <- sample[seq(from = 1, to = length(sample), by = 2)]
  x.all <- sample[seq(from = 2, to = length(sample), by = 2)]
  optim.lambda <- matrix(0, ncol = 2, nrow = n.sim)
  
  cl <- parallel::makeCluster(4)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  envir <- environment(fun.cv.lasso.alt)
  parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
  parallel::setDefaultCluster(cl = cl)
  
  for (i in 1:n.sim){
    optim.tmp <- optim(par = start.val, 
                       function(a) fun.cv.lasso.alt(y.all[[i]], x.all[[i]], kf = 10, lambda = a, pen.factor = pen.factor),
                       method = "L-BFGS-B", lower = 0.001)
    optim.lambda[i, 1] <- optim.tmp$par
    optim.lambda[i, 2] <- optim.tmp$value
  }
  return(optim.lambda)
  parallel::stopCluster(cl)
}


fun.cv.lasso.alt <- function(y, x, lambda, kf, pen.factor = TRUE){ 
  n.obs <- length(y)
  cv.index <- createTimeSlices(y, initialWindow = 0.2 * length(y), horizon = 1)
  n.trainset <- length(cv.index$train)
  n.valset <- length(cv.index$test)
  mspe.mat <- rep(1, n.valset)
  
  for (i in 1:n.trainset) {
    y.train <- y[cv.index$train[[i]]]
    y.val <- y[cv.index$test[[i]]]
    x.train <- x[cv.index$train[[i]], ]
    x.val <- matrix(x[cv.index$test[[i]], ], nrow = 1)
    if (pen.factor){
      pen.temp <- 1/abs(lm(y ~ x)$coefficients[-1]) ^ 1
    } else {
      pen.temp <- NULL
    }
    fitted <- fun.fit.lasso(y.train, x.train, lambda, pen.factor = pen.temp)
    mspe.mat[i] <- fun.pred.lasso(y.val, x.val, fitted)
  }
  cv.mspe <- mean(mspe.mat)
  return(cv.mspe)
}

fun.optim.lambda.alt <- function(sample, start.val, pen.factor = TRUE){
  n.sim <- length(sample)/2
  y.all <- sample[seq(from = 1, to = length(sample), by = 2)]
  x.all <- sample[seq(from = 2, to = length(sample), by = 2)]
  optim.lambda <- matrix(0, ncol = 2, nrow = n.sim)
  for (i in 1:n.sim){
    optim.tmp <- optim(par = start.val, 
                       function(a) fun.cv.lasso.alt(y.all[[i]], x.all[[i]], kf = 10, lambda = a, pen.factor = pen.factor),
                       method = "L-BFGS-B", lower = 0.001)
    optim.lambda[i, 1] <- optim.tmp$par
    optim.lambda[i, 2] <- optim.tmp$value
  }
  return(optim.lambda)
}

# OLD FIT FUNCTION
fun.fit.lasso.all <- function(data, c_lambda, pen.factor = TRUE) {
  n.rep <- length(data) / 2
  setting <- length(data[[1]])
  test.portion <- .1
  y.list <- data[seq(from = 1, to = length(data), by = 2)]
  x.list <- data[seq(from = 0, to = length(data), by = 2)]
  index <- fun.train.test.split(1:setting, test.portion = test.portion)
  fit.list <- list()
  beta <- Matrix(0, nrow = n.rep, ncol = ncol(x.list[[1]]), sparse = T)
  mspe <- rep(0, n.rep)
  
  # Lambda
  if (pen.factor) {
    lambda <- c_lambda * (setting * (1 - test.portion)) ^ (1 / 2) / log(log(setting * (1 - test.portion)))
  } else {
    # lambda <- c_lambda * log(setting) / log(log(setting))
    lambda <- c_lambda * (setting * (1 - test.portion)) ^ (1 / 2)
  }
  
  for (i in 1:n.rep) {
    x.train <- x.list[[i]][index$train, ]
    y.train <- y.list[[i]][index$train]
    
    x.test <- as.matrix(x.list[[i]][index$test, ])
    y.test <- y.list[[i]][index$test]
    
    # Check penalty factor
    if (pen.factor) {
      penalty <- (1 / abs(lm(y.train ~ scale(x.train, TRUE, TRUE))$coefficients[- 1])) ^ 1
    } else {
      penalty <- NULL
    }
    
    fit <- fun.fit.lasso(y.train, x.train, lambda = lambda, pen.factor = penalty)
    beta[i, ] <- fit$beta
    mspe[i] <- fun.pred.lasso(y.test, x.test, fit)
    fit.list[[i]] <- fit
  }
  return(list(mspe = mspe, beta = beta))
}

fun.fit.lasso <- function(y.train, x.train, lambda, pen.factor = NULL){
  x.train <- Matrix(x.train)
  x.train.s <- scale(x.train, TRUE, TRUE)
  if (is.null(pen.factor)){
    pen.factor <- rep(1, ncol(x.train))
  }
  fit <- glmnet(x.train.s, y.train, lambda = lambda, penalty.factor = pen.factor, standardize = FALSE)
  return(fit)
}

fun.pred.lasso <- function(y.val, x.val, fitted){ # calculate mspe for given fit, validation set, lambda
  x.val.s <- scale(x.val, TRUE, TRUE)
  y.val.hat <- predict.glmnet(fitted, x.val.s)
  mspe <- fun.mse(y.val, y.val.hat)
  return(mspe)
}

fun.lasso.predict.par <- function(cl, y, x, pen.factor = TRUE, pen.type = 'ridge', 
                                  window = 10, horizon = 1){
  ## Splitting data into rolling windows
  index.roll <- createTimeSlices(y, initialWindow = 12 * window, horizon = horizon)
  n.test <- length(index.roll$test)
  p <- ncol(x)
  mspe <- rep(NA, n.test)
  lambda <- rep(NA, n.test)
  coef <- matrix(NA, nrow = n.test, ncol = p)
  
  result <- parLapply(cl, c(1:n.test), function(i) {
    x.train <- x[index.roll$train[[i]], ]
    y.train <- y[index.roll$train[[i]]]
    x.test <- matrix(x[index.roll$test[[i]], ], nrow = horizon)
    y.test <- y[index.roll$test[[i]]]
    # First find lambda by optimizing w.r.t  cv
    lambda.seq <- sapply(seq(from = 0.001, to = 0.3, length.out = 30),
                         function(s) fun.cv.lasso(y.train, x.train, kf = 5, lambda = s,
                                                  pen.factor = pen.factor, pen.type = pen.type))
    lambda.start <- seq(from = 0.0001, to = 0.3, length.out = 30)[which.min(lambda.seq)]
    lambda.optim <- optim(par = lambda.start, function(s) fun.cv.lasso(y.train, x.train, lambda = s, kf = 5,
                                                                       pen.factor = pen.factor, pen.type = pen.type),
                          method = "L-BFGS-B", lower = 1e-4)$par
    
    
    
    # Then find the weights
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
    coef[i, ] <- coef(fit, s = lambda.optim, exact = TRUE, x = x.train, y = y.train)[- 1]
    lambda[i] <- lambda.optim
    return(list(mspe = mspe, coef = coef, lambda = lambda))
  })
  mspe <- foreach(i = 1:n.test, .combine = c) %do% result[[i]]$mspe
  coef <- foreach(i = 1:n.test, .combine = rbind) %do% result[[i]]$coef
  lambda <- foreach(i = 1:n.test, .combine = c) %do% result[[i]]$lambda
  return(list(mspe = mspe, coef = coef, lambda = lambda))
}