dgp4.test.1000 <- fun.dgp4(n.burn = 1000, n.obs = 1000)
plot(ts(dgp4.test.1000$y))

dgp4.test.fit <- fun.fit.lasso.all(dgp4.test.1000, c_lambda = 0.139, pen.factor = TRUE)
print(dgp4.test.fit$mspe)
print(dgp4.test.fit$beta)

plot(ts(dgp4.test.1000$y))
plot(ts(dgp4.test.1000$x))

adf.test(dgp4.test.1000$y, 'stationary', k = 1)

################################## Testsbereich -----

y.all <- dgp2.sim.optim.200[seq(from = 1, to = length(dgp2.sim.optim.200), by = 2)]
x.all <- dgp2.sim.optim.200[seq(from = 2, to = length(dgp2.sim.optim.200), by = 2)]

dgp.all <- c(2, 3, 4)

fun.rep <- function(dgp, n.burn, n.obs, lambda, penalty = TRUE){
  sample <- eval(parse(text = paste('fun.dgp', dgp, '(n.burn = ', n.burn, ', n.obs = ', n.obs, ')', sep = '')))
  fit <- fun.fit.lasso.all(sample, c_lambda = lambda, pen.factor = penalty)
  return(list(mspe = fit$mspe, beta = fit$beta))
}

fun.rep(4, 100, 100, dgp4.alasso.c_lambda, TRUE)

tic()
test.rep <- replicate(100, fun.rep(4, n.burn = 500, n.obs = 10000, lambda = 0.36, penalty = TRUE))
test.rep.mspe <- mean(do.call(c, test.rep[(seq(from = 1, to = length(test.rep), by = 2))]))
print(test.rep.mspe)
# test.rep.beta <- do.call(rbind, test.rep[(seq(from = 0, to = length(test.rep), by = 2))])
toc()

tic()
dgp4.alasso.mspe.new <- rep(0, length(n.settings))
dgp4.alasso.beta.new <- list()
dgp4.alasso.sr.new <- matrix(0, nrow = length(n.settings), ncol = 3)
for (i in 1:length(n.settings)) {
  dgp4.alasso.rep <- replicate(1000, fun.rep(4, n.burn = 500, n.obs = n.settings[i], 
                                      lambda = 0.1, penalty = TRUE))
  dgp4.alasso.mspe.new[i] <- mean(do.call(c, 
                                          dgp4.alasso.rep[(seq(from = 1, to = length(dgp4.alasso.rep), by = 2))]))
  dgp4.alasso.beta.new[[i]] <- do.call(rbind, 
                                  dgp4.alasso.rep[(seq(from = 0, to = length(dgp4.alasso.rep), by = 2))])
  dgp4.alasso.sr.new[i, ] <- fun.perf.lasso(beta.true = c(0.4, 0, 1 / sqrt(n.settings[i]), 0, 0.3, - 0.3, 0, 0),
                                                 beta.est = dgp4.alasso.beta.new[[i]], n.sim = 1000)
}
toc()


#####

plot(ts(dgp4.test.1000$y))
plot(ts(dgp4.test.1000$x))

plot(ts(x.all.1000[[1]][, 7]))
lines(ts(x.all.1000[[1]][, 8]), col = 'red')
plot(ts(x.all.1000[[1]][, 8]))

plot(function(n) (n - 3) / n, from = 1, to = 5e2)
plot(function(n) 1 - 1 / (n ^ 3 + 1), from = 1, to = 5e2, col = 'red', add = T)
plot(function(n) n / (n + 1), from = 1, to = 5e2, col = 'blue', add = T)

plot(function(n) 1 / (n ^ (1 / 6) * log(log(n))) + log(log(n)) / n ^ (1 / 3), from = 1, to = 1e5, col = 'blue', add = T)
plot(function(n) n^(1/2), from = 1, to = 5e3, col = 'green', add = T)
plot(function(n) n / log(n), from = 1, to = 5e3, col = 'purple', add = T)


par(mfrow = c(1, 2))
plot(function(n) n ^ (1 / 4) / log(log(n)), from = 10, to = 1e6, ylim = c(0, 80))
plot(function(n) n ^ (1 / 3) / log(log(n)), from = 10, to = 1e6, col = 'red', add = T)
plot(function(n) n ^ (1 / 2) / log(log(n)), from = 10, to = 1e6, col = 'blue', add = T)

plot(function(n) 1 / (n ^ (1 / 4) * log(log(n))) + log(log(n)) / n ^ (1 / 4), from = 10, to = 1e6, ylim = c(0, 1.75))
plot(function(n) 1 / (n ^ (1 / 6) * log(log(n))) + log(log(n)) / n ^ (1 / 3), from = 10, to = 1e6, col = 'red', add = T)
plot(function(n) 1 / log(log(n)) + log(log(n)) / n ^ (1 / 2), from = 10, to = 1e6, col = 'blue', add = T)
par(mfrow = c(1, 1))

test.lasso.fit <- fun.fit.lasso(dgp4.test.1000$y, dgp4.test.1000$x)

dim(1 / abs(lm(dgp4.test.1000$y ~ dgp4.test.1000$x)$coefficients[-1]))

as.vector(1 / abs(test.lasso.fit$beta))

fun.cv.lasso(dgp4.test.1000$y, dgp4.test.1000$x, lambda = 0.16, kf = 10, pen.factor = TRUE)
fun.cv.lasso(dgp4.test.1000$y, dgp4.test.1000$x, lambda = 0.20, kf = 10, pen.factor = TRUE)
fun.cv.lasso(dgp4.test.1000$y, dgp4.test.1000$x, lambda = 0.9, kf = 10, pen.factor = TRUE)


test.fit <- glmnet(dgp4.test.1000$x, dgp4.test.1000$y)
x.val <- matrix(dgp4.test.1000$x[1:10, ], ncol = 8)
test.pred <- coef.glmnet(test.fit, s = 0.3)
test.pred[-1]

####

dgp2.test.1000 <- fun.dgp2(n.burn = 500, n.obs = 1000)

plot(ts(dgp2.test.1000$y))
plot(ts(dgp2.test.1000$x))
adf.test(dgp2.test.1000$x[, 1])
adf.test(dgp2.test.1000$x[, 2])
adf.test(dgp2.test.1000$x[, 3])
adf.test(dgp2.test.1000$x[, 4])
adf.test(dgp2.test.1000$x[, 5])
adf.test(dgp2.test.1000$x[, 6])
adf.test(dgp2.test.1000$x[, 7])
adf.test(dgp2.test.1000$x[, 8])
adf.test(dgp2.test.1000$x[, 9])


tic()
emp.test <- fun.lasso.predict(y.1_12, x.1_12)
toc()

tic()
emp.test1 <- fun.cv.lasso(y.1_12[2:121], x.1_12[2:121, ], lambda = 0.12, kf = 10, pen.factor = TRUE, pen.type = 'ridge')
toc()


test.obj <- rep(NA, 500)
tic()
for (i in 1:(500 - 120)){
  test.obj[i] <- mean(sapply(seq(from = 0.001, to = 0.3, length.out = 5), 
                             function(s) fun.cv.lasso(y.1_12[i:(i + 120 - 1)], x.1_12[i:(i + 120 - 1), ], kf = 5, 
                                                      lambda = s, pen.factor = TRUE, pen.type = 'ridge')))
}
toc()


a <- c(1:6, NA, 7, NA, 8:10)
a[is.na(a)] <- 0
a



#### DEBUGGING NOTE:
# After reading the error trace back, I know that something is wrong with the CV function.
# Specifically, something is wrong with the scaling, since lasso is unaffected but alasso is affected.

### First I set up a test function to record the scaled matrix of predictors x

fun.cv.test <- function(y, x, lambda, kf, pen.factor = TRUE, pen.type){ # calculate k-fold CV mspe
  n.obs <- length(y)
  p <- ncol(x)
  cv.index <- fun.ts.split(y, k = kf)
  n.trainset <- length(cv.index$train)
  n.valset <- length(cv.index$val)
  mspe.mat <- rep(1, n.valset)
  x.train.all <- list()
  x.val.all <- list()
  
  
  for (i in 1:n.trainset) {
    y.train <- y[cv.index$train[[i]]]
    y.val <- y[cv.index$val[[i]]]
    x.train <- x[cv.index$train[[i]], ]
    x.val <- x[cv.index$val[[i]], ]
    # Check penalty factor
    if (pen.factor){
      if ('lm' %in% pen.type) {
        penalty <- (1 / abs(lm(y.train ~ x.train)$coefficients[-1])) ^ 1
      } else if ('ridge' %in% pen.type) {
        ridge <- glmnet(x.train, y.train, alpha = 0, nlambda = 50)
        ridge.coef <- coef(ridge, s = 1e-10, exact = TRUE, x = x.train, y = y.train)[- 1]
        penalty <- (1 / abs(ridge.coef)) ^ 1
      } else {
        stop('Error: imput either lm or ridge for penalty factor calculation')
      }
      
      penalty <- p * penalty / sum(penalty)
      x.train <- scale(x.train, FALSE, penalty)
      # capture all NA, NaN, Inf and setting them to 0
      x.train[is.na(x.train)] <- 0
      x.train[is.nan(x.train)] <- 0
      x.train[is.infinite(x.train)] <- 0
      x.val <- scale(x.val, FALSE, penalty)
      x.val[is.na(x.val)] <- 0
      x.val[is.nan(x.val)] <- 0
      x.val[is.infinite(x.val)] <- 0
    } else {
      penalty <- NULL
    }
    x.train.all[[i]] <- x.train
    x.val.all[[i]] <- x.val
  }
  return(list(train = x.train.all, val = x.val.all))
}

### Second, I check for zeroes in the matrix, as I have modified the CV lasso function to give zero at NA, NaN, and Inf values


test.mat1 <- matrix(NA, ncol = 10, nrow = 380)
test.mat2 <- matrix(NA, ncol = 10, nrow = 380)
test.list <- list()
tic()
for (i in 1:(500 - 120)){
  temp.list <- lapply(seq(from = 0.001, to = 0.3, length.out = 5), 
                             function(s) fun.cv.test(y.1_12[i:(i + 120 - 1)], x.1_12[i:(i + 120 - 1), ], kf = 10, 
                                                      lambda = s, pen.factor = TRUE, pen.type = 'ridge'))
  test.list[[i]] <- temp.list
  for (j in 1:5){
    test.mat1[i, ] <- sum(foreach(k = 1:10, .combine = 'c') %do% sum(temp.list[[j]][["train"]][[k]] == 0))
    test.mat2[i, ] <- sum(foreach(k = 1:10, .combine = 'c') %do% sum(temp.list[[j]][["val"]][[k]] == 0))
  }
}
toc()

sum(test.mat1)
sum(test.mat2)
# Apparently, something is wrong with these, as there are almost all zeros in the range of data from row 186 to 235.
# We can see it here
test.mat1[185:236]
test.mat2[185:236]

### Next, I check the weighted data (ridge weight) at the same range above to see if I am right
View(test.list[186])
# Here, some folds are all zeros, which indicate either NA, NaN, or Inf
# Take a look at the scale attribute reveals that all scales are 0, except for one NaN, which correspond to predictor tbl,
#   which is the T-bill rate.
var.all.lag[186:245, 8]

# The T-bill rate is surprisingly constant in this time. Since we use 10-fold fixed rolling window CV with train set
#   of size 120, each train set only starts from 11 obs at fold 1, to 110 obs at fold 10. The length of 'constant tbl'
#   period is 59 (months), so it is understandable that some folds consist of most to all of this period.
ridge.test <- glmnet(x[186:245, ], y[186:245], alpha = 0)
ridge.coef.test <- 1 / abs(coef.glmnet(ridge.test, s = 1e-8)[- 1])

pen.test <- ncol(x) * ridge.coef.test / sum(ridge.coef.test) 
lm.test <- lm(y[186:245] ~ x[186:245, ])
coef(lm.test)



### We then need to solve this problem
## The most apparent choice here is to omit the data altogether, starting the data from point t = 246.
##  However, this is not ideal because we effectively omit 1/5 of our data. That is about more than 20 years of data
##  for only 5 years of constant values.
## The solution is that we will try to follow here will be about replacing the troublesome parts with its
##  historical moving average of 20 years

var.all.lag[180:245, 8]
check <- rep(NA, 60)
for (i in 60:1){
  check[i] <- mean(var.all.lag[(186 - i - 1):(245 - i - 1), 8])
}
rev(check)

## Another solution is that we will add a random term \epsilon ~ N(mu, sigma^2) to the data
##  where mu and sigma^2 are the mean and variance of tbl for the 10-year period of tbl, from 186 - 30 to 245 + 30

eps <- rnorm(60, mean = mean(var.all.lag[(186 - 30):(245 + 30), 8]), sd = sqrt(var(var.all.lag[(186 - 30):(245 + 30), 8])))
var.all.lag[186:245, 8] + eps

##### Another test for parallelization
cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
parallel::setDefaultCluster(cl = cl)
# on.exit(parallel::stopCluster(cl), add = TRUE)
envir <- environment(fun.lasso.predict)
parallel::clusterExport(cl, varlist = ls(envir), envir = envir)


tic()
par.alasso.test <- parLapply(cl, horizon.setting,
                                  function(a) fun.lasso.predict.par(cl, get(paste('y.', a, sep = '')), get(paste('x.', a, sep = '')),
                                                                init.val = 0.03))
toc()
stopCluster(cl)

tic()
test1 <- fun.lasso.predict.par(cl, y.1_12, x.1_12, init.val = 0.03)
toc()

tic()
fun.lasso.predict(y.1_12, x.1_12, init.val = 0.03)
toc()

