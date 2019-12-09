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

plot(function(n) n ^ (1 / 2) / log(log(n)), from = 1, to = 1e5)
plot(function(n) 1 / log(log(n)) + log(log(n)) / n ^ (1 / 2), from = 1, to = 1e5, col = 'red', add = T)
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


require(stats)
centre <- function(x, type) {
  switch(type,
         mean = sum(x),
         median = median(x),
         trimmed = mean(x, trim = .1))
}
x <- rcauchy(10)
a <- centre(x, "sum")
centre(x, "median")
centre(x, "trimmed")


fun.test <- function(x, type = c('sum', 'mean', 'median')){
  if (type == 'sum') {
    
  }
}
