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
        ridge <- glmnet(x.train, y.train, alpha = 0, nlambda = 50, standardize = FALSE)
        ridge.coef <- coef(ridge, s = 1e-10, exact = TRUE, x = x.train, y = y.train)[- 1]
        penalty <- (1 / abs(ridge.coef)) ^ 1
      } else {
        stop('Error: input either lm or ridge for penalty factor calculation')
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

### Second, I check for zeroes in the matrix, as I have modified the CV lasso function to give zeroes at NA, NaN, and Inf values


test.mat1 <- matrix(NA, ncol = 10, nrow = 380)
test.mat2 <- matrix(NA, ncol = 10, nrow = 380)
test.list <- list()
tic()
# Check the part of data from index 1 to 380.
# Check 10 years of data (120 months/data points) forward starting from each index i from 1 to 380
# In each replication

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

tic()
for (i in 1:(500 - 120)){
  temp.list <- fun.cv.test(y.1_12[i:(i + 120 - 1)], x.1_12[i:(i + 120 - 1), ], kf = 10, 
                                              lambda = 0.01, pen.factor = TRUE, pen.type = 'ridge')
  test.list[[i]] <- temp.list
  test.mat1[i, ] <- sum(foreach(k = 1:10, .combine = 'c') %do% sum(temp.list[["train"]][[k]] == 0))
  test.mat2[i, ] <- sum(foreach(k = 1:10, .combine = 'c') %do% sum(temp.list[["val"]][[k]] == 0))
}
toc()

sum(test.mat1)
sum(test.mat2)
# Apparently, something is wrong with these, as there are almost all zeros in the range of data from row 186 to 235:
test.mat1[185:236]
test.mat2[185:236]

### Next, I check the weighted data (ridge weight) at the same range above to see if I am right
# Example: 
View(test.list[[186]])
# Here, some folds are all zeros, which indicate either all data are NA, NaN, or Inf after scaling. Something is wrong with the scale calculation!
# Take a look at the scale attribute reveals that all scales are 0, except for one NaN, which correspond to predictor tbl,
#   which is the T-bill rate.
var.all.lag[186:245, 8]

# The T-bill rate is unexpectedly constant in this time. Since we use 10-fold fixed rolling window CV with train set
#   of size 120, each train set only starts with 11 obs at fold 1, and 110 obs at fold 10. The length of 'constant tbl' period
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
## The solution that I try to follow here will be about replacing the troublesome parts with its
##  historical moving average of 20 years

var.all.lag[180:245, 8]
check <- rep(NA, 60)
for (i in 60:1){
  check[i] <- mean(var.all.lag[(186 - i - 1):(245 - i - 1), 8])
}
rev(check) # reverse elements


