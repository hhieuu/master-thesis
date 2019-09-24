source('myfun.R')

set.seed(2607)
dgp1.data <- fun.dgp2(n.burn = 1000, n.obs = 200)
var.y <- dgp1.data$y
var.x <- dgp1.data$x

#### Setting for simulation: ----
## Number of different settings: n = 40, 80, 120, 200, 500, 1000
## Number of model for evaluation: OLS, Oracle OLS, AR, and Adaptive Lasso
## Evaluation is carried out by computing Mean Squared Prediction Error
## For each setting, the study is replicated 100 times. Final MSPE is the average of all MSPE in all replications
## 






## OLS
# fun.cv.ols(var.y, var.x, kf = 10)
# split.index <- fun.ts.split(var.y, k = 10)
# fun.pred.ols(var.y[split.index$train[[2]]], var.x[split.index$train[[2]], ],
#              var.y[split.index$val[[2]]], var.x[split.index$val[[2]], ])
# fun.cv.lasso(var.y, var.x, lambda = 0.5, kf = 10)

# Now we will 

a <- matrix(c(1, 0, 2, 3, 5, 
              0, 1, 6, 9, 3),
            nrow = 2, byrow = T)

b <- matrix(c(-2, -6, 1, 0, 0,
              -3, -9, 0, 1, 0,
              -5, -3, 0, 0, 1),
            nrow = 3, byrow = T)

a%*%t(b)
b%*%t(a)
