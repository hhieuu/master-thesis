#### FIGURE 2

## Preparing variables
alasso.yhat.10 <- par.alasso.result.10[[4]]$y.hat
lasso.yhat.10 <-par.lasso.result.10[[4]]$y.hat
ols.yhat.10 <- par.ols.result.10[[4]]$y.hat
n.data.10 <- length(alasso.yhat.10)

alasso.yhat.15 <- par.alasso.result.15[[4]]$y.hat
lasso.yhat.15 <-par.lasso.result.15[[4]]$y.hat
ols.yhat.15 <- par.ols.result.15[[4]]$y.hat
n.data.15 <- length(alasso.yhat.15)

alasso.yhat.20 <- par.alasso.result.20[[4]]$y.hat
lasso.yhat.20 <-par.lasso.result.20[[4]]$y.hat
ols.yhat.20 <- par.ols.result.20[[4]]$y.hat
n.data.20 <- length(alasso.yhat.20)

horizon <-  1
h <- 12 * horizon - 1
date.ref <-  data[- length(data[, 1]), 1][- 1 : - h,  ]

n.ref <- dim(date.ref)[1]
n.y <- length(y.1)

## Preparing data
# Here I shift date and true return y to appropriately match prediction result.
# First, the last entry of date needs to be removed so it matches y (due to lag = 1, 1st entry X and last entry y is removed)
# Then, base on window-setting (10, 15, 20-year), first windows-setting * 12 entries need to be removed
# Lastly, last 12 * horizon - 1 entries also need to be removed

plot.data.10 <- data.frame(date = date.ref[1:n.data.10, 1], 
                           y = y.1[(length(y.1) - n.data.10 + 1): length(y.1)], 
                           alasso = alasso.yhat.10, 
                           lasso = lasso.yhat.10, 
                           ols = ols.yhat.10) %>%
  select(yyyymm, y, ols, lasso, alasso) %>%
  gather(key = "variable", value = "value", - yyyymm)
plot.data.10$variable <- factor(plot.data.10$variable, levels = unique(plot.data.10$variable))


plot.data.15 <- data.frame(date = date.ref[1:n.data.15, 1], 
                           y = y.1[(length(y.1) - n.data.15 + 1): length(y.1)], 
                           alasso = alasso.yhat.15, 
                           lasso = lasso.yhat.15, 
                           ols = ols.yhat.15) %>%
  select(yyyymm, y, alasso, lasso, ols) %>%
  gather(key = "variable", value = "value", - yyyymm)
plot.data.15$variable <- factor(plot.data.15$variable, levels = unique(plot.data.15$variable))

plot.data.20 <- data.frame(date = date.ref[1:n.data.20, 1], 
                           y = y.1[(length(y.1) - n.data.20 + 1): length(y.1)], 
                           alasso = alasso.yhat.20, 
                           lasso = lasso.yhat.20, 
                           ols = ols.yhat.20) %>%
  select(yyyymm, y, alasso, lasso, ols) %>%
  gather(key = "variable", value = "value", - yyyymm)
plot.data.20$variable <- factor(plot.data.20$variable, levels = unique(plot.data.20$variable))

## Plot
plot.yhat.10 <- ggplot(data = plot.data.10, aes(x = yyyymm, y = value)) +
  geom_line(aes(color = variable),  size = 0.75) +
  scale_color_manual(values = c('lightblue', 'lightgreen', 'red', 'deepskyblue2')) +
  labs(y = '10-year, h = 1', x = '') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.y = element_text(face = 'bold'))
print(plot.yhat.10)

plot.yhat.15 <- ggplot(data = plot.data.15, aes(x = yyyymm, y = value)) +
  geom_line(aes(color = variable),  size = 0.75) +
  scale_color_manual(values = c('lightblue', 'lightgreen', 'red', 'deepskyblue2')) +
  labs(y = '15-year, h = 1', x = '') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.y = element_text(face = 'bold'))
print(plot.yhat.15)

plot.yhat.20 <- ggplot(data = plot.data.20, aes(x = yyyymm, y = value)) +
  geom_line(aes(color = variable),  size = 0.75) +
  scale_color_manual(values = c('lightblue', 'lightgreen', 'red', 'deepskyblue2'),
                     labels = c('True return', 'OLS prediction', 'Lasso prediction', 'Alasso prediction')) +
  labs(y = '20-year, h = 1', x = '') +
  theme_bw() +
  theme(legend.position = c(0.5, 0.1),
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.title.y = element_text(face = 'bold')) +
  guides(shape = guide_legend(override.aes = list(size = 0.3)))
print(plot.yhat.20)

## Arrange plots into one canvas
grid.arrange(plot.yhat.10, plot.yhat.15, plot.yhat.20, nrow = 3)

