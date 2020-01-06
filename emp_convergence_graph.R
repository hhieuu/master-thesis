#### FIGURE 1

alasso.lambda.seq.1 <- function(n) n ^ (1 / 2) / log(log(n))
alasso.lambda.seq.2 <- function(n) n ^ (1 / 2) / log(n)
alasso.lambda.seq.3 <- function(n) n ^ (1 / 2) / log(n) ^ 2
alasso.lambda.seq.4 <- function(n) n ^ (1 / 2) / log(n) ^ 3

alasso.lambda.tot.1 <- function(n) 1 / log(log(n)) + log(log(n)) / (n ^ (1 / 2))
alasso.lambda.tot.2 <- function(n) 1 / log(n) + log(n) / (n ^ (1 / 2))
alasso.lambda.tot.3 <- function(n) 1 / log(n) ^ 2 + log(n) ^ 2 / (n ^ (1 / 2))
alasso.lambda.tot.4 <- function(n) 1 / log(n) ^ 3 + log(n) ^ 3 / (n ^ (1 / 2))



plot.alasso.lambda1 <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
  stat_function(fun = alasso.lambda.seq.1, aes(colour = 'lambda1')) +
  stat_function(fun = alasso.lambda.seq.2, aes(colour = 'lambda2')) +
  stat_function(fun = alasso.lambda.seq.3, aes(colour = 'lambda3')) +
  stat_function(fun = alasso.lambda.seq.4, aes(colour = 'lambda4')) +
  scale_colour_manual(values = c('lambda1' = 'black', 'lambda2' = 'red', 'lambda3' = 'blue', 'lambda4' = 'green'),
                      name = '',
                      labels = expression(lambda[n]^(al1), lambda[n]^(al2), lambda[n]^(al3), lambda[n]^(al4))) +

  xlim(40, 1e4) +
  labs(title = '(a)', x = 'n', y = expression(lambda[n]^(al))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom')

plot.alasso.lambda2 <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
  stat_function(fun = alasso.lambda.tot.1, aes(colour = 'lambda1')) +
  stat_function(fun = alasso.lambda.tot.2, aes(colour = 'lambda2')) +
  stat_function(fun = alasso.lambda.tot.3, aes(colour = 'lambda3')) +
  stat_function(fun = alasso.lambda.tot.4, aes(colour = 'lambda4')) +
  scale_colour_manual(values = c('lambda1' = 'black', 'lambda2' = 'red', 'lambda3' = 'blue', 'lambda4' = 'green'),
                      name = '',
                      labels = expression(g[n]^(al1), g[n]^(al2), g[n]^(al3), g[n]^(al4))) +
  xlim(40, 1e4) +
  labs(title = '(b)', x = 'n', y = expression(g[n]^(al))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom')

grid.arrange(plot.alasso.lambda1, plot.alasso.lambda2, ncol = 2)

# tikz('latex\\convergence_rate_1.tex', width = 6, height = 6)
# grid.arrange(plot.alasso.lambda1, plot.alasso.lambda2, ncol = 2)
# dev.off()
