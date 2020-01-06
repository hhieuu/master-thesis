
## Preparing variables/functions

# function for plotting given estimated coefficients data
fun.plot.coef <- function(date, alasso.coef, lasso.coef, ols,
                          horizon = 1, 
                          legend = TRUE, ylab = '') { 
  n.coef <- length(alasso.coef)
  h <- 12 * horizon - 1
  date.ref <-  date[- length(date)][- 1 : - h]
  n.ref <- length(date.ref)
  # Prepare data
  plot.data <- data.frame(date = date.ref[1:n.coef], 
                          alasso = alasso.coef, 
                          lasso = lasso.coef) %>%
    select(date, lasso, alasso) %>%
    gather(key = "variable", value = "value", - date)
  plot.data$variable <- factor(plot.data$variable, levels = unique(plot.data$variable))
  # Plot
  if (legend){
    plot <- ggplot(data = plot.data, aes(x = date, y = value)) +
      geom_line(aes(color = variable)) +
      scale_color_manual(values = c('red', 'deepskyblue2')) +
      labs(y = paste(ylab), x = '') +
      theme_bw() +
      theme(legend.position = c(0.5, 0.05),
            legend.direction = 'horizontal',
            legend.title = element_blank(),
            legend.background = element_blank(),
            axis.title.y = element_text(face = 'bold'))
  } else {
    plot <- ggplot(data = plot.data, aes(x = date, y = value)) +
      geom_line(aes(color = variable), size = 0.75) +
      scale_color_manual(values = c('red', 'deepskyblue2')) +
      labs(y = paste(ylab), x = '') +
      theme_bw() +
      theme(legend.position = 'none',
            axis.title.y = element_text(face = 'bold'))
  }
  return(plot)
}

## Plot
# Get variable names and date vector
var.names <- c('d/p', 'd/y', 'e/p', 'd/e', 'svar', 'bm', 'ntis',
               'tbl', 'lty', 'ltr', 'tms', 'dfy', 'dfr', 'infl')
plot.date <- data$yyyymm

#### FIGURE 4: Plot for each variable, window = 10, h = 1 ----
for (i in 1:length(var.names)) {
  plot.temp <- fun.plot.coef(date = plot.date,
                             alasso.coef = par.alasso.result.10[[4]]$coef[, i],
                             lasso.coef = par.lasso.result.10[[4]]$coef[, i],
                             horizon = 1,
                             legend = FALSE, ylab = var.names[i])
  assign(paste('plot.coef.', var.names[i], sep = ''), plot.temp)
}

# Add legend
plot.coef.infl <- plot.coef.infl +
  theme(legend.position = c(0.5, 0.2),
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.background = element_blank())

# Arrange plots into one canvas
selection_graph <- grid.arrange(plot.coef.dp, plot.coef.dy, plot.coef.ep, plot.coef.de, plot.coef.svar, plot.coef.bm, plot.coef.ntis, 
             plot.coef.tbl, plot.coef.lty, plot.coef.ltr, plot.coef.tms, plot.coef.dfy, plot.coef.dfr, plot.coef.infl,
             ncol = 2)

#### FIGURE 3: Plot for each variable, window = 10, h = 1 / 12 ----
for (i in 1:length(var.names)) {
  plot.temp <- fun.plot.coef(date = plot.date,
                             alasso.coef = par.alasso.result.10[[1]]$coef[, i],
                             lasso.coef = par.lasso.result.10[[1]]$coef[, i],
                             horizon = 1 / 12,
                             legend = FALSE, ylab = var.names[i])
  assign(paste('plot.coef.1_12.', var.names[i], sep = ''), plot.temp)
}

plot.coef.1_12.infl <- plot.coef.1_12.infl +
  theme(legend.position = c(0.5, 0.2),
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.background = element_blank())

grid.arrange(plot.coef.1_12.dp, plot.coef.1_12.dy, plot.coef.1_12.ep, plot.coef.1_12.de, plot.coef.1_12.svar, plot.coef.1_12.bm, plot.coef.1_12.ntis, 
             plot.coef.1_12.tbl, plot.coef.1_12.lty, plot.coef.1_12.ltr, plot.coef.1_12.tms, plot.coef.1_12.dfy, plot.coef.1_12.dfr, plot.coef.1_12.infl,
             ncol = 2)



