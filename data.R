### Import and clean Data
data.original <- read_excel("PredictorData2018.xlsx")

# Convert dates
data.original[[1]] <- as.Date(as.character(paste(data.original[[1]], '01', sep = '')), format = "%Y%m%d")
data <- subset(data.original, data.original$yyyymm >= as.Date('1927-01-01'))
date.start.year <- as.numeric(format(data[[1]][1], '%Y'))
date.start.month <- as.numeric(format(data[[1]][1], '%m'))
date.end.year <- as.numeric(format(data[[1]][length(data[[1]])], '%Y'))
date.end.month <- as.numeric(format(data[[1]][length(data[[1]])], '%m'))
# Convert data to time series
for (i in 2:length(data)) {
  data[[i]] <- ts(as.numeric(data[[i]]), start = c(date.start.year, date.start.month), frequency = 12)
}

### Preparing variables
# Dependent Var
var.premium <- data$CRSP_SPvw - data$Rfree
var.premium.stdrd <- scale(var.premium)
n.obs <- length(var.premium)

# Regressors
var.dp <- log(data$D12) - log(data$Index) # d/p ratio
var.dy <- log(data$D12) - log(c(data$Index[-1], data$Index[n.obs])) # d/y ratio
var.ep <- log(data$E12) - log(data$Index) # e/p ratio
var.de <- log(data$D12) - log(data$E12) # d/e ratio
var.svar <- data$svar # stock variance
var.csp <- data$csp # cross-sectional premium
var.bm <- data$`b/m` # book to market ratio
var.ntis <- data$ntis # net equity expansion
# var.eqis # percent equity issuing
var.tbl <- data$tbl # treasury bill/T-bill rate from 1920 to 1933
var.lty <- data$lty # long term government bond yield, 1919 to 1925
var.ltr <- data$ltr # long term government bond return, 1926 to 2005
var.tms <- var.lty - var.tbl # term spread
var.dfy <- data$BAA - data$AAA # default yield spread
var.dfr <- data$corpr - data$ltr # default return spread
var.infl <- data$infl # inflation from Consumer Price Index, 1919 to 2005
# var.ik # 
var.all <- cbind(var.dp, var.dy, var.ep, var.de, var.svar, var.bm, var.ntis, var.tbl,
                 var.lty, var.ltr, var.tms, var.dfy, var.dfr, var.infl) # design matrix of all regressors