# master-thesis
This project is a part of my Master thesis, topic "LASSO-based predictive regression for stock returns"

INFORMATION
The project's analysis is written in R, and the paper itself is written in LaTeX (MiKTeX)
The project includes:
  - Data Generating Processes that emulate the predictive environment concerned in the paper
  - Simulation study: using LASSO (Tibshirani, 1996) and Adaptive LASSO (Zou, 2006). The theory closely follows Knight (2008) and Lee et al. (2018). Simulation study follows Lee et al. (2018) with some modification.
  - Emprical study: comparing performance of LASSO and Adaptive LASSO to OLS and Random Walk with Drift in predicting Equity Premium on S&P 500 index using valuation ratios and macroeconomic variables as predictors. Data is taken from the update dataset used in Welch and Goyal (2008). Can be found in Goyal's website.
  - Custom functions to do the job
  
INSTRUCTION
Please install and call the following packages before running any script:
```
### Package installation
install.packages('readxl')
install.packages('zoo')
install.packages('glmnet')
install.packages('tseries')
install.packages('caret')
install.packages('Matrix')
install.packages('corrplot')
install.packages('snow')
install.packages('foreach')
install.packages('tictoc')
install.packages('VIF')
install.packages('tikzDevice')
install.packages('gridExtra')
install.packages('dplyr')
install.packages('ggplot2')

### Call packages
library('readxl')
library('zoo')
library('glmnet')
library('Matrix')
library('tseries')
library('caret')
library('Matrix')
library('corrplot')
library('parallel')
library('snow')
library('MASS')
library('tictoc')
library('foreach')
library('VIF')
library('ggplot2')
library('tikzDevice')
library('gridExtra')
library('dplyr')
```
