rm(list=ls())
# The file provides a comprehensive guide to perform unit-root and structural breaks tests in R
# I employ the weekly time-series of the FTSEMIB index between and

library(readxl)
library(tseries)
library(lubridate)
library(urca)
library(strucchange)

# The data set employed is an example taken from Grillini et al. (2019) 
# "Pricing of time-varying illiquidity within the Eurozone: Evidence using Markov switching liquidity-adjusted capital asset pricing model"
# published in the International Review of Financial Analysis
# available at https://www.sciencedirect.com/science/article/abs/pii/S1057521919300924
# The present sample collects monthly observation for market illiquidity between January 1996 and December 2018 for Belgium

df <- read.table("Belgium", header = T)
bel <- ts(data = df, start = c(1996, 1), frequency = 12)

# Let's plot to have an idea
plot.ts(bel, ylab = "Belgium Illiquidity", main = "Belgium market illiquidity")

# The series looks like there are two well distinct regimes
# denoting periods of high and low illiquidity 
# The most common test is the ADF test. The one below is from the package {tseries}
adf.test(bel, alternative = "stationary")

# The p-value does not reject the null hypothesis of a unit-root, 
# but using the ADF test in the `urca` package, we can exploit different features of this test. 
# We can also test for unit-root with drift and drift and trend. 

summary(ur.df(bel, type = "none", lags = 10, selectlags = "BIC"))
summary(ur.df(bel, type = "trend", lags = 10, selectlags = "BIC"))
summary(ur.df(bel, type = "drift", lags = 10, selectlags = "BIC"))

# All models show the presence of a Unit-root with varying significance lvels, but the ADF test has some issues.
# Therefore, we also test stationarity using Phillips-Perron and KPPS tests. 
# The difference between PP and KPPS tests are that the former has a unit root as null hypothesis,
# likewise ADF, while KPPS has stationarity under the null. 
# The Phillps-Perron test allows for a non parametric correction of the test statistic
# It also allows for a model with constant only or also with trend.

pp1 <- ur.pp(bel, type = "Z-tau", model = 'constant', lags = "long")
pp2 <- ur.pp(bel, type = "Z-tau", model = 'trend', lags = "long")
summary(pp1)
summary(pp2)

# Interestingly, we find opposite evidence using the PP test. 
# What do we do? We can test stationarity using KPPS, or dig deeper using alternative models. 

## Stationarity tests
# While the previous section dealt with unit-root test, here we present a stationarity test.
# We employ the KPSS test, which consists of an LM test under the null of stationarity.
# Kwiatkowski et al (1992) test for a random walk component in the regression.
# We implement this test from the package 'urca, although a version of this test is 
# also available in the package 'tseries'. 

kpps <- ur.kpss(bel, type = "tau", lags = "long")
summary(kpps)

# Although the null hypothesis  of stationarity seems not rejected, confirming evidence 
# of the PP test, one may be interested in testing for structural breaks in the 
# time-series.  

#### Structural breaks  
# We rely on two packages: `strucchange` and `urca`. The former tests parameter instability using two approaches. 
# One based on fluctuation tests (CUSUM and MOSUM residuals). 
# The other based on F-stats, such as the Chow test. 
# Under the first approach, there is a structural change if an appropriately chosen 
# empirical fluctuation process crosses the boundaries that the corresponding limiting process 
# crosses only with a certain probability ($\alpha$). 
# In our case, we firstly need to fit the AR(1) process and then test parameter stability. 


BelLag <- lag(bel, -1)
ar <- ts.union(bel, BelLag)
ar <- na.omit(ar)
SBbel <- efp(bel~BelLag, data=ar, type = 'OLS-CUSUM')
plot(SBbel)
sctest(SBbel)

# Both the graphical and statistical output strongly not reject the presence of a structural break. 
# In particular, from the graph, the fluctuation crosses the boundary around 2013, which corresponds to the "visual" regime change in the time-series plot.  
par(mfrow=c(1,2))
plot(SBbel)
plot.ts(bel, ylab = "Belgium Illiquidity", main = "Belgium market illiquidity")


# The second class of tests are those based on F-stats. In our example:

SBbel2 <- Fstats(bel~BelLag, data=ar)
plot(SBbel2)
sctest(SBbel2)


# Both the plot and the p-value of the Chow test confirm our findings using empirical fluctuation processes. 
# In particular, we can observe again one break point around 2013. 
# There may be however different strucural breaks that are not "visible" from the plot. 
# We test for the number of break dates using the BIC criterion and the residual sum of squares. 
# In this pahse some judjement is necessary, since results may differ

par(mfrow=c(1,1))
BP <- breakpoints(bel~BelLag, data=ar)
plot(BP)

# We observe that the BIC criterion suggests 1/2 breakpoints, while RSS up to four. 
# To simplify we assume one break date in the time-series. 
# Moreover, the graphical investigation of the time-series plot seems to suggest one evident break date.

coef(BP, breaks = 1)
plot.ts(bel, ylab = "Belgium Illiquidity", main = "Belgium")
lines(fitted(BP, breaks = 1), col = 2)
lines(confint(BP, breaks = 1))


# For instance, the break date is at 2013 By also plotting the confidence interval for the break date, we can clearly see that it corresponds to the drop in illiquidity in 2013
# All these evidence provide strong support for our claim that AR(p) processes are econometrically inappropriate if 
# unit-roots and structural breaks are not adequately accounted for. 
# However, the last step is to test stationarity also in the presence of structural breaks. 
# We do this using the `ur.za()` function in the package `urca`. 


zaBel <- ur.za(bel, model = "both", 5)
summary(zaBel)
plot(zaBel)

# Results show that we cannot confidently reject the null hypothesis of a unit root 
# at all significance levels whilst accounting for structural breaks.
# In addition, the break point is identified at 202, a coherent value as compared to all other tests. 