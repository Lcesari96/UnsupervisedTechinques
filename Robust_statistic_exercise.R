# Ex 9.4 + Ex 4 mock exam

# Ex 9.4 

# Aim : Investigate the sensitivity of the two estimators to outliers.

# (a) There is an obvious outlier in the data set. Create a new data set identical to the
# original one, but with this outlier removed.

# Run both lm and lmrob on the new data set, and compare the results to the results
# of the same regression method on the original data (i.e., new lm-result to original
# lm-result, new lmrob-result to original lmrob-result)

# by (i) commenting on how ttest results of the variables have changed qualitatively 
# (significances and positive or negative direction of influence) 

# and (ii) computing the Euclidean distance between the vectors of estimators of 
# the regression coefficients ignoring the intercept (as this value is not comparable
# to the others). 

# These are available as component coef of both the lm and lmrob.output. 
# Comment on what the results mean regarding the sensitivity of the two regressions.


library(robustbase)
setwd("C:/Users/loren/Desktop/Henning_Exam/Data")
unicef <- read.table("unicef97.dat",header=TRUE)

# LS - estimator
uniceflm <- lm(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                 Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid, data=unicef)
# Summary
summary(uniceflm)
# We have 4 significant reg: Literacy.Fem, Drinking.Water, Polio.Vacc and Urban.Pop.
# Multiple R-squared:  0.7587 --> statistical measure of how close the data are to 
# the fitted regression line  

# Residuals vs fitted (NB: y vs x/ std residuals ~ residuals)
plot(uniceflm, which = 1) # LINEARITY --> Relation still well approximated by linearity.

# Residuals vs Normal QQ(normal theoretical quantiles)
plot(uniceflm, which = 2) # NORMALITY --> Overall the residuals can be approximated with a 
# N distribution since most of the points are aligned. We have issues with SaoTP, 
# Angol, Niger and some other observations. They will have to be specifically
# analyzed later.

# Scale-Location ( sqrt std residuals vs fitted values)
plot(uniceflm, which = 3) # HOMOSKEDASTICITY --> seems that is present an issue 
# with HOMOSKEDASTICITY, however weak heteroscedasticity can often be tolerated.
# Residuals vs leverage

plot(uniceflm, which = 5) # FIND INFLUENTIAL CASES --> Not all outliers and 
# leverage points are influential in linear regression analysis. 
# For example Niger is an outlier since shows an high residual values but it doesn't shows
# a significant impact on our model. 
# We can identify some points, far to the right, as good leverage points since, even 
# if are are "outliers" they don't influence our model.
# Difference is for observation SaoTP that is bounded outside the cook distances i.e
# meaning that this obs is a bad leverage points that influence our model.
# It can alter the results if we exclude them from analysis.

# MM - estimator
unicefmm <- lmrob(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                    Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid,
                  data=unicef)
# Summary
summary(unicefmm)
# With a Robust regression 5 regressor becomes significant instead of 4 and the 
# Multiple R-squared:  0.8142 increases. This is probably due to the fact that 
# Robust regressions deal automaticcaly with normality and outliers, that are 
# the principal issues checked before.
# Std residuals vs Robust distance
plot(unicefmm,which=1)
# This plot is divided in 4 regions in which we can distinguish outliers, good and bad 
# leverage points. The issues are given by the bad leverage points that can influence 
# our model as observation SaoTP.
# Residuals vs Normal QQ 
plot(unicefmm,which=2)
# As we can see Robustness weights downweights outliers and bad leverage points
# "automatically".
# Robustness weights
plot(1:121,unicefmm$rweights,xlab="Observation number",ylab="Weight",
     main="MM-estimator, robustness weights",type="n")
text(1:121,unicefmm$rweights,rownames(unicef),cex=0.7)
which(unicefmm$rweights == 0)
# The Robustness weights sets equal to 0 : Angol, Niger and SaoTP observations, since
# are considered as "influnce" observations i.e some as outliers as Angol and Niger,
# and SaoTP as bad leverage point.
# Residuals vs. fitted
par(mfrow=c(1,1))
plot(uniceflm, which = 1)
plot(unicefmm$fitted,unicefmm$residuals,xlab="Fitted values",ylab="Residuals",
     main="MM-estimator, residuals vs. fitted")
# However a robust estimator will not fit all points well and will fit the outliers 
# even worse than LS because it weights them out for fitting.
# Thus is better run it with the L1/absolute value or other robust loss when 
# attempting to fit the "good majority" well.

# New dataset without the crucial outliers: SaoTP.
unicef1 <- unicef[-91,]
uniceflm1 <- lm(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                 Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid, data=unicef1)

a <- summary(uniceflm)
b <- summary(uniceflm1)
a$coefficients
b$coefficients
# The R squared increases, the significant coeff are the same but increases
# those significance up to 0 '***', moreover the variables Literacy.Ad chages
# sign and the the intercept decreases.

plot(uniceflm1,which=1)
# Lineraity holds.

plot(uniceflm1,which=2)
# Normality is better, we also have the same issues with previous ouliers as before.
# It seems that we have a "masking" case i.e now SierraL is considered as outlier thus
# this means that it was "masked" by SaoTP.

plot(uniceflm1,which=3) # seems that is present an issue with homoskedasticy,
# however weak heteroscedasticity can often be tolerated.

plot(uniceflm1,which=5)
# No bad leverage points detected, now good leverage points are more evident and 
# also the "masking" case.

unicefmm1 <- lmrob(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                    Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid,
                  data=unicef1)
summary(unicefmm1)
# The summary measures are the same since the Robustness weights will set 
# equal to 0 all the outliers and bad leverage points.
plot(unicefmm1,which=1)
# As we said before, SierraL now is dected as critical point
# since it was "masked" by SaoTp.
plot(unicefmm1,which=2)

plot(1:120,unicefmm1$rweights,xlab="Observation number",ylab="Weight",
     main="MM-estimator, robustness weights",type="n")
text(1:120,unicefmm1$rweights,rownames(unicef),cex=0.7)
which(unicefmm1$rweights == 0)
# Now only Angol and Niger are set = 0 from the Robustness weights since SaoTP
# is were already removed. 

# Euclidean distance between the vectors of estimators of the regression 
# coefficients ignoring the intercept (as this value is not comparable
# to the others). 
# Comment on what the results mean regarding the sensitivity of the two
# regressions.
a <- uniceflm$coeff
b <- unicefmm$coeff
a
b
d <- dist(a, b, method = "euclidean")
# The LS regression have a more "high" values for the regressors, this meas that 
# is more sensitive to changes with respect the input variables.
# This is confirmed by the distances, since the abs value is involved, we give a
# "measure" of how this sensitivity changes in both directions(signs).
# The baseline is that adding or removing outliers and "noise" changes lm much more than
# lmrob.

# Ex 4 mock exam

# (a) Which of the regression estimators do you find most trustworthy here and why (you
# can comment on known characteristics of the methods but you are also expected to
# use the data analysis for arguing your decision)?

# INTRODUCTION
# LS/MM/S
# The LS estimator i.e the "classical" regression, give us the highest asymptotic efficiency but it
# is very sensitive to the outliers.
# The basic idea underlying the robust linear model (in this case with MM and S estimators) 
# is that some of the data are distributed conditionally normal and the remaining fraction comes 
# from some arbitrary distribution i.e. the outliers. 
# The goal of a robust method is to estimate the parameters (beta and sigma^2) of 
# this conditional normal distribution without giving the outliers too much influence
# by using robustness weights. 
# Thus the main difference among the two regression is that in the robust one some observations 
# are weighted.
# The MM estimator combines automatically high efficiency and high BP ~ 1/2 i.e an high
# tolerance with respect outliers. The general idea is that we combine an M estimation for the 
# regression coeff and an S estimation for the scale, then the efficiency comes from the M and 
# the robustness from the S.ù
# It uses a bi-square redescending score function, and it returns a highly robust 
# and highly efficient estimator (with 50% breakdown point and 95% asymptotic 
# efficiency for normal errors).
# The S estimator has an high BP ~ 1/2 but low asy efficiency(24%). 
# S estimator is used to inizialize the MM one.

# (a) By looking the summaries the R^2 is higher for the robust regression with the MM-estimator.
# From the diagnostics plot of the "classical" linear regression we can see that there are
# issues in normality and homoskedasticity but, most important, a bad leverage points
# identified by its high cook's score in the Residual vs leverage diagnosis plot.
# Probably most of the issues are related to this observation (SaoTP) and other two 
# identified in the normal QQ plot (Angola and Niger). For me seems reasonable choose the 
# robust regression with the MM estimator for the following reasons:
# 1) Robust regression automatically deals with normality and outliers by weighting them
# out of fitting, and those issues are the issues that characterize our dataset.
# 2) MM over S due to its higher efficiency. Moreover if we look the plot of the robustness weights,
# the S estimator put = 0 too many observations that identify as outlier, this implies a strong loss
# of information (thus less efficiency) in our model.
# Moreover, by looking the last two plot : residual vs fitted, for the MM the residuals seems more 
# "compact" and well distributed.

# (b) Which of the covariance matrix estimators do you find most trustworthy here and
# why (you can comment on known characteristics of the methods but you are also
# expected to use the data analysis for arguing your decision)?

# INTRODUCTION
# "Classical" COV/ MCD alpha = 0.5/ MCD alpha = 0.75

# ROBUST MULTIVARIATE LOCATION AND COVARIANCE ESTIMATOR
# Minimum covariance determinant estimator MCD 
# Situation : we need to estimate the location and the covariance matrix of
# Y1, ..., Yn iid with Y1 ~ Np(mu, capitalSigma) (Capital sigma since we are in the 
# multivariate case) that are potentially violated by outliers.

# Solution : MCD i.e  a robust multivariate estimator for location(mu) and 
# covariance matrix(capitalSigma)

# 1) Given n/2 < h < n we need to find the h (> n/2) observations Y* whose classical 
# covariance matrix has the lowest possible determinant.
# 2) Then the Tmcd location estimator will be Tmcd(Y1, ..., Yn) = sum(for i=1 to h){Y*i} 
# (Sample mean estimator on the h Y*)
# 2.1) while the Smcd coovariance matrix estimator will be:
# Smcd(Y1, ..., Yn) = c * Cov(Y*1, ..., Y*h) ( c * Sample covariace estimation 
# (classical covariance matrix on the h Y*))
# c > 0 (consistency factor) will be chosen for consistency under the  ~ Np(mu, capitalSigma) 
# distribution depending on h(n of chosen obs) and p(number of variables)

# The MCD method looks for the h (> n/2) (h = h(??,n,p) = h.alpha.n(alpha,n,p)) observations 
# (out of n) whose classical covariance matrix has the lowest possible determinant.
# The raw MCD estimate of location is then the average of these h points,
# whereas the raw MCD estimate of scatter is their covariance matrix, multiplied by a 
# consistency factor (.MCDcons(p, h/n)) to make it consistent at the
# normal model and unbiased at small samples. 
# h is specified through parameter alpha with default 0.5,
# h is chosen so that h/n = ?? approx ?? with automatic rounding.
# BP is about equal to 1 - ?? (equal for n to infinity).

# The MCD estimator for loc and cov identifies far more outliers than "classical" estimators 
# (in "classical" other outliers are "masked" by the largest ones i.e a "far" outlier can make a
# "closer" outlier treated as a "normal" obs and vicerversa  --> Masking.)

# NB by setting an higher ?? = 0.75 (where h is chosen so that h/n = ??) implies to have
# an higher h i.e we are a little bit more tolerate against outliers
# (in this case ?? = 0.5 and ?? = 0.75 give almost identical results.)

# The idea is the same as the robust regression i.e some of the data are distributed
# conditionally normal and the remaining fraction comes from some arbitrary distribution 
# i.e. the outliers. --> we estimate the location and the cov on the fraction of the 
# data not consider as outlier.

# OUTLIERS IDENTIFICATION 

# In standard situations, outliers can be identified by looking at distances btw obs 
# from the mean or multivariate mean (if p >= 2) or their robust counterparts
# (in regression use residuals)

# REC dist btw obs --> is given by the sum of the differences of each variables of the obs, 
# thus is a number (in the case obs vs mean -> differences of each variables with the
# multivariate means for each variables)

# Under normal distribution, squared Mahalanobis distances of the obs from mean
# are chi-squared(p)-distributed, and a distance larger than a high quantile
# (say 0.99 - this is a subjective choice) can identify an outlier.

# It's much better to identify outliers from robust distances and
# the MCD object can be used to compute Squared robust Mahalanobis distances 
# --> High values of this robust distance (thus SQRT(mcdd$mah)) from obs wrt estimated locations may identify outliers
# $mah --> squared robust Mahalanobis distances btw observation and the final estimate of the location
# and scatter.

# (b) The MCD is used to make the robust estimations of the location(mu) and covariance matrix
# (capitalSigma). 
# Outliers can be idenfied by looking the distances between observation and multivariate 
# means. 
# Under a normal distribution, the squared mahalanobis distances between observations and means
# are a chi-squared(p = n of var) distributed and we can use those to identify outliers.
# The squared mahalanobis distances larger than an arbitrary quartile (here is used 0.99) 
# identify outliers.
# The MCD object can be used to do this since it contains the squared robust mahalanobis distances
# obtained from the estimated Smcd(Robust estimation of the covariance matrix by the MCD)
# between observations and Tmcd(Robust estimation for the locations by the MCD).
# Those are better than the standard squared mahalanobis distances obtained from the classical 
# covariance matrix between observations and the classical multivariate means since are 
# robust distances.
# The idea behind the MCD estimator is to find h (> n/2) observations that gives us the 
# minimum determinant for the covariance matrix, and then apply the classical 
# estimators on this h observations: sample mean and classical covariance.
# It is inline with the basic idea behind the general robust estimation theory i.e, we try to
# split the outliers from the observations and then we try to make the estimations by consider 
# only the observations without giving importance to the outliers.
# From this, we obtain the raw location and covariance MCD robust estimations.
# R, beyond the raw ones, give us somewhat improved versions with better efficiency 
# but same robustness.
# h is specified through parameter alpha with default value of 0.5 where it is chosen 
# so that h/n = ?? approximate ?? with automatic rounding.
# Thus by changing the value of ?? we change the number of considered observations h.
# In the appendix analysis we have 4 plots: the first 2 for the squared robust mahalanobis
# distances for ?? = 0.5 and 0.75 obtained from Smcd and Tmcd, while the 2nd two for the 
# squared robust mahalanobis distances for ?? = 0.5 and 0.75 obtained from Smcd and Tmcd
# vs the standard squared mahalanobis distances. 
# The MCD with ?? = 0.5 detects more outliers with respect the MCD with ?? = 0.75 while 
# the Squared standard Mahalanobis distance detects very few of those.
# Since by the previous analysis we have seen a strong presence of outliers i exclude 
# the classical covariance estimator. 
# By looking the plots with an MCD with ?? = 0.5 seems that we are treating too much 
# observations as outliers (also observations with not so high squared robust 
# mahalanobis distances) and since we are looking to the best trade-off between 
# asymptotic efficiency and robustness i choose the MCD estimators with ?? = 0.75.

# c) Which of the countries do you think are outliers in the sense that they seem to behave
# substantially different from the others (you can use the abbreviated country names
# as in the plots), based on which plots or results?

#(c) With the appendix analysis we have multiple way to check this.
# First by checking the first 4 diagnostic plots from the Least Squares regression
# we have different highlighted points: SaoTP, Angola and Niger from the Normal QQ-plot;
# Same by the Scale-Location; Niger, Ruand and SaoTP from the Residuals vs Leverage, where 
# Niger and Ruand seems simple outliers while SaoTP is an evident bad leverage points 
# identified by its high cook's score.
# The diagnostic plots from the Robust regression confirm what we have said before and 
# moreover from the Standardized residuals vs. Robust Distances plot we can see some
# good leverage points that, however, are not a problem for our model.
# To conclude by looking the MCD plots we can say that: with ?? = 0.5 a lot 
# of observations seems to have issues but for me this is a too "robust" approach 
# that penalize to much the efficiency; with ?? = 0.75 the points with higher Squared 
# robust Mahalanobis distance are SaoTP, Mosam, Ruand and Samb; the Squared standard 
# Mahalanobis distance detects very few outlier like Haiti, nick and SaoTP.
# To summaryze we can say that:
# 1) We need to take care when we look to diagnosis that not involve robust methods
# due to the presence of the "masking" effect;
# 2) With the MCD results basically we are looking only to the distances between 
# observation and means based on a reference covariance matrix, and thus we can't
# really distinguish among outliers, good and bad leverage points for a model fitting;
# 3) SaoTP is the most reccurent problematic observation, that is classified as a
# bad leverage points for linear regression by both plots: Standardized residuals vs. Robust Distances
# and Residuals vs Leverage.
# 4) Angol and Niger are simply outliers and we can see it from Standardized residuals vs. Robust Distances
# plot due to those high values for the residuals. Also Ruand is an outlier from
# the Residuals vs Leverage plot.

# Other points detected in the MCD plots can be "lower" outlier or god leverage points.

# We can also looking, for a more precise classification, to both MM???estimator,
# residuals vs. fitted and the MCD plots and say: Niger, Angol, Liber, SierraL, CongoDr, Guin, 
# Haiti, Nic are outliers while Mosam is god leverage point and SaoTp a bad one.

# To conclude by looking the plot of the MM???estimator, robustness weight the model 
# treat as outlier or bad leverage points by putting is weights =0 and then by 
# putting those out of fitting: Angol, Niger and SaoTP, while other problematic 
# observations are simply downweighted.

# (d) For outlier identification, two different kinds of analyses were run here, namely (i)
# regression (least squares, MM, and S), and (ii) covariance matrix estimation based
# on all variables (standard, and MCD with ?? = 0.5 and ?? = 0.75). In what sense
# are outliers identified by regression different from outliers identified from covariance
# matrix estimation?

# In case (i) the outliers are identified throught residuals i.e an high residual value
# for an observation identify it as an outlier, that is the distance from the observations 
# and our regression line. 
# Moreover in (i) we can distinguish between outliers and leverage points, where
# outliers are associated with the above definition while leverage points are 
# "outliers" with respect to the scale of the explanatory variables i.e the 
# the regression variable(X) used to predict the outcome(Y).
# In (ii) we identify outliers with respect to an estimation (not wrt to a model),
# basing our decision on robust distances(squared mahalanobis) between observations 
# and multivariate robust location estimate(Tmcd), obtained from a robust estimated 
# covariance matrix(Smcd).
# In (ii) we only identify outliers as plus or less critical observations that can influence 
# our estimation, we not distinguish between outliers and leverage points since the aim
# is different.
# To summaryze in (i) we aim to predict something and thus we identify outliers
# with respect to our model in order to reweight them out of the fitting(or 
# downweight their influence), while in (ii) we identify outliers in order to try
# to estimate our location and covariance matrix by taking into account only
# the observations not considered as outliers without giving importance to those.
# This is the basic idea about robustness.

# (e) A social scientist suggests that the observation "SaoTP" should be removed, because
# its level of foreign aid makes it essentially different from all other observations, and
# it should therefore not be used in the same analysis. What do you think of this
# suggestion? Which of the three regression estimators would in your opinion be most
# affected by such a decision?

# I think that we shouldn't remove this observation since due to the nature of the 
# data where observation represent "regions" and vast areas (are a sort "clusters" of the 
# population) by removing even just one of those it would result in a big loss of
# informations for our model. It is itself a result from which we can understand
# where we need to concentrate our attention. Moreover robust regression were
# developed in order to deal with kind problem.
# The LS regression estimator is the one that will be most affected by such a decision
# since adding or removing outliers and "noise" changes lm(LS regression ) much more than
# lmrob (MM-robust regression) due to the reweighting features of the robust regression.
  
# S-estimator
unicefmm$init.S
# Summary 
unicefmm$init.S
# Robustness weights
plot(1:121,unicefmm$init.S$rweights,xlab="Observation number",ylab="Weight",
     main="S-estimator, robustness weights",type="n")
text(1:121,unicefmm$init.S$rweights,rownames(unicef),cex=0.7)
# Residuals vs fitted
plot(unicefmm$init.S$fitted,unicefmm$init.S$residuals,xlab="Fitted values",
     ylab="Residuals",main="S-estimator, residuals vs. fitted")

