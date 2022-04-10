library(pdfCluster)
library(cluster)
library(fpc)
library(prabclus)
library(smacof)
library(tidyverse) 
library(cluster)
library(dbscan)
library(factoextra) # for silhouette plots
library(mclust)
library(robustHD) # MAD std  
# install.packages("C:/Users/loren/Desktop/EMMIXskew_1.0.3.tar.gz",repos=NULL,type="source")
library(EMMIXskew)
library(poLCA)
library(nomclust)
library(fda) # Functional data analysis
library(funFEM) 
library(robustbase)
library(mixtools)
library(tclust)
set.seed(123456)

# INTRODUCTION
par(mfrow=c(1,1))
setwd("C:/Users/loren/Desktop/Henning_Exam/Data")
dortmund <- read.table("Dortmund_msbd.dat",header=TRUE,row.names=1)
mean(dortmund$birthdeath) # -0.2316692 Sample mean estimator Tn = 1/n * sum(i=1 to n){xi}
sum(dortmund$birthdeath>-0.2317) # 165 obs are over the Sample mean!

boxplot(dortmund$birthdeath,main="dortmund$birthdeath")
boxplot(dortmund$birthdeath,main="dortmund$birthdeath",ylim=c(-0.9,0.135))
abline(mean(dortmund$birthdeath),0,col=2) 

which.min(dortmund$birthdeath) # obs 133 is an outlier
mean(dortmund$birthdeath[-133]) # -0.03481514 Sample mean estimator obtained after
# removing the outlier 

# KEY PROBLEM : outlier can drive arbitrally far away our estimation

# SOLUTION : robust methods can accommodate outliers without giving them too much 
# influence.

# ROBUSTNESS THEORY AND ROBUST ESTIMATORS

# IMPORTANT
# BASIC IDEA behind the ROBUST THEORY:
# The basic idea underlying the robust theory is to assume that some of the data
# is distributed conditionally normal while the remaining fraction comes from some
# arbitrary distribution that are the outliers. 
# The goal of a robust method is to estimate the parameters of this conditional
# normal distribution without giving to the outliers too much influence. 

# NB -> The BP determine how much can be bigger the fraction of outliers that an 
# estimator can handle!

# M - ESTIMATORS (Robust estimators) (DEEPER EXPLANATION AT THE END)

# HUBER M - ESTIMATOR: M - estimator for the location (i.e replace the sample mean or the 
# median) but require a scale s (i.e if we have a N(0, sigma^2) --> sigma^2 is a scale s
# parameter)

# It behaves like the sample mean if u <= gamma (i.e when we do not have much outlier
# influence) thus it exploit the efficiency of the sample mean estimator
# (where u = y - t i.e the diff btw real values and estimates, and gamma is a sort of 
# outlier influence treshold --> if is higher we have strong outlier influence),
# while when u > gamma it behaves like the median and exploit is bounded influence 
# function and high BP i.e., robstness.

# NB both sample mean and median are location estimators :
# mean --> best EFFICIENCY low ROBUSTNESS (BP = 0), 
# median --> high ROBUSTNESS, low EFFICIENCY.

# Rec
mean(dortmund$birthdeath) # -0.2316692
median(dortmund$birthdeath) # -0.01918741
sd(dortmund$birthdeath) # -0.01918741

# MAD : robust scale estimator in order to obtain s for apply the HUBER M - ESTIMATOR
mad(dortmund$birthdeath, constant = 1.4826)  # 0.04708085
# constant = 1.4826 (default) --> is a factor used to enforce the convergence:
# MAD(Y1, ..., Yn) ->p to sigma2 i.e the convergence of our estimator to the true 
# values

# HUBER M - ESTIMATOR
HM <- huberM(dortmund$birthdeath,k=1.5) 
# k = 1.5 (default) --> tuning parameter(optimizzator) ca.(approximately): 
# 95% efficiency with Normal models by preserving the robustness

HM$mu # -0.02259647

HM$s # 0.04708085 # Uses the MAD by default and gives it out.
# ( = mad(dortmund$birthdeath, constant = 1.4826) )

HM$it # 8 # Number of iterations of algorithm

HM$SE # NA # If se=TRUE, huberM estimates its own standard error, can be used for tests 
# and confidence intervals

boxplot(dortmund$birthdeath,main="dortmund$birthdeath")
boxplot(dortmund$birthdeath,main="dortmund$birthdeath",ylim=c(-0.9,0.135))
abline(mean(dortmund$birthdeath),0,col= 1) # Sample mean (black)
abline(HM$mu,0,col= 2) # HUBER M - ESTIMATOR (red)
abline(median(dortmund$birthdeath),0,col= 3) # it behaves like the median (green)

# ROBUST MULTIVARIATE LOCATION AND COVARIANCE ESTIMATOR MCD

pairs(dortmund) # Multivariate(p) --> before by consider only $birthdeath
                #                     univariate(p = 1)

# Minimum covariance determinant estimator MCD 
# Situation : we need to estimate the location and the covariance matrix of
# Y1, ..., Yn iid with Y1 ~ Np(mu, capitalSigma) (Capital sigma since we are in the 
# multivariate case) that are potentially violated by outliers.

# Solution : MCD i.e  a robust multivariate estimator for location(mu) and 
# covariance matrix(capitalSigma)

# 1) Given n/2 < h < n we need to find the h (> n/2) observations Y* whose classical 
# covariance matrix has the lowest possible determinant.
# 2) Then the RAW Tmcd location estimator will be :
# Tmcd(Y1, ..., Yn) = sum(for i=1 to h){Y*i} (Sample mean estimator on the h Y*)
# 2.1) While the RAW Smcd covariance matrix estimator will be:
# Smcd(Y1, ..., Yn) = c * Cov(Y*1, ..., Y*h) i.e.,
# c * Sample covariace estimation (classical covariance matrix on the h Y* where 
# c > 0 (consistency factor) will be chosen for consistency under the conditionally
# normal ~ Np(mu, capitalSigma) distribution depending on the h(n of chosen obs)
# and p(number of variables).

# IMPORTANT
# RECALL THE BASIC IDEA BEHIND ROBUST TECHNIQUES: The basic idea that underlying
# the robust methods is that some of the data are assumed to be distributed conditionally
# normal(in this case the h) and the remaining fraction(in this case the n-h) comes 
# from some arbitrary distribution (i.e., the outliers). 

cd <- cov(dortmund) # Sample covariace estimation (classical covariance matrix)

mcdd <- covMcd(dortmund) # MCD
help("covMcd")
# The MCD method looks for the h (> n/2) (h = h(alpha,n,p) = h.alpha.n(alpha,n,p))
# observations (out of n) whose classical covariance matrix has the lowest possible
# determinant.
# The raw MCD estimate of location is then the average of these h points, whereas
# the raw MCD estimate of scatter is their covariance matrix, multiplied by a consistency
# factor (.MCDcons(p, h/n)) to make it consistent at the normal model
# (~ Np(mu, capitalSigma))and unbiased at small samples. 

# h is specified through parameter alpha with default 0.5,
# h is chosen so that h/n = alpha approx alpha with automatic rounding.
# BP is about equal to 1 - alpha (for n to infinity).

# IMPORTANT
# NB by setting an higher alpha = 0.75 (where h is chosen so that h/n = alpha)
# implies to have an higher h i.e we are a little bit more tolerate against outliers,
# since we are selecting an higher number of observations to consider as "good" one. 
# NB h/n is a proportion thus by setting higher its result(alpha) implies to have 
# an higher h.

mcdd75 <- covMcd(dortmund,alpha=0.75)
# MCD(alpha=0.75) is some kind of compromise between the classical mean/cov-matrix
# and the MCD with alpha = 0.5, but results are not always in between.

colMeans(dortmund) # Sample means for each variables 
mcdd$center # Tmcd with alpha = 0.5 for each variables 
mcdd75$center # Tmcd with alpha = 0.75 for each variables 

# Estmated covariance matrix is in component cov:
str(mcdd$cov)

# The covMcd output also has components raw.cov and raw.center.
# These give the MCD how it was defined above (slides), whereas cov and center
# give somewhat improved estimation versions of covariance matrix (capitalsigma)
# and location(mu) of a Y~ Np(mu, capitalSigma) with better efficiency but same
# robustness.

# These are really not very different here:
str(mcdd$raw.cov)
mcdd$raw.center
mcdd$center

# Visualize differences on 2 variables:

# Used for showing the ellipses defined by the cov-matrices 

# Looking at variables 2 and 4.
v1 <- 2 
v2 <- 4
plot(dortmund[,c(v1,v2)])

ellipse(colMeans(dortmund[,c(v1,v2)]),cd[c(v1,v2),c(v1,v2)],alpha=0.01) # Sample
# covariace estimation (classical covariance matrix) --> alpha=0.01 here means that 
# if data are normally distributed, 99% are expected to be inside the ellipse.

ellipse(mcdd$center[c(v1,v2)],mcdd$cov[c(v1,v2),c(v1,v2)],col=2,alpha=0.01) # Smcd 
# with alpha = 0.5

ellipse(mcdd75$center[c(v1,v2)],mcdd75$cov[c(v1,v2),c(v1,v2)],col=4,alpha=0.01) 
# Smcd with alpha= 0.75

text(3,10.5,"standard cov")
text(4.2,7.5,"MCD h/n=0.75",col=4)
text(3.7,6,"MCD h/n=0.5",col=2)

#NOTES ABOUT CORRELATIONS

# Sample coorelations :
cor(dortmund[,v1],dortmund[,v2])
cor(dortmund[,v1],dortmund[,v2])

# The following computes "robust correlation" for cov-entries:
mcdd$cov[v1,v2]/sqrt(mcdd$cov[v1,v1]*mcdd$cov[v2,v2])
mcdd75$cov[v1,v2]/sqrt(mcdd75$cov[v1,v1]*mcdd75$cov[v2,v2])

# NOTICE
# NB Classical mean and cov-matrix don't differ whether computed on 2 variables,
# or on 7 variables with values for 2 variables extracted.
# This does not hold for MCD! (Use computations from 7 variables here.)

#i.e.,
mcdd751 <- covMcd(dortmund[,c(2,4)],alpha=0.75) # computed on 2 var
mcdd751$cov
mcdd75$cov[,c(2,4)] # computed on 7 var (before) + extraction of 2 var
# mcdd751$cov != mcdd75$cov[,c(2,4)]

# IMPORTANT REC SUMMARY
# NB by setting an higher alpha = 0.75 (where h is chosen so that h/n = alpha)
# implies to have an higher h i.e we are a little bit more tolerate against outliers,
# since we are selecting an higher number of observations to consider as "good" one. 
# NB h/n is a proportion thus by setting higher its result(alpha) implies to have 
# an higher h.
# Since, RECALL THE BASIC IDEA BEHIND ROBUST TECHNIQUES: The basic idea that underlying
# the robust methods is that some of the data are assumed to be distributed conditionally
# normal(in this case the h) and the remaining fraction(in this case the n-h) comes 
# from some arbitrary distribution (i.e., the outliers). 

# OUTLIERS IDENTIFICATION 

# In standard situations, outliers can be identified by looking at distances btw obs 
# from the mean or multivariate mean (if p >= 2) or their robust counterparts.
# (in regression we use residuals i.e the difference between the true points Y and 
# the regression fitted line Yhat)

# REC dist btw obs --> is given by the sum of the distances btw each variables of
# the obs, thus is a number (in the case obs vs mean -> differences of each variables
# with the  means for each variables(multivariate means) --> While for the MD 
# distance: distance btw the obs and an "overall centroid", a special "centroid" 
# i.e a central point in the multivariate space where all the means from all the 
# variables intersect, a sort of overall mean for multivariate data.)
# NB LOOK AT LINK TO THE MD DISTANCE IN THE DISTANCE_BASED_METHOD SCRIPT FOR A 
# DEEPER EXPLANATION.

# Under normal distribution, squared Mahalanobis distances(thus the "normal" MD 
# since squared is implied) of the obs from mean are chi-squared(p)-distributed,
# and a distance larger than a high quantile (say 0.99 - this is a subjective 
# choice) can identify an outlier.

# It's much better to identify outliers from robust distances, the MCD object 
# contains the squared robust Mahalanobis distances (a sort of robust
# version of the traditional Leverage obtained from MD distance (look slides))
# --> High values of this robust distances(SQRT(mcdd$mah)) btw obs and the "overall 
# centroid" obtained basing on the robust estimated locations(Tmcd) and the robust
# estimated covariance matrix Smcd, may identify outliers.

# NB The leverage is obtained just by the sqrt MD basing on the sample mean and 
# cov matrix estimation, if those are robust (i.e., Tmcd and Smcd) we obtain the 
# robust mahalonobis distances.

# FOR EXAM: Since is better identify outliers from robust distances(masking issue),
# the MCD object contains the computed squared mahalanobis distancebased on the 
# Tmcd for locations and the Smcd for the covariance matrix (instead of the sample
# estimations).
# High values of this distances for observation may identifies outliers.

# covMcd$mah --> robust mahalanobis distances btw the observations and the "overall 
# centroid" obtained by using the robust estimated locations(Tmcd) and also as input,
# together with the data, the robust estimated covariance matrix Smcd.


# PLOT FOT OUTLIERS DETECTION WITH THE MCD ESTIMATOR 

# REC : mcdd <- covMcd(dortmund), 1:170 n of obs
plot(1:170,sqrt(mcdd$mah),type="n",xlab="Observation",
     ylab="Squared robust Mahalanobis distance")
text(1:170,sqrt(mcdd$mah),rownames(dortmund),cex=0.7)
abline(sqrt(qchisq(0.99,7)),0,col=2) # qchisq(0.99,7) --> 0.99 high quantile, 7 p
                                     # n of variables

# Should look at smaller values to see more precisely what's going on:
# (i.e we need to change YLIM)
plot(1:170,sqrt(mcdd$mah),type="n",ylim=c(0,30),xlab="Observation",
     ylab="Squared robust Mahalanobis distance")
text(1:170,sqrt(mcdd$mah),rownames(dortmund),cex=0.7)
abline(sqrt(qchisq(0.99,7)),0,col=2)

# Compare the simple sqrt Mahalanobis distances with the sqrt robust Mahalanobis
# distances in the task of OUTLIER DETECTION.

# (Always look at the YLIM values)
plot(sqrt(mahalanobis(dortmund,colMeans(dortmund),cd)),sqrt(mcdd$mah),type="n",
     xlim=c(0,10),ylim=c(0,30),xlab="Squared standard Mahalanobis distance",
     ylab="Squared robust Mahalanobis distance")
text(sqrt(mahalanobis(dortmund,colMeans(dortmund),cd)),sqrt(mcdd$mah),
     rownames(dortmund),cex=0.7)
abline(sqrt(qchisq(0.99,7)),0,col=2)
abline(v=sqrt(qchisq(0.99,7)),col=2)

# HOW TO INTERPRET THOSE PLOTS? --> The orizzontal red line corresponds to the y 
# axis i.e., if the y axis = sqrt robust Mahalanobis distances then all the obs
# above the orizzontal red line are identified as outliers.
# While the vertical red line correspond to the x axis where all the obs on the 
# right of the vertical red line are identified as outliers.

# Compare with alpha=0.75
plot(sqrt(mcdd75$mah),sqrt(mcdd$mah),xlim=c(0,10),ylim=c(0,30),
     xlab="Squared robust Mahalanobis distance (alpha=0.75)",
     ylab="Squared robust Mahalanobis distance (alpha=0.5)")
abline(sqrt(qchisq(0.99,7)),0,col=2)
abline(v=sqrt(qchisq(0.99,7)),col=2)

# NB by setting an higher alpha = 0.75 (where h is chosen so that h/n = alpha)
# implies to have an higher h i.e we are a little bit more tolerate against outliers,
# since we are selecting an higher number of observations to consider as "good" 
# one(y*). 
# NB h/n is a proportion thus by setting higher its result(alpha) implies to have 
# an higher h.

# (in this case alpha = 0.5 and alpha = 0.75 give almost identical results.)

# The MCD estimator for loc and cov identifies far more outliers than "classical"
# estimators(in "classical" other outliers are "masked" by the largest ones i.e.,
# a "far" outlier can make a "closer" outlier treated as a "normal" obs and vicerversa
# --> this phenomenon is called Masking.)

#IMPORTANT
# NB Through the MCD we identify only the "outliers" deviate from an assumed Gaussian 
# distributional shape.(Thus we avoid masking)

# Since, RECALL THE BASIC IDEA BEHIND ROBUST TECHNIQUES: The basic idea that underlying
# the robust methods is that some of the data are assumed to be DISTRIBUTED CONDITIONALLY 
# NORMAL(in this case the h) and the remaining fraction(in this case the n-h) comes 
# from some arbitrary distribution (i.e., the outliers). 

# IMPORTANT
# This does not mean that they are all meaningfully different from the main core of 
# observations!
# Such outliers can be faulty measurements, correct but special,or just result from
# non-normality of the process! --> 

# In ROBUST REGRESSION we can DISTINGUISH BTW OUTLIERS,GOOD or BAD LEVERAGE POINTS
# since we have a task as reference.

# ROBUST REGRESSION SSION.

# CLASSICAL LINEAR REGRESSION
set.seed(3355111)
x1 <- runif(100)
x2 <- x1+runif(100,-0.3,0.3)
x3 <- -x1+rnorm(100,sd=0.5)
y <- 1*x1-1*x2+3*x3+rnorm(100,sd=0.5)

regdata1 <- as.data.frame(cbind(y,x1,x2,x3))
pairs(regdata1)

lmrd1 <- lm(y~x1+x2+x3,data=regdata1)
summary(lmrd1)
# NB R^2 --> tell us how well the regression model fits the observed data.

par(mfrow=c(2,2))
plot(lmrd1)

# CHECKING THE MODEL ASSUMPTIONS i.e Check if the model works well for the data:

# A) Linearty -> linear relationship among indep var(Y) and dep var(X).

# B) Normality -> residuals should distributed as a ~ N, i.e., if residuals doesn't
# shows anomalies i.e., are distributed according to a normal, it means that our 
# prediction/output is close to the true values(since the residuals are given by 
# the difference btw the true Y and the Yhat = X*Bhat fitted regression line(our 
# model)).

# NB OUTLIERS and strong skewness are problems for normality; many other violations 
# of normality are not(since with large data we can exploit the CLT).

# C) Homoskedasticty -> residuals should have the same constant variances. Strong 
# violation will affect inference, weak heteroscedasticity can often be tolerated.
# NB linked to what stated before, heteroskedastic residuals can shows anomalous 
# / more varying patterns, where for residuals is not good since means that our
# model doesn't predict well.

# D) Independence -> indep amog indep variables Xs (no cov/corr tolerate). Hard 
# to diagnose, may affect inference.

# Visible issues with skewness or heteroscedasticity may point to transformations.

help("plot.lm")
# Six plots (selectable by which) are currently available:
# 1) Plot of residuals against fitted values
# 2) Scale-Location plot of sqrt(| residuals |) against fitted values
# 3) Normal Q-Q plot,
# 4) Plot of Cook's distances versus row labels 
# 5) Plot of residuals against leverages (with cook's distance also)
# 6) Plot of Cook's distances against leverage/(1-leverage). 

# NB By default: 1,2,3 and 5 are provided and here will be our main focus.

# IMPORTANT
# The numbers next to some points(i.e., its row number) in each plot indicates 
# EXTREME VALUES that could shows POTENTIALLY PROBLEMATIC OBSERVATIONS according 
# to each different criterion for each different plots(i.e., 1,2,3 and 5 plot)
# and are identified by the row numbers in the data set.

# 1) Residuals(y.axis) against fitted values(x.axis) (for checking the Linearity
# assumption):
# This plot shows if residuals have non-linear patterns. 
# There could be a non-linear relationship between residuals and an outcome variable
# (Yhat) and the pattern could show up it in this plot if the model doesn't capture 
# the non-linear relationship between Y and Xs.
# If you find equally spread residuals around a horizontal line without distinct 
# patterns then is a good indication that we don't have to deal with non-linear 
# relationships, or that we have captured it with our fitted model if it was present.

# 2) Normal Q-Q plot: std.Residuals(y.axis) vs theoretical quantiles from the normal
# distributions(x.axis) (for check the Normality assumption):
# This plot shows if residuals are normally distributed.
# It's good if residuals are lined well on the straight dashed line.

# 3) Scale location plot: sqrt(std.Residuals)(y.axis) vs fitted values(x.axis)
# (for checking the Omosckedasticity assumption):
# This plot shows if the residuals are spread equally along the range of predictors.
# It's good if you see a horizontal line with equally (randomly) spread points.

# 4) std.Residuals(y.axis) vs Leverage(x.axis) (for checking the presence of 
# INFLUENCE-OUTLIERS):
# NB the leverage is i-th diagonal element hii derived from the mahalanobis distance.
# (see robust slides).
# This plot helps us to find influential cases.
# Not all outliers are influential in linear regression analysis, they might not
# be influential to determine our regression line i.e., our model.
# On the other hand, they could be extreme cases against a regression line and can 
# alter the results if we exclude them from the analysis.
# We watch out for influential outliers at the upper right corner or at the lower 
# right corner. Those spots are the places where cases can be influential against
# a regression line(our model).
# Specifically we look for cases outside the dashed line i.e the Cook's distance.
# When cases are outside of the Cook's distance(meaning that they have an high
# Cook's distance scores), they are influential with respect our regression results
# i.e., those observations are influential outliers.
# NB is more reliable use diagnosis plot from a ROBUST REGRESSION that are explained
# after due to the "Masking" issue.

# FROM HERE WE CAN SAY --> PROBABLY A BAD LEVERAGE POINTS but is more reliable also
# see the robust diagnostic plots (due to masking)

# NB removing an influence-outlier can have a big impact on the slopes(Betas) and
# R^2 of the model in classical linear regression lm, while in the robust regression 
# lmrob they change less (robustness property).
# NB Thus adding or removing outliers and "noise" changes lm much more than lmrob.

# Influence outliers can be simple outliers(with respect to the Y scale) or Leverage
# points(with respect to the X scale).
# Notice that the Leverage points are points that do not violate the model assumptions.
# They could be: Good leverage points or Bad leverage points.
# DIFF BTW GOOD AND BAD LEVERAGE POINTS and OUTLIERS --> see word doc. 

# FOR EXAM COMMENTS:
# The four plots shows potential problematic cases indicate through row numbers
# of the observation. 
# If some cases are identified across all four plots, you might want to take a 
# close look at them individually. 
# Is there anything special for the subject? 
# Or could it be a simply errors in data entry?
# So, what does having patterns in residuals mean to your research? 
# It's not just a go-or-stop sign. It tells you about your model and data. 
# Your current model might not be the best way to understand your data if there's 
# so much good stuff left in the data i.e., non captured informations as non-linearity
# and so on...
# In that case, you may want to go back to your theory and hypotheses.
# Is it really a linear relationship between the predictors and the outcome? 
# You may want to include a quadratic term, for example. A log transformation 
# may better represent the phenomena that you'd like to model.
# Or, is there any important variable that you left out from your model? 
# Other variables that we didn't include (e.g., age or gender) may play an important
# role in your model and data. 
# Or, maybe, our data were systematically biased when collecting and we may want
# to redesign the data collection methods.

# EXAMPLES 

# Add a regression outlier(with respect the scale of Y) --> c(-6,0.05,-0.1,0.9) 
# (4 values i.e y, x1, x2, x3):

regdata2 <- rbind(regdata1,c(-6,0.05,-0.1,0.9))
pairs(regdata2)

lmrd2 <- lm(y~x1+x2+x3,data=regdata2)
summary(lmrd2) # Different summary statistics: lower R^2 and different regression
# coefficients with respect to lmrd1.

par(mfrow=c(2,2))
plot(lmrd2)

# Add a Bad leverage point(with respect the scale of X) --> c(-4,20,-0.1,0.9) 
 
regdata3 <- rbind(regdata1,c(-4,20,-0.1,0.9))
pairs(regdata3)

lmrd3 <- lm(y~x1+x2+x3,data=regdata3)
summary(lmrd3)

par(mfrow=c(2,2))
plot(lmrd3)

# Add a Good leverage point(with respect the scale of X) --> c(22,20,-0.1,0.9)

regdata4 <- rbind(regdata1,c(22,20,-0.1,0.9))
pairs(regdata4)

lmrd4 <- lm(y~x1+x2+x3,data=regdata4)
summary(lmrd4)

par(mfrow=c(2,2))
plot(lmrd4)

# Example for "masking":

data(starsCYG,package="robustbase")
light.lm<-lm(log.light~log.Te, starsCYG)
plot(starsCYG)
abline(light.lm$coef[1],light.lm$coef[2])
par(mfrow=c(2,2))
plot(light.lm)
par(mfrow=c(1,1))
plot(light.lm, which = 5)

# Masking means that several outliers can "mask" each other, i.e., outliers can stop
# each other from appearing as outliers in standard diagnostic methods such as Cook's
# statistic or leverage/Mahalanobis distance.
# Standard diagnostic methods for checking a model assumption normally require all 
# other model assumptions to be fulfilled.
# Cook's distance can find an outlier assuming it's the only one(thus masking will
# affect this).
# If the relation between x and Y is nonlinear or heteroscedastic, the normal
# QQ-plot is worthless because it aggregates residuals over different x where in 
# fact they have different distributions.

# KEY IMPORTAND
# NB MASKING is the KEY PROBLEM with outliers in classical regression methods (but
# also in all the other non robust method of identification --> see MCD), thus, 
# as i said before, is better have a look to diagnosis plots from a robust regression.

# ROBUST REGRESSION ESTIMATORS 

# Standard robust regression still assumes i.i.d. and linearity, but under 
# heteroscedasticity and non-linearity we can sometimes fit the majority of the data 
# well by a homoscedastic linear model.

# IMPORTANT
# Diagnostics from robust regression are generally more reliable than traditional 
# diagnostics since they can deal with the masking issue.
# Traditional diagnostics are not that bad ,as robust diagnostics, they will show 
# many problems but with the classical regression these problems need to be dealt 
# with. whereas robust estimation deals with outliers(i.e., also the masking issues)
# /non-normality automatically, thus make sense also look to diagnosis with robust 
# regression estimators.

# MM - estimator

# FOR EXAM COMMENT 
# The basic idea underlying the robust linear model (lmrob) is that some 
# (1-alpha > 0.5) of the data is distributed conditionally normal and the remaining 
# fraction (alpha) comes from some arbitrary distribution (i.e., the outliers). 

# The goal of a robust method is to estimate the parameters (beta and sigma^2) of 
# this conditional normal distribution without giving to the outliers too much 
# influence. 
# If the bulk of the data (i.e, the good data) is not distributed conditionally
# normal then a linear model is not appropriate regardless of whether it is fit 
# robustly or not. Of course you can still use all of the standard linear modeling 
# tricks. For instance a log transformation of the response sometimes helps with 
# heteroskedasticity.

help("lmrob") # Computes an MM-type regression estimator.
# It uses a bi-square redescending score function, and it returns a highly robust 
# and highly efficient estimator (with 50% breakdown point and 95% asymptotic 
# efficiency for normal errors)

# REC : Add a Bad leverage point --> c(-4,20,-0.1,0.9)
# regdata3 <- rbind(regdata1,c(-4,20,-0.1,0.9))
lm3mm <- lmrob(y~x1+x2+x3,data=regdata3) 
par(mfrow=c(2,3))
help("plot.lmrob")
plot(lm3mm)
plot(1:101,lm3mm$rweights) # lm3mm$rweights --> wi (from slides) additional 
# diagnostic plot (wi downweight high observation values while the setted wi=0 
# identify outliers)

par(mfrow=c(1,1))
plot(lm3mm, which = 1) # RUN IT IN ORDER TO UNDERTSAND THE INTERPRETATION
# 1) Standardised residuals vs. robust distances : The plot divides the value range 
# in four regions marked by lines and coloured points: 
# (1) Regular --> cloud of black points (bulk of the data)
# (2) Outlier --> upper left (empty) over/under (1) --> identify outliers
# (3) Leverage --> front of (1) long rectangular part --> identify GOOD leverage
# points(empty)
# (4) Outlier and Leverage --> upper and lower of (3)--> identify BAD leverage points.

# NB Points in the fourth category are of major importance as they are outliers
# points that influence the model fit if they are removed from the estimation.

plot(lm3mm, which = 2) 
# 2) Normal QQ -> is a normal QQ plot (same interpretation)

plot(lm3mm, which = 3) 
# 3) Response vs Fitted values --> is a cross validation plot that compare y(response)
# with yhat(fitted values) 
# NB robust estimators will not fit all points well, or make all look like normal!
# It will fit the outliers even worse than LS(least square i.e, classical regression),
# because it weights them out for fitting.
# This can be a problem for cross-validation.
# obs 101 don't match the cross -validation due to the high weight assigned in 
# robust regression.

plot(lm3mm, which = 4) 
# 4) Residuals vs fitted values --> check the linearity as before (same interpretation)
# MOREOVER...
# Residuals vs fitted -> The plot is used to detect non-linearity, unequal error 
# variances, and outliers. 

# We plot it with the robustness weights in order to see which and where the 
# observation are treated as outliers by the estimator i.e., where the estimator
# set the robustness weights = 0.

plot(lm3mm, which = 5) 
# 5) Sqrt of abs Residuals vs Fitted values (scale-Location plot) --> check the
# homoskedasticity as before (same interpretation)

# REC
# Robust diagnostics and robust estimation deals with outliers/non-normality
# automatically i.e., identify that they are not violated in general and the influence
# of outliers by re-weighting them.

# The initial estimator (default, equivalent to init = "S",) uses as initial
# estimator an S-estimator. The method argument takes a string that specifies the 
# estimates to be calculated as a chain i.e default method = "MM" or "SM" -->
# S init - estimator of the scale of the M - estimator used after.
# Takes the BP from the S init - estimator used for it (rho1 = u/c1) and the 
# efficiency from (rho2  = u/c2).
# IDEA : rho2 governs data information, but rho1 S-scale governs outlier tolerance.
# Default: c1 = 1.548 yields BP~1/2, c2 = 4.685 yield 95% efficiency.

a <- summary(lmrd1) # Regression with traditional estimators (LS)
b <- summary(lm3mm) # Regression with ROBUST ESTIMATORS (MM)

# NB
# The output from robust regression is also very similar to those from linear regression.
# The only exception now is that that the observations are weighted based on how 
# deviant they are, and also robust standard errors will be computed for model 
# estimates.

a$coefficients
b$coefficients

a$r.squared # 0.9201405
b$r.squared # 0.9173593

which(lm3mm$rweights == 0) # obs 101 is an outlier
regdata3[101,]# Bad leverage point --> c(-4,20,-0.1,0.9)

# FOR EXAM COMMENTS
# NB : A robust estimator will not fit all points well, or make all look like normal!
# It will fit the outliers even worse than LS, because it weights them out for fitting.
# This can be a problem for cross-validation; squared loss is strongly dominated by
# biggest misfits. Squared loss is itself non-robust!
# Better run with L1/absolute value or other robust loss when attempting to fit the
# "good majority" well.

# EXAMPLES 
# REC : regdata1 <- as.data.frame(cbind(y,x1,x2,x3)) --> Data with no outliers/
# influencial cases
par(ask=F)
lm1mm <- lmrob(y~x1+x2+x3,data=regdata1) 
par(mfrow=c(2,3))
plot(lm1mm)
plot(1:100,lm1mm$rweights)
par(mfrow=c(1,1))
plot(lm1mm, which = 1) # LOOK THE DIFFERENCE IDENTIFICATIONS AMONG OUTLIERS AND
# LEVERAGE POINTS IN THIS PLOT
summary(lm1mm)
which(lm1mm$rweights == 0)

# REC : Add a regression outlier(with respect the scale of Y) --> c(-6,0.05,-0.1,0.9) 
# regdata2 <- rbind(regdata1,c(-6,0.05,-0.1,0.9))
lm2mm <- lmrob(y~x1+x2+x3,data=regdata2) 
par(mfrow=c(2,3))
plot(lm2mm)
plot(1:101,lm2mm$rweights)
par(mfrow=c(1,1))
plot(lm2mm, which = 1)# LOOK THE DIFFERENCE IDENTIFICATIONS AMONG OUTLIERS AND
# LEVERAGE POINTS IN THIS PLOT
summary(lm2mm)
which(lm2mm$rweights == 0)

# REC : Add a Bad leverage point(with respect the scale of x) --> c(-4,20,-0.1,0.9) 
# regdata3 <- rbind(regdata1,c(-4,20,-0.1,0.9))
lm3mm <- lmrob(y~x1+x2+x3,data=regdata3) 
par(mfrow=c(2,3))
plot(lm3mm)
plot(1:101,lm2mm$rweights)
par(mfrow=c(1,1))
plot(lm3mm, which = 1)# LOOK THE DIFFERENCE IDENTIFICATIONS AMONG OUTLIERS AND
# LEVERAGE POINTS IN THIS PLOT
summary(lm3mm)
which(lm3mm$rweights == 0)

# REC : Add a Good leverage point(with respect the scale of x) --> c(22,20,-0.1,0.9)
# regdata4 <- rbind(regdata1,c(22,20,-0.1,0.9))
lm4mm <- lmrob(y~x1+x2+x3,data=regdata4)
par(mfrow=c(2,3))
plot(lm4mm)
par(mfrow=c(1,1))
plot(lm4mm, which = 1)# LOOK THE DIFFERENCE IDENTIFICATIONS AMONG OUTLIERS AND
# LEVERAGE POINTS IN THIS PLOT
plot(1:101,lm4mm$rweights)
summary(lm4mm)
which(lm4mm$rweights == 0)  # 0 

# REC : masking example
# MM -estimator
# 4 observations c(11,20,30,34) are outliers with |weight| = 0
light.mm<-lmrob(log.light~log.Te, starsCYG)
summary(light.mm)
par(mfrow=c(2,3))
plot(light.mm)
par(mfrow=c(1,2))
plot(light.mm, which = 1, main = "MASKING WITH LMROB") # MASKING WITH LMROB
plot(light.lm, which = 5,main = "MASKING WITH LM") # MASKING WITH LM
# NB We can see that the typical behaviour of the model with masking presence is
# an UP AND DOWN FITTED LINE!
# NB Lmrob is able to detect the masking!
plot(1:47,light.mm$rweights,xlab="Obs. number",ylab="Robustness weights")
which(light.mm$rweights == 0) # --> # 4 observations c(11,20,30,34) are outliers 
# with |weight| = 0

# IMPORTANT
# NB To sum up, the MM - ROBUST REGRESSION differs from the classical basically 
# since he reweight potential influencial observation in order to obtain a good 
# fitting. He takes out of fitting by putting those weight equal to 0 only OUTLIERS
# and BAD LEVERAGE POINT. 
# Moreover it can recognize and deal with the issue of MASKING(always by playing
# with the weight of the responsible observations).
# IN this way he follows the BASIC IDEA behind the ROBUST THEORY:
# The basic idea underlying the robust theory is that we assume that some of the 
# data is distributed conditionally normal and the remaining while the remaining 
# fraction comes from some arbitrary distribution that are the outliers. 
# The goal of a robust method is to estimate the parameters of this conditional
# normal distribution without giving to the outliers too much influence. 

# S-estimator
light.s<-lmrob(log.light~log.Te, method="S",starsCYG)
par(mfrow=c(2,3))
plot(light.s)
par(mfrow=c(1,2))
plot(light.mm, which = 1, main = "LMROB - MM ESTIMATOR") 
plot(light.s, which = 1, main = "LMROB - S ESTIMATOR") 
# NB As we can see the S estimator is less tollerant with respect to the outliers 
# but also less efficient, since by excluding and down-weighting to much observations
# we lose information i.e, efficiency.
# Moreover this is also because it has the role of govern the robustness in the 
# MM estimator.
par(mfrow=c(1,1))
plot(1:47,light.s$rweights,xlab="Obs. number",ylab="Robustness weights")
which(light.s$rweights == 0)
# 7 observations c(7,9,11,18,20,30,34) are outliers with |weight| = 0
# (in the MM only 4 are treated as outliers/ bad leverage points)
# S - estimator is less tolerate against outliers (but also less efficient).

# NB
# We can derive the S - estimator directly from the lmrob object since it used
# for its initialization.
light.mm$init.S 
par(mfrow=c(1,2))

plot(1:47,light.mm$init.S$rweights,xlab="Observation number",ylab="Weight",
     main="S-estimator, robustness weights",type="n")
text(1:47,light.mm$init.S$rweights,rownames(starsCYG),cex=0.7)

plot(light.mm$init.S$fitted,light.mm$init.S$residuals,xlab="Fitted values",
     ylab="Residuals",main="S-estimator, residuals vs. fitted")
# Residuals vs fitted -> The plot is used to detect non-linearity, unequal error 
# variances, and outliers. 
# We plot it with the robustness weights in order to see which and where the 
# observation are treated as outliers by the estimator i.e., where the estimator
# set the robustness weights = 0.

# NOTICE THAT WE CAN'T OBTAIN DIAGNOSTIC AND SO ON SINCE WE ARE NOT USING THE 
# ESTIMATOR FOR FITTING A ROBUST REGRESSION AS BEFORE I.E.,
# lmrob(log.light~log.Te, method="S",starsCYG)
# WE ARE JUST FITTING THE SCALE FOR THE M ESTIMATOR USED AFTER IN THE MM ROBUST
# RERESSION FITTING I.E., light.mm<-lmrob(log.light~log.Te, starsCYG), WHERE 
# METHOD IS DEFAULT method="MM".

# COMPARING ROBUST ESTIMATORS (COMMENTS FOR EXAM)

# --> WHEN YOU MAKE COMMENTS FIRST HAVE A LOOK AT THE SUMMARY TABLE IN THE SLIDES!!!

# IN ORDER:

# M-ESTIMATOR(HUBER m) --> M-estimation method is a generalization of the maximum 
# likelihood estimation(i.e., the ML is a special case of M-estimators) in the
# context of location models(i.e., estimation of the mean). 

# The introduction was motivated by the robust statistics since now we are try to 
# find a more general t estimator able to minimize p, that is a function with 
# certain properties for each observations yi --> thus we are try to find a t that
# minimize the sum of the p for i = 1 to n.
# --> i.e, Tn = argmin(t){sum(i=1 to n)[p(yi,t)]}

# NB p measures the fit loss of estimating yi by t -> thus we are trying to find
# the best t estimator basing our search on the minimization of a sort of
# "residual function" p(yi,t).

# NB this is a generalization since p can be chosen among different functions i.e.,
# The function p, or its derivative, psi, can be chosen in such a way to provide
# the estimator desirable properties (in terms of bias and efficiency) when the 
# data are truly from the assumed distribution i.e, with no outliers.
# (RECALL THE BASIC IDEA OF ROBUST THEORY --> split of the data Y*).

# NB for e.g the traditional LS or ML-estimator is given by setting p = - log f(yi,t).

# Now recall that we are in the field of location estimation where both the mean
# and the median are estimators for it.

# Basically the M-estimator behave like the mean when the proportion of outlier
# influence is under a certain threshold(in the slide |u| = |yi - t| <= gamma) 
# exploiting the efficiency properties instead of robustness of the mean, viceversa
# when |u| > gamma it behaves like the median exploiting the robustness properties
# instead of an high efficiency(as the mean) of the median.

# NB Threshold --> since we use those estimator when we suspect that we have outliers
# in the data thus the difference u = yi - t is due to outliers! This is because 
# we call gamma "outlier influence threshold".

# Also, this is because we have introducing a generalization of the maximum likelihood 
# estimation(more flexible), in order to "switch" the functions.

# The M estimator is nearly as efficient as OLS. 
# The weakness of M-estimation method is the lack of consideration on the data
# distribution since is not a function of the overall data because only using the 
# median as the weighted value for removing or lower the outliers influence.

# NB In most practical cases, the M-estimators are of psi-type where the psi
# function is the derivative of the p function which sobstitute the p one.
# Because often is simpler to differentiate with respect to Tn and solve for the 
# root of the derivative. When this differentiation is possible, the M-estimator
# is said to be of psi-type. Otherwise, the M-estimator is said to be of p-type.
#(LOOK SLIDE FOR MORE DEEPER EXPLANATION)

# Exist also an M-estimator for the scale that follows the same principle stated
# before for the location.

# S - estimator -->  the S - estimator minimize an M-estimator of scale since it 
# is based on residual scale of M-estimation methodology.
# It has an efficiency of about the 25% if p is chosen to give BP~1/2.
# Have the same asymptotic properties corresponding to M-estimators but it can
# also handle up to the 50% of the outliers appearing in the data. 
# This method uses the residual standard deviation to overcome the weaknesses of
# median used in M one.
# S-estimators are more robustly than the M-estimator, because S-estimators have
# smaller asymptotic bias and smaller asymptotic variance in the case contaminated
# data.

# MM - estimator -->  MM estimation method is a combination of high breakdown value 
# estimation method and efficient estimation method, which was the first estimation 
# with an high breakdown point and high efficiency under normal error.
# The key difference with the Huber estimator lies in the definition of the weight
# function now based on a redescending score from the S - estimator.

# The main IDEA is that we combine an M-regression estimator with an S estimator 
# for the scale (thus we have an S initialization), where the efficiency comes
# from the M-regression estimator and the robustness from the S estimator for the
# scale.
# MM-estimator has the same BP as S-estimator used for it.

# COMPARISONS + COMMENTS  (LOOK AT THE SUMMARY TABLE)

# MM is the more efficiency compared with S method and it is better in outlier
# generating error distribution. 

# The M.Huber-estimator is not robust with respect to high leverage points, so it
# should be used in situations where high leverage points do not occur.

# SUMMARYZE --> MM best trade-off btw Asymptotic efficiency and  robustness
# (S too robust less efficient, M  to efficient less robust)

# M doesn't put =0  bad leverage points (only outliers)
# S put = 0 with psi (look slide) bad leverage points and ouliers but it has 
# 25% asy efficient

# MM often best

# Ex 9.4 + Ex 4 mock exam (OTHER PART IN ANOTHER SCRIPT)

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
unicef <- read.table("unicef97.dat",header=TRUE)

# LS - estimator
uniceflm <- lm(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
               Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid, data=unicef)
# Summary
summary(uniceflm)
# We have 4 significant reg: Literacy.Fem, Drinking.Water, Polio.Vacc and Urban.Pop.
# Multiple R-squared:  0.7587 --> statistical measure of how close the data are to 
# the fitted regression line

# Residuals vs fitted (NB: y vs x/ std residuals = residuals)
plot(uniceflm, which = 1) # LINEARITY --> Relation still well approximated by linearity.
# Residuals vs Normal QQ(normal theoretical quantiles)
plot(uniceflm, which = 2) # NORMALITY --> Overall the residuals can be approximated with a 
# N distribution since most of the points are aligned. We have issues with SaoTP, 
# Angol and Niger observations. They will have to be specifically analyzed later.
# Scale-Location ( sqrt std residuals vs fitted values)
plot(uniceflm, which = 3) # HOMOSKEDASTICITY --> seems that is present an issue 
# with HOMOSKEDASTICITY, however weak heteroscedasticity can often be tolerated.
# Residuals vs leverage
plot(uniceflm, which = 5) # FIND INFLUENTIAL CASES --> Not all outliers and 
# leverage points are influential in linear regression analysis. 
# For example Niger is an outlier since shows an high residual values but he doesn't show
# a significant impact. we can identify two points, far to right o as good leverage points with respect that not influnce, we look for cases outside
# the dashed line (Cook's distance) that can alter the results if we exclude them from
# analysis.
# This plot helps us to find influential cases.
# Not all outliers are influential in linear regression analysis, they might not
# be influential to determine a regression line.
# On the other hand, they could be extreme cases against a regression line and 
# can alter the results if we exclude them from analysis.
# We watch out for outlying values at the upper right corner or at the lower right 
# corner. Those spots are the places where cases can be influential against a 
# regression line, we look for cases outside the dashed line i.e the Cook's distance.
# When cases are outside of the Cook's distance (meaning they have high Cook's distance 
# scores), the cases are influential to the regression results.


# MM - estimator
unicefmm <- lmrob(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                    Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid,
                  data=unicef)
# Summary
summary(unicefmm)
# Std residuals vs Robust distance

plot(unicefmm,which=1)
# Residuals vs Normal QQ 
plot(unicefmm,which=2)
# Robustness weights
plot(1:121,unicefmm$rweights,xlab="Observation number",ylab="Weight",
     main="MM-estimator, robustness weights",type="n")
text(1:121,unicefmm$rweights,rownames(unicef),cex=0.7)
# Residuals vs. fitted
plot(unicefmm$fitted,unicefmm$residuals,xlab="Fitted values",ylab="Residuals",
     main="MM-estimator, residuals vs. fitted")

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



