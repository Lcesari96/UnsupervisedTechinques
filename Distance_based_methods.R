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

# DISTANCE BASED METHODS
# Partitioning -> Kmeans, PAM
# Hierarchical(agglomerative) -> single, average, complete, Ward

set.seed() <- # Remember
setwd() <- # Remember

round(apply(data[,"variablese of interest e.g 3:12"],2,sd),0) # if variances are 
# too diff, is better to standardized the variables in order to avoid to obtaining 
# clusters that are dominated by variables having the largest amount of variation. 

data <- # Dataset or scaled dataset: scale(data)

plot(data) # Understand the main features of the data (pairs plot since usually
           # multivariate ndim > 2), to go deeper below look the MDS with the 
           # chosen distances (ndim = 2)
  
################################################################################

# Distances (are dissimilarity, thus greater are those values for pairs of 
# observations and more they are dissimilar)
  
# Binary/categorical
jac <- dist(data, method = "binary") # Only joint presences are meaningful.

# dL1/p is dSM for binary variables
p <- ncol(data) # n of var
sm <- dist(data,method="manhattan")/p # Both joint presences and absences are 
                                      # meaningful
# general dSM for cat var
# sm in R-package nomclust computes dSM
sm<-sm(data)

# Continuos 
eu<-dist(data,method="euclidean") # Reference distance

man<-dist(data,method="manhattan") #I decide to use the  Manhattan L1 distance since
# is preferred over the Euclidean as the dimension of the data increases. 
# If we consider the extreme case where we use the l(inf) Minkowski distance 
# (exponent = infinity), then the distance is the highest difference between any 
# two dimensions of our vectors.
# This doesn't make sense when we have many dimensions as we would be ignoring most
# of those, ending by measuring distances basing on a single attribute(the most far).
# This example is for explain that reducing the exponent makes all the dimensions
# play a bigger role in the distance calculation, this is because the L1 distance 
# is preferable when we deal with high_dimensional data sets.
# To summarize lower is the exponent and less relevant an high difference in some
# given dimensions will be 
# In practice, the euclidean distance is preferable but when it doesn't make sense 
# in our data we use the L1 one -->
# Since there is no evidence against the usage of the Euclidean distance for this 
# data (such as high dimensions or categorical variables)

# NB for MD(mahalanobis distance) use data NON SCALED
# "squared" è sottointeso --> Squared MD
mahalm <- matrix(0,ncol=nrow(data),nrow=nrow(data))
cov <- cov(data)
for (i in 1:nrow(data))
  mahalm[i,] <- mahalanobis(data,as.numeric(data[i,]),cov)
mal<-as.dist(mahalm)
# For uncorrelated variables, the EU(euclidean) is equal to the MD.
# In a normal Euclidean space the variables are represented by axes with right
# angles i.e 90 degrees to each other.
# If two or more variables are correlated, the axes are no longer at right angles,
# and the measurements become impossible with a simple ruler. Moreover if we have
# p > 3 var, we can't plot them in regular 3 dimensional space precisely due to 
# the axes issue.
# The MD solves this measurement problem, since measures the distances btw the obs
# relative to a special "centroid" i.e a central point in the multivariate space 
# where all the means from all the variables intersect, a sort of overall mean for
# multivariate data.
# The larger is the MD the further away from this "centroid" the obs is.
# This is because MD is often related to outlier detection.
# NB for the Mahalanobis distance it doesn't make a difference whether the data
# set is scaled or not.

# Mixed
setwd("~/Henning/Dataset")
data<- read.table("Boston.dat",header=TRUE)
gow <- daisy(data, metric="gower", type=list(asymm=4))
# Gower dissimilarity is used when we have mixed data i.e both quantitative and
# categorical dimensions.
# The distances to use are specified in the attribute type:
# symm -> SM
# asymm -> Jaccard
# unspecifyed ->  L1 manhattan distance
# e.g type=list(asymm=c(2,4),symm=c(3,6)) if vars 2,3,4,6 are binary.
# It returns a dist object.

# NB if we have categorical variables(i.e NON-BINARY)
# type=list(asymm=c(2,4),symm=c(3,6)) doesn't work.
# We need to specify apriori the variables as factor or ordered variables.
# e.g
data$variable <- as.ordered(data$variable)
data$variable <- as.factor(data$variable)
# NB this is a not a distance if NA are present since it doesn't fulfill the 
# triangle inequality.

# CORRELATION DISSIMILARITY (they are not distances)
# Here we treat the variables as observation, thus with this kind of distance we 
# are trying to cluster the variables.

setwd("C:/Users/loren/Desktop/Henning_Exam/Data")
covid2021 <- read.table("covid2021.dat")
covid2021cl <- as.matrix(covid2021[,5:559])
# This selects the variables for clustering
# read.table produces data.frame but matrix is needed here.

corcovid <- cor(t(covid2021cl))
# Need to transpose data matrix because cor is asked to treat
# the observations as variables.
cordistcovid <- 0.5-corcovid/2 # Correlation dissimilarity
cordistcovid <- as.dist(cordistcovid)
# For using R-functions that expect a dist object.
# We get 0 for r = 1 i.e no dissimilarity if there is maximum correlation
# and max 1 for r = ???1 i.e maximum dissimilarity if there is maximum negative 
# correlation 

# A 2nd way (identical in terms of aims)
# cordistcovid2 <- 1 - abs(corcovid)
# cordistcovid <- as.dist(cordistcovid2)
# has largest dissimilarity for r = 0 and small if abs(r) large.

# NB: Correlation measures the degree of linear dependence, i.e., how well one 
# variable can be expressed as linear function of the other.
# This could be of interest also for observations... Do they follow the same 
# "relative pattern" over the different variables?
# Can use correlation dissimilarity between observations, computed as if they
# were variables.

mds <- mds("chosen distance", ndim=2)
plot(mds$conf) # conf: Matrix of fitted configurations i.e preserved eu distances 
# after the dimensional reduction. It helps in understand the main features of the
# data. (NB it preserves the distances as Euclidean one)
mds$stress # Amount of info. lost with the dim-reduction (bounded 0-1) (lower it
# is,lower loss we had)

################################################################################

# Run distance-based clustering methods (Hierarchicals(single, average, complete
# ward) + K - medoids (PAM)) + Explanation of the "choosing" criterion (ASW)
# average silhouette width 

# IDEA: Compute the ASW(average silhouette width) of each partition obtained
# from the chosen clustering methods with different number of cluster k,
# then chose the best according to the maximum ASW.

# NB different ASW values can be compared ONLY IF they are calculated through the same 
# distance.
# NB FIRST THINGS TO DO --> LOOK THE DENDROGRAMS!! (If the data are not so 
# meaningful this can help us a lot in order to make a choiche)

# ASW is the mean of the s(i)(silhouette coeff for obs i) 

# s(i) tell us how well each observation fit into a cluster for a given partition
# s(i) = [b(i) - a(i)] / [max{a(i), b(i)}] 
# -1 < s(i) < 1 
# a(i) = mean distance between xi and all other data points in the same cluster.
# Measure of how well xi fit to its cluster.(smaller it is better is the fit)
# b(i) = smallest mean distance between xi and data points in another cluster.
# The cluster from which xi shows the smallest b(i) is the "neighboring cluster"
# (is the next best fit cluster for xi) (higher it is, better is the fit).
# An s(i) close to 1 means that xi fits well to his cluster. 
# An s(i) close to -1, means that xi should fit in its "neighbouring cluster". 
# An s(i) near to zero means that xi is on the border of two clusters.

# LOOP INIT

tsingle <- tave <- tcom <- tward <- tpam <- list()
# These lists will collect all the clustering labels for the observation
# tsingle = cluster labels from hierachical method hclust with single linkage
# ...

nc <- 2:15 # Number of k clusters that we want to try

ssil <- asil <- csil <- wsil <- psil <- list()
# These lists will collect the s(i) values for each observation
# ssil =  s(i)s obtained with hclust + single linkage 
# ...

sasw <- aasw <- casw <- wasw <- pasw <- NA
# Vectors of ASW values: 1 for each clustering method with length equal to k number
# of partitions that we try.

# Euclidean e.g (for other distance we just need to replace "eu")

# I run hierarchical methods first to produce the dendrograms, so that in the loop
# only cutree has to be applied.
single <- hclust(taidist,method="single")
average <- hclust(taidist,method="average")
complete <- hclust(taidist,method="complete")
ward <- hclust(taidist,method="ward.D2")

# Run all distance-based clustering methods with different k partitions
for(k in nc){
  print(k)
  tsingle[[k]] <- cutree(single,k)
  tave[[k]] <- cutree(average,k)
  tcom[[k]] <- cutree(complete,k)
  tward[[k]] <- cutree(ward,k)
  tpam[[k]] <- pam(eu,k)
  
  ssil[[k]] <- silhouette(tsingle[[k]],dist=eu)
  asil[[k]] <- silhouette(tave[[k]],dist=eu)
  csil[[k]] <- silhouette(tcom[[k]],dist=eu)
  wsil[[k]] <- silhouette(tward[[k]],dist=eu)
  psil[[k]] <- silhouette(tpam[[k]],dist=eu)
  
  sasw[k] <- summary(ssil[[k]],dist=eu)$avg.width
  aasw[k] <- summary(asil[[k]],dist=eu)$avg.width
  casw[k] <- summary(csil[[k]],dist=eu)$avg.width
  wasw[k] <- summary(wsil[[k]],dist=eu)$avg.width
  pasw[k] <- summary(psil[[k]],dist=eu)$avg.width
}

################################################################################

# Extraction and interpretations of the results

# FIRST THINGS TO DO --> LOOK THE DENDROGRAMS!!
plot(single)
plot(average)
plot(complete)
plot(ward)

# Extracting information examples
table(tsingle[[3]]) # hclust + single clustering labels for k = 3 partition 
# (cutree = 3)
# ...

table(tpam[[5]]$clustering) # clustering labels by using PAM for k = 5
# ... 

ssil[[3]] # s(i)s hclust + single for k = 3 (sil object)
# ...

psil[[5]] # s(i)s using PAM for k = 5
# ...

sasw[3] # ASW hclust + single for k = 3
# ...

pasw[5] # ASW PAM for k = 5
# ...

# Ko optimal number of clusters for each methods according to the ASW 
which.max(sasw) 
which.max(aasw) 
which.max(casw)
which.max(wasw)
which.max(pasw) 

# Correspondent ASW val
max(sasw,na.rm=TRUE) 
max(aasw,na.rm=TRUE) 
max(casw,na.rm=TRUE) 
max(wasw,na.rm=TRUE) 
max(pasw,na.rm=TRUE)

# Summary plot of all the ASW curves 
plot(1:15,sasw,ylim=c(0,0.5),type="l",xlab="Number of clusters",ylab="ASW")
points(1:15,aasw,ylim=c(0,0.5),type="l",col=2,lty=2)
points(1:15,casw,ylim=c(0,0.5),type="l",col=3,lty=3)
points(1:15,wasw,ylim=c(0,0.5),type="l",col=4,lty=4)
points(1:15,pasw,ylim=c(0,0.5),type="l",col=5,lty=5)
legend(8,0.2,legend=c("single","average","complete","ward","pam"),lty=1:5,col=1:5)

# If it present is better to look at local optimums for a k > 2 -->
# In some applications there's tendency of ASW to favour K = 2; local optima for 
# K > 2 may also be worth exploring.

################################################################################

# Clustering visualization

# NB FOR CLUSTERING VISUALIZATION + INTERPRETATION IS BETTER TO SAVE THE PLOT
# IN WORD AND THEN MAKE THE COMPARISON (MDS, silhouette and dendrograms)

# MDS -> The MDS is a dimensional reduction technique used to display the information 
# contained in a distance matrix. It Takes as input a distance matrix between each 
# pairs of observations and a chosen number of dimensions, then an MDS algorithm 
# places each object into a lower dimensional space such that the distances are
# preserved as well as possible.
# In this way we can visualize our clusters since the distances represents a measures
# of similarities.

# From here we can make an idea about the number of k clusters.

# NB the MDS preserve the distances as euclidean one.

mds1 <- mds(jac, ndim=2)
plot(mds1$conf) # conf: Matrix of fitted configurations i.e preserved eu distances 
# after dimension reduction
mds1$stress # Amount of info. lost with the dim-reduction (bounded 0-1) (lower it
# is,lower loss we had)

mds2 <- mds(sm, ndim=2)
plot(mds2$conf)
mds2$stress

mds3 <- mds(eu, ndim=2)
plot(mds3$conf)
mds3$stress 

mds4 <- mds(man, ndim=2)
plot(mds4$conf)
mds4$stress

mds5 <- mds(mal, ndim=2)
plot(mds5$conf)
mds5$stress

# Visualizations of the best results

# NB BEFORE WRITING PLAY WITH THE SILHOUTTE PLOTS AND THE RESULTS!!!
# From here we can check the fit of the obs to the clusters and if some clusters
# are disconnected (negative s(i)s and "neighbourings") ...
# Maybe is better go higher with k in order to split clusters with negative s(i)s?
# Or go lower to join?

# NB FOR PLOT THE SILHOUTTER WE CAN USE DIRECTLY THE COMMAND PLOT BUT NOT ALWAYS
# WORKS ... HOWEVER, IN BOTH, THERE IS A CORRESPONDECE BTW SILHOUTTE PLOTS MDS 
# AND DENDROGRAMS!
# Try to identify it for the interpretations!

# MDS
# par(mfcol=c(,)) <- for comparisons
plot(mds1$conf,pch=clusym[tsingle[[6]]],col=tsingle[[6]],
     main="Single Linkage, 6 clusters")
# ...
# NB Notice that the presence of "spread" clusters can be related to the MDS 
# visualization, since we are reducing the dimensions can happen that some obs 
# seems "distant" but in reality (i.e if we go higher with the dimension) those
# are closer. This can be verified by looking to the correspondent silhouette 
# plot i.e if the s(i) are all high and positive for the "spread clusters".

# Silhouette plots
fviz_silhouette(ssil[[6]],main="Single Linkage, 6 clusters")
# ...
# From here we can check the fit of the obs to the clusters and if some clusters
# are disconnected (negative s(i)s and "neighbourings") ...
# Maybe is better go higher with k in order to split clusters with negative s(i)s?
# Or go lower to join?

# Note: Points with negative silhouette width aren't necessarily "misclassified".
# Putting them in neighbouring cluster will change silhouette values of other
# points and may lead to worse ASW.

# Dendrograms 
plot(single)
# ...
# Dendrograms help us to understand the structure of the hierarchical partitions.
# Often for single linkage the number of clusters would need to be vastly increased
# to get at some structure in the big set of observations, thus it will be excluded.

# IF THE DATA ARE NOT MEANINGFUL AND UNCLEAR (VERY SPREAD) HIERARCHIAL METHOD ARE 
# GOOD CHOICE!

# Common comment for reject a dendrogram:
# "the number of clusters would need to be vastly increased to get at some structure
# in the big set of observations"

################################################################################

# Example of clusters interpretation (see benchmark script)

plot(mdstetra$conf,pch=clusym[tave[[10]]],col=tave[[10]],
     main="Average Linkage, 10 clusters")
# Referring to the average linkage clustering with K = 10, 
plot(mdstetra$conf,pch=clusym[tave[[10]]],col=tave[[10]],
     main="Average Linkage, 10 clusters")
# Clusters 4, 5, and 9 are rather small, and a biologist would probably not feel
# sure about whether these constitute species. 
# Cluster 3 seems to be somewhat heterogeneous, and biologists may look into that and 
# maybe give explanations.
# Cluster 10 ("0") also has some more within-cluster variation with worse silhouette 
# widths than the others, 
fviz_silhouette(asil[[10]],main="Average Linkage, 10 clusters")
# but it looks homogeneous in the MDS plots, and like a distinctive pattern of the data
# set. 
# It is split into two by PAM and complete linkage, but this doesn't seem to improve
# matters.
# The other clusters seems to be well distinguished groups of observations that may well
# qualify as different bee species.

################################################################################

# K-means

# Gap -statistic
# Since we are dealing with unsupervised case in which the number of classifications,
# i.e the k clusters, are not known apriori we need to estimate it. 
# As the ASW comparisons adopted before, clusgap for kmeans is used to estimate 
# ko.

gap <- clusGap(data, FUN = kmeans, nstart = 100,K.max = 10,B = 100)
# The idea of the Gap Statistics is to choose the number of K, where the biggest
# jump(gap) in within-cluster distance occurred, based on the overall behavior of
# uniformly drawn samples.
# Basically those uniformly drawn samples are a reference simulated and averaged 
# distributions of our within-cluster distances for each k till K.max (a "linear/ 
# non-linear" decreasing sequence of points (1 for each k), since the within
# cluster distance is monotonically decreasing as the number of clusters k increase) 
# thus our gap(k)s are computed by subtracting the within cluster distance of our 
# data with the simulated one for k step.

# The samples drawn from a uniform distribution are the worst data in order to 
# obtain clusters (due to the kind of distribution) thus we doesn't expect straight 
# decreasing or increasing(i.e,local min or max) in within cluster distances as number
# of k increase.
# While, since we known that as k increases the within-cluster distance tends to
# decrease in any case, by computing the differences i.e, the gap(k)s and if they
# are significant, those are differences are due to a STEEP(ripido) DECREASING
# of the within cluster distance for a given k.
# That, if they exceed the thershold of the gap(k+1) - SE.factor * S'(K+1), inidicates
# a local minima in the within cluster distances of our data.

# ko is given by the k that minimize (i.e argim(k)) 
# gap(k) >= gap(k+1) - SE.factor * S'(K+1).
# SE.factor * S'(K+1) acts as a threshold i.e it removes small changes (the sampling 
# noise from the data). 
# Once we know ko we plug it into kmeans.
# B indicates the number of bootstrap simulated samples in order to obtain our 
# reference simulated and averaged distributions of our within-cluster distances.

plot(gap) # Look for local maxima for k (i.e the maximum difference)
gap # --> Number of clusters (method 'firstSEmax', SE.factor=1): ko (here is where
    # we can read the optimal n of clusters)

km <- kmeans(data, k0)

# MDS on the Kmeans object
plot(mds$conf,col=km$cluster,
     pch=clusym[km$cluster], main="kmeans with K = ko")

# APPLICATION EXAMPLE 

# Gap stat

data(oliveoil)
olive_gap <- clusGap(scale(oliveoil[,3:10]), FUN = kmeans, nstart = 100,
                     K.max = 10,B = 100)
plot(olive_gap) # 5
print(olive_gap,method="firstSEmax",SE.factor=2) # Here we specify the method and 
# the SE.factor(threshold) in order to make a choice

# NB BUT ALWAYS LOOK AT THE PLOT IN ORDER TO MAKE OUR OWN CHOICE!

# kmeans 
km <- kmeans(scale(oliveoil[,3:10]), 5)

eu <- dist(scale(oliveoil[,3:10]), method = "euclidean")

mds <- mds(eu, ndim=2)
plot(mds$conf) 
mds$stress 

plot(mds$conf,pch=clusym[km$cluster],col=km$cluster,
     main="km, 5 clusters")

sil <- silhouette(km$cluster,dist=eu)
fviz_silhouette(sil, main="km, 5 clusters")
summary(sil,dist=eu)$avg.width # ASW

# GAPNC FUNCTION do both clusgap and kmeans in one run.

require(cluster)

gapnc <- function(data,FUNcluster=kmeans,
                  K.max=10, B = 100, d.power = 2,
                  spaceH0 ="scaledPCA",
                  method ="globalSEmax", SE.factor = 2,...){
  # As in original clusGap function the ... arguments are passed on
  # to the clustering method FUNcluster (kmeans).
  # Run clusGap
  gap1 <- clusGap(data,kmeans,K.max, B, d.power,spaceH0,...)
  # Find optimal number of clusters; note that the method for
  # finding the optimum and the SE.factor q need to be specified here.
  nc <- maxSE(gap1$Tab[,3],gap1$Tab[,4],method, SE.factor)
  # Re-run kmeans with optimal nc.
  kmopt <- kmeans(data,nc,...)
  out <- list()
  out$gapout <- gap1
  out$nc <- nc
  out$kmopt <- kmopt
  out
}

a <- gapnc(data, K.max = 15, nstart = 100)
a$gapout # contains all the calls of clusgap 
a$nc # contains ko from clusgap
a$kmopt # contains the kmeans obj run with nc as number of cluster

print(a$gapout,method="globalSEmax",SE.factor=2) # k0
plot(a$gapout) # gap plot

# MDS on the Kmeans object
plot(mds$conf,col=a$kmopt$cluster,
     pch=clusym[a$kmopt$cluster], main="kmeans with K = nc/ko")
     
# We can also check with the silhouette the fit of the obs to the cluster produced
# by Kmeans 
sil <- silhouette(a$kmopt$cluster,dist=eu) # s(i)s
summary(sil,dist=eu)$avg.width # ASW

# NB KMEANS WORKS ONLY WITH THE EUCLIDEAN DISTANCE DUE TO A CONVERGENCE ISSUES.

################################################################################

# ARI (Adjusted rand index) -> measure of similarity/matching between two partitions.

# Applications:
# -> Compare clustering with external grouping (for interpretation and validation)
#    Could use "true" external grouping to compare clustering methods.
# -> Investigate effect of variables or outliers removed.
# -> Investigate stability of clustering by comparing clusterings on resampled 
#    data sets. (May want to cluster clusterings.)
# -> Doesn't require partition K1 =  partition K2 and cluster matching.
# -> Typical values depend on K1, K2, cluster sizes.

#e.g 
adjustedRandIndex(wpam[[2]]$clustering,wdbcdiag)
# adjustedRandIndex(clusteringVectorK1, clusteringVectorK2)

# Where clusteringVectorsK are objects that contains the cluster labels for each 
# obs as:
# [1] 1 1 1 1 1 1 1 1 1 1 2 1 1 1 

# ARI takes values between 0 and 1:
# -> 1 perfect agreement between two partitions.
# -> 0 if are two random partitions.

# ARI LOOP e.g for comparisons
arival <- list()
for (i in 1:8) {
  for (j in 1:8) {
    if(j> i){
      arival[[j]] <- adjustedRandIndex(ari[[i]],ari[[j]])
      print(c(arival[[j]], "ARI:", i, j))
    }
  }
}

################################################################################

# Clustering method explanations and comparisons 

################################################################################

# PAM and Kmeans

# k - medoids problem is a clustering problem similar to the k-means.

# Both k-means and k-medoids algorithms are partitional (breaking the dataset up
# into groups) and attempt to minimize the distance between points labeled to be
# in a cluster and a point designated as the center of that cluster.

# In contrast to the k-means algorithm, k-medoids chooses actual data points as
# centers (medoids), and thereby allows for greater interpretability of the cluster 
# centers than in k-means, where the center of a cluster is not necessarily one of 
# the input data points (it is the average distance between the points in the cluster).

# k-medoids can be used with arbitrary dissimilarity measures, whereas k-means 
# generally requires Euclidean distance for efficient solutions. Because k-medoids
# minimizes a sum of pairwise dissimilarities instead of a sum of squared Euclidean
# distances, it is more robust to noise and outliers than k-means(NB it has a low
# robustness).

# The number k of clusters is assumed known a-priori (which implies that the programmer 
# must specify k before the execution of a k-medoids algorithm). 

# The "goodness" of the given value of k can be assessed with methods such as the
# silhouette method.

# The medoid of a cluster is defined as the object in the cluster whose average 
# dissimilarity to all the objects in the cluster is minimal, that is, it is a most
# centrally located point in the cluster. 

# WHEN USING KMEANS INSTEAD OF PAM(K-MEDOIS) AND VICEVERSA

# K-means minimizes within-cluster variance since if you look at the definition of
# variance, it is identical to the sum of squared Euclidean distances from the center.

# It is not correct to use arbitary distances: because k-means may stop converging 
# with other distance functions.

# If you want arbitrary distance functions, we need to use k-medoids PAM. The medoid 
# minimizes arbitrary distances (because it is defined as the minimum), and there 
# only exist a finite number of possible medoids, too (for convergence). 

# In practice PAM is recommended to use whenever Euclidean Distance does not make 
# sense in your data (such as high dimensions for convergence).

################################################################################

# HCLUST EXPLANATION (ALL LINKAGES (WARD INCLUDED))

# The idea behind hierarchical agglomerative clustering methods starts by treating 
# each observation as a singleton cluster and then merge one by one each observation
# according to the lower values of those dissimilarity that we give as input, until
# we reach an unique cluster.
# The "methods" i.e single, average, complete and ward tell us how to evaluate the 
# dissimilarity between clusters in order to find the minimal one for the merging.
# Essentially :

# Single --> take the dissimilarity that we give as input between the closest 
#    observations of the clusters as benchmarks. 
#    (Emphasizes between-cluster separation disregarding homogeneity)
# EXTREME

# Average --> take the average of the dissimilarities  that we give as input between
#    the observations of the clusters as benchmarks. 
#    (Trade-off btw single and complete linkage)
# PREFERABLE

# Complete --> take the dissimilarity that we give as input between between the 
#    farthest observations of the clusters as benchmarks.
#    (Emphasizes within-cluster homogeneity disregarding separation)
# EXTREME

# Ward --> Instead of measuring the distance directly, it analyzes the variance 
#    of clusters. Ward's method says that the distance between two clusters is
#    how much the sum of squares will increase when we merge them.
#    Thus we merge clusters according to the lower increase of their variances.
# WARD IS SAID TO BE THE MOST SUITABLE METHOD FOR QUANTITAVE VARIABLES
# NB  also produce good dendrograms when the data are not so meaningful

################################################################################

# Hierarchical clustering (hclust) vs partitional clustering (PAM,kmeans)

# Partitional:
# -> Faster
# -> Requires stronger assumptions: N of cluster K and the initial centers
#                                   (Initializzation).
# -> sensible to the initialization of the centers.(different runs can leads to 
#                                                         different results)
# -> Partitional clustering results in exactly k clusters.(Hierarchical in a single
# one, in the sense that not always is possible find a suitable k0 for hclust 
# methods)

# Hierarchical: 
# -> Slower 
# -> Requires only a similarity measure (does not require any input parameters)
# -> Returns a much more meaningful and subjective division of clusters.
# -> more suitable for categorical data e.g species, class, ... whenever an 
#    "hierarchy" can be built

# NB: we can use the dendrograms in order to "approximate" the number of cluster k
# to put as input in a partitional algorithm.

# NB --> In order to have a more objective way for the selection of clustering 
# methods and distances i prefer to run the clustering with different distances
# and then make a choice based on a "quality" index such as the ASW.

################################################################################



