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

setwd("C:/Users/loren/Desktop/Henning_Exam/Data")
tombdata <- read.table("tombdataX00.dat",header=TRUE)
# Rows : 77 tombs (observation)
# Col : 145 types of pottery artifacts(manufatti in ceramica)(variables)
# cells : presence or not(0/1) of a specific type of pottery artifact

# Aim : tombs with similar artifacts are likely to come from the same period.
# The purpose of clustering is to find clusters of tombs that have been created 
# in the same period(similar) or vicerversa(dissimilar).
# Produce two distance-based clusterings of the tombs, explain why you chose the
# specific clustering methods and motivate all the methodological decisions you made,
# including your choice of a distance. Produce a visualisation of each clustering.
# Compare the clusterings and discuss to what extent each of them may be helpful 
# for tomb dating so that the discussion can be understood by an archaeologist.

set.seed(123456)

# Since we are dealing with p  = 145 binary variables they have the same scale,
# thus we don't need to standardize the variables.

# CHOICHE OF THE DISTANCE 
# Since tombs with similar artifacts are likely to come from the same period and 
# since we need to find clusters of tombs that have been created in the same period
# only joint presences of pottery artifacts are meaningful thus i use the jaccard 
# distance. If we use the simple matching one we can end up with similar tomb that 
# are similar only for the joint absence of pottery artifacts, that is not what we 
# want. 

jtomb <- dist(x = tombdata, method = "binary")

# A QUICK LOOK TO THE DATA 
# Due to the nature and the dimensions of the data I decide to visualize those 
# with a dimensional reduction technique.
# The MDS is a dimensional reduction technique used to display the information 
# contained in a distance matrix. It Takes as input a distance matrix between each 
# pair of observations and a chosen number of dimensions, then an MDS algorithm 
# places each object into a lower dimensional-dimensional space such that the 
# distances are preserved as well as possible. 
# I use this technique in order to visualize our data basing on distances in two
# dimensions for making a quick idea on how the observations looks like.
mds <- mds(jtomb, ndim=2)
plot(mds$conf) 
mds$stress
# However, with this dimensional reduction technique  we may lose some informations
# from the data, in this case we lose about the 37% of the informations.
# Probably is better go higher with the dimensions but since we aim to have a clear
# visualizzation of those, we keep ndim = 2.
# As we can see no clear patterns and agglomerations are present and the data seems
# very spread.
# Maybe in the two extremes left and right we have more points and also some small
# clouds in the "center" that seems a little bit separated.
# From this i think that an hierarchical method seems more appropriate due to the 
# spread nature of this data, since it will produce more meaningful clusters with 
# respect to a partitional method that we can interpret(due to the hierarchy construction) 
# and doesn't require any input parameters as the number of cluster, since is very
# hard say something about it in this case.
# However, for having a more "objective" way of choosing, I decide to run all different
# distance-based methods with a different number of partition k in a big loop,
# and then choosing the "best" according to the ASW(average silhouette width).
# To sum up the "criterion" of choice we can say that: the ASW is the mean of the
# s(i) (silhouette values), where s(i) measures the fit of each observation(i) to
# its cluster for a given partition.
# Thus we are going to choose the method and the number of k according to the highest 
# value of the ASW.
# For me, this is a very natural way in order to choose which partition taking into 
# account.

# LOOP INITIALIZZATION
tsingle <- tave <- tcom <- tward <- tpam <- list()
# These lists will collect all clustering labels for the observations
# tsingle = cluster labels from hierachical method hclust with single linkage
# ... 

nc <- 2:15 # Number of k clusters that we want to try

ssil <- asil <- csil <- wsil <- psil <- list()
# These lists will collect the s(i) values for each observation
# ssil =  s(i)s obtained with hclust + single linkage 
# ...

sasw <- aasw <- casw <- wasw <- pasw <- NA
# Vectors of ASW values: 1 for each clustering method with length equal to k number
# of partitions.

# RUN OF THE ALL DISTANCE BASED CLUSTERING METHODS WITH DIFFERENT K CLUSTERS

# I run hierarchical methods first to produce the dendrograms, so that in the loop
# only cutree has to be applied.
single <- hclust(jtomb,method="single")
average <- hclust(jtomb,method="average")
complete <- hclust(jtomb,method="complete")
ward <- hclust(jtomb,method="ward.D2")

for(k in nc){
  print(k)
  tsingle[[k]] <- cutree(single,k)
  tave[[k]] <- cutree(average,k)
  tcom[[k]] <- cutree(complete,k)
  tward[[k]] <- cutree(ward,k)
  tpam[[k]] <- pam(jtomb,k)
  
  ssil[[k]] <- silhouette(tsingle[[k]],dist=jtomb)
  asil[[k]] <- silhouette(tave[[k]],dist=jtomb)
  csil[[k]] <- silhouette(tcom[[k]],dist=jtomb)
  wsil[[k]] <- silhouette(tward[[k]],dist=jtomb)
  psil[[k]] <- silhouette(tpam[[k]],dist=jtomb)
  
  sasw[k] <- summary(ssil[[k]],dist=jtomb)$avg.width
  aasw[k] <- summary(asil[[k]],dist=jtomb)$avg.width
  casw[k] <- summary(csil[[k]],dist=jtomb)$avg.width
  wasw[k] <- summary(wsil[[k]],dist=jtomb)$avg.width
  pasw[k] <- summary(psil[[k]],dist=jtomb)$avg.width
}

# Ko optimal number of clusters for each methods according to the ASW 
which.max(sasw) # [1] 2
which.max(aasw) # [1] 2
which.max(casw) # [1] 15
which.max(wasw) # [1] 15 # seems the best 
which.max(pasw) # [1] 14

# Correspondent ASW val
max(sasw,na.rm=TRUE) # 0.04971287
max(aasw,na.rm=TRUE) # 0.04971287
max(casw,na.rm=TRUE) # 0.04017532
max(wasw,na.rm=TRUE) # 0.06027678 # max
max(pasw,na.rm=TRUE) # 0.05734671

# Summary plot of all the ASW curves 
plot(1:15,sasw,ylim=c(0,0.3),type="l",xlab="Number of clusters",ylab="ASW")
points(1:15,aasw,ylim=c(0,0.3),type="l",col=2,lty=2)
points(1:15,casw,ylim=c(0,0.3),type="l",col=3,lty=3)
points(1:15,wasw,ylim=c(0,0.3),type="l",col=4,lty=4)
points(1:15,pasw,ylim=c(0,0.3),type="l",col=5,lty=5)
legend(8,0.3,legend=c("single","average","complete","ward","pam"),lty=1:5,col=1:5)

# CHOICHE OF THE CLUSTERING METHODS (PAM and WARD)

# Lets first have a look to the dendrograms.
plot(single)
plot(average)
plot(complete)
# For single, average and complete the number of clusters would need to be vastly
# increased to get at some structure in the big set of observations.
plot(ward)
# With ward linkage seems that we reach a certain hierarchy among the observations,
# thus i decide to focus on this method. 
# The ASW overall plot seems no suggest a clear number of K for ward method, thus 
# i decide to choose the k according to the dendrogram.
# Moreover, we can notice that we have the highest ASW values both for WARD and for 
# PAM instead of single, average and complete, thus i decide to choose the 
# hierarchical clustering method with ward linkage and the partitioning clustering
# method k-medoids(PAM).

# I choose the hierarchical clustering method with ward linkage for the reasons 
# explained above, since we can reach some structure in the big set of observations
# with a reasonable number of clusters.

# While i choose the partitioning clustering method k-medoids(PAM) due to: the high
# ASW values; since, instead of k-means, it can deal with any dissimilarity measure
# because the euclidean distance doesn't make sense in our data; and due to his better 
# interpretability, since as centers(medoids) it takes observations from the data.

# CLUSTERING VISUALIZZATION

plot(ward)
# By looking the dendrogram i decide to take 6 as optimal number of clusters for ward.
plot(mds$conf,pch=clusym[tward[[6]]],col=tward[[6]], main="ward, 6 clusters")

# By looking the overall ASW plot i decide to take 9 as optimal number of clusters
# for PAM (seems that for k=9 we a sort of local maximum for the ASW)
plot(mds$conf,pch=clusym[tpam[[9]]$clustering],col=tpam[[9]]$clustering, main="pam, 9 clusters")

# silhouette plots
plot(wsil[[6]],main="ward, 6 clusters")
plot(psil[[9]],main="pam, 9 clusters")

#ADJ
adjustedRandIndex(tward[[6]],tpam[[9]]$clustering) # 0.3117985


