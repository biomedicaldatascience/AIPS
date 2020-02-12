#To Infer ancestry memberships among subpopulations in Europe.

load("pca.euro.RData")

#For ancestry unknown samples, replace FID in .fam file with 0 before running AIPS-AI.R.
replace(score0[,1], (score0[,1]==0) , 9999) -> score0[,1]

#Generate a index table for ancestry known samples (see "eurosubpop7.txt" for 7 clusters or "intraeuro.952samples.22idx.txt" for 22 subpopulations)
n.subpop <- read.table("eurosubpop7.txt",header=T);#for 7 subclusters
#n.subpop <- read.table("intraeuro.4376samples.22idx.txt",header=T);#for 22 subpopulations
n.pop=nrow(n.subpop)
n <- nrow(score0)
unk.n <- n-n.subpop[n.pop,2]

# Calculate the Centroids from each ancestry known subpopulation on the top five scores (PCs).
n.score   <- 5                         # input : number of scores / number of eigenvalues which will be used on analysis
n.mixture <- 3                         # assumtion: 3= max number of admixtures
score1     <- score0[,3:(2+n.score)]

# for example, we define 7 European subpopulations. 
#If you have less or more subpopulations, need to modify the number of subpopulations.
known.n <- n.subpop

# To compute averages from different number of subpopulations and/or more/less scores, please modify the below commands.
#### this is based on reclassified 7 subpopulatioin
#### require to sort data by FID
score <- score1[order(score0[,1]),]

#nrow=n.subpopulations and ncol=n.top scores
mean_g <- matrix(rep(NA, (n.pop*n.score)), nrow=n.pop, ncol=n.score,byrow = TRUE);

for(i in 1:n.pop){
  for(j in 1:n.score) {
    mean_g[i,j] <- mean(score[known.n[i,1]:known.n[i,2], j], na.rm=TRUE)
  }
}

pop.mean  <- mean_g[,1:n.score]

#eigenval  <- as.vector(eigenvalues[1:n.score,])
eigenval  <- eigenvalues[1:n.score]

########## Calculate distance from  each score to the center points of known ancestry
centers=t(pop.mean);scores=t(score1)

#calculate distance
n.center  <- length(pop.mean[,1])

dist.ancestry = apply(centers,2,function(center){sqrt(colSums((scores-center)^2))});
dist.order <- t(apply(dist.ancestry,1,order))

eigenval.sum <- sum(eigenval)
wdist.ancestry = apply(centers,2,function(center){sqrt(t((scores-center)^2) %*% eigenval/eigenval.sum)});
wdist.order <- t(apply(wdist.ancestry,1,order))

#dsistance among centers
dist.center  = apply(centers,2,function(center){sqrt(colSums((centers-center)^2))});
wdist.center = apply(centers,2,function(center){sqrt(t((centers-center)^2) %*% eigenval/eigenval.sum)});

#compute  proportion of line segments to centroids for defining group menbership
powerD.prob    <- matrix(0, nrow=n , ncol=n.center)
power.dist     <- matrix(NA, nrow=n , ncol=n.mixture)
powerD.sum     <- rep(NA, n)

exponD.prob    <- matrix(0, nrow=n , ncol=n.center)
expon.dist     <- matrix(NA, nrow=n , ncol=n.mixture)
exponD.sum     <- rep(NA, n)

eigenvD.prob    <- matrix(0, nrow=n , ncol=n.center)
eigenv.dist     <- matrix(NA, nrow=n , ncol=n.mixture)
eigenvD.sum     <- rep(NA, n)

alpha=1

for(i in 1:n){
  #Power Distance weights
  power.dist[i,]    <- 1/{dist.ancestry[i,dist.order[i,1:n.mixture]]}^alpha
  power.dist[i, 1+ which(dist.center[dist.order[i,1],dist.order[i,2:n.mixture]] < dist.ancestry[i,dist.order[i,2:n.mixture]])] <- 0
  powerD.sum[i] <- sum(power.dist[i,], na.rm=TRUE)
  powerD.prob[i,dist.order[i,1:n.mixture]] <- power.dist[i,]/powerD.sum[i]
  
  #Exponential Distance weights
  expon.dist[i,]    <- 1/exp(alpha*dist.ancestry[i,dist.order[i,1:n.mixture]])
  expon.dist[i, 1+ which(dist.center[dist.order[i,1],dist.order[i,2:n.mixture]] < dist.ancestry[i,dist.order[i,2:n.mixture]])] <- 0
  exponD.sum[i] <- sum(expon.dist[i,], na.rm=TRUE)
  exponD.prob[i,dist.order[i,1:n.mixture]] <- expon.dist[i,]/exponD.sum[i]
  
  #EVD
  eigenv.dist[i,]    <- 1/{wdist.ancestry[i,wdist.order[i,1:n.mixture]]}^alpha
  eigenv.dist[i, 1+ which(wdist.center[wdist.order[i,1],wdist.order[i,2:n.mixture]] < wdist.ancestry[i,wdist.order[i,2:n.mixture]])] <- 0
  eigenvD.sum[i] <- sum(eigenv.dist[i,], na.rm=TRUE)
  eigenvD.prob[i,dist.order[i,1:n.mixture]] <- eigenv.dist[i,]/eigenvD.sum[i]
}


final.pw <- cbind(score0[,1:7],power.dist,powerD.prob) 
final.pweval <- cbind(score0[,1:7],eigenv.dist,eigenvD.prob) 
final.pwexp <- cbind(score0[,1:7],expon.dist,exponD.prob)  

write.table(final.pw, "./test1/final.pw.csv",col.names=T,row.names=T, sep=",")
write.table(final.pweval,"./test1/final.pweval.csv",col.names=TRUE,row.names=T, sep=",")
write.table(final.pwexp, "./test1/final.pwexp.csv",col.names=TRUE,row.names=T, sep=",")

center.order <- t(apply(dist.center,1,order))
write.table(dist.center, "./test1/distance_among_centers.csv",col.names=TRUE,row.names=TRUE, sep=",")
write.table(center.order,"./test1/rank_distance_among_centers.csv",col.names=TRUE,row.names=TRUE, sep=",")

wcenter.order <- t(apply(wdist.center,1,order))
write.table(wdist.center, "./test1/wdistance2.csv",col.names=TRUE,row.names=TRUE, sep=",")
write.table(wcenter.order,"./test1/wrank2.csv",col.names=TRUE,row.names=TRUE, sep=",")
