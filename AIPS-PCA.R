#AIPS.pc computes PC scores using function either "eigen" or "svd".

AIPS.pc <- function(infile, K=NULL, method="eigen",outplot) {
  #read merged data coded with Additive Components
  #When merging data in Plink, the data with Ancestry Informtive Markers(AIMs) should be in the first part of the merged data.
  #If the data with AIMs is different part, please modify the script by ancestry index.
  plink.dta <- read.table(infile, header =TRUE)
  X<-plink.dta[,7:length(plink.dta)];
  p<-length(plink.dta)-6;p;#number of markers
  n<-nrow(plink.dta);n;#number of sample size
  
  x_mean <- rep(NA, p)
  x_sd <- as.vector(rep(NA, p) )
  
  for(j in 1:p){
    x_mean[j] <- mean(X[,j], na.rm=TRUE)
    x_sd[j] <- sd(X[,j], na.rm=TRUE)
    replace(X[,j], is.na(X[,j]) , x_mean[j]) -> X[,j]  # replace mean if there is missing genotype value.
  }
  
  adj.std <- cbind(x_mean,x_sd)
  
  vec1_sample <- as.vector(rep(1,n))
  vec1_snps  <- as.vector(rep(1,p))
  X <- as.matrix(X)
  
  # Standardize the data
  X_sig <- vec1_sample%*% t(x_sd)
  X_std <- (X - vec1_sample %*% t(x_mean))/X_sig;
  
  ## Compare n and p (n<p or n>=p) and compute the eigenvectors and scores 
  ####on larger number of samples than markers;
  #Choose the method from either "eigen" or "svd"
  if(n >= p & method=="eigen"){
     r <-eigen(cov(X_std))
     evec <- r$vectors
     eval <- r$values
     sdev <- sqrt(eval)
    if (is.null(K)) {
      scores <- X_std%*%evec
    }
    else {
      scores <- X_std%*%evec[,1:K]  
    }
  }
  
  if(n >= p & method=="svd"){
     r <- svd(X_std)
     evec <- r$v
     sdev <- r$d/sqrt(p-1)
     eval <- sdev*sdev
    if (is.null(K)) {
      scores <- X_std%*%evec
    }
    else {
      scores <- X_std%*%evec[,1:K]  
    }
  }
   
  #on larger number of markers(n<p);    
  #Choose the method from either "eigen" or "svd"
  if(n < p & method=="eigen"){
     r <- eigen(cov(t(X_std)))
     eval <- r$values
     evec <- t(X_std)%*%r$vectors
    if (is.null(K)) {
      scores <- X_std%*%evec
    }
    else {
      scores <- X_std%*%evec[,1:K]  
    }
  }
     
  if(n < p & method=="svd"){
     r <- svd(t(X_std))
     evec <- t(X_std)%*%r$v
     sdev <- r$d/sqrt(p-1)
     eval <- sdev*sdev
    if (is.null(K)) {
      scores <- X_std%*%evec
    }
    else {
      scores <- X_std%*%evec[,1:K]  
    }
  }
  
  score0 <- cbind(plink.dta[,1:2],scores) 
  # to generate plot to check eigenvalues
  png(outplot)
  plot(c(1:10),eval[1:10],xaxt="n",type="o",xlab="Order of Eigenvalues",ylab="Eigenvalues")
  title(main="Plot of Top Eigenvalues", col.main="black", font.main=1)
  dev.off()
  #Return the results as a list and save
  list (snp.weights=evec, eval=eval,pcs.discovery=score0,adj.discovery=adj.std)
}

W <- AIPS.pc(infile="outputA.raw", K=10, method="svd",outplot="output.png")

score0 <- W$pcs.discovery
eigenvalues    <- W$eval
save(score0,eigenvalues,file="pca.euro.RData")
