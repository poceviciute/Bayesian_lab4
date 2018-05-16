# Lab 4

# Question 1

#a)
data <- read.table("eBayNumberOfBidderdata.dat", header = TRUE)
cols <- c(1,3:ncol(data))
data_noC <- data[,cols]
glm_fit <- glm(nBids~., family=poisson, data_noC)
summary(glm_fit)
#Intercept, VerifyID, Sealed, MajBlem, LogBook, MinBidShare

#b)
library(mvtnorm)
#xtx <- t(as.matrix(data))%*%as.matrix(data)

log_poisson <- function(beta, y, x){
  x <- as.matrix(x)
  xtx <- solve(t(x)%*%x)
  y <- as.matrix(y)
  #print("log poisson")
  #print(dim(x))
  #print(dim(xtx))
  lin_pred <- x%*%beta
  loglik <- sum(y*lin_pred-exp(lin_pred))/sum(log(factorial(y)))
  logprior <- dmvnorm(beta, mean=as.vector(rep(0,ncol(x))), sigma=100*xtx, log = TRUE)
  logpost <- loglik+logprior
  logpost
}
#log_poisson(beta = c(rep(1,ncol(data)-1)), y,x)
betas_init <- as.vector(rep(0,ncol(data)-1))
y <- data$nBids
x<- data[,2:ncol(data)]
optim_fit <- optim(betas_init, log_poisson, gr=NULL, y, x, method=c("BFGS"),control=list(fnscale=-1), hessian=TRUE)

beta_tilde <- optim_fit$par
hessian <- -1*optim_fit$hessian
inv_hessian <- solve(hessian)

post_approx <- rmvnorm(n=1000, mean=beta_tilde, sigma=inv_hessian)
colnames(post_approx) <- colnames(x)
#hist(post_approx[,"VerifyID"], freq = FALSE)


#c)

metropolis <- function(n, c, sigma, logpostfun,theta,...){
  thetas <- matrix(nrow=n+1, ncol=length(theta))
  thetas[1,]<-theta
  #logpostfun(theta,...)
  temp1 <-logpostfun(thetas[1,],...)
  #print(thetas[1,])
  acc_prob <- vector(length=n)
  acc_prob[1] <- 0
  for(i in 1:n){
    temp_theta<-rmvnorm(n=1, thetas[i,], c*sigma)
    temp2 <-logpostfun(as.vector(temp_theta),...)
    #print(paste("iter",i,"temp2",temp2))
    acc_prob[i+1] <- min(1,exp(temp2-temp1))
    #print(acc_prob)
    u <- runif(n = 1, 0,1)
    #print(u)
    if(u>acc_prob[i+1]){
      thetas[i+1,]<-thetas[i,]
    }
    else{
      thetas[i+1,] <- temp_theta
      temp1 <- temp2
    }
    
  }
  data.frame(thetas, "acc.prob"=acc_prob)
}
metro2 <- metropolis(1000, c=1, sigma=inv_hessian, log_poisson, theta=betas_init, y=y, x=x)
metro <- metropolis(1000, c=0.5, sigma=inv_hessian, log_poisson, theta=betas_init, y=y, x=x)

acc_prob <- as.vector(metro$acc.prob)
avg_acc <- mean(acc_prob)

beta_metro <- as.matrix(metro[,1:length(betas_init)])
phi <- exp(beta_metro)
plot(phi[,1], type="l")
for(i in 1:ncol(phi)){
  plot(phi[,i], type="l", xlab="Iteration", ylab=paste0("phi_", i), main=paste0("phi_", i))
}
phi_means <- colMeans(phi)

for(i in 1:ncol(phi)){
  hist(phi[,i], freq=FALSE, xlab=paste0("phi_", i), main=paste0("phi_", i))
  lines(density(phi[,i]), col="red")
}


#d)

x_vec <- c(1,1,1,1,0,0,0,1,0.5)
lambda <- exp(x_vec%*%t(beta_metro))
y_mat <- apply(lambda, 1, function(i){rpois(1000, i)})
hist(y_mat, freq = FALSE, xlab="y_pred")
lines(density(y_mat))
hist(lambda)
prob_y0 <- length(y_mat[y_mat==0])/length(y_mat)

#data[which(data$nBids==0),]
#data[,x_vec]
