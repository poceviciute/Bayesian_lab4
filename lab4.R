# Lab 4

# Question 1

#a)
data <- read.table("eBayNumberOfBidderdata.dat", header = TRUE)
c <- c(1,3:ncol(data))
data_noC <- data[,c]
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

metropolis <- function(n, c, sigma, logpostfun,theta,...){
  thetas <- matrix(nrow=n+1, ncol=length(theta))
  thetas[1,]<-theta
  #logpostfun(theta,...)
  temp1 <-logpostfun(thetas[1,],...)
  #print(thetas[1,])
  for(i in 1:n){
    temp_theta<-rmvnorm(n=1, thetas[i,], c*sigma)
    temp2 <-logpostfun(as.vector(temp_theta),...)
    #print(paste("iter",i,"temp2",temp2))
    acc_prob <- min(1,exp(temp2-temp1))
    #print(acc_prob)
    u <- runif(n = 1, 0,1)
    #print(u)
    if(u>acc_prob){
      thetas[i+1,]<-thetas[i,]
    }
    else{
      thetas[i+1,] <- temp_theta
      temp1 <- temp2
    }
    
  }
  thetas
}
beta_metro <- metropolis(1000, c=1, sigma=inv_hessian, log_poisson, theta=betas_init, y=y, x=x)
phi <- exp(beta_metro)
plot(phi[,1], type="l")
