library(smoothmest) # double exponential distribution
library(MASS) # MV norm distribution
library(MCMCpack)
library(mcmcse) 
library(coda)

i35e <- read.csv("./i35e.csv")[,-1]

create.w <- function(data,order=1){
  n <- nrow(data)
  w <- matrix(0,nrow=n,ncol=n)
  if(order==1){
    for(i in 1:(n-1)){
      w[i+1,i] <-  w[i,i+1] <- w[1,2] <- w[n-1,n] <- 1
    }
  }
  
  if(order==2){
    for(i in 1:(n-2)){
      w[i+1,i] <-  w[i+2,i] <- w[i,i+1] <- w[i,i+2] <- w[1,2] <- w[1,3] <- w[n-1,n] <- w[n,n-1] <- 1
    }
  }
  
  return(w)
}

AMCMC.update <- function(draw,cur.mn,cur.var,cur.it){
  if(cur.it >0){
    mn <- ((cur.it-1)*cur.mn+draw)/cur.it
    if(cur.it==1){
      v <- matrix(0,nrow=length(draw),ncol=length(draw))
    } else {
      v <- (cur.it-2)*cur.var+(cur.it-1)*(cur.mn%*%t(cur.mn))+draw%*%t(draw)
      v <- (v-cur.it*(mn%*%t(mn)))/(cur.it-1)
    }
  } else {
    mn <- matrix(0,nrow=nrow(cur.mn),ncol=1)
    v <- matrix(0,nrow=nrow(draw),ncol=nrow(draw))
  }
  return(list(mn=mn,var=v))
}

likelihood <- function(theta,x.star,y,expected.crash,log=TRUE){
  l <- sum( dpois(y, lambda = expected.crash * exp((x.star)%*%theta), log=log) )
  return(l)
}

prior.beta <- function(s,beta,log=TRUE){
  p <- sum( log(ddoublex(beta, mu=0, lambda=s)) )
  return(p)
}

prior.phi.star <- function(phi.star,tau,Q.s,q){
  (q/2)*log(tau)-(tau/2)*t(phi.star)%*%Q.s%*%phi.star
}

# markov-chain monte carlo function (obtaining draws from posterior distribution)
mcmc <- function(data,nreps=10000,q){
  data <- data[order(data$BMP,decreasing=FALSE),] # order the data from smallest to largest beginning mileposts, and keep only variables of interest
  y <- data$TOTALCRASHES
  
  X <- model.matrix(TOTALCRASHES ~ ., family="poisson",offset=log(data$expected.crash),data=data[,c(4,16:21,24,25)])
  W <- create.w(data) # make sure that the data is is order, smallest to largest
  
  ##  Create new X matrix ##
  P.orthog <- diag(nrow(X)) - X %*% solve(t(X)%*%X) %*% t(X)
  M.eigen = eigen(P.orthog %*% W %*% P.orthog)
  M <- M.eigen$vectors[,order(M.eigen$values,decreasing=TRUE)[1:q]]
  X.star <- cbind(X,M)
  Q = diag(rowSums(W)) - W
  Q.s = t(M)%*%Q%*%M
  
  ## Expected crashes ##
  data$length = (data$EMP-data$BMP)
  data$expected.crash = sum(data$TOTALCRASHES) / (sum(data$AADT*data$length)) * data$AADT*data$length
  ## ##
  
  ## Center and scale ##
  data[,c(18,20,21,23,24,25)] = scale(data[,c(18,20,21,23,24,25)],center=TRUE,scale=TRUE)
  ## ##
  
  theta.draws <- matrix(0,nrow=nreps,ncol=ncol(X)+q)
  theta.draws[1,] <- (glm(y ~ X.star[,-1], family="poisson",offset=log(data$expected.crash),data=data[,c(4,16:21,24,25)])$coef)
  theta.draws[1,] = ifelse(is.na(theta.draws[1,]),0,theta.draws[1,])
  beta = theta.draws[,1:ncol(X)]
  colnames(beta) = names(glm(TOTALCRASHES ~ ., family="poisson",offset=log(data$expected.crash),data=data[,c(4,16:21,24,25)])$coef)
  phi = matrix(0,nrow=nreps,ncol=length(y))
  phi.star = matrix(0,nrow=nreps,ncol=q)
  phi.star[1,] = theta.draws[1,(ncol(X)+1):ncol(theta.draws)]
  s = numeric(nreps)
  s[1] = 1
  tau = numeric(nreps)
  tau[1] = 1
  d.theta = numeric(nreps)
  d.theta.2 = matrix(0,nrow=nreps,ncol=ncol(X)+q)
  mu = matrix(0,nrow=nreps,ncol=nrow(X))
  prediction = matrix(0,nrow=nreps,ncol=nrow(X))
  
  amcmc <- list(mn=theta.draws[1,],var=matrix(0,nrow=length(theta.draws[1,]),ncol=length(theta.draws[1,])))
  
  c <- 0
  for(i in 2:nreps){
    theta.current <- theta.draws[i-1,]
    if(i<250){
      
      ## Update theta ##
      theta.proposal <- mvrnorm(1, theta.current, Sigma=(0.01^2) * diag(ncol(theta.draws))) 
    } else {
      theta.proposal <- mvrnorm(1, theta.current, Sigma= 2.4^2/ncol(theta.draws) * (amcmc$var + (0.01^2) * diag(ncol(theta.draws))) )
    }
    
    mh.ratio <- likelihood(as.matrix(theta.proposal),X.star,y,data$expected.crash) + 
      prior.beta(s=s[i-1],theta.proposal[1:ncol(X)]) + 
      prior.phi.star(theta.proposal[(ncol(X)+1):length(theta.proposal)],tau=tau[i-1],Q.s,q) - 
      likelihood(as.matrix(theta.current),X.star,y,data$expected.crash) - 
      prior.beta(s=s[i-1],theta.current[1:ncol(X)]) - 
      prior.phi.star(theta.current[(ncol(X)+1):length(theta.current)],tau=tau[i-1],Q.s,q)
    
    if(log(runif(1)) < mh.ratio) {
      theta.draws[i,] <- theta.proposal
      c <- c + 1
    } else { 
      theta.draws[i,] <- theta.current
    }
    
    amcmc <- AMCMC.update(theta.draws[i,],amcmc$mn,amcmc$var,i)
    
    beta[i,] = theta.draws[i,(1:ncol(X))]
    phi.star[i,] = cbind(theta.draws[i,(ncol(X)+1):ncol(theta.draws)])
    phi[i,] = M %*% cbind(theta.draws[i,(ncol(X)+1):ncol(theta.draws)]) 
    ## ##
    
    ## Update s ##
    a = 2.01
    b = 1
    s[i] = rinvgamma(1,a + ncol(beta), sum(abs(beta[i,]) + b))
    ## ##
    
    ## Update tau ##
    a.tau = 2.01
    b.tau = 1
    tau[i] = rgamma(1,q/2 + a.tau, 1/2*t(phi.star[i,])%*%Q.s%*%phi.star[i,] + 1/b.tau ) 
    ## ##
    
    
    ## DIC ##
    d.theta[i] =  -2 * likelihood(as.matrix(theta.draws[i,]),X.star,y,data$expected.crash)
    ## ##
    
    ## Prediction ##
    mu[i-1,] = exp(X%*%beta[i,] + phi[i,])
    prediction[i-1,] = rpois(ncol(mu),data$expected.crash * mu[i-1,])
    ## ##
  }
  
  ## DIC ##
  theta.mean = amcmc$mn
  dbar = mean(d.theta)
  pd = dbar - -2 * likelihood(as.matrix(theta.mean),X.star,y,data$expected.crash)
  DIC = pd + dbar
  ## ##
  
  print(c/nreps)
  return(list("beta"=beta,"phi"=phi,"s"=s,"tau"=tau,"dic"=DIC,"pd"=pd,"prediction"=prediction,"mu"=mu))
}


# run monte carlo markov chain #
test <- mcmc(i35e,nreps=120000,q=70)
burn <- 15000

beta <- test$beta[-c(1:burn),]
phi <- test$phi[-c(1:burn),]
prediction <- test$prediction[-c(1:burn),]
# #

# Beta credible intervals #
beta.quantiles <- apply(beta,2,function(x) quantile(x,probs=c(.025,.975)))
beta.quantiles <- rbind(beta.quantiles,apply(beta,2,mean)) 
row.names(beta.quantiles)[3] <- "Mean"
stargazer(t(beta.quantiles))

# #

# Prediction intervales #
pred.quantiles <- apply(prediction,2,function(x) quantile(x,probs=c(.025,.975)))
# #

# plot intervales #
plot(1:ncol(prediction),pred.quantiles[2,],pch="-",xlab="Road segment",ylab="Number of crashes",main="Prediction Intervals (I35E)")
points(1:ncol(prediction),pred.quantiles[1,],pch="-")
segments(1:ncol(prediction),pred.quantiles[1,],1:ncol(prediction),pred.quantiles[2,])
points(1:ncol(prediction),i35e$TOTALCRASHES,pch=20,cex=.8,col='blue')

for(i in 1:ncol(prediction)){
  if(i35e$TOTALCRASHES[i] < pred.quantiles[1,i] | i35e$TOTALCRASHES[i] > pred.quantiles[2,i])
    points(i,i35e$TOTALCRASHES[i],col='red',pch=18,cex=1.3)
}
# #

# plot phi #
phi.quantiles <- apply(phi,2,function(x) quantile(x,probs=c(.025,.5,.975)))

plot(1:ncol(prediction),phi.quantiles[3,],pch="-",xlab="Road segment",ylab=expression(phi),main="Spatial effects")
points(1:ncol(prediction),phi.quantiles[1,],pch="-")
segments(1:ncol(prediction),phi.quantiles[1,],1:ncol(prediction),phi.quantiles[3,])
abline(h=0,lty=2,col='dark green')

for(i in 1:ncol(prediction)){
  if(0 < phi.quantiles[1,i]){
    print(i)
    segments(i,phi.quantiles[1,i],i,phi.quantiles[3,i],col='red')
  }  
}

for(i in 1:ncol(prediction)){
  if(0 > phi.quantiles[3,i]){
    segments(i,phi.quantiles[1,i],i,phi.quantiles[3,i],col='blue')
    
  }
}
# #