# import libraries
library(smoothmest) # double exponential distribution
library(MASS) # mv-normal distribution
library(MCMCpack) 
library(mcmcse) 
library(coda)

# import functions
source('functions.R')

# import data, I35
i35e = read.csv("./i35e.csv")[,-1]

# markov-chain monte carlo function (obtaining draws from posterior distribution)
mcmc = function(data, nreps=100000, q){
  # order the data from smallest to largest beginning mileposts (for spatial estimation)
  data = data[order(data$BMP, decreasing=FALSE),] 
  
  # create proximity matrix
  W = create.w(data, order=1)
  
  # response variable
  y = data$TOTALCRASHES
  
  # X matrix
  covariates = c("TOTALCRASHES","MEDWID","MED_TYPE","LSHLDWID","LSHL_TYP","SURF_WID","RSHLDWID","NO_LANES","LANEWID")
  X = model.matrix(TOTALCRASHES~., data=data[,covariates])
  
  # orthogonalize spatial effects
  P.orthog = diag(nrow(X)) - X%*%solve(t(X)%*%X) %*% t(X)
  M.eigen = eigen(P.orthog %*% W %*% P.orthog)
  M = M.eigen$vectors[,order(M.eigen$values,decreasing=TRUE)[1:q]]
  X.star = cbind(X,M)
  Q = diag(rowSums(W)) - W
  Q.s = t(M)%*%Q%*%M
  
  # calculate expected crashes
  data$length = (data$EMP-data$BMP)
  data$expected.crash = sum(data$TOTALCRASHES) / (sum(data$AADT*data$length)) * data$AADT*data$length

  # center and scale covariates (for better convergence)
  data[,c(18,20,21,23,24,25)] = scale(data[,c(18,20,21,23,24,25)],center=TRUE,scale=TRUE)

  # initialize matrix for storing draws
  theta.draws <- matrix(0,nrow=nreps,ncol=ncol(X)+q)
  
  # begin draws at maximum likelihood estimates
  mod = glm(y ~ X.star[,-1], family="poisson", offset=log(data$expected.crash), data=data[,c(4,16:21,24,25)])
  theta.draws[1,] <- mod$coef
  
  # replace NA's with 0's
  theta.draws[1,] = ifelse(is.na(theta.draws[1,]),0,theta.draws[1,])
  
  # initialize beta vector
  beta = theta.draws[,1:ncol(X)]
  colnames(beta) = names(glm(TOTALCRASHES ~ ., family="poisson", offset=log(data$expected.crash), data=data[,c(4,16:21,24,25)])$coef)
  
  # initialize phi matrix
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
      # update theta
      theta.proposal <- mvrnorm(1, theta.current, Sigma=(0.01^2) * diag(ncol(theta.draws))) 
    } else {
      theta.proposal <- mvrnorm(1, theta.current, Sigma= 2.4^2/ncol(theta.draws) * (amcmc$var + (0.01^2) * diag(ncol(theta.draws))) )
    }
    
    # calculate metropolis hastings ratio
    mh.ratio <- likelihood(as.matrix(theta.proposal),X.star,y,data$expected.crash) + 
      prior.beta(s=s[i-1],theta.proposal[1:ncol(X)]) + 
      prior.phi.star(theta.proposal[(ncol(X)+1):length(theta.proposal)],tau=tau[i-1],Q.s,q) - 
      likelihood(as.matrix(theta.current),X.star,y,data$expected.crash) - 
      prior.beta(s=s[i-1],theta.current[1:ncol(X)]) - 
      prior.phi.star(theta.current[(ncol(X)+1):length(theta.current)],tau=tau[i-1],Q.s,q)
    
    
    # update draws
    if(log(runif(1)) < mh.ratio) {
      theta.draws[i,] <- theta.proposal
      c <- c + 1
    } else { 
      theta.draws[i,] <- theta.current
    }
    
    # adaptive mcm update
    amcmc <- AMCMC.update(theta.draws[i,],amcmc$mn,amcmc$var,i)
    
    beta[i,] = theta.draws[i,(1:ncol(X))]
    phi.star[i,] = cbind(theta.draws[i,(ncol(X)+1):ncol(theta.draws)])
    phi[i,] = M %*% cbind(theta.draws[i,(ncol(X)+1):ncol(theta.draws)]) 
    
    
    # update s
    a = 2.01
    b = 1
    s[i] = rinvgamma(1,a + ncol(beta), sum(abs(beta[i,]) + b))
    
    # update tau
    a.tau = 2.01
    b.tau = 1
    tau[i] = rgamma(1,q/2 + a.tau, 1/2*t(phi.star[i,])%*%Q.s%*%phi.star[i,] + 1/b.tau ) 
    
    # store -2 log likelihood (for DIC)
    d.theta[i] =  -2 * likelihood(as.matrix(theta.draws[i,]),X.star,y,data$expected.crash)
    
    # prediction
    mu[i-1,] = exp(X%*%beta[i,] + phi[i,])
    prediction[i-1,] = rpois(ncol(mu),data$expected.crash * mu[i-1,])
    print(i)
  }
  
  # calculate DIC
  theta.mean = amcmc$mn
  dbar = mean(d.theta)
  pd = dbar - -2 * likelihood(as.matrix(theta.mean),X.star,y,data$expected.crash)
  DIC = pd + dbar
  
  # 
  print(c/nreps)
  return(list("beta"=beta,"phi"=phi,"s"=s,"tau"=tau,"dic"=DIC,"pd"=pd,"prediction"=prediction,"mu"=mu))
}


# run monte carlo markov chain
draws <- mcmc(i35e,nreps=1000,q=70)


# remove burn-in
burn <- 150
beta <- draws$beta[-c(1:burn),]
phi <- draws$phi[-c(1:burn),]
prediction <- draws$prediction[-c(1:burn),]


# Beta credible intervals #
beta.quantiles <- apply(beta,2,function(x) quantile(x,probs=c(.025,.975)))
beta.quantiles <- rbind(beta.quantiles,apply(beta,2,mean)) 
row.names(beta.quantiles)[3] <- "Mean"
# #

# Prediction intervales #
pred.quantiles <- apply(prediction,2,function(x) quantile(x,probs=c(.025,.975)))
# #

# plot prediction intervals #
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
