# defines proximity matrix with either 1st or 2nd orders (any proximity > 2 in linear road data is likely noise)
create.w <- function(data,order=1){
  n <- nrow(data)
  W <- matrix(0,nrow=n,ncol=n)
  if(order==1){
    for(i in 1:(n-1)){
      W[i+1,i] <-  W[i,i+1] <- W[1,2] <- W[n-1,n] <- 1
    }
  }
  if(order==2){
    for(i in 1:(n-2)){
      W[i+1,i] <-  W[i+2,i] <- W[i,i+1] <- W[i,i+2] <- W[1,2] <- W[1,3] <- W[n-1,n] <- W[n,n-1] <- 1
    }
  }
  return(W)
}

# adaptive markov chain monte carlo, specified in section 2.4 
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

# likelihood function (note on the log scale for computational ease)
likelihood <- function(theta,x.star,y,expected.crash,log=TRUE){
  sum(dpois(y, lambda=expected.crash*exp((x.star)%*%theta), log=log))
}

# laplace (double exponential) distributional functions for LASSO penalty 
prior.beta <- function(s,beta,log=TRUE){
  sum(log(ddoublex(beta, mu=0, lambda=s)))
}

# prior on phi-star as specified by 3.5
prior.phi.star <- function(phi.star,tau,Q.s,q){
  (q/2)*log(tau)-(tau/2)*t(phi.star)%*%Q.s%*%phi.star
}