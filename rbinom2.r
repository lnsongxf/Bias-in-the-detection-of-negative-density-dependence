rbinom2 <- function(n,S,p){
  
  # ranadom number generator of S trials given R number success with Bernully
  # probability p
  z = runif(n)
  R=numeric(length(S))
  
  
  use=S==0
  R[use] <- 0
  
  
  S0=unique(S[S>0])
  glS1 <- lgamma(S0+1)
  if (p>0 & p<1){
  for (i in 1:length(S0)){
    
    use2=S==S0[i]
    k<-seq.int(0,S0[i])
    C <- cumsum(exp(glS1[i] - lgamma(k+1) -lgamma(S0[i]-k+1) + k*log(p) + (S0[i]-k)*log(1-p)))
    R[use2] <- .bincode(z[use2], c(0,C))-1
    
  }
  } else {
    R = S*p
  }
  
  return(R)
}



