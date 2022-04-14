ln_Beta_Binomial<-function(x,N,M,a){
  a[a<=1e-3]=1e-3
  n=N
  b=a*(N-M)/M
  L= (lgamma(n + 1)-lgamma(x + 1)-lgamma(n - x + 1)) + lbeta((a + x),(b + n - x)) - lbeta(a,b)

  L[(n*M==0) & (x==0)]=0

  return(L)
}


My_BBin_Dist2<-function(x, N, p, a){
  # Calculates
  #
  #  -log( 0.5 BBin(x,N,p,a) + sum_{k=x+1}^N BBin(y,N,p,a) )
  #
  #   where BBin(x,N,p,a) is the probaility of getting x from a beta binomial
  #   with N events, a mean p*N and a skewness parameter a.
  statement = (N*p==0) & (x==0)
  if (any(is.na(statement))){
    Z = NaN
  }else if (all(statement)){
    Z=-log(0.5)
  }else{
    if (is.null(a)){
      # this value works well for England.
      P4=c(2.3703, 3.9310, -0.9484, 4.0321)
      a=P4[1]+P4[2]/(1+exp(P4[3]*(N*p-P4[4])))
    }
    #sum over k, work out the log Beta Binomials and then add up their
    #exponentials, this seems to work better if there are extreme values.
    k=x[1]:N[1]
    L=ln_Beta_Binomial(k,N,p*N,a)
    L[1]=L[1] + log(0.5)

    mL=max(L)
    L=L-mL
    Z= - mL - log(sum(exp(L)))
  }
  return(Z)
}
