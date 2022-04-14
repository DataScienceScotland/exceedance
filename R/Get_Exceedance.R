Get_Exceedance<-function(D, N, EP){
# This gets the Cumulative Exceedance going back over time !

  SHORT_TIME_FRAME=7;
  worstT = vector(mode='numeric', length=ncol(D))
  worstP = vector(mode='numeric', length=ncol(D))
  for (L in 1:ncol(D)){
    PP=rep(1, nrow(D))
    end = nrow(D)
    #% First see if going back over 7 days gives anything useful!
    for (T.idx in (nrow(D):(nrow(D) - SHORT_TIME_FRAME))){
      DD = sum(D[T.idx:end, L])
      NN = sum(N[T.idx:end, L])
      EE = sum(N[T.idx:end, L] * EP[T.idx:end, L]) / NN
      PP[T.idx] = My_BBin_Dist2(DD, NN, EE, c())
    }
    # If there is a high exceedance go back further.
    if (max(PP, na.rm=T)>5){
      for (T.idx in (dim(D)[1]-SHORT_TIME_FRAME):1){
        DD = sum(D[T.idx:end, L])
        NN = sum(N[T.idx:end, L])
        EE = sum(N[T.idx:end, L] * EP[T.idx:end, L]) / NN
        PP[T.idx] = My_BBin_Dist2(DD, NN, EE, c())
      }

      # ... but truncate if it drops below the expected percentage 4
      # times.
      m = which(EP[1:end,L] * N[1:end,L] > D[1:end,L])
      l_m = length(m)
      if (l_m>4){
        PP[1:m[l_m-4]] = 0
      }
    }
    # find the highest exceedance and the duration back in time.
    i = which.max(PP)
    worstT[L]=i
    worstP[L]=PP[i];
  }
return (data.frame(worstT=worstT, worstP=worstP))
}
