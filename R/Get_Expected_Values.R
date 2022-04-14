seq2<-function(start, end){
  if (end >= start){
    return(start:end)
  }else{
    return(c())
  }
}

Get_Expected_Values<-function(D, N){
  #
  # Expect D to be array (Location, Time)

  # Set an upper threshold TH which gets ignored.
  TH=3

  # loop round to get rid of high points !
  for (loop in 1:5){
    E=D
    PP=matrix(1, nrow = nrow(D), ncol = ncol(D))
    for (L in 1:dim(D)[2]){
      # define a population trend
      Q=D
      Q = Q[,-L]
      n=N[,-L]

      # mean over all other locations across all times.
      Mn=rowMeans(Q, na.rm = TRUE)/rowMeans(n, na.rm = TRUE)

      # loop over all time points
      for (t.idx in 1:dim(D)[1]){
        # go back 6 weeks to get an average.
        # Get rid of early time.
        TT=(t.idx-42):t.idx
        TT = TT[TT>=5]

        # Get rid of current time point and anywhere with high exceedance

        Q=D[TT, L]
        idx=(TT == t.idx) | (PP[TT, L] > TH)
        Q = Q[!idx]
        n=N[TT, L]
        idx=(TT==t.idx) | (PP[TT, L]>TH)
        n =  n[!idx]

        # calculate the mean
        MM=mean(Q, na.rm = TRUE)/mean(n, na.rm = TRUE)

        if ((MM==0) | (is.na(MM))){
          # if that fails same again but across all time with low enough exceedane
          TT=seq2(5, t.idx)
          TT = TT[!(PP[TT, L]>TH)]
          MM=mean(D[TT, L], na.rm = TRUE)/mean(N[TT, L], na.rm = TRUE)

          if ((MM==0) | is.na(MM)){
            # and if that fails over all time no-matter what
            TT=seq2(5, t.idx)
            MM=mean(D[TT, L], na.rm = TRUE) / mean(N[TT, L], na.rm = TRUE)
          }
        }

        # Calculate expected value at that time.
        E[t.idx, L] = Mn[t.idx] * MM / mean(Mn, na.rm = TRUE)
      }
    }

    # Set expected value to zero if there are no swabs.
    E[N==0] = 0;

    # Calculate the log(exceedance probabilities) so that we know if times
    # need to be dicounted.

    PP = matrix(
      apply(expand.grid(1:nrow(N), 1:ncol(N)),
            1,
            function(x){
              My_BBin_Dist2(D[x['Var1'], x['Var2']], # x
                            N[x['Var1'], x['Var2']], # N
                            min(E[x['Var1'], x['Var2']], 0.9, na.rm = TRUE), # P
                            c()) # a
            }),
      ncol=ncol(N))

    colnames(PP) = colnames(D)
    rownames(PP) = rownames(D)
  }
  return(E)
}
