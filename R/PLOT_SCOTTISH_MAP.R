rgb2hex<-function(x){
  if (any(is.na(x))){
    NA
  }else{
    r = x[1]
    g = x[2]
    b = x[3]
    rgb(r, g, b)
  }
}


PLOT_SCOTTISH_MAP<-function(tmp_V, tmp_Colors, limits){
  V=tmp_V[Map2PlaceName]
  COL=tmp_Colors[Map2PlaceName]

  m=1:32
  # White borders unless Exceedance Value > limits(1), in which case its
  # black.
  EC = matrix(0.999, ncol = 3, nrow=32)
  EC[V>limits[1],]=0
  EC1 = apply(EC, 1, rgb2hex)

  plot(M$geometry,  col=COL[m], border=EC, lwd = 0.2)

  # Find regions greater than limits(2) and thicken their borders
  m.idx = which(V > limits[2])
  for (m in m.idx){
      plot(M[m,]$geometry, add=T, lwd=0.4)
  }
}

