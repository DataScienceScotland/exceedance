#!/usr/bin/env Rscript
#setwd("C:/Work/COVID/Exceedance")
setwd("H://Exceedance_modelling//Exceedance")
rm(list = ls())

inputfile = commandArgs(trailingOnly=TRUE)
if (length(inputfile) == 0){
  inputfile = "LocalAuth Time Series 17092020.xlsx"
}

source("R/Read_In_Data.R")
read(inputfile)
load("Scottish_Data.RData")

source("./R/functions.R")
source("./R/Get_Exceedance.R")
source("./R/Get_Expected_Values.R")

####

Numbers=Positives+Negatives

T.index=160:(dim(Positives)[1] - 1) # Remove the last day again, give it a whirl...

D=as.matrix(Positives[T.index,])
N=as.matrix(Numbers[T.index,])

EP=Get_Expected_Values(D, N)

PP = matrix(
      apply(expand.grid(1:nrow(D), 1:ncol(D)),
           1,
           function(x){
              My_BBin_Dist2(D[x['Var1'], x['Var2']], N[x['Var1'], x['Var2']], EP[x['Var1'], x['Var2']], c())           }),
      ncol=ncol(D))

colnames(PP) = colnames(D)
rownames(PP) = rownames(D)

Todays_Exceedance=PP
end = dim(PP)[1]
Todays_Order = y = sort(PP[end,], decreasing = TRUE)
i = order(PP[end,], decreasing = TRUE)

# highest 5 exceedances
Q=sort(i[1:5])
TodaysQ=Q;

# highest 5 expected values
i = order(EP[end,], decreasing = TRUE)
H=sort(i[1:5])
TodaysH=H
#


#CREATE COLORMAP
Cb=array(dim=c(100,50,3))
COL=matrix(data=c(0,0,0,0,0,0.0256,0,0,0.0513,0,0,0.0769,0,0,0.1026,0,0,0.1282,0,0,0.1538,0,0,0.1795,0,0,0.2051,0,0,0.2308,0,0,0.2564,0,0,0.2821,0,0,0.3077,0,0,0.3333,0,0,0.3590,0,0,0.3846,0,0,0.4103,0,0,0.4359,0,0,0.4615,0,0,0.4872,0,0,0.5128,0,0,0.5385,0,0,0.5641,0,0,0.5897,0,0,0.6154,0,0,0.6410,0,0,0.6667,0,0,0.6923,0,0,0.7179,0,0,0.7436,0,0,0.7692,0,0,0.7949,0,0,0.8205,0,0,0.8462,0,0,0.8718,0,0,0.8974,0,0,0.9231,0,0,0.9487,0,0,0.9744,0,0,1.0000,0.0167,0,0.9833,0.0333,0,0.9667,0.0500,0,0.9500,0.0667,0,0.9333,0.0833,0,0.9167,0.1000,0,0.9000,0.1167,0,0.8833,0.1333,0,0.8667,0.1500,0,0.8500,0.1667,0,0.8333,0.1833,0,0.8167,0.2000,0,0.8000,0.2167,0,0.7833,0.2333,0,0.7667,0.2500,0,0.7500,0.2667,0,0.7333,0.2833,0,0.7167,0.3000,0,0.7000,0.3167,0,0.6833,0.3333,0,0.6667,0.3500,0,0.6500,0.3667,0,0.6333,0.3833,0,0.6167,0.4000,0,0.6000,0.4167,0,0.5833,0.4333,0,0.5667,0.4500,0,0.5500,0.4667,0,0.5333,0.4833,0,0.5167,0.5000,0,0.5000,0.5167,0,0.4833,0.5333,0,0.4667,0.5500,0,0.4500,0.5667,0,0.4333,0.5833,0,0.4167,0.6000,0,0.4000,0.6167,0,0.3833,0.6333,0,0.3667,0.6500,0,0.3500,0.6667,0,0.3333,0.6833,0,0.3167,0.7000,0,0.3000,0.7167,0,0.2833,0.7333,0,0.2667,0.7500,0,0.2500,0.7667,0,0.2333,0.7833,0,0.2167,0.8000,0,0.2000,0.8167,0,0.1833,0.8333,0,0.1667,0.8500,0,0.1500,0.8667,0,0.1333,0.8833,0,0.1167,0.9000,0,0.1000,0.9167,0,0.0833,0.9333,0,0.0667,0.9500,0,0.0500,0.9667,0,0.0333,0.9833,0,0.0167,1.0000,0,0),
           ncol=3,
           byrow = TRUE)
MAX_Col=6 # maximum exceedance
MAX_bright=0.005 # maximum expected proportion positive

for (v in 1:100){
  for (b in 1:50){
    B = 0.3+0.7*b/50
    C = COL[v,]
    C = 1-B*(1-C)
    Cb[v,b,1:3]= round(C*255)
  }
}


XT = c(0, 25, 50)
XTL = sapply(XT, function(x){paste0(x*MAX_bright*100/50,'%')})

YT= seq(0,MAX_Col,1)*100/MAX_Col
YTL=sapply(YT, function(x){paste0(x*MAX_Col/100) })

# Generate colours for each local authority
end = dim(Todays_Exceedance)[1]
V=round(1+99*Todays_Exceedance[end,] / MAX_Col)
V[V>100]=100
C=COL[V,]
B=0.3+0.7*EP[end,]/MAX_bright
B[B>1]=1;

C=1-(B %*% t(rep(1,3))) * (1-C)


source("R/Get_Map.R")
source("R/PLOT_SCOTTISH_MAP.R")

date_from_PP = rownames(PP)
cat("Plot maps of todays exceendance.\n")
for (q in 25:length(T.index)){
  date_string = date_from_PP[q]
  DD = paste0('Todays_Exceedance_Map_', date_string, '.png')
  if (!file.exists(DD)){ #if map does not exist
    png(DD, width = 3.25,
        height    = 3.25,
        units     = "in",
        res       = 1200,
        pointsize = 4)
    par(fig=c(0,1,0,1))
    V=round(1+99*PP[q,]/MAX_Col)
    V[V>100]=100
    C=COL[V,]

    B=0.3+0.7*EP[q,] / MAX_bright
    B[B>1]=1

    C=1-(B %*% t(rep(1,3))) * (1-C)
    C1 = apply(C / max(C), 1, rgb2hex)

    PLOT_SCOTTISH_MAP(PP[q,], C1, -log(c(0.05, 0.01)))

    text(-6e5, 8.5e6, "Today's Exceedance", cex=2)
    text(-6e5,8.4e6,
         date_string,
         cex=1.5)
    par(fig=c(0.8, 0.88, 0.65, 0.9), new=T)
    par(mar=c(0, 0, 0, 0))

    plot(c(1,dim(Cb)[2]), c(1, dim(Cb)[1]), yaxt="n", xaxt='n', col='white', frame.plot=F)
    box(lwd=0.1)
    # rasterImage( Cb / max(Cb), 1, 1, 50, 100)
    rasterImage( Cb[seq(100,1,-1),,] / max(Cb), 1, 1, 50, 100)
    axis(1, at = XT, labels = XTL, las = 1, cex.axis=0.9, lwd=0.1)
    axis(4, at = YT, labels = YTL, las = 1, cex.axis=0.8, lwd=0.1)
    mtext(expression("Expected percent.\n  of positives"), 1, cex=1, padj=4.1)
    mtext("Exceedance", 4, padj=2.1)

    dev.off()
  }
}


#  Get the Historical Cumulative Exceedance - just update if the data is there.
status=tryCatch({
	load("BBnin_Data_to_Date.Rdata")
  	last_time=nrow(Historical_worstP)
 	},
  error=function(e){
    writeLines("creating data for BBnin_Data_to_Date.Rdata")
    return(0)
  }
)
if (status == 0){
  Historical_worstP = matrix(0, nrow = nrow(D), ncol = ncol(D))
  last_time=14;
}
if (nrow(Historical_worstP) < nrow(D)){
  Historical_worstP_tmp = matrix(0, nrow = nrow(D), ncol = ncol(D))
  Historical_worstP_tmp[1:nrow(Historical_worstP),] = Historical_worstP
  Historical_worstP = Historical_worstP_tmp
  Historical_worstP_tmp = NULL

}

for (i in seq2((last_time+1), nrow(D))) {
  writeLines(paste('Updating up-to time', as.character(i)))

  EP=Get_Expected_Values(D[1:i, ], N[1:i,])

  ep=EP
  d=D[1:i,]
  n=N[1:i,]

  worst = Get_Exceedance(d, n, ep)
  worstT = worst$worstT
  worstP = worst$worstP
  Historical_worstP[i,]=worstP;
}

save("Historical_worstP", file = "BBnin_Data_to_Date.RData")


# Make sure we have today's values loaded in.
EP=Get_Expected_Values(D, N)
worst = Get_Exceedance(D, N, EP);
worstT = worst$worstT
worstP = worst$worstP

# worst last week and this week.
end = nrow(Historical_worstP)

idx=(end-13):(end-7)
M1=apply(Historical_worstP[idx,], 2, max)
idx=(end-6):end
M2=apply(Historical_worstP[idx,], 2, max)

# Note Exceedance of 3 is approximately p=0.05, Exceedance of 4.6 is
# approximately p=0.01;

cat("Plot summary of recent increasers.\n")

# Select some "bad" locations:
m=which(
  (M2>3) & (M2>M1) & (Historical_worstP[end,]>4.6) | (M2>4.6) & (Historical_worstP[end,]>3))

o=order(Historical_worstP[end,m],decreasing = TRUE)
m=m[o]
Y=y=Historical_worstP[end,m][o]
To_Watch=m

png("Recent_increasers.png",
    width     = 4.0,
    height    = 3.25,
    units     = "in",
    res       = 1000,
    pointsize = 4)

plot(c(15, nrow(Historical_worstP[,m]) + 20),
     range(Historical_worstP[,m]) + 0.99,
     type = 'n',
     log='y',  cex.axis=1.5,
     ylab='Recent cumulative Exceedance', cex.lab=1.5,
     xlab='time', xaxt="n",
     main='Scottish Local Authorities with High Exceedance')
axis(1, at = seq(10,80,10), date_from_PP[seq(10,80,10)])
palette=rep(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999",
	"#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),10)

for (i in 1:length(m)){
  if(i<=10){
    lwd=1.9
    font=2
    lwd2 = 1
  }else{
    lwd=0.9
    font=1
    lwd2 = 0.5
  }
  lines(Historical_worstP[,m[i]], lwd=lwd, col=palette[i])

  x0 = length(Historical_worstP[,m[i]])
  y0 = Historical_worstP[x0,m[i]]
  x1 = x0 + 7
  y1 = ifelse(i==1, y0, min(y0, y1 * 0.88))
  arrows(x0, y0, x1, y1, code = 2, col = palette[i], lwd=lwd2, length=0.01)
  text(x1+0.5, y1, Place_names[m[i]], col=palette[i], cex=1, adj = 0, font = font)
}

dev.off()

# NOW DRAW MAP OF RECENT CUMULATIVE EXCEEDANCE.
cat("Draw maps of recent cumulative exceedances.\n")



MAX_Col=10
Max_bright=0.005

for (v in 1:100){
  for (b in 1:50){
    B=0.3+0.7*b/50
    C=COL[v,]
    C=1 - B * (1-C)
    Cb[v,b,1:3] = round(C * 255)
  }
}



XT = c(0, 25, 50)
XTL = sapply(XT, function(x){paste0(x*MAX_bright*100/50,'%')})


YT= seq(0,MAX_Col,2)*100/MAX_Col
YTL=sapply(YT, function(x){paste0(x*MAX_Col/100) })



for (q in 25:length(T.index)){
  date_string = date_from_PP[q]
  DD = paste0('Recent_Exceedance_Map_', date_string, '.png')
  #
  if (!file.exists(DD)){ #if map does not exist
    png(DD, width     = 3.25,
        height    = 3.25,
        units     = "in",
        res       = 1200,
        pointsize = 4)
    par(fig=c(0,1,0,1))

    V=round(1+99*Historical_worstP[q,]/MAX_Col)
    V[V>100]=100
    C=COL[V,]

    B=0.3+0.7*EP[q,] / MAX_bright
    B[B>1]=1

    C=1-(B %*% t(rep(1,3))) * (1-C)
    C1 = apply(C / max(C), 1, rgb2hex)

    # higher exceedance thresholds as we find maximum across time.
    PLOT_SCOTTISH_MAP(Historical_worstP[q,], C1, c(6,8))

    text(-6e5, 8.5e6, "Recent Cumulative Exceedance", cex=2)
    text(-6e5,8.4e6,
         date_string,
         cex=1.5)
    par(fig=c(0.8, 0.88, 0.65, 0.9), new=T)
    par(mar=c(0, 0, 0, 0))

    plot(c(1,dim(Cb)[2]), c(1, dim(Cb)[1]), yaxt="n", xaxt='n', col='white', frame.plot=F)
    box(lwd=0.1)
    # rasterImage(Cb / max(Cb), 1, 1, 50, 100)
    rasterImage(Cb[seq(100,1,-1),,] / max(Cb), 1, 1, 50, 100)
    axis(1, at = XT, labels = XTL, las = 1, cex.axis=0.9, lwd=0.1)
    axis(4, at = YT, labels = YTL, las = 1, cex.axis=0.8, lwd=0.1)
    mtext(expression("Expected percent.\n  of positives"), 1, cex=1, padj=4.1)
    mtext("Exceedance", 4, padj=2.1)

    dev.off()
  }
}


## Worst_LAs_
cat("Plot worst LAs.\n")

T.index = 1:nrow(EP)
XT=length(T.index) + seq(-42,0,7)
XTL = rownames(EP)[XT]
rangeXT = range(XT)
for (Z in 0:0){
  i = order(worstP, decreasing = TRUE)
  y = worstP[i]
  O = i[1:10+Z*8]
  DD = paste0('Worst_LAs_',  as.character(10*Z+1),'_', as.character(10*Z+10), '.png')

  png(DD, width     = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)

  par(mfrow=c(4, 2))

  B = matrix(NA, ncol=2, nrow=nrow(N))

  for (i in 1:8){
      main = paste0(
        Place_names[O[i]],
        "(Exceedance=",
        as.character(0.01*round(100*y[i+Z*10])),
        ", Expected Proportion Positive=",
        as.character(0.01*round(100*EP[nrow(PP),O[i]]*100)),
        "%)")

      X=D[,O[i]] / N[,O[i]]
      for (j in 1:nrow(N)){
        PPtmp = pbinom( 0:N[j,O[i]], N[j,O[i]], EP[j,O[i]])

        B[j, 1]= ( which(PPtmp>0.025)[1] - 1) / N[j,O[i]]
        B[j, 2]= ( which(PPtmp>0.975)[1] - 1) / N[j,O[i]]
      }

      plot(rangeXT, c(0, max(B[rangeXT[1]:rangeXT[2],2], na.rm=T)) * 100,
           ylab = '% Pillar2 positive', xlab = NA, type='n',
           yaxt="n", xaxt="n",  frame.plot=F,
           main=main)
      box(lwd=0.1)

      polygon( c(T.index, rev(T.index)),   100 * c(B[,1], rev(B[,2])),
               col =rgb(0.8, 0.8, 0.8), border=NA)

      lines(T.index, 100* B[,1], lwd=0.1)
      lines(T.index, 100* B[,2], lwd=0.1)

      lines(T.index, 100*X, col='red', lwd=0.4)
      points(T.index[worstT[O[i]]:length(T.index)],
            100 * X[worstT[O[i]]:length(T.index)], col='red', cex=1.2, pch=16)

      lines(T.index, 100*EP[,O[i]], lwd=0.1)
      points(T.index, 100*EP[,O[i]], pch=16, cex=0.9)

      axis(1, at = XT, labels = XTL, las = 1, cex.axis=0.9, lwd=0.1)
      axis(2, las = 1, cex.axis=0.8, lwd=0.1)
  }
  dev.off()
}


# Worst today LA PP
cat("Plot worst today LA PPs.\n")

T.index = 1:nrow(PP)
XT=length(T.index) + seq(-21,0,3)
XTL = rownames(PP)[XT]
rangeXT = range(XT)

for (Z in 0:0){
  i = order(PP[nrow(PP),], decreasing = TRUE)
  y = PP[nrow(PP),][i]
  O = i[1:10+Z*8]
  DD = paste0('Worst_Today_LAs_PP', as.character(10*Z+1), '_', as.character(10*Z+10), '.png')

  png(DD, width = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)

  par(mfrow=c(4, 2))

  B = matrix(NA, ncol=2, nrow=nrow(N))

  for (i in 1:8){
    if (y[i]>2){
      main = paste0(
        Place_names[O[i]],
        "(Exceedance=",
        as.character(0.01*round(100*y[i+Z*10])),
        ", Expected Proportion Positive=",
        as.character(0.01*round(100*EP[nrow(PP),O[i]]*100)),
        "%)")

      X=D[,O[i]] / N[,O[i]]
      for (j in 1:nrow(N)){
        PPtmp = pbinom( 0:N[j,O[i]], N[j,O[i]], EP[j,O[i]])

        B[j, 1]= ( which(PPtmp>0.025)[1] - 1) / N[j,O[i]]
        B[j, 2]= ( which(PPtmp>0.975)[1] - 1) / N[j,O[i]]
      }

      plot(rangeXT, c(0, max(B[rangeXT[1]:rangeXT[2],2], na.rm=T)) * 100,
           ylab = '% Pillar2 positive', xlab = NA, type='n',
           yaxt="n", xaxt="n",  frame.plot=F,
           main=main)
      box(lwd=0.1)

      polygon( c(T.index, rev(T.index)),   100 * c(B[,1], rev(B[,2])),
               col =rgb(0.8, 0.8, 0.8), border=NA)

      lines(T.index, 100* B[,1], lwd=0.1)
      lines(T.index, 100* B[,2], lwd=0.1)

      lines(T.index, 100*X, col='blue', lwd=0.4)
      points(T.index[worstT[O[i]]:length(T.index)],
             100 * X[worstT[O[i]]:length(T.index)], col='blue', cex=1.2, pch=16)

      lines(T.index, 100*EP[,O[i]], lwd=0.1)
      points(T.index, 100*EP[,O[i]], pch=16, cex=0.9)



      axis(1, at = XT, labels = XTL, las = 1, cex.axis=0.9, lwd=0.1)
      axis(2, las = 1, cex.axis=0.8, lwd=0.1)

    }
  }
  dev.off()
}

# Worst today LA EP
cat("Plot worst today LA EP.\n")

T.index = 1:nrow(EP)
XT=length(T.index) + seq(-21,0,3)
XTL = rownames(EP)[XT]
rangeXT = range(XT)

for (Z in 0:0){
  i = order(EP[nrow(EP),], decreasing = TRUE)
  y = EP[nrow(EP),][i]
  O = i[1:10+Z*8]
  DD = paste0('Worst_Today_LAs_EP', as.character(10*Z+1), '_', as.character(10*Z+10), '.png')

  png(DD, width = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)

  par(mfrow=c(4, 2))

  B = matrix(NA, ncol=2, nrow=nrow(N))

  for (i in 1:8){
    if (y[i]>0.0001){
      main = paste0(
        Place_names[O[i]],
        "(Exceedance=",
        as.character(0.01*round(100*y[i+Z*10])),
        ", Expected Proportion Positive=",
        as.character(0.01*round(100*EP[nrow(PP),O[i]]*100)),
        "%)")

      X=D[,O[i]] / N[,O[i]]
      for (j in 1:nrow(N)){
        PPtmp = pbinom( 0:N[j,O[i]], N[j,O[i]], EP[j,O[i]])
        # PPtmp=cdf('Bin',[0:N(O(i),j)],N(O(i),j),EP(O(i),j))

        B[j, 1]= ( which(PPtmp>0.025)[1] - 1) / N[j,O[i]]
        B[j, 2]= ( which(PPtmp>0.975)[1] - 1) / N[j,O[i]]
      }

      plot(rangeXT, c(0, max(B[rangeXT[1]:rangeXT[2],2], na.rm=T)) * 100,
           ylab = '% Pillar2 positive', xlab = NA, type='n',
           yaxt="n", xaxt="n",  frame.plot=F,
           main=main)
      box(lwd=0.1)

      polygon( c(T.index, rev(T.index)),   100 * c(B[,1], rev(B[,2])),
               col =rgb(0.8, 0.8, 0.8), border=NA)

      lines(T.index, 100* B[,1], lwd=0.1)
      lines(T.index, 100* B[,2], lwd=0.1)

      lines(T.index, 100*X, col='blue', lwd=0.4)
      points(T.index[worstT[O[i]]:length(T.index)],
             100 * X[worstT[O[i]]:length(T.index)], col='blue', cex=1.2, pch=16)

      lines(T.index, 100*EP[,O[i]], lwd=0.1)
      points(T.index, 100*EP[,O[i]], pch=16, cex=0.9)



      axis(1, at = XT, labels = XTL, las = 1, cex.axis=0.9, lwd=0.1)
      axis(2, las = 1, cex.axis=0.8, lwd=0.1)

    }
  }
  dev.off()
}



fprintf <- function(file, ...){
  # fprintf(file, "Hello World %d", 10)
  writeLines(sprintf(...), file)
}


## Recent List tex
cat("Printing recent lists in `.tex`\n")

if (file.exists("Recent_List.tex")){
  status = file.rename("Recent_List.tex", "Recent_List_old.tex");
}
fp<-file("Recent_List.tex", open='a')


end = nrow(Historical_worstP)
idx=(end-13):(end-7)
M1=apply(Historical_worstP[idx,], 2, max)
idx=(end-6):end
M2=apply(Historical_worstP[idx,], 2, max)


# Note Exceedance of 3 is approximately p=0.05, Exceedance of 4.6 is
# approximately p=0.01;

# Select some "bad" locations:
m=which(
  (M2>3) & (M2>M1) & (Historical_worstP[end,]>4.6) | (M2>4.6) & (Historical_worstP[end,]>3))

o=order(Historical_worstP[end,m],decreasing = TRUE)
m=m[o]
Y=y=Historical_worstP[end,m][o]

if (length(m)>30){
  m = m[1:30]
}
if (length(m)>1){
  fprintf(fp, 'Recent cumulative exceedance values highlight the following LAs as being high risk:\\\\')
}
if (length(m)==1){
  fprintf(fp,'Recent cumulative exceedance values highlight the following LA as being high risk:\\\\')
}
if (length(m)>=1){
  for (i in (1:length(m))){
    if (i<=10){
      fprintf(fp,'$\\bullet$ {\\bf %s} (exceedance=%g) -', Place_names[m[i]], 0.01*round(100*worstP[m[i]]))
    }else{
      fprintf(fp,'$\\bullet$ %s (exceedance=%g) -', Place_names[m[i]], 0.01*round(100*worstP[m[i]]))
    }

    # 2 week trend
    n=m[i]
    HwP=Todays_Exceedance[nrow(Todays_Exceedance)-13:0,n]
    tt=0:13
    mdl = lm(HwP~tt)
    ci=confint(mdl)

    flag=1
    if (ci[2,1] > 1){
       fprintf(fp,' increasing dramatically'); flag=0;
    }
    if (flag & ci[2,1] > 0.5){
       fprintf(fp,' increasing quickly'); flag=0;
    }
    if (flag & ci[2,1] > 0){
       fprintf(fp,' increasing slowly'); flag=0;
    }
    if (flag & ci[2,2] < -0.5){
       fprintf(fp,' decreasing'); flag=0;
    }
    if (flag & ci[2,2]<0){
       fprintf(fp,' decreasing slowly'); flag=0
    }
    if (flag){
      fprintf(fp,' no significant recent trend')
    }

    if (i == length(m)){
      fprintf(fp,'.\\\\ \n')
    }else{
      fprintf(fp,';\\\\ \n')
    }
  }
}else{
  fprintf(fp,'Recent cumulative exceedance values highlight that no LAs are currently high risk.\\');
}
close(fp);

Keepm=m


if (file.exists("Today_List.tex")){
  status = file.rename("Today_List.tex", "Today_List_old.tex")
}
fp<-file("Today_List.tex", open='a')

end = nrow(PP)
o=order(PP[end, ], decreasing = TRUE)
y=PP[end,][o]
m=o[1:10]
idx = y[1:10] < 3
m = m[!idx]

if (length(m)>1){
  fprintf(fp,"Today's exceedance values highlight the following LAs as being high risk:\\\\");
}
if (length(m)==1){
  fprintf(fp,"Today''s exceedance values highlight the following LA as being high risk:\\\\");
}
if (length(m) >= 1){
  for (i in 1:(length(m)-1)){
    fprintf(fp,"$\\bullet$ %s (exceedance=%g);\\\\ \n",Place_names[m[i]], 0.01*round(100*PP[end, m[i]]))
  }
  i=length(m);
  fprintf(fp,"$\\bullet$ %s (exceedance=%g).\\\\ \n", Place_names[m[i]], 0.01*round(100*PP[end, m[i]]))
}else{
  fprintf(fp,"Today's exceedance values highlight that no LAs are currently high risk.\\");
}
close(fp);




## Combined Plots
cat("Plot combined summaries.\n")

m=sort(unique(c(m, Keepm)))
end = nrow(Historical_worstP)
o=order(Historical_worstP[end, m], decreasing = TRUE)
y=Historical_worstP[end, m][o]
m = m[o]

XT=length(T.index) + seq(-55,0,3) # plot the last 55 days
XTL = rownames(EP)[XT]

for (Z in 1:ceiling(length(m)/6)){
  DD = paste0('Combined_Plots', as.character(Z), '.png')
  png(DD, width = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)

  par(mfrow=c(6, 1), mar=c(2, 4, 2, 5) + 0.3)
  n=1:6+(Z-1)*6
  idx = n>length(m)
  n = n[!idx]
  for (i in 1:length(n)){
    z=m[n[i]]

    t1=rep(NaN, length(T.index))
    idx = (Todays_Exceedance[,z]>3) & (Todays_Exceedance[,z]<4.6)
    t1[idx]=1
    T1=rep(NaN, length(T.index))
    idx = Todays_Exceedance[,z]>4.6
    T1[idx]=1

    t2=rep(NaN, length(T.index))
    idx = (Historical_worstP[,z] > 6) & (Historical_worstP[,z] < 8)
    t2[idx]=1
    T2=rep(NaN, length(T.index))
    idx = (Historical_worstP[,z]>8)
    T2[idx]=1
    main = paste0(Place_names[z], '(Daily Exceedance=', as.character(0.01*round(100*Todays_Exceedance[end,z])),
                  ', Recent Cumulative Exceedance=', as.character(0.01*round(100*Historical_worstP[end,z])),
                  ', Expected Proportion Positive=', as.character(0.01*round(100*EP[end, z]*100)), '%)')

    plot(range(XT), c(0, max(c(Historical_worstP, Todays_Exceedance), na.rm=T)),
         type='n',
         main=main, xlab=NA, ylab=NA,
         yaxt="n", xaxt="n",  frame.plot=F, xlim=c(max(XT)-21, max(XT)+0.2))
    box(lwd=0.4)
    lines(T.index, Todays_Exceedance[,z], col='blue', lwd=0.8)
    points(t1 * T.index, t1*Todays_Exceedance[,z], col='blue', pch=16, cex=1.0)
    points(T1 * T.index, T1*Todays_Exceedance[,z], col='blue', pch=17, lwd=0.3, cex=1.3)
    lines(x=c(max(XT)-28, max(XT)+2), y=c(3,3), lty=2, col = 4, lwd = 0.4)
    lines(x=c(max(XT)-28, max(XT)+2), y=c(6,6), lty=2, col = 2, lwd=0.4)
    par(new=TRUE)

    lines(T.index, Historical_worstP[,z], col='red', lwd=0.8)
    points(t2 * T.index, t2*Historical_worstP[,z], col='red', pch=16, cex=1.0)
    points(T2 * T.index, T2*Historical_worstP[,z],  col='red', pch=17, lwd=0.3, cex=1.3)

    axis(1, at = XT, labels = XTL, las = 1, cex.axis=0.9, lwd=0.4)
    axis(2, las = 1, cex.axis=0.8, lwd=0.1, col='blue')
    axis(4, las = 1, cex.axis=0.8, lwd=0.1, col='red')

    mtext(expression("Recent Cumulative\n  Exceedance"), 4, padj=3.6, cex=0.9, col='red')
    mtext(expression("  Daily\nExceedance"), 2, padj=-1.5, cex=0.9, col='blue')
  }
  dev.off()
}

############################
##  END OF NORMAL CODE    ##
##  START OF EXTRA PLOTS  ##
############################




#CREATE COLORMAP
#Cb=array(dim=c(100,50,3))
COL=matrix(data=c(0,0,0,0,0,0.0256,0,0,0.0513,0,0,0.0769,0,0,0.1026,0,0,0.1282,0,0,0.1538,0,0,0.1795,0,0,0.2051,0,0,0.2308,0,0,0.2564,0,0,0.2821,0,0,0.3077,0,0,0.3333,0,0,0.3590,0,0,0.3846,0,0,0.4103,0,0,0.4359,0,0,0.4615,0,0,0.4872,0,0,0.5128,0,0,0.5385,0,0,0.5641,0,0,0.5897,0,0,0.6154,0,0,0.6410,0,0,0.6667,0,0,0.6923,0,0,0.7179,0,0,0.7436,0,0,0.7692,0,0,0.7949,0,0,0.8205,0,0,0.8462,0,0,0.8718,0,0,0.8974,0,0,0.9231,0,0,0.9487,0,0,0.9744,0,0,1.0000,0.0167,0,0.9833,0.0333,0,0.9667,0.0500,0,0.9500,0.0667,0,0.9333,0.0833,0,0.9167,0.1000,0,0.9000,0.1167,0,0.8833,0.1333,0,0.8667,0.1500,0,0.8500,0.1667,0,0.8333,0.1833,0,0.8167,0.2000,0,0.8000,0.2167,0,0.7833,0.2333,0,0.7667,0.2500,0,0.7500,0.2667,0,0.7333,0.2833,0,0.7167,0.3000,0,0.7000,0.3167,0,0.6833,0.3333,0,0.6667,0.3500,0,0.6500,0.3667,0,0.6333,0.3833,0,0.6167,0.4000,0,0.6000,0.4167,0,0.5833,0.4333,0,0.5667,0.4500,0,0.5500,0.4667,0,0.5333,0.4833,0,0.5167,0.5000,0,0.5000,0.5167,0,0.4833,0.5333,0,0.4667,0.5500,0,0.4500,0.5667,0,0.4333,0.5833,0,0.4167,0.6000,0,0.4000,0.6167,0,0.3833,0.6333,0,0.3667,0.6500,0,0.3500,0.6667,0,0.3333,0.6833,0,0.3167,0.7000,0,0.3000,0.7167,0,0.2833,0.7333,0,0.2667,0.7500,0,0.2500,0.7667,0,0.2333,0.7833,0,0.2167,0.8000,0,0.2000,0.8167,0,0.1833,0.8333,0,0.1667,0.8500,0,0.1500,0.8667,0,0.1333,0.8833,0,0.1167,0.9000,0,0.1000,0.9167,0,0.0833,0.9333,0,0.0667,0.9500,0,0.0500,0.9667,0,0.0333,0.9833,0,0.0167,1.0000,0,0),
           ncol=3,
           byrow = TRUE)

COL <- COL[c(rep(16,17), rep(33,17), rep(50,16), rep(66,17), rep(82,16), rep(99,17)),]
MAX_Col=6 # maximum exceedance
MAX_bright=0.005 # maximum expected proportion positive

for (v in 1:100){
  for (b in 1:50){
    B = 0.3+0.7*b/50
    C = COL[v,]
    C = 1-B*(1-C)
    Cb[v,b,1:3]= round(C*255)
  }
}


XT = c(0, 25, 50)
XTL = sapply(XT, function(x){paste0(x*MAX_bright*100/50,'%')})

YT= seq(0,MAX_Col,1)*100/MAX_Col
YTL=sapply(YT, function(x){paste0(x*MAX_Col/100) })

# Generate colours for each local authority
end = dim(Todays_Exceedance)[1]
V=round(1+99*Todays_Exceedance[end,] / MAX_Col)
V[V>100]=100
C=COL[V,]
B=0.3+0.7*EP[end,]/MAX_bright
B[B>1]=1;

C=1-(B %*% t(rep(1,3))) * (1-C)


source("R/Get_Map.R")
source("R/PLOT_SCOTTISH_MAP.R")
q <- dim(PP)[1]
date_from_PP = rownames(PP)
cat("Plot maps of todays exceendance.\n")
#for (q in 25:length(T.index)){
date_string = date_from_PP[q]
DD = paste0('TEST_Todays_Exceedance_Map_', date_string, '.png')
#  if (!file.exists(DD)){ #if map does not exist
png(DD, width = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 1200,
    pointsize = 4)
par(fig=c(0,1,0,1))
V=round(1+99*PP[q,]/MAX_Col)
V[V>100]=100
C=COL[V,]

B=0.3+0.7*EP[q,] / MAX_bright
B[B>1]=1

#    C=1-(B %*% t(rep(1,3))) * (1-C)
C1 = apply(C , 1, rgb2hex)

PLOT_SCOTTISH_MAP(PP[q,], C1, -log(c(0.05, 0.01)))

#    text(-6e5, 8.5e6, "Today's Exceedance", cex=2)
#    text(-6e5,8.4e6,
#         date_string,
#         cex=1.5)
par(fig=c(0.8, 0.88, 0.65, 0.9), new=T)
par(mar=c(0, 0, 0, 0))

plot(x=c(0,50), y=c(0,100), type="l", yaxt="n", xaxt='n', col='white', frame.plot=F)
box(lwd=0.1)

for(i in (1:99)){ 
  lines(x=c(0,50), y=c(i,i), lwd=2, col=rgb2hex(COL[i,]), lend = 1)
}
axis(4, at = YT, labels = YTL, las = 1, cex.axis=0.8, lwd=0.1)
mtext("Exceedance", 4, padj=2.1)

dev.off()
#  }
#}




# NOW DRAW MAP OF RECENT CUMULATIVE EXCEEDANCE.
cat("Draw maps of recent cumulative exceedances.\n")

COL=matrix(data=c(0,0,0,0,0,0.0256,0,0,0.0513,0,0,0.0769,0,0,0.1026,0,0,0.1282,0,0,0.1538,0,0,0.1795,0,0,0.2051,0,0,0.2308,0,0,0.2564,0,0,0.2821,0,0,0.3077,0,0,0.3333,0,0,0.3590,0,0,0.3846,0,0,0.4103,0,0,0.4359,0,0,0.4615,0,0,0.4872,0,0,0.5128,0,0,0.5385,0,0,0.5641,0,0,0.5897,0,0,0.6154,0,0,0.6410,0,0,0.6667,0,0,0.6923,0,0,0.7179,0,0,0.7436,0,0,0.7692,0,0,0.7949,0,0,0.8205,0,0,0.8462,0,0,0.8718,0,0,0.8974,0,0,0.9231,0,0,0.9487,0,0,0.9744,0,0,1.0000,0.0167,0,0.9833,0.0333,0,0.9667,0.0500,0,0.9500,0.0667,0,0.9333,0.0833,0,0.9167,0.1000,0,0.9000,0.1167,0,0.8833,0.1333,0,0.8667,0.1500,0,0.8500,0.1667,0,0.8333,0.1833,0,0.8167,0.2000,0,0.8000,0.2167,0,0.7833,0.2333,0,0.7667,0.2500,0,0.7500,0.2667,0,0.7333,0.2833,0,0.7167,0.3000,0,0.7000,0.3167,0,0.6833,0.3333,0,0.6667,0.3500,0,0.6500,0.3667,0,0.6333,0.3833,0,0.6167,0.4000,0,0.6000,0.4167,0,0.5833,0.4333,0,0.5667,0.4500,0,0.5500,0.4667,0,0.5333,0.4833,0,0.5167,0.5000,0,0.5000,0.5167,0,0.4833,0.5333,0,0.4667,0.5500,0,0.4500,0.5667,0,0.4333,0.5833,0,0.4167,0.6000,0,0.4000,0.6167,0,0.3833,0.6333,0,0.3667,0.6500,0,0.3500,0.6667,0,0.3333,0.6833,0,0.3167,0.7000,0,0.3000,0.7167,0,0.2833,0.7333,0,0.2667,0.7500,0,0.2500,0.7667,0,0.2333,0.7833,0,0.2167,0.8000,0,0.2000,0.8167,0,0.1833,0.8333,0,0.1667,0.8500,0,0.1500,0.8667,0,0.1333,0.8833,0,0.1167,0.9000,0,0.1000,0.9167,0,0.0833,0.9333,0,0.0667,0.9500,0,0.0500,0.9667,0,0.0333,0.9833,0,0.0167,1.0000,0,0),
           ncol=3,
           byrow = TRUE)

COL <- COL[c(rep(10,10), rep(20,10), rep(30,10), rep(40,10), rep(50,10), rep(60,11),
             rep(70,10), rep(80,10), rep(90,10), rep(99,10)),]


MAX_Col=10
Max_bright=0.005

for (v in 1:100){
  for (b in 1:50){
    B=0.3+0.7*b/50
    C=COL[v,]
    C=1 - B * (1-C)
    Cb[v,b,1:3] = round(C * 255)
  }
}



XT = c(0, 25, 50)
XTL = sapply(XT, function(x){paste0(x*MAX_bright*100/50,'%')})


YT= seq(0,MAX_Col,2)*100/MAX_Col
YTL=sapply(YT, function(x){paste0(x*MAX_Col/100) })


q<- dim(PP)[1]
#for (q in 25:length(T.index)){
date_string = date_from_PP[q]
DD = paste0('TEST_Recent_Exceedance_Map_', date_string, '.png')
#
if (!file.exists(DD)){ #if map does not exist
  png(DD, width     = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)
  par(fig=c(0,1,0,1))
  
  V=round(99*(Historical_worstP[q,]/MAX_Col))
  V[V>100]=100
  C=COL[V,]
  
  B=0.3+0.7*EP[q,] / MAX_bright
  B[B>1]=1
  
  #  C=1-(B %*% t(rep(1,3))) * (1-C)
  C1 = apply(C , 1, rgb2hex)
  
  # higher exceedance thresholds as we find maximum across time.
  PLOT_SCOTTISH_MAP(Historical_worstP[q,], C1, c(6,8))
  
  #   text(-6e5, 8.5e6, "Recent Cumulative Exceedance", cex=2)
  #    text(-6e5,8.4e6,
  #         date_string,
  #         cex=1.5)
  par(fig=c(0.8, 0.88, 0.65, 0.9), new=T)
  par(mar=c(0, 0, 0, 0))
  
  plot(c(1,dim(Cb)[2]), c(1, dim(Cb)[1]), yaxt="n", xaxt='n', col='white', frame.plot=F)
  for(i in (1:99)){ 
    lines(x=c(0,50), y=c(i,i), lwd=2, col=rgb2hex(COL[i,]), lend = 1)
  }
  axis(4, at = YT, labels = YTL, las = 1, cex.axis=0.8, lwd=0.1)
  mtext("Exceedance", 4, padj=2.1)
  
  
  dev.off()
}


### for blue palette

temp.COL <- colorRampPalette(c("lightskyblue", "blue", "navy"))
COL <- temp.COL(100)

### for grey and blue only

COL <- c(rep("#808080", 50), rep("#0000ff", 50))


MAX_Col=6 # maximum exceedance
MAX_bright=0.005 # maximum expected proportion positive


XT = c(0, 25, 50)
XTL = sapply(XT, function(x){paste0(x*MAX_bright*100/50,'%')})

YT= seq(0,MAX_Col,1)*100/MAX_Col
YTL=sapply(YT, function(x){paste0(x*MAX_Col/100) })

# Generate colours for each local authority
end = dim(Todays_Exceedance)[1]
V=round(1+99*Todays_Exceedance[end,] / MAX_Col)
V[V>100]=100
C=COL[V]


source("R/Get_Map.R")
source("R/PLOT_SCOTTISH_MAP.R")
q <- dim(PP)[1]
date_from_PP = rownames(PP)
cat("Plot maps of todays exceendance.\n")
#for (q in 25:length(T.index)){
date_string = date_from_PP[q]
DD = paste0('TEST_Todays_Exceedance_Map_', date_string, '_(blue_pallete2).png')
#  if (!file.exists(DD)){ #if map does not exist
png(DD, width = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 1200,
    pointsize = 4)
par(fig=c(0,1,0,1))
V=round(1+99*PP[q,]/MAX_Col)
V[V>100]=100
C=COL[V]

C1 = C

PLOT_SCOTTISH_MAP(PP[q,], C1, -log(c(0.05, 0.01)))

#    text(-6e5, 8.5e6, "Today's Exceedance", cex=2)
#    text(-6e5,8.4e6,
#         date_string,
#         cex=1.5)
par(fig=c(0.8, 0.88, 0.65, 0.9), new=T)
par(mar=c(0, 0, 0, 0))

plot(x=c(0,50), y=c(0,100), type="l", yaxt="n", xaxt='n', col='white', frame.plot=F)
box(lwd=0.1)

for(i in (1:99)){ 
  lines(x=c(0,50), y=c(i,i), lwd=2, col=COL[i], lend = 1)
}
axis(4, at = YT, labels = YTL, las = 1, cex.axis=0.8, lwd=0.1)
mtext("Exceedance", 4, padj=2.1)

dev.off()
#  }
#}




# NOW DRAW MAP OF RECENT CUMULATIVE EXCEEDANCE.
cat("Draw maps of recent cumulative exceedances.\n")

#COL=matrix(data=c(0,0,0,0,0,0.0256,0,0,0.0513,0,0,0.0769,0,0,0.1026,0,0,0.1282,0,0,0.1538,0,0,0.1795,0,0,0.2051,0,0,0.2308,0,0,0.2564,0,0,0.2821,0,0,0.3077,0,0,0.3333,0,0,0.3590,0,0,0.3846,0,0,0.4103,0,0,0.4359,0,0,0.4615,0,0,0.4872,0,0,0.5128,0,0,0.5385,0,0,0.5641,0,0,0.5897,0,0,0.6154,0,0,0.6410,0,0,0.6667,0,0,0.6923,0,0,0.7179,0,0,0.7436,0,0,0.7692,0,0,0.7949,0,0,0.8205,0,0,0.8462,0,0,0.8718,0,0,0.8974,0,0,0.9231,0,0,0.9487,0,0,0.9744,0,0,1.0000,0.0167,0,0.9833,0.0333,0,0.9667,0.0500,0,0.9500,0.0667,0,0.9333,0.0833,0,0.9167,0.1000,0,0.9000,0.1167,0,0.8833,0.1333,0,0.8667,0.1500,0,0.8500,0.1667,0,0.8333,0.1833,0,0.8167,0.2000,0,0.8000,0.2167,0,0.7833,0.2333,0,0.7667,0.2500,0,0.7500,0.2667,0,0.7333,0.2833,0,0.7167,0.3000,0,0.7000,0.3167,0,0.6833,0.3333,0,0.6667,0.3500,0,0.6500,0.3667,0,0.6333,0.3833,0,0.6167,0.4000,0,0.6000,0.4167,0,0.5833,0.4333,0,0.5667,0.4500,0,0.5500,0.4667,0,0.5333,0.4833,0,0.5167,0.5000,0,0.5000,0.5167,0,0.4833,0.5333,0,0.4667,0.5500,0,0.4500,0.5667,0,0.4333,0.5833,0,0.4167,0.6000,0,0.4000,0.6167,0,0.3833,0.6333,0,0.3667,0.6500,0,0.3500,0.6667,0,0.3333,0.6833,0,0.3167,0.7000,0,0.3000,0.7167,0,0.2833,0.7333,0,0.2667,0.7500,0,0.2500,0.7667,0,0.2333,0.7833,0,0.2167,0.8000,0,0.2000,0.8167,0,0.1833,0.8333,0,0.1667,0.8500,0,0.1500,0.8667,0,0.1333,0.8833,0,0.1167,0.9000,0,0.1000,0.9167,0,0.0833,0.9333,0,0.0667,0.9500,0,0.0500,0.9667,0,0.0333,0.9833,0,0.0167,1.0000,0,0),
#           ncol=3,
#           byrow = TRUE)

#COL <- COL[c(rep(10,10), rep(20,10), rep(30,10), rep(40,10), rep(50,10), rep(60,11),
#             rep(70,10), rep(80,10), rep(90,10), rep(99,10)),]


### for blue palette

temp.COL <- colorRampPalette(c("lightskyblue", "blue", "navy"))
COL <- temp.COL(100)

### for grey and blue only

COL <- c(rep("#808080", 60), rep("#0000ff", 40))



MAX_Col=10
Max_bright=0.005


XT = c(0, 25, 50)
XTL = sapply(XT, function(x){paste0(x*MAX_bright*100/50,'%')})


YT= seq(0,MAX_Col,2)*100/MAX_Col
YTL=sapply(YT, function(x){paste0(x*MAX_Col/100) })


q<- dim(PP)[1]

date_string = date_from_PP[q]
DD = paste0('TEST_Recent_Exceedance_Map_', date_string, '_(blue_palette2).png')
#
if (!file.exists(DD)){ #if map does not exist
  png(DD, width     = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)
  par(fig=c(0,1,0,1))
  
  V=round(99*(Historical_worstP[q,]/MAX_Col))
  V[V>100]=100
  C=COL[V]
  
  #  C=1-(B %*% t(rep(1,3))) * (1-C)
  C1 = C
  
  # higher exceedance thresholds as we find maximum across time.
  PLOT_SCOTTISH_MAP(Historical_worstP[q,], C1, c(6,8))
  
  #   text(-6e5, 8.5e6, "Recent Cumulative Exceedance", cex=2)
  #    text(-6e5,8.4e6,
  #         date_string,
  #         cex=1.5)
  par(fig=c(0.8, 0.88, 0.65, 0.9), new=T)
  par(mar=c(0, 0, 0, 0))
  
  plot(c(1,dim(Cb)[2]), c(1, dim(Cb)[1]), yaxt="n", xaxt='n', col='white', frame.plot=F)
  for(i in (1:99)){ 
    lines(x=c(0,50), y=c(i,i), lwd=2, col=COL[i], lend = 1)
  }
  axis(4, at = YT, labels = YTL, las = 1, cex.axis=0.8, lwd=0.1)
  mtext("Exceedance", 4, padj=2.1)
  
  
  dev.off()
}


