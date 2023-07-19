#This file contains functions:
##########################
#gradient: Helper Function: Take the gradient of the likelihood function: Input w: theta, A: comparison matrix, Y1,Y2,Y3: responses denoting whether the winner is 1,2,or 3; Output: Gradient of the Function
#MLE: Helper Function: Compute the MLE of the Loss function. Input: A: comparison Matrix, Y1,Y2,Y3: responses denoting whether the winner is 1,2,or 3, Output: MLE estimator
#variance: Helper Function: compute the variance of the estimator: Input: m-th entry, w:estimated MLE, A: comparison matrix; output: values of g^{(m)} function at parameter w
#F: Helper Function, compute the function $f^{(m)}$ in paper. Input: m:mth entry, w estimator, A: comparison matrix, YY1,YY2,YY3:response variables, L number of comparisons; #output values of f^{(m)} at parameter w.
#GA.CI: Helper Function, compute the confidence intervals for item m. Input: m:mth entry, w estimator, A: comparison matrix, Y1,Y2,Y3,YY1,YY2,YY3:response variables, L number of comparisons,B:How many times of bootstrap samples, W: true theta, seq.C: c_0 values on Page 26 of the main text; Output, confidence intervals----#
#GA: Helper Function, function of Gaussian approximation. Input: m:mth entry, w estimator, A: comparison graph, Y1,Y2,Y3,YY1,YY2,YY3:response variables, L number of comparisons,B:How many times of bootstrap samples, W: true theta, para: theoretical p values; Output, indicator on whether exceeding the GA approximation threshold or not----
##########################
#-------------gradient: Helper Function: Take the gradient of the likelihood function: Input w: theta, A: comparison matrix, Y1,Y2,Y3: responses denoting whether the winner is 1,2,or 3; Output: Gradient of the Function#--------------
gradient = function(w,A,Y1,Y2,Y3){
  n = length(w)
  R1 = array(exp(w),c(n,n,n))
  R2 = aperm(R1,c(2,3,1))
  R3 = aperm(R1,c(3,1,2))
  RR1 = R1/(R1+R2+R3)
  RR2 = R3/(R1+R2+R3)
  RR3 = R2/(R1+R2+R3)
  G1 = A*(Y1-RR1); G2 = A*(Y2-RR2); G3 = A*(Y3-RR3)
  G1.tmp = aperm(G1,c(2,3,1))
  G1.gradient = colSums(colSums(G1.tmp))
  G2.tmp = aperm(G2,c(3,1,2))
  G2.gradient = colSums(colSums(G2.tmp))
  G3.tmp = G3
  G3.gradient = colSums(colSums(G3.tmp))
  G = G1.gradient+G2.gradient+G3.gradient 
  return(G)    #output gradient at value w.
}
##--------------MLE: Helper Function: Compute the MLE of the Loss function. Input: A: comparison Matrix, Y1,Y2,Y3: responses denoting whether the winner is 1,2,or 3, Output: MLE estimator-------------------------------------------------------------------##
MLE = function(A,Y1,Y2,Y3){
  w0 = rnorm(n); w0 = w0 - mean(w0)
  Gw = gradient(w0,A,Y1,Y2,Y3)
  V = max(abs(Gw))
  w = w0
  k0 = 1
  while(V > 0.0000001 && k0 <= 10000){
    w = w + 0.01*Gw
    Gw = gradient(w,A,Y1,Y2,Y3)
    V = max(abs(Gw))
    k0 = k0 + 1
  }
  return(w)   #output MLE estimator
}
#---------------variance: Helper Function: compute the variance of the estimator: Input: m-th entry, w:estimated MLE, A: comparison matrix; output: values of g^{(m)} function------------------------------------------------------------------##
variance = function(m,w,A){
  n = length(w)
  R1 = array(exp(w),c(n,n,n))
  R2 = aperm(R1,c(2,3,1))
  R3 = aperm(R1,c(3,1,2))
  RR1 = R1/(R1+R2+R3)
  RR2 = R3/(R1+R2+R3)
  RR3 = R2/(R1+R2+R3)
  V1.tmp = A*RR1*(1-RR1); V2.tmp = A*RR2*(1-RR2); V3.tmp = A*RR3*(1-RR3)
  VV = sum(V1.tmp[m,,]+V2.tmp[,m,]+V3.tmp[,,m])
  return(VV)   #output values of g^{(m)} function
}
##------------------F: Helper Function, compute the function $f^{(m)}$ in paper. Input: m:mth entry, w estimator, A: comparison matrix, YY1,YY2,YY3:response variables, L number of comparisons---------------------------------------------------------------##
F = function(m,w,A,YY1,YY2,YY3,L){
  FF = numeric(L)
  n = length(w)
  R1 = array(exp(w),c(n,n,n))
  R2 = aperm(R1,c(2,3,1))
  R3 = aperm(R1,c(3,1,2))
  RR1 = R1/(R1+R2+R3)
  RR2 = R3/(R1+R2+R3)
  RR3 = R2/(R1+R2+R3)
  for(ffff in 1:L){
    GG1 = A*(YY1[,,,ffff]-RR1); GG2 = A*(YY2[,,,ffff]-RR2); GG3 = A*(YY3[,,,ffff]-RR3)
    FF[ffff] = sum(GG3[,,m])+sum(GG2[,m,])+sum(GG1[m,,])
  }
  return(FF)    #output values of f^{(m)}
}
##-----------------------------GA.CI: Helper Function, compute the confidence intervals for item m. Input: m:mth entry, w estimator, A: comparison matrix, Y1,Y2,Y3,YY1,YY2,YY3:response variables, L number of comparisons,B:How many times of bootstrap samples, W: true theta, seq.C: c_0 values on Page 26 of the main text; Output, confidence intervals----#
GA.CI = function(n,m,A,Y1,Y2,Y3,YY1,YY2,YY3,L,B,W,seq.C){
  w.est = MLE(A,Y1,Y2,Y3) 
  var.V = numeric(n)
  F.V = matrix(0,n,L)
  for(kk in 1:n){
    var.V[kk] = variance(kk,w.est,A)
    F.V[kk,] = F(kk,w.est,A,YY1,YY2,YY3,L)
  }
  FG.V = F.V/var.V
  F.M = FG.V[m,] - t(FG.V[-m,])
  sd.M = sqrt(1/var.V[m]+1/var.V[-m])
  F.M.N = t(F.M)/sd.M
  B.M = matrix(rnorm(L*B),L,B)
  B.M.P = F.M.N%*%B.M
  cut.tmp = apply(abs(B.M.P),2,max)
  cut.v = quantile(cut.tmp,0.95)/sqrt(L)
  w.V = w.est[m]-w.est[-m]
  w.V.N = w.V*sqrt(L)/sd.M
  R.left = 1 + sum(1*(-w.V.N>cut.v))
  R.right = n - sum(1*(w.V.N>cut.v))
  R.length = R.right-R.left
  R.CI = sum(1*(m<R.left))+sum(1*(m>R.right))
  ww.V = w.est[m]-w.est[-m]-(W[m]-W[-m])
  stat.V = max(abs(ww.V/sd.M))*sqrt(L)
  vv = 1*(stat.V>cut.v)
  cutZ = qnorm(0.975)
  C.len = length(seq.C)
  cut.C = numeric(C.len)
  R.left.C = numeric(C.len); R.right.C = numeric(C.len)
  R.CI.C = numeric(C.len)
  R.length.C = numeric(C.len)
  R.C.I = numeric(C.len)
  R.left.N = numeric(C.len); R.right.N = numeric(C.len)
  R.CI.N = numeric(C.len)
  R.length.N = numeric(C.len)
  R.N.I = numeric(C.len)
  for(uu in 1:C.len){
    tmp.C = seq.C[uu]
    sd.M.C = cutZ/sqrt(var.V[m])/(1+tmp.C)/sqrt(2*log(n))+1/sqrt(var.V[-m])
    F.M.N.C = t(F.M)/sd.M.C
    B.M.P.C = F.M.N.C%*%B.M
    cut.tmp.C = apply(abs(B.M.P.C),2,max)
    cut.C[uu] = quantile(cut.tmp.C,0.95)/sqrt(L)
    w.V.N.C = w.V*sqrt(L)/sd.M.C
    R.left.C[uu] = 1 + sum(1*(-w.V.N.C>cut.C[uu]))
    R.right.C[uu] = n - sum(1*(w.V.N.C>cut.C[uu]))
    R.length.C[uu] = R.right.C[uu]-R.left.C[uu]
    tmp.N.cut = (1+tmp.C)*sqrt(2*log(n))
    R.left.N[uu] = 1 + sum(1*(-w.V.N.C>tmp.N.cut))
    R.right.N[uu] = n - sum(1*(w.V.N.C>tmp.N.cut))
    R.length.N[uu] = R.right.N[uu]-R.left.N[uu]
    R.CI.C[uu] = sum(1*(m<R.left.C[uu]))+sum(1*(m>R.right.C[uu]))
    R.CI.N[uu] = sum(1*(m<R.left.N[uu]))+sum(1*(m>R.right.N[uu]))
    stat.V.C = max(abs(ww.V/sd.M.C))*sqrt(L)
    R.C.I[uu] = 1*(stat.V.C>cut.C[uu])
    R.N.I[uu] = 1*(stat.V.C>tmp.N.cut)
  }
  #####################################################################################
  ##Output
  #####################################################################################
  result = matrix(0,3,(1+2*C.len))
  result[1,] = c(R.CI,R.CI.C,R.CI.N)  #Record length of confidence intervals R.CI is our confidence interval using variance of normalization factor
  #R.CI.C and R.CI.N are confidence intervals using our method and Chao's method via \eta_m as normalization parameter
  result[2,] = c(R.length,R.length.C,R.length.N) #length of confidence intervals 
  result[3,] = c(vv,R.C.I,R.N.I) #one-sided confidence intervals.
  return(result)
}

#-----------------------------GA: Helper Function, function of gaussian approximation. Input: m:mth entry, w estimator, A: comparsion graph, Y1,Y2,Y3,YY1,YY2,YY3:response variables, L number of comparisons,B:How many times of bootstrap samples, W: true theta, para: theoretical p values; Output, indicator on whether exceeding the GA approxmation threshold or not----
GA = function(m,A,Y1,Y2,Y3,YY1,YY2,YY3,L,B,W,para){
  w.est = MLE(A,Y1,Y2,Y3)
  var.V = numeric(n)
  F.V = matrix(0,n,L)
  for(kk in 1:n){
    var.V[kk] = variance(kk,w.est,A)
    F.V[kk,] = F(kk,w.est,A,YY1,YY2,YY3,L)
  }
  FG.V = F.V/var.V
  F.M = FG.V[m,] - t(FG.V[-m,])
  sd.M = sqrt(1/var.V[m]+1/var.V[-m])
  F.M.N = t(F.M)/sd.M
  B.M = matrix(rnorm(L*B),L,B)
  B.M.P = F.M.N%*%B.M
  cut.tmp = apply(abs(B.M.P),2,max)
  cut.v = quantile(cut.tmp,para)/sqrt(L)
  w.V = w.est[m]-w.est[-m]-(W[m]-W[-m])
  stat.V = max(abs(w.V/sd.M))*sqrt(L)
  vv = 1*(stat.V>cut.v)
  return(vv)  #output whether the test statistic exceeds the GA approxmation threshold or not
}