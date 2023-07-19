#This file contains functions:
###################################
#gradient.new: Helper Function: Take the gradient of the likelihood function: Input w: theta, A: comparison matrix, Y: response y; Output: gradient of loss function at parameter $w$#
#MLE.new: Helper Function: Compute the MLE of the Loss function. Input: A: comparison matrix, Y: Response varialble, Output: MLE estimator
#variance.new: Helper Function: compute the variance of the estimator: Input: m-th entry (m is the entry of interest), w:estimated MLE, A: comparison matrix; output: variance of the estimator
#F.new: Helper Function, compute the function $f^{(m)}$ in paper. Input: m:mth entry, w estimator, comparison matrix, Y.A:response variable: y in the paper, L number of comparisons;#Output value of function $f^{(m)}$ at $w$.
####################################
#-------------gradient.new: Helper Function: Take the gradient of the likelihood function: Input w: theta, A: comparison matrix, Y: response y; Output: Gradient of the Function#--------------
gradient.new = function(w,A,Y){ 
  n = length(w)
  rowA = nrow(A)
  R0 = matrix(exp(w),n,rowA)
  R0 = t(R0)
  RR0 = rowSums(A*R0)
  RR1 = R0/RR0
  RR1 = RR1*A
  G1 = (YY0-RR1)*A
  G = colSums(G1)
  return(G)  #Outout gradient of loss function at parameter $w$.
}
####--------------MLE.new: Helper Function: Compute the MLE of the Loss function. Input: A: comparison matrix, Y: Response varialble, Output: MLE estimator-------------------------------------------------------------------##
MLE.new = function(A,Y){
  w0 = rnorm(n); w0 = w0 - mean(w0)
  Gw = gradient.new(w0,A,Y)
  V = max(abs(Gw))
  w = w0
  k0 = 1
  while(V > 0.00001 && k0 <= 10000){
    w = w + 0.01*Gw
    print(mean(w))
    ##w = w - mean(w)
    Gw = gradient.new(w,A,Y)
    V = max(abs(Gw))
    print(V)
    k0 = k0 + 1
  }
  return(w)     #Output MLE Estimator
}
##---------------variance.new: Helper Function: compute the variance of the estimator: Input: m-th entry (m is the entry of interest), w:estimated MLE, A: comparison matrix; output: variance of the estimator------------------------------------------------------------------##
variance.new = function(m,w,A){
  n = length(w)
  rowA = nrow(A)
  R0 = matrix(exp(w),n,rowA)
  R0 = t(R0)
  RR0 = rowSums(A*R0)
  RR1 = R0/RR0
  RR1 = RR1*A
  Idx.m = (A[,m]==1)
  tmp.sub.A.m = A[Idx.m,m]
  tmp.sub.R.m = RR1[Idx.m,m]
  V1.tmp = tmp.sub.R.m*(1-tmp.sub.R.m)*tmp.sub.A.m
  VV = sum(V1.tmp)
  return(VV) #Output: value of function $g^{(m)}$ at $w$.
}
##------------------F.new: Helper Function, compute the function $f^{(m)}$ in paper. Input: m:mth entry, w estimator, comparison matrix, Y.A:response variable: y in the paper, L number of comparisons---------------------------------------------------------------##
F.new = function(m,w,A,Y.A,L){
  n = length(w)
  rowA = nrow(A)
  R0 = matrix(exp(w),n,rowA)
  R0 = t(R0)
  RR0 = rowSums(A*R0)
  RR1 = R0/RR0
  RR1 = RR1*A
  Idx.m = (A[,m]==1)
  tmp.sub.R.m = RR1[Idx.m,m]
  Y.A.tmp = Y.A[Idx.m,m,]-tmp.sub.R.m
  FF = colSums(Y.A.tmp)
  return(FF)  #Output value of function $f^{(m)}$ at $w$.
}