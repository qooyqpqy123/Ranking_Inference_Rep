source("Source_Real_Data.R") #load some helper functions.
##------------------Read In Data-----------------------##
xx1 = read.csv('jester_1.csv')
xx1 = xx1[,-1]
xx2 = read.csv('jester_2.csv')
xx2 = xx2[,-1]
xxyy = rbind(xx1,xx2)   ## dim 14116*100
##---------------End Section------------------------------------------------------------------##



###------------   Section: Instruction for the input ----------------------------------------###############
set.seed(2022) 
MM0 = 6; p = 0.0000063;      #Input MM0: Tuple Size, we can set MM0=2,...,6; p: sampling probability, the value can be found in section 5.3, 
                             #By setting M,p as suggested in Table 5, or Table 6, we can reproduce the results.
n = 100; L = 80; B = 1000    #n: number of items to be compared, L sampled comparisons, B: number of bootstraps.
seq.C = 1                    #constant c in normalization parameter
#-------End Section --------------------------------

##------------Section :Initializarion of Comparison Graphs---------------------------------------------------------------------##
edge.NN = choose(n,MM0)          #all M tuples from n items 
graph.NN = rbinom(1,edge.NN,p)   #we sample observed edges with probability p
AB = matrix(0,graph.NN+1000,n)
AB[,1:MM0] = 1
AB = apply(AB,1,sample)
AB = t(AB)
A00.tmp = unique(AB)
A = A00.tmp[1:graph.NN,]
YY0 = matrix(0,graph.NN,n)
YYY0 = array(0,c(graph.NN,n,L))
for(o in 1:graph.NN){            #Construct comparison graph
  tmp.A.v = (A[o,]==1)
  IdxL = sample(1:14116,L)
  for(oo in 1:L){                # Construct comparisons Y via L independent comparisons
    Y.value = as.numeric(xxyy[IdxL[oo],tmp.A.v])
    Y.max = which.max(Y.value)
    YY0.tmp = numeric(MM0)
    YY0.tmp[Y.max] = 1
    YYY0[o,tmp.A.v,oo] = YY0.tmp
    YY0[o,tmp.A.v] = YY0[o,tmp.A.v] + YY0.tmp/L
  }
}
#--------------------End Section------------------------##


##------------------Section: Compute MLE----------------------------------------------------##
w.est = MLE.new(A,YY0)                              ## MLE estimator
w.sort = sort(w.est)
MM = 15                                             ## MM
w.tmp = tail(w.sort,MM)
w.idx = numeric(MM)                                  ## w.idx is the index of largest MM values of theta
for(pp in 1:MM){
  w.idx[pp] = which(w.est == w.tmp[pp])
}
w.idx_new=c(10,30,50,70,90)                         #The index we are interested in
#write.csv(w.idx,"w.idx.csv")
##-----------------------End Section----------------------------------------------------------##



##---------------------------This section aims at computing f and g functions, the definitions are given in section 3.2------------------------------------------------------##
var.V = numeric(n)                                                             
F.V = matrix(0,n,L)
for(kk in c(1:n)){
  var.V[kk] = variance.new(kk,w.est,A)
  F.V[kk,] = F.new(kk,w.est,A,YYY0,L)
}                                         ## compute f and gcompute f and g, the definitions are given in section 3.2
FG.V = F.V/var.V
##--------------------------End Section-------------------------------------------------------##



##-------------------Section: Record the confidence intervals-----------------------------------##
B.M = matrix(rnorm(L*B),L,B)                
R.left_vec<-c()                             ## [R.left,R.right] is our two-sided CI for every p and M    
R.right_vec<-c()
R.one_vec<-c()                               ## [R.one,n] is our one-sided CI (column OC in Table 4)
for (i in c(1:length(w.idx_new))){
  m=w.idx_new[i]                             #construct confidence intervals for items 10,30,50,70,90
  F.M = FG.V[m,] - t(FG.V[-m,])              #For every entry m construct test statistics
  sd.M = sqrt(1/var.V[m]+1/var.V[-m])
  F.M.N = t(F.M)/sd.M
  B.M.P = F.M.N%*%B.M
  cut.tmp = apply(abs(B.M.P),2,max)                 #maximum of test statistics
  cut.v = quantile(cut.tmp,0.95)/sqrt(L)            #cut-off value for two sided confidence intervals
  cut.tmp.one = apply(B.M.P,2,max)                 #maximum of test statistics
  cut.v.one = quantile(cut.tmp.one,0.95)/sqrt(L)    #cut-off value for one sided confidence intervals
  w.V = w.est[m]-w.est[-m]
  w.V.N = w.V*sqrt(L)/sd.M
  R.left = 1 + sum(1*(-w.V.N>cut.v))                        ## [R.left,R.right] CI for item m
  R.right = n - sum(1*(w.V.N>cut.v))
  R.length = R.right-R.left
  R.CI = sum(1*(m<R.left))+sum(1*(m>R.right))
  R.one = 1 + sum(1*(-w.V.N>cut.v.one))               ## [R.one,n] one-sided CI for item m
  cutZ = qnorm(0.975)
  #####################################################################################
  R.left_vec<-c(R.left_vec, R.left)                                   ## [R.left,R.right] is our two-sided CI      
  R.right_vec<-c(R.right_vec,R.right) 
  R.one_vec<-c(R.one_vec,R.one)                                       ## [R.left,n] is our one-sided CI
  }
  #####################################################################################
  print(R.left_vec)                                                   ## [R.left,R.right] is our two-sided CI (column TC1 in Table 5 or 6)  for different choice of $MM0$ and p.
  print(R.right_vec) 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  