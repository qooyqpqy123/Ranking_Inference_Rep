source("Source_Real_Data.R") #load some helper functions.
##------------------Read In Data-----------------------##
xx1 = read.csv('jester_1.csv')
xx1 = xx1[,-1]
xx2 = read.csv('jester_2.csv')
xx2 = xx2[,-1]
xxyy = rbind(xx1,xx2)   ## dim 14116*100
##---------------- End Section-----------------------------------------------------------------##


##-------------This section: Basic Setting in paper, as suggested in section 5.3--------------------------------------------------------------------##
set.seed(2022)                       
MM0 = 3; n = 100; p = 0.05; L = 80; B = 1000 #MM0: 3 tuples, n: 100 items, p:sampling probability 0.05, L:number Of comparisions, B: times of bootstrap.
seq.C = 1                         #sequence of c when comparing with the benchmark
##--------------End Section-------------------------------------------------------------------##


##------------This section: initialization of Comparison Graphs A, comparisons Y---------------------------------------------------------------------##
edge.NN = choose(n,MM0)           #all M tuples from n items 
graph.NN = rbinom(1,edge.NN,p)    #we sample observed edges with probability p
AB = matrix(0,graph.NN+1000,n)
AB[,1:MM0] = 1
AB = apply(AB,1,sample)
AB = t(AB)
A00.tmp = unique(AB)
A = A00.tmp[1:graph.NN,]
YY0 = matrix(0,graph.NN,n) 
YYY0 = array(0,c(graph.NN,n,L)) 
for(o in 1:graph.NN){          #Construct comparison graph
  tmp.A.v = (A[o,]==1)
  IdxL = sample(1:14116,L)
  for(oo in 1:L){              #Construct comparisons Y
    Y.value = as.numeric(xxyy[IdxL[oo],tmp.A.v])
    Y.max = which.max(Y.value)
    YY0.tmp = numeric(MM0)
    YY0.tmp[Y.max] = 1
    YYY0[o,tmp.A.v,oo] = YY0.tmp
    YY0[o,tmp.A.v] = YY0[o,tmp.A.v] + YY0.tmp/L
  }
}
##------------------End Section---------------------------------------------------------------##


##------------------This section: Compute MLE----------------------------------------------------##
w.est = MLE.new(A,YY0)                            ## MLE estimator
#write.csv(w.est,"w_est.csv")                     #store the MLE Estimator
w.sort = sort(w.est)
MM = 15                                           ## Pick the top 15 elements
w.tmp = tail(w.sort,MM)
w.idx = numeric(MM)                               ## w.idx is the index of largest MM values of theta
for(pp in 1:MM){
  w.idx[pp] = which(w.est == w.tmp[pp])
}
w.idx                                             #Index of top 15 index Column (Joke ID) in Table 4 from low to high (score). 
w.est[w.idx]                                      #correponding Column Scores (top 15 items)
#write.csv(w.idx,"w.idx.csv")
##-----------------------End Section----------------------------------------------------------##


##---------------------This section aims at computing f and g functions, the definitions are given in section 3.2------------------------------------------------------------##
var.V = numeric(n)                                                             
F.V = matrix(0,n,L)
for(kk in c(1:n)){
  var.V[kk] = variance.new(kk,w.est,A)
  F.V[kk,] = F.new(kk,w.est,A,YYY0,L)
}                                                                             
FG.V = F.V/var.V #output values of f/g
##-----------------------End Section----------------------------------------------------------##


##----------------In the following section, we compute uniform one-sided CIs for all variables-----------------------------------------------------------------##
B.M = matrix(rnorm(L*B),L,B)               #Gaussian variables that help on Gaussian multiplier bootstrap
FF = FG.V[1,] - t(FG.V[-1,])               #construct test statistics for the first entry
SDD = sqrt(1/var.V[1]+1/var.V[-1])         #standard error
FF.N = t(FF)/SDD                           #normalized statistics
for(mm in 2:n){                            #construct test statistics for the 2:n entries
  FF.tmp = FG.V[mm,] - t(FG.V[-mm,])
  SDD.tmp = sqrt(1/var.V[mm]+1/var.V[-mm])
  FF.N.tmp = t(FF.tmp)/SDD.tmp
  FF.N = rbind(FF.N,FF.N.tmp)
}
FF.B  = FF.N%*%B.M                         #bootstrap statistics  
cut.FF = apply(FF.B,2,max)                 #maximum of bootstrap statistics                            
cut.v.FF = quantile(cut.FF,0.95)/sqrt(L)   #cut-off value                                     
R.one.uniform = numeric(n)                 #record uniform one-sided confidence values
for(mmm in 1:n){                           #for every item in [n], build their normalized uniform one-sided CIs.
  w.V.uniform = w.est[mmm]-w.est[-mmm]
  SDD.tmp.uniform = sqrt(1/var.V[mm]+1/var.V[-mm])
  w.V.N.uniform = w.V.uniform*sqrt(L)/SDD.tmp.uniform
  R.one.uniform[mmm] = 1 + sum(1*(-w.V.N.uniform>cut.v.FF))
}
R.one.uniform[w.idx]                       #Output the uniform One-Sided CI (normalized)
##-----------------------End Section----------------------------------------------------------##


##----------------------In the following section, we compute uniform one-sided CIs for all variables without normalize by the variance-----------------------------------------------------------##
FFF = FG.V[1,] - t(FG.V[-1,])             #construct test statistics for the first entry
FFF.N = t(FFF)
for(mm in 2:n){                           #construct test statistics for the 2:n entries
  FFF.tmp = FG.V[mm,] - t(FG.V[-mm,])
  FFF.N.tmp = t(FFF.tmp)
  FFF.N = rbind(FFF.N,FFF.N.tmp)
}
FFF.B  = FFF.N%*%B.M                      #bootstraped statistics
cut.FFF = apply(FFF.B,2,max)              #maximum of bootstraped statistics                                              
cut.v.FFF = quantile(cut.FFF,0.95)/sqrt(L)  #cut-off value                                          
RR.one.uniform = numeric(n)
for(mmm in 1:n){                         #for every item in [n], build their normalized uniform one-sided CIs.
  w.V.u = w.est[mmm]-w.est[-mmm]
  w.V.N.u = w.V.u*sqrt(L)
  RR.one.uniform[mmm] = 1 + sum(1*(-w.V.N.u>cut.v.FFF))   
}
RR.one.uniform[w.idx]                    #Output the uniform One-Sided CI Line UOC of Table 4
#write.csv(RR.one.uniform,"RR.one.uniform.csv")
##-----------------------End Section----------------------------------------------------------##


##-------------------In the following section, we record the confidence intervals-----------------------------------##
R.left_vec<-c()            ## [R.left,R.right] is our two-sided CI (column TC1 in Table 4) 
R.right_vec<-c()
R.one_vec<-c()             ## [R.one,n] is our one-sided CI (column OC in Table 4)
R.left.C_vec<-c()          ## [R.left.C,R.right.C] is our two-sided CI (column TC2 in Table 4)
R.right.C_vec<-c()
R.left.N_vec<-c()          ## [R.left.N,R.right.N] is Chao Gao CI   (Column TC3 in Table 4)
R.right.N_vec<-c()
for (i in c(1:length(w.idx))){                       #For any index in top-15, compute its confidence interval
  m=w.idx[i]
  F.M = FG.V[m,] - t(FG.V[-m,])                      #For every entry m construct test statistics
  sd.M = sqrt(1/var.V[m]+1/var.V[-m])
  F.M.N = t(F.M)/sd.M
  B.M.P = F.M.N%*%B.M
  cut.tmp = apply(abs(B.M.P),2,max)                  #maximum of test statistics
  cut.v = quantile(cut.tmp,0.95)/sqrt(L)             #cut-off value for two sided confidence intervals
  cut.tmp.one = apply(B.M.P,2,max)                   #maximum of test statistics
  cut.v.one = quantile(cut.tmp.one,0.95)/sqrt(L)     #cut-off value for one sided confidence intervals
  w.V = w.est[m]-w.est[-m]
  w.V.N = w.V*sqrt(L)/sd.M
  R.left = 1 + sum(1*(-w.V.N>cut.v))                    ## [R.left,R.right] CI for item m
  R.right = n - sum(1*(w.V.N>cut.v))
  R.length = R.right-R.left
  R.CI = sum(1*(m<R.left))+sum(1*(m>R.right))
  R.one = 1 + sum(1*(-w.V.N>cut.v.one))                 ## [R.one,n] one-sided CI for item m
  cutZ = qnorm(0.975)
  C.len = length(seq.C)
  cut.C = numeric(C.len)
  R.left.C = numeric(C.len); R.right.C = numeric(C.len)   ## record [R.left.C,R.right.C] is our two-sided CI (column TC2 in Table 4)  
  R.CI.C = numeric(C.len)
  R.length.C = numeric(C.len)
  R.C.I = numeric(C.len)
  R.left.N = numeric(C.len); R.right.N = numeric(C.len)  ## record [R.left.N,R.right.N] is Chao Gao CI   (Column TC3 in Table 4)
  R.CI.N = numeric(C.len)
  R.length.N = numeric(C.len)
  R.N.I = numeric(C.len)
  for(uu in 1:C.len){             #record these confidence intervals for every different constant C in normalization factor (see section 5.3 of papers for more details)
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
    R.CI.C[uu] = sum(1*(m<R.left.C[uu]))+sum(1*(m>R.right.C[uu]))     #coverage rate of Our confidence intervals (rank)
    R.CI.N[uu] = sum(1*(m<R.left.N[uu]))+sum(1*(m>R.right.N[uu]))     #coverage rate of Confidence intervals of Chao's method (rank)
  }
  R.left_vec<-c(R.left_vec, R.left)                                    ## [R.left,R.right] is our two-sided CI      
  #write.csv(R.left_vec,"R.left_vec.csv")
  R.right_vec<-c(R.right_vec,R.right)   
  #write.csv(R.right_vec,"R.right_vec.csv")
  R.one_vec<-c(R.one_vec,R.one)                                       ## [R.left,n] is our one-sided CI
  #write.csv(R.one_vec,"R.one_vec.csv.csv")
  R.left.C_vec<-c(R.left.C_vec,R.left.C)  
  #write.csv(R.left.C_vec,"R.left.C_vec.csv")
  R.right.C_vec<-c(R.right.C_vec,R.right.C) 
  #write.csv(R.right.C_vec,"R.right.C_vec.csv")
  R.left.N_vec<-c(R.left.N_vec,R.left.N) 
  #write.csv(R.left.N_vec,"R.left.N_vec.csv")                         ## [R.left.N,R.right.N] is Chao Gao CI       
  R.right.N_vec<-c(R.right.N_vec,R.right.N)} 
  #write.csv(R.right.N_vec,"R.right.N_vec.csv")}
  #####################################################################################
  
  print(R.left_vec)    ## [R.left,R.right] is our two-sided CI (column TC1 in Table 4) 
  print(R.right_vec)  
  print(R.one_vec)     ## [R.one,n] is our one-sided CI (column OC in Table 4)
  print(R.left.C_vec)  ## [R.left.C,R.right.C] is our two-sided CI (column TC2 in Table 4)
  print(R.right.C_vec)
  print(R.left.N_vec)  ## [R.left.N,R.right.N] is Chao Gao CI   (Column TC3 in Table 4)
  print(R.right.N_vec)
   #--------------------------End Section------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
