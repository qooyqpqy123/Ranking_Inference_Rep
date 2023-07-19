#####################################################################################
source("Source_Simulation.R")
n = 80; p = 0.05; L = 60; B = 1000; NN = 500; m = 10                                 ## item m  of interest
seq.C = seq(0,1,length=10)                                                           ## sequence of C (used in normalization parameter \eta_m) (see section 5.2 for more details)
C.len = length(seq.C)
###########################Section: Initialization storing variables##########################################################
W.est.M = matrix(0,n,NN)                          ## MLE estimator 
R.one.uniform.M = matrix(0,n,NN)                  ## normalized one sided CI for all items 
R.one.uniform.V = numeric(NN)                     ## normalized one sided CI length for all items    
RR.one.uniform.M = matrix(0,n,NN)                 ## not normalized 
RR.one.uniform.V = numeric(NN)
R.two.V = numeric(NN)                             ## empirical size of two sided CI (theta)
R.two.length.V = numeric(NN)                      ## two sided CI length 
R.CI.V = numeric(NN)                              ## empirical size of two sided CI (rank)
R.one.V = numeric(NN)                             ## empirical size of one sided CI (theta)   
##################################################################################### ## compare with Chao Gao
R.CI.C.M = matrix(0,C.len,NN)                     ## empirical size of two sided CI (rank) based on GA 
R.length.C.M = matrix(0,C.len,NN)                 ## length CI based on GA 
R.C.I.M = matrix(0,C.len,NN)                      ## empirical size CI (theta)
R.CI.N.M = matrix(0,C.len,NN)                     ## empirical size of two sided CI (rank) Chao Gao  
R.length.N.M = matrix(0,C.len,NN)                 ## length CI Chao Gao 
R.N.I.M = matrix(0,C.len,NN)                      ## empirical size CI (theta) Chao gao 
#####################################################################################
#####################################################################################
K = 10
seq.m = (K-5):(K+5)                            ## item index 
#####################################################################################
len.m = length(seq.m)
K.test.M = matrix(0,len.m,NN)                     ## empirical size of one-sided CI (theta)
K.test.rank.M = matrix(0,len.m,NN)                ## empirical size of one-sided CI (rank)
#####################################################################################
#####################################################################################
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
KK.seq = c(5,10,15)                                                                ##--------## sequence of KK, number of ground truth that need to be included.
KK.len = length(KK.seq)                                                          
KK.set.R = matrix(0,KK.len,NN)                                                     ##--------## set cardinality uniform normalized CI 
KK.power.R = matrix(0,KK.len,NN)                                                   ##--------## uniform normalized CI empirical power 
KK.set.RR = matrix(0,KK.len,NN)                                                    ##--------## set cardinality uniform CI 
KK.power.RR = matrix(0,KK.len,NN)                                                  ##--------## uniform CI empirical power
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------##
#####################################################################################
#####################################################################################
##########################Simulation  Begin###########################################################
for(hh in 1:NN){       #We repeat NN times for every simulation 
  #----------Section: Storation Variables Initialization ----------------------
A = array(0,c(n,n,n))    #Matrix to store comparison matrix
Y1 = array(0,c(n,n,n))   #store comparison results Y1, Y2, Y3
Y2 = array(0,c(n,n,n))
Y3 = array(0,c(n,n,n))
YY1 = array(0,c(n,n,n,L))
YY2 = array(0,c(n,n,n,L))
YY3 = array(0,c(n,n,n,L))
W = seq(4,2,length=n); W = W - mean(W)                                               ## true theta
for(o in 1:(n-2)){  #Construct comparison graph and response variables (A, YY1,YY2,YY3,Y1,Y2,Y3) where $YY1,YY2, YY3$ 
                    # are variables for all L comparisons, and Y1,Y2,Y3 are those averaged values.
	for(oo in (o+1):(n-1)){
		for(ooo in (oo+1):n){
			A[o,oo,ooo] = rbinom(1,1,p)
			p1.tmp = exp(W[o])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p2.tmp = exp(W[oo])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p3.tmp = exp(W[ooo])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p.tmp = c(p1.tmp,p2.tmp,p3.tmp)
			for(oooo in 1:L){		
				Y.tmp = rmultinom(1,1,p.tmp)
				YY1[o,oo,ooo,oooo] = Y.tmp[1]
				YY2[o,oo,ooo,oooo] = Y.tmp[2]
				YY3[o,oo,ooo,oooo] = Y.tmp[3]
				Y1[o,oo,ooo] = Y1[o,oo,ooo] + YY1[o,oo,ooo,oooo]/L
				Y2[o,oo,ooo] = Y2[o,oo,ooo] + YY2[o,oo,ooo,oooo]/L        
				Y3[o,oo,ooo] = Y3[o,oo,ooo] + YY3[o,oo,ooo,oooo]/L
			}
		}
	}
}
#####################################################################################
####################Section: Compute MLE #################################################################
	w.est = MLE(A,Y1,Y2,Y3)                                                        ## MLE estimator
	W.est.M[,hh] = w.est
####################Section: Compute values of f^{(m)} and g^{(m)} functions at estimated MLE#################################################################
	var.V = numeric(n)                                                             
	F.V = matrix(0,n,L)
	for(kk in 1:n){
		var.V[kk] = variance(kk,w.est,A)
		F.V[kk,] = F(kk,w.est,A,YY1,YY2,YY3,L)
	}                                                                              ## compute f and g
	FG.V = F.V/var.V
#####################################################################################

########################Section: Construct normalized uniform confidence intervals using Gaussian multiplier bootstrap################################################################
	B.M = matrix(rnorm(L*B),L,B)                          ## Bootstrap matrix of Gaussian variables to conduct Gaussian multiplier bootstrap
	FF = FG.V[1,] - t(FG.V[-1,])
	SDD = sqrt(1/var.V[1]+1/var.V[-1])          
	FF.N = t(FF)/SDD
	ww.V.R = w.est[1]-w.est[-1]-(W[1]-W[-1])   #MLE score differences
	ww.V.R.N = ww.V.R/SDD                      #test statistics of the first entry
	for(mm in 2:n){                            #test statistics for 2:n entires
		FF.tmp = FG.V[mm,] - t(FG.V[-mm,])
		SDD.tmp = sqrt(1/var.V[mm]+1/var.V[-mm])
		FF.N.tmp = t(FF.tmp)/SDD.tmp
		FF.N = rbind(FF.N,FF.N.tmp)
		ww.V.R.tmp = w.est[mm]-w.est[-mm]-(W[mm]-W[-mm])
		ww.V.R.N.tmp = ww.V.R.tmp/SDD.tmp
		ww.V.R.N = c(ww.V.R.N,ww.V.R.N.tmp)
	}
	FF.B  = FF.N%*%B.M*(-1)                  #bootstrap values
	cut.FF = apply(FF.B,2,max)               #take the maximum                             ## uniform normalized one-sided CI
	cut.v.FF = quantile(cut.FF,0.95)/sqrt(L) #compute cutoff value                                              ## R.one.uniform 
	stat.V.R = max(-ww.V.R.N)*sqrt(L)        #final test statistics
	vv.V.R = 1*(stat.V.R>cut.v.FF)            #whether exceed threshold
	R.one.uniform.V[hh] = vv.V.R             
	R.one.uniform = numeric(n)               #record uniform one-sided (normalized) intervals
	for(mmm in 1:n){
		w.V.uniform = w.est[mmm]-w.est[-mmm]
		SDD.tmp.uniform = sqrt(1/var.V[mmm]+1/var.V[-mmm])
		w.V.N.uniform = w.V.uniform*sqrt(L)/SDD.tmp.uniform
		R.one.uniform[mmm] = 1 + sum(1*(-w.V.N.uniform>cut.v.FF))
	}
	R.one.uniform.M[,hh] = R.one.uniform  #uniform one-sided (normalized) intervals
#############################Section: Construct (un-normalized) uniform confidence intervals using Gaussian multiplier bootstrap, comments are similar with the last section except that here we do not normalize test statistics by its variance########################################################
	FFF = FG.V[1,] - t(FG.V[-1,])
	FFF.N = t(FFF)
	ww.V.RR = w.est[1]-w.est[-1]-(W[1]-W[-1])
	for(mm in 2:n){
		FFF.tmp = FG.V[mm,] - t(FG.V[-mm,])
		FFF.N.tmp = t(FFF.tmp)
		FFF.N = rbind(FFF.N,FFF.N.tmp)
		ww.V.RR.tmp = w.est[mm]-w.est[-mm]-(W[mm]-W[-mm])
		ww.V.RR = c(ww.V.RR,ww.V.RR.tmp)
	}
	FFF.B  = FFF.N%*%B.M*(-1)
	cut.FFF = apply(FFF.B,2,max)                                                      ## uniform one-sided CI
	cut.v.FFF = quantile(cut.FFF,0.95)/sqrt(L)                                        ## RR.one.uniform
	stat.V.RR = max(-ww.V.RR)*sqrt(L)
	vv.V.RR = 1*(stat.V.RR>cut.v.FFF)
	RR.one.uniform.V[hh] = vv.V.RR
	RR.one.uniform = numeric(n)
	for(mmm in 1:n){
		w.V.u = w.est[mmm]-w.est[-mmm]
		w.V.N.u = w.V.u*sqrt(L)
		RR.one.uniform[mmm] = 1 + sum(1*(-w.V.N.u>cut.v.FFF))
	}
	RR.one.uniform.M[,hh] = RR.one.uniform
#####################################################################################
##---------------------------------------------------------------------------------##
##-----------------------Section: Compute the power of uniform one-sided confidence interval----------------------------------------------------------##
for(kkkk in 1:KK.len){
	KK.tmp = KK.seq[kkkk]
	KK.vector = 1:KK.tmp
	KK.set.R[kkkk,hh] = n - sum(1*(R.one.uniform>KK.tmp)) #set cardinality of sure screened set using normalized estimator
	KK.set.RR[kkkk,hh] = n - sum(1*(RR.one.uniform>KK.tmp)) #set cardinality of sure screened set using non-normalized estimator
	KK.seq.tmp.R = c(1*(R.one.uniform[KK.vector]>KK.tmp)) #record uncovered rate
	KK.seq.tmp.RR = c(1*(RR.one.uniform[KK.vector]>KK.tmp))
	KK.power.tmp.R = 1*(sum(KK.seq.tmp.R)>0)         #indicator on whether there are uncovered elements
	KK.power.tmp.RR = 1*(sum(KK.seq.tmp.RR)>0)
	KK.power.R[kkkk,hh] = KK.power.tmp.R             #record the uncovered indicator
	KK.power.RR[kkkk,hh] = KK.power.tmp.RR
}
##---------------------------------------------------------------------------------##
####################Section: Compute two-sided and one-sided Test Statistics#################################################################
	F.M = FG.V[m,] - t(FG.V[-m,])           #difference between f^{m} and other f^{(i)}, i neq m.
	sd.M = sqrt(1/var.V[m]+1/var.V[-m])     #variance
	F.M.N = t(F.M)/sd.M                     #normalize test statistics 
	B.M.P = F.M.N%*%B.M                     #bootstrap value 
	cut.tmp = apply(abs(B.M.P),2,max)       #take maximum of bootstrap statistics
	cut.v = quantile(cut.tmp,0.95)/sqrt(L)  #two-sided cut-off value
	cut.tmp.one = apply(-B.M.P,2,max)
	cut.v.one = quantile(cut.tmp.one,0.95)/sqrt(L)  #one-sided cut-off value
	w.V = w.est[m]-w.est[-m]                        #difference of estimators between the m-th item and the others
	w.V.N = w.V*sqrt(L)/sd.M  
	ww.V = w.est[m]-w.est[-m]-(W[m]-W[-m])   #numeritor of the test statistics  (two-sided)
	stat.V = max(abs(ww.V/sd.M))*sqrt(L)     #final test statistics             (two-sided)
	vv = 1*(stat.V>cut.v)                   #indicate whether the test statistic (two-sided) exceed the threshold
	R.two.V[hh] = vv
	stat.V.one = max(-ww.V/sd.M)*sqrt(L)      #test statistic for one-sided interval
	vv.one = 1*(stat.V.one>cut.v.one)        #indicate whether the test statistic exceed the threshold (one-sided)
	R.one.V[hh] = vv.one                   
#####################################################################################
################# Section: Record two-sided and one-sided Confidence Intervals####################################################################
	R.left = 1 + sum(1*(-w.V.N>cut.v))                                             ## [R.left,R.right] CI for item m
	R.right = n - sum(1*(w.V.N>cut.v))
	R.length = R.right-R.left
	R.two.length.V[hh] = R.length
	R.CI = sum(1*(m<R.left))+sum(1*(m>R.right))
	R.CI.V[hh] = R.CI
	R.one = 1 + sum(1*(-w.V.N>cut.v.one))                                          ## [R.one,n] one-sided CI for item m
#####################################################################################
###################Section: Compute two-sided confidence intervals using \eta_m as normalization parameter, compare our method with Chao's method##################################################################
	cutZ = qnorm(0.975)
	C.len = length(seq.C)   #choice of c
	cut.C = numeric(C.len)
	R.left.C = numeric(C.len); R.right.C = numeric(C.len)     #left and right boundary of our two-sided intervals 
	R.CI.C = numeric(C.len)                                   #coverage level of CI based on rank (our method)
	R.length.C = numeric(C.len)                               #length of confidence intervals
	R.C.I = numeric(C.len)                                    #coverage level of CI based on rank (our method)
	R.left.N = numeric(C.len); R.right.N = numeric(C.len)   #left and right boundary of our two-sided intervals  by Chao's method
	R.CI.N = numeric(C.len)                                   #coverage level of CI based on rank (Chao's method)
	R.length.N = numeric(C.len)                               #length of confidence intervals
	R.N.I = numeric(C.len)                                 #coverage level of CI based on rank (Chao's method)
	for(uu in 1:C.len){                                    #For every c, we compute these two-sided confidence intervals
		tmp.C = seq.C[uu]
		sd.M.C = cutZ/sqrt(var.V[m])/(1+tmp.C)/sqrt(2*log(n))+1/sqrt(var.V[-m])
		F.M.N.C = t(F.M)/sd.M.C
		B.M.P.C = F.M.N.C%*%B.M
		cut.tmp.C = apply(abs(B.M.P.C),2,max)
		cut.C[uu] = quantile(cut.tmp.C,0.95)/sqrt(L)
		w.V.N.C = w.V*sqrt(L)/sd.M.C
		R.left.C[uu] = 1 + sum(1*(-w.V.N.C>cut.C[uu]))      #our method, left boundary, right boundary, length
		R.right.C[uu] = n - sum(1*(w.V.N.C>cut.C[uu]))
		R.length.C[uu] = R.right.C[uu]-R.left.C[uu]
		tmp.N.cut = (1+tmp.C)*sqrt(2*log(n))              
		R.left.N[uu] = 1 + sum(1*(-w.V.N.C>tmp.N.cut))     #Chao's method, left boundary, right boundary, length
		R.right.N[uu] = n - sum(1*(w.V.N.C>tmp.N.cut))
		R.length.N[uu] = R.right.N[uu]-R.left.N[uu]
		R.CI.C[uu] = sum(1*(m<R.left.C[uu]))+sum(1*(m>R.right.C[uu]))   #coverage rate of our confidence interval based on theta
		R.CI.N[uu] = sum(1*(m<R.left.N[uu]))+sum(1*(m>R.right.N[uu]))  #coverage rate of Chao's method based based on rank
		stat.V.C = max(abs(ww.V/sd.M.C))*sqrt(L)
		print(stat.V.C)
		R.C.I[uu] = 1*(stat.V.C>cut.C[uu])   #coverage rate of our confidence interval based on theta
		R.N.I[uu] = 1*(stat.V.C>tmp.N.cut)   #coverage rate of Chao's method based based on theta
	}
	R.CI.C.M[,hh] = R.CI.C                #record all hh repetitions
	R.length.C.M[,hh] = R.length.C
	R.C.I.M[,hh] = R.C.I
	R.CI.N.M[,hh] = R.CI.N
	R.length.N.M[,hh] = R.length.N
	R.N.I.M[,hh] = R.N.I
########################Section Compute the power of one-sided intervals#############################################################
	for(hhhh in 1:len.m){  
		tmp.mm = seq.m[hhhh]
		F.M = FG.V[tmp.mm,] - t(FG.V[-tmp.mm,])
		sd.M = sqrt(1/var.V[tmp.mm]+1/var.V[-tmp.mm])
		F.M.N = t(F.M)/sd.M
		B.M.P = F.M.N%*%B.M                  #bootstrap value
		cut.tmp.one = apply(-B.M.P,2,max)    #take the maximum of bootstraped values
		cut.v.one = quantile(cut.tmp.one,0.95)/sqrt(L) #cutoff value 
		w.V = w.est[tmp.mm]-w.est[-tmp.mm] 
		w.V.N = w.V*sqrt(L)/sd.M            
		ww.V = w.est[tmp.mm]-w.est[-tmp.mm]-(W[tmp.mm]-W[-tmp.mm])
		stat.V.one = max(-ww.V/sd.M)*sqrt(L)  #maximum of test statistics
		vv.one.K = 1*(stat.V.one>cut.v.one)   #indicator on whether exceeding thresholds
		K.test.M[hhhh,hh] = vv.one.K   
		R.one = 1 + sum(1*(-w.V.N>cut.v.one))      
		K.test.rank.M[hhhh,hh] = 1*(R.one>K)  #record testing power 
	}
}

write.csv(W.est.M,"Westm.csv")   ## MLE estimator 
               
write.csv(R.one.uniform.M,"Roneuniform.csv")   ## normalized one sided CI for all items 

write.csv(R.one.uniform.V,"Roneuniformv.csv")   ## coverage level of normalized one sided CI length for all items  
             
write.csv(RR.one.uniform.M,"RRoneuniform.csv")  ## not normalized uniform one sided interval

write.csv(RR.one.uniform.V,"RRoneuniform.csv")  ## coverage level of not normalized uniform one sided interval
#---------------------Two sided confidence interval using \sigma as normalization Column 1-3 of Table 1----------------------------
write.csv(R.two.V,"Rtwov.csv")  ## empirical size of two sided CI (theta) 

write.csv(R.two.length.V,"Rtwolengthv.csv")   ## two sided CI length 
                          
write.csv(R.CI.V,"RCIV.csv")    ## empirical size of two sided CI (rank)

write.csv(R.one.V,"Ronev.csv")  ## empirical size of one sided CI (theta) 
######################Two sided confidence interval using \eta as normalization factor, comparison with benchmark: Column 4-6 and Column 7-9 in Table 1############################################################### 
## compare with bentchmark
write.csv(R.CI.C.M,"Rcicm.csv")            ## empirical size of two sided CI (rank) based on GA

write.csv(R.length.C.M,"Rlengthcm.csv")    ## length CI based on GA

write.csv(R.C.I.M,"Rcim.csv")              ## empirical size CI (theta) 
  
write.csv(R.CI.N.M,"Rcinm.csv")           ## empirical size of two sided CI (rank) benchmark  

write.csv(R.length.N.M,"Rlengthnm.csv")   ## length CI benckmark

write.csv(R.N.I.M,"rnim.csv")           ## empirical size CI (theta) benchmark
#---------------------------------------Empirical Size of one sided intervals for ranks in Table 2 -------------
write.csv(K.test.M,"Ktestm.csv")         ## empirical size of one-sided CI (theta)

write.csv(K.test.rank.M,"Ktestrankm.csv") ## empirical size of one-sided CI (rank)

#---------------------------- set cardinality uniform normalized CI in Table 3 -------------------------
write.csv(KK.set.R,"KKsetr.csv")     ##--------## set cardinality uniform normalized CI 
                                                
write.csv(KK.power.R,"KKpowerR.csv")        ##--------## uniform normalized CI empirical power 
                                                
write.csv(KK.set.RR,"KKsetRR.csv")          ##--------## set cardinality uniform CI 
                                                  
write.csv(KK.power.RR,"KKpowerRR.csv")     ##--------## uniform CI empirical power



















































































































