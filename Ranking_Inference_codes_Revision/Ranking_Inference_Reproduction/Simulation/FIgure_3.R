source("Source_Simulation.R")
paraseq=seq(0.05,0.9,0.05) #p-values ranges
n = 60; p = 0.05; L = 80; B = 500; NN = 500  #Setting of Figure 4 in section 5.1  
m = 1 #we are interested in the first entry.
vvvv = matrix(NN,length(paraseq)) #record the empirical p-value
#---------#Section Matrix Initialization for generating comparison graph A, Single comparisons YY1, YY2, YY3, responses Y1,Y2,Y3 (averages of YY1,YY2,YY3 through L independent comparisons)#----------
for(hh in 1:NN){
A = array(0,c(n,n,n))
Y1 = array(0,c(n,n,n))
Y2 = array(0,c(n,n,n))
Y3 = array(0,c(n,n,n))
YY1 = array(0,c(n,n,n,L))
YY2 = array(0,c(n,n,n,L))
YY3 = array(0,c(n,n,n,L))
W = runif(n,2,4); W = W - mean(W)  ## true theta
for(o in 1:(n-2)){
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
#for all paraseq, record the empirical size of the test and compare with the theoretical ones.
for (i in length(paraseq)){
vvvv[hh,i] = GA(m,A,Y1,Y2,Y3,YY1,YY2,YY3,L,B,W,paraseq[i])  #empirical p-values computed via Gaussian Approximation
}
#write.csv(vvvv,'gauss_app.csv') #store the empirical p-values
}


