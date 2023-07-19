source("Source_Simulation.R")
n = 60; L = 20; NN = 200 #Setting of section 5.1 Here we let L is fixed and let p vary. 
W = runif(n,2,4)         ## true theta 
W = W - mean(W)          #centered theta
r.v = seq(1.4*log(60)/sqrt(n*0.2*L),1.4*log(60)/sqrt(n*0.01*L),length=8)        ## ratio sequence  
EEE = matrix(0,NN,8)    #initialize matrix to record the results
EEE2=matrix(0,NN,8) 
EEE3=matrix(0,NN,8)
for(ee in 1:8){
	p = 1/n/n/L/r.v[ee]/r.v[ee]## p from 0.01 to 0.2 
	#---------#Section Matrix Initialization for generating comparison graph A, responses Y1,Y2,Y3.#----------
for(eee in 1:NN){   
A = array(0,c(n,n,n))
Y1 = array(0,c(n,n,n))
Y2 = array(0,c(n,n,n))
Y3 = array(0,c(n,n,n))
for(o in 1:(n-2)){
	for(oo in (o+1):(n-1)){
		for(ooo in (oo+1):n){
			A[o,oo,ooo] = rbinom(1,1,p)
			p1.tmp = exp(W[o])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p2.tmp = exp(W[oo])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p3.tmp = exp(W[ooo])/(exp(W[o])+exp(W[oo])+exp(W[ooo]))
			p.tmp = c(p1.tmp,p2.tmp,p3.tmp)
			Y.tmp = rmultinom(L,1,p.tmp)
			Y.tmp.mean = rowMeans(Y.tmp)
			Y1[o,oo,ooo] = Y.tmp.mean[1]
			Y2[o,oo,ooo] = Y.tmp.mean[2]
			Y3[o,oo,ooo] = Y.tmp.mean[3]
		}
	}
}
#-----------------Section: Compute MLE statistical rates---------------------
w.est = MLE(A,Y1,Y2,Y3)      ## MLE estimate
EEE2[eee,ee]=w.est-W       ## difference between the estimated MLE and the true parameter
EEE[eee,ee] = max(abs(w.est-W)) ## record the L infty norm difference
EEE3[eee,ee]=sum((w.est-W)^2)   #record the L_2 norm difference
}}


