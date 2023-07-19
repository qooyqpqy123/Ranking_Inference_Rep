source("Source_Simulation.R")
n = 60; L = 10; NN = 500   #Repeat 500 times We can take L=5/10/20. Here we use L=10 as an example
W = runif(n,2,4)          ## true theta 
W = W - mean(W)           #Normarlize theta
p_list= c(0.008,0.015,0.03)      #p: sampling probablity
qq=matrix(0,NN,length(p_list))   #Store the output
for (index in c(1:length(p_list))){  #Compute the histogram.
  #---------#Section Matrix Initialization for generating comparison graph A, responses Y1,Y2,Y3.#----------
for(eee in 1:NN){
A = array(0,c(n,n,n))
Y1 = array(0,c(n,n,n))
Y2 = array(0,c(n,n,n))
Y3 = array(0,c(n,n,n))
for(o in 1:(n-2)){
	for(oo in (o+1):(n-1)){
		for(ooo in (oo+1):n){
			A[o,oo,ooo] = rbinom(1,1,p_list[index])
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
#-----------------Section: Compute normalized MLE in order to plot Histograms---------------------
w.est = MLE(A,Y1,Y2,Y3)# # MLE estimate
var=variance(1,w.est,A)  #estimate variance
qq[eee,index]=(w.est[1]-W[1])*sqrt(L*var)  #normalize
write.csv(qq,"qq.csv")   #store the value
}}
	



