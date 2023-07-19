data1<-read.csv("jester-data-1.csv",header=FALSE)
Data1<-data1[which(data1$V1==100),]
Data1<-Data1[,-1]
rownames(Data1)<-NULL
write.csv(Data1,"jester_1.csv")

data2<-read.csv("jester-data-2.csv",header=FALSE)
Data2<-data2[which(data2$V1==100),]
Data2<-Data2[,-1]
rownames(Data2)<-NULL
write.csv(Data2,"jester_2.csv")

