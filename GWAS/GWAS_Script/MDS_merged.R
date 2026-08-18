# data<- read.table(file="MDS_merge2.mds",header=TRUE)
data<- read.table(file="MDS_merge2.eigenvec",header=FALSE)
colnames(data) <- c("FID", "IID", "PC1","PC2",	"PC3",	"PC4",	"PC5",	"PC6",	"PC7",	"PC8",	"PC9",	"PC10")
print(colnames(data))
race<- read.table(file="racefile.txt",header=TRUE)
print(colnames(race))
datafile<- merge(data,race,by=c("IID","FID"))
head(datafile)

pdf("MDS_2.pdf",width=7,height=7)
for (i in 1:nrow(datafile))
{
if (datafile[i,13]=="EUR") {plot(datafile[i,3],datafile[i,4],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="green")}
par(new=T)
if (datafile[i,13]=="ASN") {plot(datafile[i,3],datafile[i,4],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="red")}
par(new=T)
if (datafile[i,13]=="AMR") {plot(datafile[i,3],datafile[i,4],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col=470)}
par(new=T)
if (datafile[i,13]=="AFR") {plot(datafile[i,3],datafile[i,4],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="blue")}
par(new=T)
if (datafile[i,13]=="OWN") {plot(datafile[i,3],datafile[i,4],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=3,cex=0.7,col="black")}
par(new=T)
}

abline(v=-0.003,lty=3)
abline(h=0.01,lty=3)
legend("topright", pch=c(1,1,1,1,3),c("EUR","ASN","AMR","AFR","OWN"),col=c("green","red",470,"blue","black"),bty="o",cex=1)
