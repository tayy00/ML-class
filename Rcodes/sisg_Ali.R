setwd('/Users/ashojaie/Dropbox/ClassesEtc/SISG')
library(MASS)
library(ISLR)
library(mclust)
library(cluster)
library(genefilter)
library(igraph)
library(glasso)
library(bnlearn)

##
## example of multiple hypothesis testing
##
p <- 1000
niter <- 100
n <- 50
nrej <- numeric(niter)
for(i in 1:niter){
  x <- matrix(rnorm(p*n),p,n)
  pvals <- rowttests(x, fac=factor(rep(1:2,each=n/2)))$p.value
  nrej[i] <- sum(pvals < 0.05)
}

pdf('multipletest.pdf', height=4, width=7)
par(mfrow=c(1,2))

hist(nrej, freq=F, main='no of rejected hypotheses', xlab='')
abline(v=mean(nrej), col='red')

plot((1:p)/p,sort(ps), main='QQ plot of observed p-valuse', 
	xlab='U(0,1) quantiles', ylab='Observed quantiles')
abline(0,1,col='red')
dev.off()

#library(samr)
#SAM(x,y=rep(1:2,each=n/2))


##
##testing corr
##
x <- seq(-2,2,length.out=1000)
y <- x^2 + runif(1000)

#plot(x,y)
cor(x,y)
cor(x,y,method="kendall")
cor(x,y,method="spearman")

##
##
##
n1 <- 300
n2 <- 200

s1 <- 0.4
s2 <- 0.5

m1 <- 0.7
m2 <- -1.9

rho1 = 0.2
rho2 = 0.9
Ip <- diag(c(1,1))
Sig1 <- matrix(c(1, rho1, rho1, 1), 2, 2)
Sig2 <- matrix(c(1, -rho2, -rho2, 1), 2, 2)

d1 <- mvrnorm(n = n1, mu=c(m1,1), Sigma=s1*Sig1)
d2 <- mvrnorm(n = n2, mu=c(m2,-0.5), Sigma=s2*Sig2)

d <- dist(rbind(d1,d2))
d.s <- dist(rbind(sqr1,sqr2))

hc.s <- hclust(d=d, method='single')
g.s <- cutree(hc.s, k=2)

hc.c <- hclust(d=d, method='complete')
g.c <- cutree(hc.c, k=2)

hc.a <- hclust(d=d, method='average')
g.a <- cutree(hc.a, k=2)

pdf('linkcomp_s.pdf')
plot(-4:4, -4:4, type='n', xlab='', ylab='', main='single linkage')
points(rbind(d1,d2), cex=0.7, col=g.s+2)
dev.off()

pdf('linkcomp_c.pdf')
plot(-4:4, -4:4, type='n', xlab='', ylab='', main='complete linkage')
points(rbind(d1,d2), cex=0.7, col=g.c+2)
dev.off()

pdf('linkcomp_a.pdf')
plot(-4:4, -4:4, type='n', xlab='', ylab='', main='average linkage')
points(rbind(d1,d2), cex=0.7, col=g.a+2)
dev.off()


side.1 <- seq(-2, 2, by=0.01)
side.2 <- seq(-1, 1, by=0.01)

nside.1 <- length(side.1)
nside.2 <- length(side.2)

sqr1 <- matrix(0, nrow=4*nside.1, ncol=2)
sqr2 <- matrix(0, nrow=4*nside.2, ncol=2)

uniwidth <- 0.4

sqr1[1:nside.1,1] <- side.1 + runif(nside.1, -uniwidth, uniwidth)
sqr1[1:nside.1,2] <- -2 + runif(nside.1, -uniwidth, uniwidth)
sqr1[(nside.1+1):(2*nside.1),1] <- -2 + runif(nside.1, -uniwidth, uniwidth)
sqr1[(nside.1+1):(2*nside.1),2] <- side.1 + runif(nside.1, -uniwidth, uniwidth)
sqr1[(2*nside.1+1):(3*nside.1),1] <- side.1 + runif(nside.1, -uniwidth, uniwidth)
sqr1[(2*nside.1+1):(3*nside.1),2] <- 2 + runif(nside.1, -uniwidth, uniwidth)
sqr1[(3*nside.1+1):(4*nside.1),1] <- 2 + runif(nside.1, -uniwidth, uniwidth)
sqr1[(3*nside.1+1):(4*nside.1),2] <- side.1 + runif(nside.1, -uniwidth, uniwidth)

sqr2[1:nside.2,1] <- side.2 + runif(nside.2, -uniwidth, uniwidth)
sqr2[1:nside.2,2] <- -1 + runif(nside.2, -uniwidth, uniwidth)
sqr2[(nside.2+1):(2*nside.2),1] <- -1 + runif(nside.2, -uniwidth, uniwidth)
sqr2[(nside.2+1):(2*nside.2),2] <- side.2 + runif(nside.2, -uniwidth, uniwidth)
sqr2[(2*nside.2+1):(3*nside.2),1] <- side.2 + runif(nside.2, -uniwidth, uniwidth)
sqr2[(2*nside.2+1):(3*nside.2),2] <- 1 + runif(nside.2, -uniwidth, uniwidth)
sqr2[(3*nside.2+1):(4*nside.2),1] <- 1 + runif(nside.2, -uniwidth, uniwidth)
sqr2[(3*nside.2+1):(4*nside.2),2] <- side.2 + runif(nside.2, -uniwidth, uniwidth)

plot(-4:4, -4:4, type='n', xlab='', ylab='', main='single linkage')
points(sqr1, cex=0.7, col=3)
points(sqr2, cex=0.7, col=4)

ds <- dist(rbind(sqr1,sqr2))

hc.s <- hclust(d=ds, method='single')
g.s <- cutree(hc.s, k=2)

hc.c <- hclust(d=ds, method='complete')
g.c <- cutree(hc.c, k=2)

hc.a <- hclust(d=ds, method='average')
g.a <- cutree(hc.a, k=2)

pdf('linkcomp_s.pdf')
plot(-4:4, -4:4, type='n', xlab='', ylab='', main='single linkage')
points(rbind(sqr1,sqr2), cex=0.7, col=g.s+2)
dev.off()

pdf('linkcomp_c.pdf')
plot(-4:4, -4:4, type='n', xlab='', ylab='', main='complete linkage')
points(rbind(sqr1,sqr2), cex=0.7, col=g.c+2)
dev.off()

pdf('linkcomp_a.pdf')
plot(-4:4, -4:4, type='n', xlab='', ylab='', main='average linkage')
points(rbind(sqr1,sqr2), cex=0.7, col=g.a+2)
dev.off()

##
## Analysis NCI60 data
##
dat <- read.table('nci.data')
dat <- t(dat)
dim(dat)

labs <- read.table('nci.info', skip=14)
labs <- as.vector(as.matrix(labs))

table(labs)

dat <- scale(dat, center=F, scale=T)

##define the distance
d <- dist(dat)

par(mfrow=c(1,3))
plot(hclust(d), labels=labs, col="green", main="Complete Linkage", hang=-1,
	xlab="", sub="", ylab="", cex.main=1.5, cex.lab=0.6)
plot(hclust(d, method="average"), labels=labs, col="orange", hang=-1,
	main="Average Linkage", xlab="", sub="", ylab="", cex.main=1.5, cex.lab=0.8)
plot(hclust(d, method="single"), labels=labs, col="blue", hang=-1,
	main="Single Linkage", xlab="", sub="", ylab="", cex.main=1.5, cex.lab=0.8)

par(mfrow=c(1,3))
plot(hclust(d), labels=NULL, col="green", main="Complete Linkage", hang=-.1,
	xlab="", sub="", ylab="", cex.main=1.5, cex.lab=0.6)
plot(hclust(d, method="average"), labels=NULL, col="orange", hang=-.1,
	main="Average Linkage", xlab="", sub="", ylab="", cex.main=1.5, cex.lab=0.8)
plot(hclust(d, method="single"), labels=NULL, col="blue", hang=-.1,
	main="Single Linkage", xlab="", sub="", ylab="", cex.main=1.5, cex.lab=0.8)



##
##K-means clustering
##
# a 2-dimensional example
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(x) <- c("x", "y")
cl <- kmeans(x, 2)

pdf('kmeans-ex.pdf')
par(mar=c(4.5,4.5,1,1))
plot(x, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)
dev.off()

K <- 10
wss <- numeric(K)
bss <- numeric(K)
tss <- numeric(K)
for(k in 1:K){
	cl <- kmeans(x, k, nstart = 25)
	tss[k] <- cl$totss
	wss[k] <- cl$tot.withinss
	bss[k] <- cl$betweenss
}

pdf('kmeans-SS-1.pdf', width=6, height=4)
par(mar=c(4.5,4.5,1,1))
plot(1:K, wss, type='b', lty=2, lwd=2, col=2, ylab='Sum of Squares', 
	xlab='Number of clusters', ylim=c(min(c(bss, wss)), max(c(bss, wss))))
#lines(1:K, bss, lty=4, lwd=2, col=4)
legend('right', c('WSS'), col=c(2), lwd=2, lty=c(2), bty='n', pch=1)
arrows(x0=3,y0=30,x1=2.1,y1=19, col=3, lwd=2)
text(x=4,y=33,'Optimal Number of Clusters', col=3)
dev.off()

pdf('kmeans-SS-2.pdf', width=6, height=4)
par(mar=c(4.5,4.5,1,1))
plot(1:K, wss, type='b', lty=2, lwd=2, col=2, ylab='Sum of Squares', 
	xlab='Number of clusters', ylim=c(min(c(bss, wss)), max(c(bss, wss))))
lines(1:K, bss, type='b', lty=4, lwd=2, col=4, pch=2)
legend('right', c('WSS','BSS'), col=c(2,4), lwd=2, lty=c(2,4), pch=c(1,2), bty='n')
dev.off()

##
## K-means clustering for the NCI data
##
dat <- read.table('nci.data')
dim(dat)
dat <- t(dat)
dim(dat)

labs <- read.table('nci.info', skip=14)
labs <- as.vector(as.matrix(labs))
#table(labs)

K <- 10
wss <- numeric(K)
bss <- numeric(K)
tss <- numeric(K)
for(k in 1:K){
	cl <- kmeans(dat, k, nstart = 25)
	tss[k] <- cl$totss
	wss[k] <- cl$tot.withinss
	bss[k] <- cl$betweenss
}

pdf('NCIkmeans-1.pdf')
par(mar=c(4.5,4.5,1,1), width=6, height=4)
plot(1:K, wss, type='b', lty=1, lwd=2, col=2, ylab='Sum of Squares', 
	xlab='Number of clusters')
dev.off()

##remove subtypes with only one member
indx <- names(table(labs))[table(labs) == 1]
dat2 <- dat[!(labs %in% indx),]
#dim(dat2)

##select subset of data and order the samples
labs2 <- labs[!(labs %in% indx)]
ord <- order(labs2)
dat2 <- dat2[ord,]
labs2 <- labs2[ord]

K <- 20
wss <- numeric(K)
bss <- numeric(K)
tss <- numeric(K)
for(k in 1:K){
	cl <- kmeans(dat2, k, nstart = 25)
	tss[k] <- cl$totss
	wss[k] <- cl$tot.withinss
	bss[k] <- cl$betweenss
}

pdf('NCIkmeans-2.pdf', width=6, height=4)
par(mar=c(4.5,4.5,2,1))
plot(1:K, wss, type='b', lty=1, lwd=2, col=2, ylab='Sum of Squares', 
	xlab='Number of clusters')
dev.off()

##
##Model Based Clustering
##

library(mclust)
?Mclust

irisMclust <- Mclust(iris[,-5])
plot(irisMclust)


##
##Spectral Clustering
##
p <- 50
r <- 1
s1 <- matrix(r,p,p)
diag(s1) <- 0
s <- matrix(0,2*p,2*p)
s[1:p,1:p] <- -s1; s[(p+1):(2*p),(p+1):(2*p)] <- -s1
diag(s) <- apply(s1,1,sum)[1]

ed <- eigen(s)
evals <- ed$values; evals[abs(evals) < 1E-8] <- 0
indx <- which(evals == 0)
evecs <- ed$vectors[,indx]

pdf('spectclust.pdf', height=4, width=8)

par(mfrow=c(1,2), mar=c(4.5,4.5,2.5,.5))
plot(evals, pch=4, col=3, main='eigenvalues')

plot(evecs[,1], pch=1, col=2, , main='eigenvectors', xlab='Index', 
	ylab='evecs with 0 eval')
points(evecs[,2], pch=16, col=4)
legend('center', c(expression(v[1]),expression(v[2])), col=c(2,4), 
	pch=c(1,16), bty='n')

dev.off()


##
##PCA
##
n <- 500
s <- 0.5
rho = 0.9

Ip <- diag(c(1,1))
Sig <- matrix(c(1, rho, rho, 1), 2, 2)

X <- mvrnorm(n=n, mu=c(0,0), Sigma=s*Sig)

pc.1 <- prcomp(X, center=TRUE, retx=TRUE)
pc.2 <- princomp(X, center=TRUE)
Z <- X %*% pc.1$rot

b1 <- pc.1$rot[1,1]/pc.1$rot[1,1]
b2 <- pc.1$rot[1,2]/pc.1$rot[2,2]	#not quite perp, prob rounding err
#b2 <- -1/b1

pdf('pcaplot.pdf', width=4, height=4)
par(mar=c(4.5,4.5,1,1))
plot(X, xlab=expression(X[1]), ylab=expression(X[2]), cex=0.6, col=4)

abline(a=0, b=b1, lty=1, lwd=2, col=2)
#abline(a=0, b=b2, lty=1, lwd=2, col=2)
segments(x0=-d,y0=-d*b2,x1=d,y1=d*b2, lwd=2, col=2)
dev.off()


pdf('pcaplot_wrot.pdf', width=7.5, height=4)
par(mfrow=c(1,2), mar=c(4.5,4.5,1,1))

plot(X, xlab=expression(X[1]), ylab=expression(X[2]), cex=0.6, col=4)
abline(a=0, b=b1, lty=1, lwd=2, col=2)
#abline(a=0, b=b2, lty=1, lwd=2, col=2)
segments(x0=-d,y0=-d*b2,x1=d,y1=d*b2, lwd=2, col=2)

plot(pc.1$x, ylim=c(-2.5,2.5), cex=0.6, col=4, xlab=expression(PC[1]), ylab=expression(PC[2]))

dev.off()


##
## PCA for USArrests
##
data(USArrests)

pdf('arrests-1.pdf')
pairs(USArrests, col=4, cex=0.7, upper.panel=panel.smooth)
dev.off()

states <- row.names(USArrests)
states[1:5]

apply(USArrests, 2, mean)
apply(USArrests, 2, var)

pc.out <- prcomp(USArrests, scale=TRUE)
print(pc.out$rot)

pc.out <- prcomp(USArrests, scale=TRUE, retx=TRUE)
names(pc.out)

pdf('arrests-2.pdf', height=5, width=7)
par(mar=c(4.5,4.5,.5,.5))
plot(pc.out$x[,1], pc.out$x[,2], type="n", xlab="1st PC", ylab="2nd PC")
text(pc.out$x[,1], pc.out$x[,2], labels=states, cex=0.7, col=4)
dev.off()

pdf('arrests-3.pdf', height=5, width=7)
par(mar=c(4.5,4.5,.5,.5))
plot(pc.out$rot, type="n", xlim=c(-0.6, -0.2))
text(pc.out$rot, names(USArrests), col=4)
dev.off()

screeplot(pc.out)

pdf('arrests-4.pdf', height=5, width=7)
par(mar=c(2.5,4.5,.5,.5))
plot(pc.out, col='orange', main='', type='l')
axis(1, at=c(0.7, 1.8, 3.2, 4.5), labels=paste('PC', 1:4))
dev.off()

print(pc.out$sdev^2)
print((pc.out$sdev^2)/sum(pc.out$sdev^2))

pdf('arrests-5.pdf', height=4, width=8)
par(mfrow=c(1,2), mar=c(4.5,4.5,.5,.5))
plot((pc.out$sdev^2)/sum(pc.out$sdev^2), xlab="PC", ylab="PVE", 
	col=4, ylim=c(0, 1), type='b')
plot(cumsum(pc.out$sdev^2)/sum(pc.out$sdev^2), xlab="PC", col=4, 
	ylab="Cumulative PVE", ylim=c(0, 1), type='b')
dev.off()

biplot(pc.out)

##
##MDS on the arrests data
##
library(cluster)

dist.1 <- daisy(USArrests)
dist.2 <- daisy(USArrests, metric="manhattan")

mds.1 <- cmdscale(dist.1)
mds.2 <- cmdscale(dist.2)
mds.3 <- isoMDS(dist.1)

pdf('arrests-mds.pdf', height=4, width=8)
par(mfrow=c(1,2), mar=c(4.5,4.5,2.5,.5))
x <- mds.1[,1]
y <- mds.1[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "n", main='Euclidean Dist')
text(x, y, labels = rownames(USArrests), col='blue', cex=0.7)

x <- mds.2[,1]
y <- mds.2[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "n", main='Manhattan Dist')
text(x, y, labels = rownames(USArrests), col='blue', cex=0.7)
dev.off()

x <- mds.3$points[,1]
y <- mds.3$points[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "n", main='Manhattan Dist')
text(x, y, labels = rownames(USArrests), col='blue', cex=0.7)


##
## PCA for NCI60 data
##
dat <- read.table('nci.data')
dim(dat)
dat <- t(dat)
dim(dat)

labs <- read.table('nci.info', skip=14)
labs <- as.vector(as.matrix(labs))
#table(labs)

##remove subtypes with only one member
indx <- names(table(labs))[table(labs) == 1]
dat2 <- dat[!(labs %in% indx),]
#dim(dat2)

##select subset of data and order the samples
labs2 <- labs[!(labs %in% indx)]
ord <- order(labs2)
dat2 <- dat2[ord,]
labs2 <- labs2[ord]

pc.nci <- prcomp(dat2, scale=TRUE, retx=TRUE)
#summary(pc.nci)

pdf('nci60pca-1.pdf', width=8, height=4)
par(mfrow=c(1,2), mar=c(4.5, 4.5, 1, .5))
plot((pc.nci$sdev^2)/sum(pc.nci$sdev^2), xlab="PC", ylab="PVE", 
	col=4, ylim=c(0, 1), type='b', cex=0.5)
#plot(cumsum(pc.nci$sdev^2)/sum(pc.nci$sdev^2), xlab="PC", col=4, 
#	ylab="Cumulative PVE", ylim=c(0, 1), type='b')
screeplot(pc.nci, type='l', col='blue', main='') #main='NCI60 data', 
dev.off()

mycol <- as.integer(as.factor(labs2))
pdf('nci60pca-2.pdf', height=3, width=9)
par(mfrow=c(1,3), mar=c(4.5,4.5,.5,.5))
plot(pc.nci$x[,1], pc.nci$x[,2], type="p", xlab="1st PC", ylab="2nd PC",
    cex=1.2, col=mycol, pch=mycol)
legend('topright', unique(labs2), col=unique(mycol), pch=unique(mycol), cex=0.6)
plot(pc.nci$x[,1], pc.nci$x[,3], type="p", xlab="1st PC", ylab="3rd PC",
    cex=1.2, col=mycol, pch=mycol)
plot(pc.nci$x[,2], pc.nci$x[,3], type="p", xlab="2nd PC", ylab="3rd PC",
    cex=1.2, col=mycol, pch=mycol)
dev.off()

##original data
tmp <- scale(dat2, center=F, scale=T)
d1 <- dist(tmp)

##project the data into a 4D space and perform k-menas clustering
d2 <- dist(pc.nci$x[,1:4])

par(mfrow=c(2,1), mar=c(0.5,3,3,0.5))
plot(hclust(d1, method="average"), labels=labs2, hang=0.2, col="blue",
	main="Original Data (Average Linkage)", xlab="", sub="", ylab="", 
	cex.main=1.5, cex.lab=0.6)
plot(hclust(d2, method="average"), labels=labs2, hang=0.2, 
	main="Data Projected to First 10 PCs (Average Linkage)", xlab="", 
	sub="", ylab="", cex.main=1.5, cex.lab=0.6, col="orange")

library(lattice)

tmp <- data.frame(pc.nci$x[,1], pc.nci$x[,2], pc.nci$x[,3])
names(tmp) <- c('PC1', 'PC2', 'PC3')
print(cloud(pc.nci$x[,1] ~ pc.nci$x[,2] * pc.nci$x[,3]), #data=tmp, 
print(cloud(PC1 ~ PC2 * PC3, data=tmp,
	groups=labs2, screen=list(z=20, x=-70), perspective = FALSE))
	key=list(title="NCI60 Data", x=.15, y=.85, corner = c(0,1),
	border=TRUE, points=Rows(trellis.par.get("superpose.symbol"), 1:3),
          text=list(levels(iris$Species)))))

pdf('nci60pca-4.pdf')
#print(cloud(pc.nci$x[,1] ~ pc.nci$x[,2] * pc.nci$x[,3], #data = tmp,
print(cloud(PC1 ~ PC2 * PC3, data=tmp,
	groups=labs2, screen=list(z=20, x=-70),
	perspective=FALSE, key=list(title="NCI60 Data", x=.15, y=.85, 
	corner=c(-2.95,0.45), border=TRUE,
	points=Rows(trellis.par.get("superpose.symbol"), 
	1:length(unique(labs2))), text = list(unique(labs2))) 
	))
dev.off()


##
##net inference
##
p <- 5
s <- 0.2
A <- 1*matrix(runif(p*p) > (1-s), p, p)
A <- 1*((A + t(A)) > 0); diag(A) <- 0

A
library(igraph)
g <- graph.adjacency(A, mode='undirected')

rho <- 0.6
Sinv <- rho*A; diag(Sinv) <- 1
Sig <- solve(Sinv)

n <- 1000
X <- chol(Sig) %*% matrix(rnorm(p*n), p, n)

dim(X)
S <- (1/n) * X %*% t(X)
dim(S)

Ahat <- solve(S)
Ahat <- 1*(Ahat > 1E-3); diag(Ahat) <- 0
ghat <- graph.adjacency(Ahat, mode='undirected')

pdf('gmex-2.pdf', width=6, height=4)
par(mfrow=c(1,2), mar=c(0.5, 0.5, 2.5, 0.5))
plot(g, vertex.label=c(1:5), layout=layout.circle(g), main='True Graph')
plot(ghat, vertex.label=c(1:5), layout=layout.circle(ghat), main='Est Graph')
dev.off()


n <- 7
X <- chol(Sig) %*% matrix(rnorm(p*n), p, n)

dim(X)
S <- (1/n) * X %*% t(X)
dim(S)

Ahat <- solve(S)
Ahat <- 1*(Ahat > 1E-3); diag(Ahat) <- 0
ghat <- graph.adjacency(Ahat, mode='undirected')

pdf('gmex-3.pdf', width=6, height=4)
par(mfrow=c(1,2), mar=c(0.5, 0.5, 2.5, 0.5))
plot(g, vertex.label=c(1:5), layout=layout.circle(g), main='True Graph')
plot(ghat, vertex.label=c(1:5), layout=layout.circle(ghat), main='Est Graph')
dev.off()

library(glasso)

sachs <- as.matrix(read.table("sachscov.txt"))
dim(sachs)

##glasso
est.1 <- glasso(s=sachs, rho=5, approx=FALSE, penalize.diagonal=FALSE)

##neighborhood selection
est.2 <- glasso(s=sachs, rho=5, approx=TRUE, penalize.diagonal=FALSE)



