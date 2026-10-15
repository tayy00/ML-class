library(MASS)
library(xtable)

setwd('/Users/ashojaie/Dropbox/ClassesEtc/UW/578Cspring14_ML4Omics/exam/')

data(crabs)

pairs(crabs)

dat <- crabs
# the first two are factors; turn to integers
dat[,1] <- as.integer(crabs[,1])
dat[,2] <- as.integer(crabs[,2])
dat <- dat[,-3]

pca.1 <- prcomp(dat, scale=T, retx=T)
round(summary(pca.1,3))
pcs <- round(pca.1$rot,3)
pcs[pcs < 0.1] = ''
print(pcs)
print(round(pca.1$rot,3))
print(round(summary(pca.1,3)$importance,4))
xtable(round(summary(pca.1,3)$importance,4))

plot(pca.1$x[,1], pca.1$x[,2], type="p", xlab="1st PC", ylab="2nd PC", col=dat$sex+1)
#plot(pca.1$x[,1], pca.1$x[,2], type="n", xlab="1st PC", ylab="2nd PC", col=dat$sex+1)
#text(pca.1$x[,1], pca.1$x[,2], labels=states, cex=0.7, col=4)

pdf('crabspca.pdf', width=12, height=6)
par(mfrow=c(1,2), mar=c(4.5,4.5,4.5,4.5))
biplot(pca.1, main='crabs data, scaled')
biplot(pca.2, main='crabs data, unscaled')
dev.off()

pdf('crabspca_orig.pdf', width=6, height=6)
pca.2 <- prcomp(dat, retx=T)
summary(pca.2)
print(round(pca.2$rot,3))
plot(pca.2$x[,1], pca.2$x[,2], type="p", xlab="1st PC", ylab="2nd PC", col=dat$sex+1)
dev.off()


d <- dist(dat)
mds <- cmdscale(d)
plot(mds, col=dat$sex+1)

set.seed(1234)
x <- rbind(matrix(rnorm(2*3),3,2), matrix(rnorm(2*3,2,1),3,2), matrix(rnorm(2*3,-2,1),3,2))
x[6,2] <- x[6,2] - 0.2
d <- dist(x)
hc <- hclust(d)

labs <- c(1:9)
labs[2:3] <- c('A','B')
labs[8:9] <- c('C','D')

pdf('hcaex.pdf', width=8, height=4)
par(mfrow=c(1,2))
plot(hc, xlab='', sub='')
plot(x, type='n', main='Locations of Points', xlab='x', ylab='y')
text(x=x[,1],y=x[,2],labels=labs)
dev.off()


