setwd('/Users/ashojaie/Dropbox/Teaching/UW/BIOST546_ML4BBD/Spring2016/Rcodes/')
dat <- read.csv('Ch10Ex11.csv', header=F)

d1 <- dist(1-cor(dat))
s1 <- hclust(d1,method='single')
c1 <- hclust(d1,method='complete')
a1 <- hclust(d1,method='average')


d2 <- as.dist(1-cor(dat))
s2 <- hclust(d2,method='single')
c2 <- hclust(d2,method='complete')
a2 <- hclust(d2,method='average')


par(mfrow=c(2,3))
plot(s2, main='single, as.dist()')
plot(c2, main='complete, as.dist()')
plot(a2, main='average, as.dsit()')

plot(s2, main='single, as.dist()', hang=-1)
plot(c2, main='complete, as.dist()', hang=-1)
plot(a2, main='average, as.dsit()', hang=-1)


plot(s1, main='single, dist()')
plot(c1, main='complete, dist()')
plot(a1, main='average, dist()')

