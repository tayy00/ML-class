library(MASS)

n = 500
s1 = matrix(c(1,0,0,1),ncol=2)
s2 = matrix(c(1,0.7,0.7,1),ncol=2)
s3 = matrix(c(1,-0.7,-0.7,1),ncol=2)
ss = list(s1,s2,s3)

set.seed(1)
#par(mfrow=c(1,3))
for(i in 1:3){
	x = mvrnorm(n=n, mu=c(0,0), Sigma=ss[[i]])
	evals = eigen(as.matrix(ss[[i]]), symmetric=TRUE)
	evecs = evals$vectors; evals = evals$values
	pdf(paste('~/Dropbox/Teaching/UW/578C_ML4BBD/Spring2015/Lectures/Figures/covmat-',i,'.pdf',sep=''), height=4, width=4)
#	par(mar=c(4,4,1,1))
	par(mar=c(1,1,1,1))
	plot(x, xlim=c(-3,3), ylim=c(-4,4), col='gray', axes=F, pch=16, cex=.7)#,xlab=expression(X[1]), ylab=expression(X[2]))
	abline(h=0)
	abline(v=0)
	arrows(x0=0,y0=0,x1=evecs[1,1]*evals[1]*2,y1=evecs[2,1]*evals[1]*2,
		length=0.10, col='blue', lwd=2)
	arrows(x0=0,y0=0,x1=evecs[1,2]*evals[2]*2,y1=evecs[2,2]*evals[2]*2,
		length=0.10, col='magenta', lwd=2)
	dev.off()
}
