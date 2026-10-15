pdfplot = TRUE
par(mar=c(4,4,2,2))
if !(pdfplot) par(mfrow=c(2,3), mar=c(4,4,2,2))

## prcomp
pc.out <- prcomp(USArrests, scale=TRUE, retx=TRUE)

if(pdfplot){
	pdf('/Users/ashojaie/Dropbox/Teaching/UW/578C_ML4BBD/Spring2015/Lectures/Figures/arrests-2.pdf')
	plot(pc.out$x[,1], pc.out$x[,2], type="n", xlab="PC1", ylab="PC2")	#, xlim=c(-3,3), ylim=c(-3,3)
	text(pc.out$x[,1], pc.out$x[,2], labels=row.names(USArrests), col='deepskyblue3')
	dev.off()
}


if(pdfplot){
	pdf('/Users/ashojaie/Dropbox/Teaching/UW/578C_ML4BBD/Spring2015/Lectures/Figures/arrests-3.pdf')
	plot(pc.out$rot, type="n", xlim=c(-1,1), ylim=c(-1,1))
	text(pc.out$rot, names(USArrests), col='orange')
	dev.off()
}

if(pdfplot){
	pdf('/Users/ashojaie/Dropbox/Teaching/UW/578C_ML4BBD/Spring2015/Lectures/Figures/arrestbiplot.pdf')
	biplot(pc.out, pc.biplot=TRUE, col=c('deepskyblue3','orange'))	#, xlim=c(-3,3), ylim=c(-3,3)
	dev.off()
}



## princomp
pc.out <- princomp(USArrests, cor=TRUE, retx=TRUE)

plot(pc.out$scores[,1], pc.out$scores[,2], type="n", xlab="PC1", ylab="PC2")	#, xlim=c(-3,3), ylim=c(-3,3)
text(pc.out$scores[,1], pc.out$scores[,2], labels=row.names(USArrests), col='deepskyblue3')

plot(pc.out$loadings, type="n", xlim=c(-1,1), ylim=c(-1,1), xlab="PC1", ylab="PC2")
text(pc.out$loadings, names(USArrests), col='orange')

biplot(pc.out, pc.biplot=TRUE, col=c('deepskyblue3','orange'), xlab="PC1", ylab="PC2")	#, xlim=c(-3,3), ylim=c(-3,3)
