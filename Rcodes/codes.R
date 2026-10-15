# EXAMPLE 1
xtr <- matrix(rnorm(100*100),ncol=100)
xte <- matrix(rnorm(100000*100),ncol=100)
beta <- c(rep(1,10),rep(0,90))
ytr <- xtr%*%beta + rnorm(100)
yte <- xte%*%beta + rnorm(100000)
rsq <- trainerr <- testerr <-  NULL
for(i in 2:100){
	mod <- lm(ytr~xtr[,1:i])
	rsq <- c(rsq,summary(mod)$r.squared)
	beta <- mod$coef[-1]
	intercept <- mod$coef[1]
	trainerr <- c(trainerr, mean((xtr[,1:i]%*%beta+intercept - ytr)^2))
	testerr <- c(testerr, mean((xte[,1:i]%*%beta+intercept - yte)^2))
}
pdf("ex1.pdf",width=10,height=4)
par(mfrow=c(1,3))
plot(2:100,rsq, xlab='Number of Variables', ylab="R Squared", log="y")
abline(v=10,col="red")
plot(2:100,trainerr, xlab='Number of Variables', ylab="Training Error",log="y")
abline(v=10,col="red")
plot(2:100,testerr, xlab='Number of Variables', ylab="Test Error",log="y")
abline(v=10,col="red")
dev.off()

# EXAMPLE 2
xtr <- matrix(rnorm(100*100),ncol=100)
beta <- c(rep(1,10),rep(0,90))
ytr <- xtr%*%beta + rnorm(100)
cv.err <- NULL
for(i in 2:50){
	dat <- data.frame(x=xtr[,1:i],y=ytr)
	mod <- glm(y~.,data=dat)
	cv.err <- c(cv.err, cv.glm(dat,mod,K=6)$delta[1])
}
pdf("ex2.pdf")
plot(2:50, cv.err, xlab="Number of Variables", ylab="6-Fold CV Error", log="y",cex.lab=1.3)
abline(v=10, col="red")
dev.off()

# EXAMPLE 3
xtr <- matrix(rnorm(100*100),ncol=100)
beta <- c(rep(1,10),rep(0,90))
ytr <- xtr%*%beta + rnorm(100)
cors <- cor(xtr,ytr)
whichers <- which(abs(cors)>.2)
mod <- lm(ytr~xtr[,whichers])
print(summary(mod))


# EXAMPLE 4
xtr <- matrix(rnorm(100*100),ncol=100)
beta <- c(rep(1,10),rep(0,90))
ytr <- xtr%*%beta + rnorm(100)
library(leaps)
out <- regsubsets(xtr,ytr,nvmax=30,method="forward")
print(summary(out))
print(coef(out,1:10))

# EXAMPLE 4
xtr <- matrix(rnorm(100*100),ncol=100)
beta <- c(rep(1,10),rep(0,90))
ytr <- xtr%*%beta + rnorm(100)
library(glmnet)
cv.out <- cv.glmnet(xtr,ytr,alpha=0,nfolds=5)
print(cv.out$cvm)
pdf("exridge.pdf")
plot(cv.out)
dev.off()
cat("CV Errors", cv.out$cvm,fill=TRUE)
cat("Lambda with smallest CV Error", cv.out$lambda[which.min(cv.out$cvm)],fill=TRUE)
cat("Coefficients", as.numeric(coef(cv.out)),fill=TRUE)
cat("Number of Zero Coefficients",sum(abs(coef(cv.out))<1e-8),fill=TRUE)


# EXAMPLE 5
xtr <- matrix(rnorm(100*100),ncol=100)
beta <- c(rep(1,10),rep(0,90))
ytr <- xtr%*%beta + rnorm(100)
library(glmnet)
cv.out <- cv.glmnet(xtr,ytr,alpha=1,nfolds=5)
print(cv.out$cvm)
pdf("exlasso.pdf")
plot(cv.out)
dev.off()
cat("CV Errors", cv.out$cvm,fill=TRUE)
cat("Lambda with smallest CV Error", cv.out$lambda[which.min(cv.out$cvm)],fill=TRUE)
cat("Coefficients", as.numeric(coef(cv.out)),fill=TRUE)
cat("Number of Zero Coefficients",sum(abs(coef(cv.out))<1e-8),fill=TRUE)

# EXAMPLE 6
xtr <- matrix(rnorm(100*100),ncol=100)
beta <- c(rep(1,10),rep(0,90))
ytr <- xtr%*%beta + rnorm(100)
library(pls)
out <- pcr(ytr~xtr,scale=TRUE,validation="CV")
summary(out)
validationplot(out,val.type="MSEP")

pdf("expcr.pdf")
validationplot(out,val.type="MSEP")
dev.off()

# EXAMPLE 7
xtr <- matrix(rnorm(1000*20),ncol=20)
beta <- c(rep(1,10),rep(0,10))
ytr <- 1*((xtr%*%beta + .2*rnorm(1000)) >= 0)
mod <- glm(ytr~xtr,family="binomial")
print(summary(mod))

# EXAMPLE 8
xtr <- matrix(rnorm(1000*20),ncol=20)
beta <- c(rep(1,5),rep(0,15))
ytr <- 1*((xtr%*%beta + .5*rnorm(1000)) >= 0)
cv.out <- cv.glmnet(xtr, ytr, family="binomial", alpha=1)
plot(cv.out)

# EXAMPLE 9
x <- matrix(rnorm(1000*50),ncol=50)
y <- sample(c(0,1),50,rep=TRUE)
ps <- NULL
for(i in 1:1000) ps <- c(ps, t.test(x[i,y==0],x[i,y==1])$p.value)
cat("Around 5% of p-values are below 0.05:", 
mean(ps<.05),fill=TRUE)
fdrs.bh <- p.adjust(ps, method="BH")
plot(ps,fdrs.bh)
plot(fdrs.bh)

# EXAMPLE 10
x <- matrix(rnorm(1000*50),ncol=50)
y <- sample(c(0,1),50,rep=TRUE)
x[1:100,y==0] <- x[1:100,y==0] + 1
ps <- NULL
for(i in 1:1000) ps <- c(ps, t.test(x[i,y==0],x[i,y==1])$p.value)
cat("Way more than 5% of p-values are below 0.05:", 
mean(ps<.05),fill=TRUE)
fdrs.bh <- p.adjust(ps, method="BH")
plot(ps,fdrs.bh)
plot(fdrs.bh)
cat("Number of Tests with FDR below 0.4:", sum(fdrs.bh<0.4), fill=TRUE)
cat("Compute the BH FDR Directly:", max(which(sort(ps,decreasing=FALSE) < .4*(1:1000)/1000)), fill=TRUE)
pdf("bh.pdf")
plot(sort(ps,decreasing=FALSE),ylab="P-Values")
abline(a=0, b=0.4/1000,col="red")
dev.off()


##
## Multiple Hypothesis Testing
##
setwd('/Users/ashojaie/Dropbox/ClassesEtc/UW/578Cspring14_ML4Omics/Lectures/Figures')
M = 10000
n = 100

##
## First, no effects at all
##
y = as.factor(c(rep("normal",n/2),rep("cancer",n/2)))
ttst.p <- function(x){
	t.test(x[y == "cancer"], x[y == "normal"])$p.value
}
ttst.s <- function(x){
	t.test(x[y == "cancer"], x[y == "normal"])$statistic
}

set.seed(1)
expdat <- matrix(rnorm(n*M), nrow=M)		#NOTE: nothing's significant
pvals <- apply(expdat, 1, ttst.p)
stats <- apply(expdat, 1, ttst.s)

min(pvals)

pdf("HDtesting-1.pdf")#,width=6,height=4)
hist(pvals, freq=F, main="Histogram of p-values", xlab='p-values', breaks=100)
abline(h=1, col=2, lwd=2, lty=2)
dev.off()

pdf("HDtesting-2.pdf",width=6,height=4)
hist(stats, freq=F, breaks=100, main='t-statistic Histogram', xlab='t-statistic')
curve(dnorm(x, mean=0, sd=1), add=TRUE, col=2, lty=2, lwd=2)
dev.off()

pdf("HDtesting-3.pdf")#,width=5,height=4)
plot(c(1:M)/M,sort(pvals),xlab='expected',ylab='observed',
	main='U[0,1] QQ-plot')
abline(0,1,col=2,lty=2,lwd=2)
dev.off()
		
indx <- which.min(pvals)
t.test(expdat[indx,y=='cancer'], expdat[indx,y=='normal'])

maxstat <- function(m){
	expdat <- matrix(rnorm(n*m),nrow=m)
	tststat <- apply(expdat, 1, ttst.s)
	return(max(tststat))
}

set.seed(2)
max.stats = replicate(100, maxstat(1000))
mean(max.stats)
hist(max.stats)

##
## What happens if there are some effects?
##
M1 = 1000
set.seed(3)
shift = c(runif(M1/2,0.1,0.5), runif(M1/2,-0.5,-0.1))
shift = matrix(shift, nrow=M1, ncol=n/2, byrow=F)
expdat2 = expdat
expdat2[1:M1,1:(n/2)] = expdat2[1:M1,1:(n/2)] + shift

pvals2 <- apply(expdat2, 1, ttst.p)
stats2 <- apply(expdat2, 1, ttst.s)

min(pvals2)

pdf("HDtesting-4.pdf")#,width=6,height=4)
hist(pvals2, freq=F, main="Histogram of p-values", xlab='p-values', breaks=100)
abline(h=1, col=2, lwd=2, lty=2)
dev.off()

pdf("HDtesting-5.pdf")#,width=6,height=4)
hist(stats2, freq=F, ylim=c(0,.4), breaks=100, main='t-statistic Histogram', xlab='t-statistic')
curve(dnorm(x, mean=0, sd=1), add=TRUE, col=2, lty=2, lwd=2)
dev.off()

pdf("HDtesting-6.pdf")#,width=5,height=4)
plot(c(1:M)/M,sort(pvals2),xlab='expected',ylab='observed',
	main='U[0,1] QQ-plot')
abline(0,1,col=2,lty=2,lwd=2)
dev.off()

##
## FDR in R
##
?p.adjust

# everything null
x <- matrix(rnorm(1000*50),ncol=50)
y <- sample(c(0,1),50,rep=TRUE)
ps <- NULL
for(i in 1:1000) ps <- c(ps, t.test(x[i,y==0],x[i,y==1])$p.value)
cat("Around 5% of p-values are below 0.05:", 
mean(ps<.05),fill=TRUE)
fdrs.bh <- p.adjust(ps, method="BH")
plot(ps,fdrs.bh)
plot(fdrs.bh)

# some non-null
x <- matrix(rnorm(1000*50),ncol=50)
y <- sample(c(0,1),50,rep=TRUE)
x[1:100,y==0] <- x[1:100,y==0] + 1
ps <- NULL
for(i in 1:1000) ps <- c(ps, t.test(x[i,y==0],x[i,y==1])$p.value)
cat("Way more than 5% of p-values are below 0.05:", 
mean(ps<.05),fill=TRUE)
fdrs.bh <- p.adjust(ps, method="BH")
plot(ps,fdrs.bh)
plot(fdrs.bh)
cat("Number of Tests with FDR below 0.4:", sum(fdrs.bh<0.4), fill=TRUE)
cat("Compute the BH FDR Directly:", max(which(sort(ps,decreasing=FALSE) < .4*(1:1000)/1000)), fill=TRUE)
pdf("bh.pdf", height=4, width=5)
plot(sort(ps,decreasing=FALSE),ylab="P-Values",cex=0.7)
abline(a=0, b=0.4/1000,col=2, lwd=2, lty=2)
abline(h=0.1/1000, col=4, lwd=2, lty=4)
legend('topleft',c('Bonferroni','FDR'),lty=c(4,2),col=c(4,2),lwd=2,bty='n')
dev.off()






