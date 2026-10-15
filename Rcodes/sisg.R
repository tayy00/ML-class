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





