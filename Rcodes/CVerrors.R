## 
## this code compares the bias and variance of estiamtes of test error
##

library(boot)

set.seed(1)
B <- 100
p <- 10
n <- 100

xtr <- matrix(rnorm(p*n),n,p)
betas <- c(rep(1,5),rep(0,p-5))
ytr <- xtr %*% betas + rnorm(n)

dat <- data.frame(ytr,xtr)
lm.fit <- glm(ytr~., data=dat)

err.LOOCV <- numeric(B)
err.valid <- numeric(B)

for(b in 1:B){
  xtr <- matrix(rnorm(p*n),n,p)
  ytr <- xtr %*% betas + rnorm(n)
  
  dat <- data.frame(ytr,xtr)
  lm.fit <- glm(ytr~., data=dat)

  dat <- data.frame(ytr,xtr)
  lm.fit <- glm(ytr~., data=dat)
    
  # LOOCV
  err.LOOCV[b] <- cv.glm(data=dat, glmfit=lm.fit, K=n)$delta[1]
  
  # LOOCV
  err.valid[b] <- cv.glm(data=dat, glmfit=lm.fit, K=10)$delta[1]
}

boxplot(cbind(err.LOOCV, err.valid))
