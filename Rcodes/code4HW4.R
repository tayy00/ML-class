library(lars)
library(tree)
library(randomForest)
library(gbm)
data(diabetes)
attach(diabetes)
set.seed(1234)
train=sample(1:nrow(diabetes), 300)
bag.diabetes=randomForest(x=x2, y=y, data=diabetes, subset=train, n.trees=500, interaction.depth=4)

boost.diabetes=gbm(x=x2,y=y,data=diabetes, subset=train,,n.trees=5000,interaction.depth=4)


set.seed(1)
boost.boston=gbm(medv~.,data=Boston,distribution="gaussian",n.trees=5000,interaction.depth=4)
summary(boost.boston)

xx <- matrix(x2,nrow(x2),ncol(x2))
newdiab <- data.frame(y,xx)
boost.boston=gbm(y~.,data=newdiab,distribution="gaussian",n.trees=500,interaction.depth=2)

boost.boston=gbm(y~x2,data=diabetes,distribution="gaussian",n.trees=500,interaction.depth=2)



