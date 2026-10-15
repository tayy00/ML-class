#### Used the liposarcoma samples from:
#### http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14533

#### The data is cancer samples with one of 2 different mechanisms for replacing Telomeres:
#### Without replacing Telomeres, fast replicating cells will quickly die off
#### Telo --> Telomerase is the standard mechanism
#### Alt ---> Is the alternative mechanism (poorly understood)

#### Can we find genes which are over/under expressed in Telo vs Alt? (ideally over expressed in Alt)

#### Grabbing the data (GEO didn't work, so I did it by hand)

expression.list <- list()
dataset <- c(363613 + 0:17)
 id <- 26946 + 0:17
for(i in 1:18){
  site <- paste("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSM",dataset[i]
                          ,"&id=",id[i], "&db=GeoDb_blob35", sep = "")
  expression.list[[i]] <- read.table(site, skip = 22, sep = "\t",
                          na.strings = "", fill = TRUE, as.is = TRUE,
                           col.names = c("probe_id","value"))
cat(i)
}

###### These next two steps may be unnecessary (the expression platforms were the same) but better safe than sorry!

#### Taking the intersection of the probe-names across samples
use.names <- expression.list[[1]][,1]
for(i in length(expression.list)){
  use.names <- intersect(use.names, expression.list[[i]][,1])
}

#### Matching up the expression values by probe across samples
exp.matrix <- matrix(0,nrow = length(use.names), ncol = length(expression.list))
for(i in 1:length(expression.list)){
  match.ind <- match(use.names, expression.list[[i]][,1])
  exp.matrix[,i] <- expression.list[[i]][match.ind,2]
}

### Storing the response (you can see this from the website)
immort <- c(rep("alt",10),rep("telo",8))

## need to impute the missing values, for now drop them
indx <- (apply(is.na(exp.matrix), 1, sum) > 0)
sum(indx)
X <- exp.matrix[!indx,]

### Running our tests!
p.vals <- apply(X,1,function(x){t.test(x[immort=="alt"], x[immort=="telo"])$p.val})
t.stats <- apply(X,1,function(x){t.test(x[immort=="alt"], x[immort=="telo"])$statistic})

hist(p.vals, freq = FALSE, xlab = "p-value", main = "p-value Histogram", breaks = 100)
abline(h = 1, col = 2, lty = 2, lwd = 5)

hist(t.stats, freq = FALSE, breaks = 100, xlab="t-statistic", main = "t-statistic Histogram")
lines(seq(-4,4,by = 0.01), dt(seq(-4,4,by = 0.01),15), col = 2, lty = 2, lwd = 5)

plot(1:2000,sort(p.adjust(p.vals, method = "BH"), decreasing = FALSE)[1:2000], xlab = "number rejected",ylab ="FDR control", type = 'l', lwd = 3)
