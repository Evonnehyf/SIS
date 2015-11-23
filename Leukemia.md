# Iterative Sure Independence Screening Procedures on the Leukemia Data

```
# The following code reproduces the training and test error results for the Leukemia data set (Golub et al. 1999)
# in Table 6 of our paper "SIS: An R Package for Sure Independence Screening in Ultrahigh Dimensional Statistical Models".

# REPLICATION CODE TABLE 6
### Read Leukemia Data Sets ###
library(SIS)
data("leukemia.train", package = "SIS")
data("leukemia.test", package = "SIS")
y1 = leukemia.train[, dim(leukemia.train)[2]]
x1 = as.matrix(leukemia.train[, -dim(leukemia.train)[2]])
y2 = leukemia.test[, dim(leukemia.test)[2]]
x2 = as.matrix(leukemia.test[, -dim(leukemia.test)[2]])
### Combine Datasets ###
x = rbind(x1, x2)
y = c(y1, y2)

### Global Parameters ###
penalty = "lasso"
tune = "cv"
nsis = 100
q = 0.95
st = FALSE
tot.sim = 100
results.leukemia.1 = matrix(0, nrow=tot.sim, ncol=3)
results.leukemia.2 = matrix(0, nrow=tot.sim, ncol=3)
results.leukemia.3 = matrix(0, nrow=tot.sim, ncol=3)
results.leukemia.4 = matrix(0, nrow=tot.sim, ncol=3)
results.leukemia.5 = matrix(0, nrow=tot.sim, ncol=3)
results.leukemia.6 = matrix(0, nrow=tot.sim, ncol=3)
results.leukemia.7 = matrix(0, nrow=tot.sim, ncol=3)
results.leukemia.8 = matrix(0, nrow=tot.sim, ncol=3)

for(randSeed in 1:tot.sim){
set.seed(randSeed)

### Split Datasets ###
n = dim(x)[1]; aux = 1:n
ind.train1 = sample(aux[y == 0], 23, replace = FALSE)
ind.train2 = sample(aux[y == 1], 13, replace = FALSE)
ind.train = c(ind.train1, ind.train2)
x.train = scale(x[ind.train,])
y.train = y[ind.train]
ind.test1 = setdiff(aux[y == 0], ind.train1)
ind.test2 = setdiff(aux[y == 1], ind.train2)
ind.test = c(ind.test1, ind.test2)
x.test = scale(x[ind.test,])
y.test = y[ind.test]

### VanSIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, iter=FALSE, standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.leukemia.1[randSeed,] = c(train.error,test.error,length(r$ix))

### VanISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.leukemia.2[randSeed,] = c(train.error,test.error,length(r$ix))

### Var1ISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, varISIS="aggr", seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.leukemia.3[randSeed,] = c(train.error,test.error,length(r$ix))

### Var2ISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, varISIS="cons", seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.leukemia.4[randSeed,] = c(train.error,test.error,length(r$ix))

### permISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, perm=TRUE, q=q, seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.leukemia.5[randSeed,] = c(train.error,test.error,length(r$ix))

### permgISIS###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, perm=TRUE, q=q, greedy=TRUE, seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.leukemia.6[randSeed,] = c(train.error,test.error,length(r$ix))

### permVarISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, varISIS="cons", perm=TRUE, q=q, seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.leukemia.7[randSeed,] = c(train.error,test.error,length(r$ix))

### permVargISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, varISIS="cons", perm=TRUE, q=q, greedy=TRUE, seed=as.numeric(randSeed), standardize=st) 
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.leukemia.8[randSeed,] = c(train.error,test.error,length(r$ix))
}

results.leukemia = matrix(0, nrow=8, ncol=6)
results.leukemia[1,c(1,3,5)] = apply(results.leukemia.1,2,median); results.leukemia[1,c(2,4,6)] = apply(results.leukemia.1,2,IQR)/1.34
results.leukemia[2,c(1,3,5)] = apply(results.leukemia.2,2,median); results.leukemia[2,c(2,4,6)] = apply(results.leukemia.2,2,IQR)/1.34
results.leukemia[3,c(1,3,5)] = apply(results.leukemia.3,2,median); results.leukemia[3,c(2,4,6)] = apply(results.leukemia.3,2,IQR)/1.34
results.leukemia[4,c(1,3,5)] = apply(results.leukemia.4,2,median); results.leukemia[4,c(2,4,6)] = apply(results.leukemia.4,2,IQR)/1.34
results.leukemia[5,c(1,3,5)] = apply(results.leukemia.5,2,median); results.leukemia[5,c(2,4,6)] = apply(results.leukemia.5,2,IQR)/1.34
results.leukemia[6,c(1,3,5)] = apply(results.leukemia.6,2,median); results.leukemia[6,c(2,4,6)] = apply(results.leukemia.6,2,IQR)/1.34
results.leukemia[7,c(1,3,5)] = apply(results.leukemia.7,2,median); results.leukemia[7,c(2,4,6)] = apply(results.leukemia.7,2,IQR)/1.34
results.leukemia[8,c(1,3,5)] = apply(results.leukemia.8,2,median); results.leukemia[8,c(2,4,6)] = apply(results.leukemia.8,2,IQR)/1.34
results.leukemia
```
