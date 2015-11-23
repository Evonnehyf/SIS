# Iterative Sure Independence Screening Procedures on the Prostate Data

```
# The following code reproduces the training and test error results for the Prostate data set (Singh et al. 2002)
# in Table 6 of our paper "SIS: An R Package for Sure Independence Screening in Ultrahigh Dimensional Statistical Models".

# REPLICATION CODE TABLE 6
### Read Prostate Cancer Data Sets ###
library(SIS)
data("prostate.train", package = "SIS")
data("prostate.test", package = "SIS")
y1 = prostate.train[, dim(prostate.train)[2]]
x1 = as.matrix(prostate.train[, -dim(prostate.train)[2]])
y2 = prostate.test[, dim(prostate.test)[2]]
x2 = as.matrix(prostate.test[, -dim(prostate.test)[2]])
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
results.prostate.1 = matrix(0, nrow=tot.sim, ncol=3)
results.prostate.2 = matrix(0, nrow=tot.sim, ncol=3)
results.prostate.3 = matrix(0, nrow=tot.sim, ncol=3)
results.prostate.4 = matrix(0, nrow=tot.sim, ncol=3)
results.prostate.5 = matrix(0, nrow=tot.sim, ncol=3)
results.prostate.6 = matrix(0, nrow=tot.sim, ncol=3)
results.prostate.7 = matrix(0, nrow=tot.sim, ncol=3)
results.prostate.8 = matrix(0, nrow=tot.sim, ncol=3)

for(randSeed in 1:tot.sim){
set.seed(randSeed)

### Split Datasets ###
n = dim(x)[1]; aux = 1:n
ind.train1 = sample(aux[y == 0], 38, replace = FALSE)
ind.train2 = sample(aux[y == 1], 30, replace = FALSE)
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
results.prostate.1[randSeed,] = c(train.error,test.error,length(r$ix))

### VanISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.prostate.2[randSeed,] = c(train.error,test.error,length(r$ix))

### Var1ISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, varISIS="aggr", seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.prostate.3[randSeed,] = c(train.error,test.error,length(r$ix))

### Var2ISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, varISIS="cons", seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.prostate.4[randSeed,] = c(train.error,test.error,length(r$ix))

### permISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, perm=TRUE, q=q, seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.prostate.5[randSeed,] = c(train.error,test.error,length(r$ix))

### permgISIS###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, perm=TRUE, q=q, greedy=TRUE, seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.prostate.6[randSeed,] = c(train.error,test.error,length(r$ix))

### permVarISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, varISIS="cons", perm=TRUE, q=q, seed=as.numeric(randSeed), standardize=st)
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.prostate.7[randSeed,] = c(train.error,test.error,length(r$ix))

### permVargISIS ###
r = SIS(x.train, y.train , family="binomial", penalty=penalty, tune=tune, nsis=nsis, varISIS="cons", perm=TRUE, q=q, greedy=TRUE, seed=as.numeric(randSeed), standardize=st) 
train.error = length(which(y.train!=predict.SIS(r, x.train, type="class")))
test.error = length(which(y.test!=predict.SIS(r, x.test, type="class")))
results.prostate.8[randSeed,] = c(train.error,test.error,length(r$ix))
}

results.prostate = matrix(0, nrow=8, ncol=6)
results.prostate[1,c(1,3,5)] = apply(results.prostate.1,2,median); results.prostate[1,c(2,4,6)] = apply(results.prostate.1,2,IQR)/1.34
results.prostate[2,c(1,3,5)] = apply(results.prostate.2,2,median); results.prostate[2,c(2,4,6)] = apply(results.prostate.2,2,IQR)/1.34
results.prostate[3,c(1,3,5)] = apply(results.prostate.3,2,median); results.prostate[3,c(2,4,6)] = apply(results.prostate.3,2,IQR)/1.34
results.prostate[4,c(1,3,5)] = apply(results.prostate.4,2,median); results.prostate[4,c(2,4,6)] = apply(results.prostate.4,2,IQR)/1.34
results.prostate[5,c(1,3,5)] = apply(results.prostate.5,2,median); results.prostate[5,c(2,4,6)] = apply(results.prostate.5,2,IQR)/1.34
results.prostate[6,c(1,3,5)] = apply(results.prostate.6,2,median); results.prostate[6,c(2,4,6)] = apply(results.prostate.6,2,IQR)/1.34
results.prostate[7,c(1,3,5)] = apply(results.prostate.7,2,median); results.prostate[7,c(2,4,6)] = apply(results.prostate.7,2,IQR)/1.34
results.prostate[8,c(1,3,5)] = apply(results.prostate.8,2,median); results.prostate[8,c(2,4,6)] = apply(results.prostate.8,2,IQR)/1.34
results.prostate
```
