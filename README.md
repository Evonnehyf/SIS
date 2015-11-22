# SIS
(Iterative) Sure Independence Screening for Generalized Linear Models and Cox's Proportional Hazards Models:

Variable selection techniques are essential tools for model selection and estimation in high-dimensional statistical models. Through this publicly available ```R``` package, we provide a unified environment to carry out variable selection using iterative sure independence screening (SIS) and all of its variants in generalized linear models and the Cox proportional hazards model.

To install the current source version of the package, you should first have installed the ```glmnet```,  ```ncvreg``` and ```survival``` packages available from the <a href="https://cran.r-project.org/" target="_blank">CRAN</a> repository. Users with a ```Windows 10``` operating system simply have to type the follwing instruction within the current ```R``` session:

```
install.packages("/your_file_path/SIS_0.7-6.tar.gz", lib="/your_R_packages_library", repos=NULL, type="source")
```
A series of worked examples follow below. More documentation to come.

```
library(SIS)
set.seed(0)
n = 400; p = 50; rho = 0.5
corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
corrmat[,4] = sqrt(rho)
corrmat[4, ] = sqrt(rho)
corrmat[4,4] = 1
corrmat[,5] = 0
corrmat[5, ] = 0
corrmat[5,5] = 1
cholmat = chol(corrmat)
x = matrix(rnorm(n*p, mean=0, sd=1), n, p)
x = x\%*\%cholmat

# gaussian response 
set.seed(1)
b = c(4,4,4,-6*sqrt(2),4/3)
y=x[, 1:5]\%*\%b + rnorm(n)
model11=SIS(x, y, family="gaussian", tune="bic")
model12=SIS(x, y, family="gaussian", tune="bic", varISIS="aggr", seed=11)
model11$ix
model12$ix

# binary response 
set.seed(2)
feta = x[, 1:5]\%*\%b; fprob = exp(feta)/(1+exp(feta))
y = rbinom(n, 1, fprob)
model21=SIS(x, y, family="binomial", tune="bic")
model22=SIS(x, y, family="binomial", tune="bic", varISIS="aggr", seed=21)
model21$ix
model22$ix

# poisson response
set.seed(3)
b = c(0.6,0.6,0.6,-0.9*sqrt(2))
myrates = exp(x[, 1:4]\%*\%b)
y = rpois(n, myrates)
model31=SIS(x, y, family="poisson", tune="bic", perm=TRUE, q=0.9, 
            greedy=TRUE, seed=31)
model32=SIS(x, y, family="poisson", tune="bic", varISIS="aggr", 
            perm=TRUE, q=0.9, seed=32)
model31$ix
model32$ix

# Cox model
set.seed(4)
b = c(4,4,4,-6*sqrt(2),4/3)
myrates = exp(x[, 1:5]\%*\%b)
Sur = rexp(n,myrates); CT = rexp(n,0.1)
Z = pmin(Sur,CT); ind = as.numeric(Sur<=CT)
y = Surv(Z,ind)
model41=SIS(x, y, family="cox", penalty="lasso", tune="bic", 
            varISIS="aggr", seed=41)
model42=SIS(x, y, family="cox", penalty="lasso", tune="bic", 
            varISIS="cons", seed=41)
model41$ix
model42$ix
```
