# SIS
(Iterative) Sure Independence Screening for Generalized Linear Models and Cox's Proportional Hazards Models:

Variable selection techniques are essential tools for model selection and estimation in high-dimensional statistical models. Through this publicly available ```R``` package, we provide a unified environment to carry out variable selection using iterative sure independence screening (SIS) and all of its variants in generalized linear models and the Cox proportional hazards model.

To install the current source version of the package, you should first have installed the ```glmnet```,  ```ncvreg``` and ```survival``` packages available from the <a href="https://cran.r-project.org/" target="_blank">CRAN</a> repository. Users with a ```Windows 10``` operating system simply have to type the follwing instruction within the current ```R``` session:

```
install.packages("/your_file_path/SIS_0.7-6.tar.gz", lib="/your_R_packages_library", repos=NULL, type="source")
```
A series of worked examples follow below. More documentation to come.
