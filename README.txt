1) simudata.R generates simulated data based on the number of features, 
sample size, the number of external variables, the number of neighbors of
an external variable

2) simu100AIC implements logistic regression Lasso to generate networks

3) simu100AICGGM implements linear regression Lasso to generate networks

4) In every file in this format, n#1-s#2-c#3-d#4-nwc0-corr60vary.R, 
#1 is the number of features, #2 is sample size, #3 is the number of 
external variables and #4 is the number of neighbors of an external variable 
This file generates simulated data, implements Lasso and produces estimated 
precision matrix outcomes.   

5) tptn.R computes accuracy in detecting latent and original interactions
of each case.  