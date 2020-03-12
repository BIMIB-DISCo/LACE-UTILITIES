# set the working directory
setwd("/data_directory/BimiB/share/LACE/Breast_Cancer/real_data/p494")

# load required scripts
library("parallel")
library("Rfast")
source("R/frontend.R")
source("R/inference.R")

# load data
data = readRDS("final_LACE_input_data/D_breast.rds")
D = list()
D[[1]] = data[[1]]
D[[2]] = data[[2]]

# define likelihood weights
lik_w1 = (1-(dim(D[[1]])[1]/(dim(D[[1]])[1]+dim(D[[2]])[1])))
lik_w2 = (1-(dim(D[[2]])[1]/(dim(D[[1]])[1]+dim(D[[2]])[1])))
lik_w = c(lik_w1,lik_w2)
lik_w = lik_w/sum(lik_w)

# define parameters for grid search

alpha1 = list(c(0.001,0.001))
alpha2 = list(c(0.002,0.001))
alpha3 = list(c(0.001,0.002))
alpha4 = list(c(0.001,0.001))
alpha5 = list(c(0.002,0.001))
alpha6 = list(c(0.001,0.002))
alpha7 = list(c(0.010,0.010))
alpha8 = list(c(0.020,0.010))
alpha9 = list(c(0.010,0.020))
alpha10 = list(c(0.010,0.010))
alpha11 = list(c(0.020,0.010))
alpha12 = list(c(0.010,0.020))
alpha13 = list(c(0.050,0.050))
alpha14 = list(c(0.100,0.050))
alpha15 = list(c(0.050,0.100))
alpha16 = list(c(0.050,0.050))
alpha17 = list(c(0.100,0.050))
alpha18 = list(c(0.050,0.100))

beta1 = list(c(0.001,0.001))
beta2 = list(c(0.002,0.001))
beta3 = list(c(0.001,0.002))
beta4 = list(c(0.010,0.010))
beta5 = list(c(0.020,0.010))
beta6 = list(c(0.010,0.020))
beta7 = list(c(0.050,0.050))
beta8 = list(c(0.100,0.050))
beta9 = list(c(0.050,0.100))
beta10 = list(c(0.001,0.001))
beta11 = list(c(0.002,0.001))
beta12 = list(c(0.001,0.002))
beta13 = list(c(0.010,0.010))
beta14 = list(c(0.020,0.010))
beta15 = list(c(0.010,0.020))
beta16 = list(c(0.050,0.050))
beta17 = list(c(0.100,0.050))
beta18 = list(c(0.050,0.100))

alpha = c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,alpha11,alpha12,alpha13,alpha14,alpha15,alpha16,alpha17,alpha18)
beta = c(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8,beta9,beta10,beta11,beta12,beta13,beta14,beta15,beta16,beta17,beta18)

# perform inference
inference = LACE(D=D,lik_w=lik_w,alpha=alpha,beta=beta,num_rs=50,num_iter=10000,n_try_bs=500,marginalize=FALSE,num_processes=18,seed=5876033,verbose=TRUE)

# save the results
save(inference,file="results/inference.RData")
