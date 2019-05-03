library(corpcor)
library(MASS)
library(bindata)
simudataBIN <- function(nodes, samples, covars,degrees, covhubcor, normalcor, binary, starshape){
#nodes indicate the number of features
# samples are number of samples
#covars number of covariates
# covhubcor is the correlation strength between covariates and some nodes
# normalcor is the correaltion among all other nodes
# degrees is the number of neighbors for covariates 


sigma= matrix(normalcor,(nodes+covars),(nodes+covars))

for (i in 1:covars){
   
sigma[i,(covars+1+(i-1)*degrees):((covars+1+(i-1)*degrees)+degrees-1)]=covhubcor
sigma[(covars+1+(i-1)*degrees):((covars+1+(i-1)*degrees)+degrees-1),i]=covhubcor
}

#intra-interaction starting index
intrastart= covars*degrees+covars
newcovars=3
newdegrees=3
newcor=0.4
if (starshape){

for (j in 1:newcovars){

sigma[intrastart+j,(intrastart+newcovars+1+(j-1)*newdegrees):((intrastart+newcovars+1+(j-1)*newdegrees)+newdegrees-1)]=newcor
sigma[(intrastart+newcovars+1+(j-1)*newdegrees):((intrastart+newcovars+1+(j-1)*newdegrees)+ newdegrees-1),intrastart+j]=newcor
}

}

diag(sigma)=1
sigma=make.positive.definite(sigma)

out = rmvbin(samples,commonprob=sigma)
#
#out <- mvrnorm(samples, mu = rep(0,(nodes+covars)), Sigma = sigma,
#               empirical = TRUE)


#if(binary){
#for (i in 1:(nodes+covars)){
#  out[,i]= (out[,i]> median(out[,i])+0)
#}
#}

return(out)
}