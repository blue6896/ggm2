set.seed(505)
library(corpcor)
library(corrplot)
library(MASS)
# numweakcor= nwc => degrees * nwc / 5   nwc=2 => 40% level, nwc=4 => 80% level 

simudata <- function(nodes, samples, covars,degrees, covhubcor, normalcor, binary, starshape,randomcorr, numweakcor,dich,valmet){
#nodes indicate the number of features
# samples are number of samples
#covars number of covariates
# covhubcor is the correlation strength between covariates and some nodes
# normalcor is the correaltion among all other nodes
# degrees is the number of neighbors for covariates 

#numweakcor= degrees * numweakcor / 5 

sigma= matrix(normalcor,(nodes+covars),(nodes+covars))
sigma= matrix(runif((nodes+covars)*(nodes+covars),0,0.1),(nodes+covars),(nodes+covars))


for (i in 1:covars){

if (randomcorr | (covhubcor > 0.6)){
  # for (iii in 1:degrees){
  #     sigma[i,(covars+1+(i-1)*degrees +iii-1)]=(0.1+ (iii-1)*(covhubcor/degrees))
  #     sigma[(covars+1+(i-1)*degrees+iii-1),i]= (0.1+(iii-1)*(covhubcor/degrees))
  # }
  if(numweakcor>0){
    for (iii in 1:numweakcor){
       sigma[i,(covars+1+(i-1)*degrees +iii-1)]=runif(1,0.10,0.20)        #0.13#(0.1+ (iii-1)*(covhubcor/degrees))
       sigma[(covars+1+(i-1)*degrees+iii-1),i]=runif(1,0.10,0.20)  # 0.13#(0.1+(iii-1)*(covhubcor/degrees))
   }
  }
    for (iii in (numweakcor+1):degrees){
       sigma[i,(covars+1+(i-1)*degrees +iii-1)]=0.6#(0.1+ (iii-1)*(covhubcor/degrees))
       sigma[(covars+1+(i-1)*degrees+iii-1),i]= 0.6#(0.1+(iii-1)*(covhubcor/degrees))
   }
}
else {
sigma[i,(covars+1+(i-1)*degrees):((covars+1+(i-1)*degrees)+degrees-1)]=covhubcor
sigma[(covars+1+(i-1)*degrees):((covars+1+(i-1)*degrees)+degrees-1),i]=covhubcor
}

}

#intra-interaction starting index
intrastart=50+covars# covars*degrees+covars
newcovars=3
newdegrees=3
newcor=0.1
#if (starshape){
#for (j in 1:newcovars){
#sigma[intrastart+j,(intrastart+newcovars+1+(j-1)*newdegrees):((intrastart+newcovars+1+(j-1)*newdegrees)+newdegrees-1)]=newcor
#sigma[(intrastart+newcovars+1+(j-1)*newdegrees):((intrastart+newcovars+1+(j-1)*newdegrees)+ newdegrees-1),intrastart+j]=newcor
#}
#}
#a submatrix 51-55
sigma[(intrastart+1):(intrastart+5),(intrastart+1):(intrastart+5)]=matrix(runif(25,0.05,0.8),5,5)
# a submatrix 56-65
sigma[(intrastart+1+5):(intrastart+5+10),(intrastart+1+5):(intrastart+5+10)]=matrix(runif(100,0.05,0.8),10,10)
# a hub 66-70
sigma[(intrastart+1+5+10),(intrastart+1+5+10):(intrastart+5+5+10)]=runif(5,0.05,0.8)
sigma[(intrastart+1+5+10):(intrastart+5+5+10),(intrastart+1+5+10)]=runif(5,0.05,0.8)
# a hub 71-80
sigma[(intrastart+1+5+10+5),(intrastart+1+5+10+5):(intrastart+10+5+10+5)]=runif(10,0.05,0.8)
sigma[(intrastart+1+5+10+5):(intrastart+10+5+10+5),(intrastart+1+5+10+5)]=runif(10,0.05,0.8)


diag(sigma)=1
sigma=make.positive.definite(sigma)
out <- mvrnorm(samples, mu = rep(0,(nodes+covars)), Sigma = sigma,
               empirical = TRUE)
#print(var(out))
#print(max(var(out)))
#print(mean(var(out)))

if(binary){
if(dich){
for (i in 1:dim(out)[2]){out[,i]=(out[,i]>median(out[,i])+0)}
}
else{
m=sigma
diag(m)=1
out=rmvbin(samples, margprob = c(rep(0.5,(nodes+covars))), bincorr = m)
}
}

sigma=sigma/max(sigma)

print(max(sigma))
print(min(sigma))
# nodes, samples, covars,degrees, covhubcor
if(binary){
png(height=1200, width=1200, file=paste("binaryoverlap-n",nodes,"s",samples,"c",covars,"d",degrees,"cor",covhubcor,".png",sep=''))
write.table(sigma,paste("bincorr-n",nodes,"s",samples,"c",covars,"d",degrees,"cor",covhubcor,'nwc',numweakcor,".txt",sep=''))
corrplot(sigma,is.corr=F)
}
else{
png(height=1200, width=1200, file=paste("contioverlap-n",nodes,"s",samples,"c",covars,"d",degrees,"cor",covhubcor,".png",sep=''))
write.table(sigma,paste("gaucorr-n",nodes,"s",samples,"c",covars,"d",degrees,"cor",covhubcor,'nwc',numweakcor,".txt",sep=''))
corrplot(sigma,is.corr=F)
}

dev.off()
write.table(out, paste("rawdata",valmet,nodes,"s",samples,"c",covars,"d",degrees,"cor",covhubcor,'nwc',numweakcor,".txt",sep=''))
return(out)
}