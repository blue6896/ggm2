library(corpcor)
source('simudata.R')
source('simudataBIN.R')
source('simu100AIC.R')
source('simu100AICGGM.R')
source('tptn.R')
library(SparseM)
#nodes, samples, covars,degrees, covhubcor, normalcor, binary
newcovars=3
newdegrees=5

nodes=100
samples=110
covars=3
degrees=5
covhubcor=0.65
normalcor=0.01
alpha=1
valmet='AIC'
repits=100
nwc=0

# no-covariate-adjustment binary 
Xbin=simudata(nodes, samples, covars,degrees ,covhubcor,normalcor, T,T,F,nwc,F,valmet)
#print(cor(X))
Y=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AIC( Xbin, covars,nodes,samples,alpha,valmet,F)
if (counting ==1){
write.table(Ytemp,paste('NONCOVBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60vary.txt',sep=''))
}
print(dim(Y))
print(dim(Ytemp))
print(class(Y))
print(class(Ytemp))
Y= Y+Ytemp 
}
write.table((Y/repits),paste('NONCOVBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))

#  covariate- Binary   DICHO
#print(cor(X))
Ybin=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AIC( Xbin, covars,nodes,samples,alpha,valmet,T)
if (counting ==1){
write.table(Ytemp,paste('COVBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60vary.txt',sep=''))
}
Ybin=Ybin+Ytemp
}
write.table((Ybin/repits),paste('COVBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))

# DICH BIN AIC
Xbin=simudata(nodes, samples, covars,degrees ,covhubcor,normalcor, T,T,F,nwc,T,valmet)
#print(cor(X))
Y=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AIC( Xbin, covars,nodes,samples,alpha,valmet,F)
if (counting ==1){
write.table(Ytemp,paste('NONCOVDICHBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60vary.txt',sep=''))
}
print(dim(Y))
print(dim(Ytemp))
print(class(Y))
print(class(Ytemp))
Y= Y+Ytemp
}
write.table((Y/repits),paste('NONCOVDICHBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))

#  covariate- Binary   DICHO
#print(cor(X))
Ybin=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AIC( Xbin, covars,nodes,samples,alpha,valmet,T)
if (counting ==1){
write.table(Ytemp,paste('COVDICHBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60vary.txt',sep=''))
}
Ybin=Ybin+Ytemp
}
write.table((Ybin/repits),paste('COVDICHBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))


#non-covariate Gauss

XGGM=simudata(nodes, samples, covars,degrees ,covhubcor,normalcor, F,T,F,nwc,F,valmet)
#print(cor(X))
YGGM=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AICGGM( XGGM, covars,nodes,samples,alpha,valmet,F)
if(counting==1){
write.table(Ytemp,paste('NONCOVGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60vary.txt',sep=''))
}
YGGM= YGGM+Ytemp
}
write.table((YGGM/repits),paste('NONCOVGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))


# covariate-GGM 

#print(cor(X))
YGGM=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AICGGM( XGGM, covars,nodes,samples,alpha,valmet,T)
if(counting==1){
write.table(Ytemp,paste('COVGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60vary.txt',sep=''))
}
YGGM=YGGM+Ytemp
}
write.table((YGGM/repits),paste('COVGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))

# GGM without covariates at all 
XGGM= XGGM[,(covars+1):(covars+nodes)]
covars=0
YGGM=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AICGGM( XGGM, covars,nodes,samples,alpha,valmet,T)
if(counting==1){
write.table(Ytemp,paste('ONLYGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60vary.txt',sep=''))
}
YGGM=YGGM+Ytemp
}
write.table((YGGM/repits),paste('ONLYGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))

covars=3
valmet='BIC'

# no-covariate-adjustment binary  DICHO
Xbin=simudata(nodes, samples, covars,degrees ,covhubcor,normalcor, T,T,F,nwc,F,valmet)
#print(cor(X))
Y=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AIC( Xbin, covars,nodes,samples,alpha,valmet,F)
if (counting ==1){
write.table(Ytemp,paste('NONCOVBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyBIC.txt',sep=''))
}
print(dim(Y))
print(dim(Ytemp))
print(class(Y))
print(class(Ytemp))
Y= Y+Ytemp
}
write.table((Y/repits),paste('NONCOVBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20BIC.txt',sep=''))

#  covariate- Binary   DICHO
#print(cor(X))
Ybin=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AIC( Xbin, covars,nodes,samples,alpha,valmet,T)
if (counting ==1){
write.table(Ytemp,paste('COVBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyBIC.txt',sep=''))
}
Ybin=Ybin+Ytemp
}
write.table((Ybin/repits),paste('COVBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20BIC.txt',sep=''))


Xbin=simudata(nodes, samples, covars,degrees ,covhubcor,normalcor, T,T,F,nwc,T,valmet)
#print(cor(X))
Y=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AIC( Xbin, covars,nodes,samples,alpha,valmet,F)
if (counting ==1){
write.table(Ytemp,paste('NONCOVDICHBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyBIC.txt',sep=''))
}
print(dim(Y))
print(dim(Ytemp))
print(class(Y))
print(class(Ytemp))
Y= Y+Ytemp
}
write.table((Y/repits),paste('NONCOVDICHBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20BIC.txt',sep=''))

#  covariate- Binary   DICHO
#print(cor(X))
Ybin=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AIC( Xbin, covars,nodes,samples,alpha,valmet,T)
if (counting ==1){
write.table(Ytemp,paste('COVDICHBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyBIC.txt',sep=''))
}
Ybin=Ybin+Ytemp
}
write.table((Ybin/repits),paste('COVDICHBINn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20BIC.txt',sep=''))



#non-covariate Gauss

XGGM=simudata(nodes, samples, covars,degrees ,covhubcor,normalcor, F,T,F,nwc,F,valmet)
#print(cor(X))
YGGM=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AICGGM( XGGM, covars,nodes,samples,alpha,valmet,F)
if(counting==1){
write.table(Ytemp,paste('NONCOVGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyBIC.txt',sep=''))
}
YGGM= YGGM+Ytemp
}
write.table((YGGM/repits),paste('NONCOVGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20BIC.txt',sep=''))

# covariate-GGM

#print(cor(X))
YGGM=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AICGGM( XGGM, covars,nodes,samples,alpha,valmet,T)
if(counting==1){
write.table(Ytemp,paste('COVGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyBIC.txt',sep=''))
}
YGGM=YGGM+Ytemp
}
write.table((YGGM/repits),paste('COVGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20BIC.txt',sep=''))

# GGM without covariates at all
XGGM= XGGM[,(covars+1):(covars+nodes)]
covars=0
YGGM=matrix(0,nodes+covars, nodes)
Ytemp=matrix(0,nodes+covars,nodes)
for (counting in 1:repits) {
Ytemp= simu100AICGGM( XGGM, covars,nodes,samples,alpha,valmet,T)
if(counting==1){
write.table(Ytemp,paste('ONLYGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60vary.txt',sep=''))
}
YGGM=YGGM+Ytemp
}
write.table((YGGM/repits),paste('ONLYGGMn',nodes,'-s',samples,'-c',covars,'-d',degrees,'-nwc',nwc,'-corr60varyREP20.txt',sep=''))

result=c()
result=c(result,unlist(tptn(Y,covars, degrees,newcovars,newdegrees)))
print(result)
result=c(result, unlist(tptn(Ybin,covars, degrees,newcovars,newdegrees)))
print(result)
result=c(result,unlist(tptn(YGGM,covars, degrees,newcovars,newdegrees)))
print(result)
write(result,paste('TEST-n',nodes,'-s',samples,'-c',covars,'-d',degrees,'-corr60varyREP20.txt',sep=''))